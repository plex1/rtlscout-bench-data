//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Schedule (identical to the original, required for cycle accuracy):
//   product A[m][i]*B[i][j] is formed at compute edge E(m)+i+j (row skew i,
//   horizontal east flow j), partial sums flow south, column j's result is
//   de-skewed by N-1-j cycles so the whole result row leaves together,
//   2N-1 = 15 compute cycles after its A row was accepted.
//
// Datapath (v5):
//   * Partial sums are exact at PW = 2*DW + clog2(N) = 19 signed bits
//     (|sum of 8 products| <= 8*2^14), results sign-extended to 32b at the
//     outputs only.
//   * Each PE accumulates in REDUNDANT carry-save form (sum, carry) - the
//     per-cycle path is a 4-level 3:2 compressor tree with NO carry-propagate
//     adder in the loop; the single CPA per column sits on the output side.
//   * The Booth-4 recoding of the stationary weight depends only on the
//     weight, so each PE stores the 12 recoding bits (x1/x2/ng per digit)
//     and the LOAD CHAIN shifts recodings, not raw weights.  The per-cycle
//     magnitude select is a plain 2-AND + OR per bit:
//         m[b] = (x1 & a[b]) | (x2 & a[b-1])
//     followed by the ng inversion mask, then the compressor tree.
//#############################################################################

module tpu_pe_cs2
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // accumulator width (signed, exact)
   (
    input wire		 clk,
    input wire		 rst,	// synchronous, active high
    input wire		 ld,	// weight-load shift enable
    input wire		 en,	// compute enable
    input wire [11:0]	 enc_north, // Booth recoding in from PE above
    input wire [DW-1:0]	 a_west,
    input wire [PW-1:0]	 ps_north_s, // partial sum (carry-save sum row)
    input wire [PW-2:0]	 ps_north_c, // partial sum (carry row, bit0 = 0)
    output wire [11:0]	 enc_south,
    output reg  [DW-1:0] a_east,
    output reg  [PW-1:0] ps_south_s,
    output reg  [PW-2:0] ps_south_c
    );

   //#########################################################################
   // Stored radix-4 Booth recoding of the stationary weight w:
   //   digit k = -2*w[2k+1] + w[2k] + w[2k-1]  (k = 0..3, w[-1] = 0)
   // x1_k = |d|=1, x2_k = |d|=2 (mutually exclusive), ng_k = d negative.
   //#########################################################################
   reg [11:0] enc;
   assign enc_south = enc;

   wire [3:0] x1 = enc[11:8];
   wire [3:0] x2 = enc[7:4];
   wire [3:0] ng = enc[3:0];

   // sign-extended a and 2a (9-bit rows)
   wire [8:0] ae  = {a_west[DW-1], a_west};
   wire [8:0] a2e = {a_west[DW-1], a_west[DW-2:0], 1'b0};

   // per-digit magnitude row: m = (x1 & a) | (x2 & 2a);  p = m ^ ng mask
   wire [8:0] p0 = (({9{x1[0]}} & ae) | ({9{x2[0]}} & a2e)) ^ {9{ng[0]}};
   wire [8:0] p1 = (({9{x1[1]}} & ae) | ({9{x2[1]}} & a2e)) ^ {9{ng[1]}};
   wire [8:0] p2 = (({9{x1[2]}} & ae) | ({9{x2[2]}} & a2e)) ^ {9{ng[2]}};
   wire [8:0] p3 = (({9{x1[3]}} & ae) | ({9{x2[3]}} & a2e)) ^ {9{ng[3]}};

   // partial-product rows, sign-extended to PW bits, offsets 0/2/4/6
   wire [PW-1:0] r0 = {{(PW-9){p0[8]}}, p0};
   wire [PW-1:0] r1 = {{(PW-11){p1[8]}}, p1, 2'b00};
   wire [PW-1:0] r2 = {{(PW-13){p2[8]}}, p2, 4'b0000};
   wire [PW-1:0] r3 = {{(PW-15){p3[8]}}, p3, 6'b000000};
   // incoming partial sum, carry-save form
   wire [PW-1:0] r4 = ps_north_s;
   wire [PW-1:0] r5 = {ps_north_c, 1'b0};
   // negative-digit carry-in bits (positions 0,2,4,6 - mutually disjoint)
   wire [PW-1:0] r6 = {{(PW-7){1'b0}}, ng[3], 1'b0, ng[2], 1'b0, ng[1], 1'b0, ng[0]};

   //#########################################################################
   // 4-level 3:2 compressor tree: 7 rows -> 2 rows (mod 2^PW)
   //#########################################################################
   wire [PW-1:0] g1 = (r0 & r1) | (r2 & (r0 | r1));
   wire [PW-1:0] s1 = r0 ^ r1 ^ r2;
   wire [PW-1:0] c1 = {g1[PW-2:0], 1'b0};

   wire [PW-1:0] g2 = (r3 & r4) | (r5 & (r3 | r4));
   wire [PW-1:0] s2 = r3 ^ r4 ^ r5;
   wire [PW-1:0] c2 = {g2[PW-2:0], 1'b0};

   wire [PW-1:0] g3 = (s1 & c1) | (s2 & (s1 | c1));
   wire [PW-1:0] s3 = s1 ^ c1 ^ s2;
   wire [PW-1:0] c3 = {g3[PW-2:0], 1'b0};

   wire [PW-1:0] g4 = (s3 & c3) | (c2 & (s3 | c3));
   wire [PW-1:0] s4 = s3 ^ c3 ^ c2;
   wire [PW-1:0] c4 = {g4[PW-2:0], 1'b0};

   wire [PW-1:0] g5 = (s4 & c4) | (r6 & (s4 | c4));
   wire [PW-1:0] s5 = s4 ^ c4 ^ r6;
   wire [PW-1:0] c5 = {g5[PW-2:0], 1'b0};

   always @(posedge clk) begin
      if (rst) begin
	 enc	   <= 12'b0;
	 a_east	   <= {DW{1'b0}};
	 ps_south_s <= {PW{1'b0}};
	 ps_south_c <= {(PW-1){1'b0}};
      end
      else begin
	 if (ld)
	   enc <= enc_north;
	 if (en) begin
	    a_east	   <= a_west;
	    ps_south_s <= s5;
	    ps_south_c <= c5[PW-1:1];
	 end
      end
   end

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input wire		 clk,
    input wire		 rst,	 // synchronous, active high
    input wire		 w_valid,// weight-row valid (load phase)
    input wire [N*DW-1:0]	w_data, // weight row, packed { col N-1 ... col 0 }
    input wire		 a_valid,// activation-row valid (compute phase)
    input wire [N*DW-1:0]	a_data, // activation row, packed { col N-1 ... col 0 }
    output wire		 c_valid,// result-row valid
    output wire [N*ACCW-1:0]	c_data  // result row, packed { col N-1 ... col 0 }
    );

   localparam PW  = 2*DW + $clog2(N);
   localparam LAT = 2*N - 1;

   wire en = ~w_valid;

   //##########################################################################
   // Booth recoding of a raw weight byte (used only at the column tops).
   //##########################################################################
   function [11:0] booth;
      input [DW-1:0] w;
      reg [3:0] b0, b1, b2;
      integer   k;
      begin
	 for (k = 0; k < 4; k = k + 1) begin
	    b0 = (k == 0) ? 4'd0 : {4{w[2*k-1]}};
	    b1 = {4{w[2*k]}};
	    b2 = {4{w[2*k+1]}};
	    booth[8+k]  = b0[0] ^ b1[0];	 // x1_k = |d| = 1
	    booth[4+k]  = ~(b0[0] ^ b1[0]) & (b1[0] ^ b2[0]); // x2_k = |d| = 2
	    booth[k]    = b2[0];		 // ng_k = d negative
	 end
      end
   endfunction

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles.  Zero is
   // injected while a_valid is low so idle/drain cycles add nothing.
   //##########################################################################
   wire [N*DW-1:0] a_left;

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer	    k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < i; k = k + 1)
		    sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*DW +: DW] = sr[i-1];
	 end
      end
   endgenerate

   //##########################################################################
   // The systolic array: carry-save partial sums flow south, Booth recodings
   // shift south during load.
   //##########################################################################
   wire [PW-1:0]  ps_s [0:N][0:N-1];
   wire [PW-2:0]  ps_c [0:N][0:N-1];
   wire [11:0]	  enc_v [0:N][0:N-1];
   wire [DW-1:0]  a_h  [0:N-1][0:N];
   wire [PW-1:0]  res  [0:N-1];

   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign ps_s[0][j]  = {PW{1'b0}};
	 assign ps_c[0][j]  = {(PW-1){1'b0}};
	 assign enc_v[0][j] = booth(w_data[j*DW +: DW]);
	 // the only carry-propagate add per column: on the output side
	 assign res[j] = ps_s[N][j] + {ps_c[N][j], 1'b0};
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe_cs2 #(.DW(DW), .PW(PW)) u_pe
	      (.clk(clk), .rst(rst), .ld(w_valid), .en(en),
	       .enc_north(enc_v[i][j]), .a_west(a_h[i][j]),
	       .ps_north_s(ps_s[i][j]), .ps_north_c(ps_c[i][j]),
	       .enc_south(enc_v[i+1][j]), .a_east(a_h[i][j+1]),
	       .ps_south_s(ps_s[i+1][j]), .ps_south_c(ps_c[i+1][j]));
	 end
      end
   endgenerate

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles so all columns
   // of a result row line up.  Values are sign-extended to ACCW here.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*ACCW +: ACCW] =
	      {{(ACCW-PW){res[j][PW-1]}}, res[j]};
	 end
	 else begin : g_deskewn
	    reg [PW-1:0] dr [0:D-1];
	    integer	     k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1)
		    dr[k] <= {PW{1'b0}};
	       end
	       else if (en) begin
		  dr[0] <= res[j];
		  for (k = 1; k < D; k = k + 1)
		    dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*ACCW +: ACCW] =
	      {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
	 end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

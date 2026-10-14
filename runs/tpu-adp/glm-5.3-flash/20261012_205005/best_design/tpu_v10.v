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
// Datapath (v10):
//   * Partial sums are exact at PW = 2*DW + clog2(N) = 19 signed bits
//     (extreme vectors reach |C| = 8*128^2 = 131072); results are
//     sign-extended to 32 bits at the outputs only.
//   * Each PE accumulates in REDUNDANT carry-save form (sum, carry).
//   * The Booth-4 recoding of the stationary weight depends only on the
//     weight, so each PE stores the 12 recoding bits (x1/x2/ng per digit)
//     and the LOAD CHAIN shifts recodings, not raw weights.  The partial
//     product is one AOI22-style select plus the ng inversion:
//         m[b] = (x1 & a[b]) | (x2 & a[b-1]),   p = m ^ ng
//     (2 gate levels), and the +1 carry-ins of negative digits ride a sparse
//     row that enters the tree at level 2.
//   * 7 rows -> 2 rows in a 3-level 3:2 compressor tree (no carry-propagate
//     adder in the accumulate loop).
//   * The single carry-propagate add per column (output side) is a 19-bit
//     Kogge-Stone prefix adder, so it does not dominate the per-cycle path.
//   * Input skew and eastward activation hops are one shared per-row shift
//     line, tapped once per PE (no per-PE activation register).
//#############################################################################

module add19ks
  #(parameter W = 19)
   (
    input wire [W-1:0] a,
    input wire [W-1:0] b,
    output wire [W-1:0] sum
    );

   // Kogge-Stone prefix adder: 5 levels of span 1/2/4/8/16 (W = 19).
   wire [W-1:0] p0 = a ^ b;
   wire [W-1:0] g0 = a & b;

   wire [W-1:0] p1 = p0 & {p0[W-2:0], 1'b0};
   wire [W-1:0] g1 = g0 | (p0 & {g0[W-2:0], 1'b0});

   wire [W-1:0] p2 = p1 & {{2{1'b0}}, p1[W-1:2]};
   wire [W-1:0] g2 = g1 | (p1 & {{2{1'b0}}, g1[W-1:2]});

   wire [W-1:0] p3 = p2 & {{4{1'b0}}, p2[W-1:4]};
   wire [W-1:0] g3 = g2 | (p2 & {{4{1'b0}}, g2[W-1:4]});

   wire [W-1:0] p4 = p3 & {{8{1'b0}}, p3[W-1:8]};
   wire [W-1:0] g4 = g3 | (p3 & {{8{1'b0}}, g3[W-1:8]});

   wire [W-1:0] p5 = p4 & {{16{1'b0}}, p4[W-1:16]};
   wire [W-1:0] g5 = g4 | (p4 & {{16{1'b0}}, g4[W-1:16]});

   // carry into bit k is the span generate of bits [k-1:0]
   assign sum = p0 ^ {g5[W-1:1], 1'b0};

endmodule

module tpu_pe_cs3
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
   // 3-level 3:2 compressor tree: 7 rows -> 2 rows (mod 2^PW)
   //#########################################################################
   wire [PW-1:0] g1 = (r0 & r1) | (r2 & (r0 | r1));
   wire [PW-1:0] s1 = r0 ^ r1 ^ r2;
   wire [PW-1:0] c1 = {g1[PW-2:0], 1'b0};

   wire [PW-1:0] g2 = (r3 & r4) | (r5 & (r3 | r4));
   wire [PW-1:0] s2 = r3 ^ r4 ^ r5;
   wire [PW-1:0] c2 = {g2[PW-2:0], 1'b0};

   wire [PW-1:0] g3 = (s1 & s2) | (r6 & (s1 | s2));
   wire [PW-1:0] s3 = s1 ^ s2 ^ r6;
   wire [PW-1:0] c3 = {g3[PW-2:0], 1'b0};

   wire [PW-1:0] g4 = (c1 & c2) | (c3 & (c1 | c2));
   wire [PW-1:0] s4 = c1 ^ c2 ^ c3;
   wire [PW-1:0] c4 = {g4[PW-2:0], 1'b0};

   always @(posedge clk) begin
      if (rst) begin
	 enc	   <= 12'b0;
	 ps_south_s <= {PW{1'b0}};
	 ps_south_c <= {(PW-1){1'b0}};
      end
      else begin
	 if (ld)
	   enc <= enc_north;
	 if (en) begin
	    ps_south_s <= s4;
	    ps_south_c <= c4[PW-1:1];
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
   // Activation distribution: one shift line per row; PE(i,j) taps the value
   // delayed by i+j compute cycles (input skew i + eastward hops j).  Zero is
   // injected while a_valid is low so idle/drain cycles add nothing.
   //##########################################################################
   wire [DW-1:0] atap [0:N-1][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_aline
	 localparam DEPTH = i + N - 1;
	 reg [DW-1:0] sh [0:DEPTH-1];
	 integer	    k;
	 wire [DW-1:0]	    din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 always @(posedge clk) begin
	    if (rst) begin
	       for (k = 0; k < DEPTH; k = k + 1)
		 sh[k] <= {DW{1'b0}};
	    end
	    else if (en) begin
	       sh[0] <= din;
	       for (k = 1; k < DEPTH; k = k + 1)
		 sh[k] <= sh[k-1];
	    end
	 end
	 for (j = 0; j < N; j = j + 1) begin : g_tap
	    if (i + j == 0) begin : g_tap0
	       assign atap[i][0] = din;	 // only row 0, col 0 is undelayed
	    end
	    else begin : g_tapn
	       assign atap[i][j] = sh[i+j-1]; // delay i+j compute cycles
	    end
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
   wire [PW-1:0]  res  [0:N-1];

   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign ps_s[0][j]  = {PW{1'b0}};
	 assign ps_c[0][j]  = {(PW-1){1'b0}};
	 assign enc_v[0][j] = booth(w_data[j*DW +: DW]);
	 // the only carry-propagate add per column: on the output side
	 add19ks #(.W(PW)) u_cpa
	   (.a(ps_s[N][j]), .b({ps_c[N][j], 1'b0}), .sum(res[j]));
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe_cs3 #(.DW(DW), .PW(PW)) u_pe
	      (.clk(clk), .rst(rst), .ld(w_valid), .en(en),
	       .enc_north(enc_v[i][j]), .a_west(atap[i][j]),
	       .ps_north_s(ps_s[i][j]), .ps_north_c(ps_c[i][j]),
	       .enc_south(enc_v[i+1][j]),
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

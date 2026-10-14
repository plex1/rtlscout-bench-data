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
// Datapath (v3):
//   * Partial sums are exact at PW = 2*DW + clog2(N) = 19 signed bits
//     (|sum of 8 products| <= 8*2^14), results sign-extended to 32b at the
//     outputs only.
//   * Each PE accumulates in REDUNDANT carry-save form: two registers
//     (sum, carry).  The per-cycle path is Booth-4 encoding of the stationary
//     weight (4 partial products) compressed together with the incoming
//     (sum, carry) pair by a 4-level 3:2 compressor tree - no carry-propagate
//     adder in the loop.  The single CPA per column happens once, on the
//     output side (bottom row -> de-skew), off the accumulate critical path.
//#############################################################################

module tpu_pe_cs
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // accumulator width (signed, exact)
   (
    input wire		 clk,
    input wire		 rst,	// synchronous, active high
    input wire		 ld,	// weight-load shift enable
    input wire		 en,	// compute enable
    input wire [DW-1:0]	 w_north,
    input wire [DW-1:0]	 a_west,
    input wire [PW-1:0]	 ps_north_s, // partial sum (carry-save sum row)
    input wire [PW-2:0]	 ps_north_c, // partial sum (carry row, bit0 = 0)
    output wire [DW-1:0] w_south,
    output reg  [DW-1:0] a_east,
    output reg  [PW-1:0] ps_south_s,
    output reg  [PW-2:0] ps_south_c
    );

   reg [DW-1:0] w;
   assign w_south = w;

   //#########################################################################
   // Radix-4 Booth recoding of the stationary weight w:
   //   digit k = -2*w[2k+1] + w[2k] + w[2k-1]  (k = 0..3, w[-1] = 0)
   // with B0 = w[2k-1], B1 = w[2k], B2 = w[2k+1]:
   //   |d| = 1 iff B0^B1, |d| = 2 iff ~（B0^B1) & (B1^B2), d = 0 otherwise,
   //   sign(d) = B2.
   //#########################################################################
   wire x1_0 = w[0];
   wire x2_0 = w[0] ^ w[1];
   wire ng_0 = w[1];
   wire x1_1 = w[1] ^ w[2];
   wire x2_1 = w[2] ^ w[3];
   wire ng_1 = w[3];
   wire x1_2 = w[3] ^ w[4];
   wire x2_2 = w[4] ^ w[5];
   wire ng_2 = w[5];
   wire x1_3 = w[5] ^ w[6];
   wire x2_3 = w[6] ^ w[7];
   wire ng_3 = w[7];

   wire [8:0] mag0 = x1_0 ? {a_west[7], a_west} :
		     (x2_0 ? {a_west, 1'b0} : 9'b0);
   wire [8:0] mag1 = x1_1 ? {a_west[7], a_west} :
		     (x2_1 ? {a_west, 1'b0} : 9'b0);
   wire [8:0] mag2 = x1_2 ? {a_west[7], a_west} :
		     (x2_2 ? {a_west, 1'b0} : 9'b0);
   wire [8:0] mag3 = x1_3 ? {a_west[7], a_west} :
		     (x2_3 ? {a_west, 1'b0} : 9'b0);

   wire [8:0] pp0 = ng_0 ? ~mag0 : mag0;
   wire [8:0] pp1 = ng_1 ? ~mag1 : mag1;
   wire [8:0] pp2 = ng_2 ? ~mag2 : mag2;
   wire [8:0] pp3 = ng_3 ? ~mag3 : mag3;

   // partial-product rows, sign-extended to PW bits, offsets 0/2/4/6
   wire [PW-1:0] r0 = {{(PW-9){pp0[8]}}, pp0};
   wire [PW-1:0] r1 = {{(PW-11){pp1[8]}}, pp1, 2'b00};
   wire [PW-1:0] r2 = {{(PW-13){pp2[8]}}, pp2, 4'b0000};
   wire [PW-1:0] r3 = {{(PW-15){pp3[8]}}, pp3, 6'b000000};
   // incoming partial sum, carry-save form
   wire [PW-1:0] r4 = ps_north_s;
   wire [PW-1:0] r5 = {ps_north_c, 1'b0};
   // negative-digit carry-in bits (positions 0,2,4,6 - mutually disjoint)
   wire [PW-1:0] r6 = {{(PW-7){1'b0}}, ng_3, 1'b0, ng_2, 1'b0, ng_1, 1'b0, ng_0};

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
	 w	   <= {DW{1'b0}};
	 a_east	   <= {DW{1'b0}};
	 ps_south_s <= {PW{1'b0}};
	 ps_south_c <= {(PW-1){1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
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
   // The systolic array: carry-save partial sums flow south.
   //##########################################################################
   wire [PW-1:0]  ps_s [0:N][0:N-1];
   wire [PW-2:0]  ps_c [0:N][0:N-1];
   wire [DW-1:0]  a_h  [0:N-1][0:N];
   wire [DW-1:0]  w_v  [0:N][0:N-1];
   wire [PW-1:0]  res  [0:N-1];

   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign ps_s[0][j] = {PW{1'b0}};
	 assign ps_c[0][j] = {(PW-1){1'b0}};
	 assign w_v[0][j]  = w_data[j*DW +: DW];
	 // the only carry-propagate add per column: on the output side
	 assign res[j] = ps_s[N][j] + {ps_c[N][j], 1'b0};
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe_cs #(.DW(DW), .PW(PW)) u_pe
	      (.clk(clk), .rst(rst), .ld(w_valid), .en(en),
	       .w_north(w_v[i][j]), .a_west(a_h[i][j]),
	       .ps_north_s(ps_s[i][j]), .ps_north_c(ps_c[i][j]),
	       .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
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

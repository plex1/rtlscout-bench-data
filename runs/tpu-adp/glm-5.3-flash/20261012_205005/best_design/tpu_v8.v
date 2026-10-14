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
// Datapath (v8) - non-Booth signed multiply, carry-save accumulation:
//   * Partial sums are exact at PW = 2*DW + clog2(N) = 19 signed bits
//     (extreme vectors reach |C| = 8*128^2 = 131072); results are
//     sign-extended to 32 bits at the outputs only.
//   * Signed product via sign-bit correction (verified on corner cases):
//        a*w = Al*Wl - 2^7*a7*Wl - 2^7*w7*Al + 2^14*a7*w7
//     with Al/Wl the low 7 bits (unsigned) and a7/w7 the sign bits.
//     Rows: 7 AND rows (a_i & w[6:0] << i), two negated sign rows
//     ({5{x7}, x7&~y[6:0]} << 7), one 3-bit correction row
//     {[14]: a7w7, [8]: a7w7, [7]: a7^w7}, plus the incoming carry-save pair:
//     12 rows in total -> 10 CSAs in 5 levels.  The per-cycle path is ONE
//     and/nand level feeding the compressor tree - no magnitude select, no
//     Booth encode, no carry-propagate adder in the loop; the single CPA per
//     column sits on the output side.
//   * Input skew and eastward activation hops are one shared per-row shift
//     line, tapped once per PE (no per-PE activation register).
//#############################################################################

module tpu_pe_nb
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
    output reg  [PW-1:0] ps_south_s,
    output reg  [PW-2:0] ps_south_c
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire	a7 = a_west[DW-1];
   wire	w7 = w[DW-1];
   wire [DW-2:0] al = a_west[DW-2:0];
   wire [DW-2:0] wl = w[DW-2:0];

   // 7 unsigned product rows: a_i * w[6:0] placed at bit i
   wire [PW-1:0] q0 = {{(PW-7){1'b0}},  {7{a_west[0]}} & wl};
   wire [PW-1:0] q1 = {{(PW-8){1'b0}},  {7{a_west[1]}} & wl, 1'b0};
   wire [PW-1:0] q2 = {{(PW-9){1'b0}},  {7{a_west[2]}} & wl, 2'b0};
   wire [PW-1:0] q3 = {{(PW-10){1'b0}}, {7{a_west[3]}} & wl, 3'b0};
   wire [PW-1:0] q4 = {{(PW-11){1'b0}}, {7{a_west[4]}} & wl, 4'b0};
   wire [PW-1:0] q5 = {{(PW-12){1'b0}}, {7{a_west[5]}} & wl, 5'b0};
   wire [PW-1:0] q6 = {{(PW-13){1'b0}}, {7{a_west[6]}} & wl, 6'b0};

   // negated sign-bit rows: -2^7*a7*w[6:0] and -2^7*w7*a[6:0] (mod 2^PW)
   // payload at bits [14:8], sign-extended to bit 18: -X (X in [0,2^7-1])
   // at position 2^7 equals (2^7 - X) at 2^7 - 2^7 = (2^14 - 2^7*X) mod 2^19
   wire [PW-1:0] rna = {{4{a7}}, a7 & ~wl, 8'b0};
   wire [PW-1:0] rnw = {{4{w7}}, w7 & ~al, 8'b0};
   // +1 carry-in of both complements (a7+w7 at bits 7/8) + the 2^14*a7*w7 term
   wire [PW-1:0] rmx = {{(PW-15){1'b0}}, a7 & w7, 5'b0, a7 & w7, a7 ^ w7, 7'b0};

   // incoming partial sum, carry-save form
   wire [PW-1:0] rps = ps_north_s;
   wire [PW-1:0] rpc = {ps_north_c, 1'b0};

   //#########################################################################
   // 5-level 3:2 compressor tree: 12 rows -> 2 rows (mod 2^PW)
   //#########################################################################
   wire [PW-1:0] g1 = (q0 & q1) | (q2 & (q0 | q1));
   wire [PW-1:0] s1 = q0 ^ q1 ^ q2;
   wire [PW-1:0] c1 = {g1[PW-2:0], 1'b0};

   wire [PW-1:0] g2 = (q3 & q4) | (q5 & (q3 | q4));
   wire [PW-1:0] s2 = q3 ^ q4 ^ q5;
   wire [PW-1:0] c2 = {g2[PW-2:0], 1'b0};

   wire [PW-1:0] g3 = (rna & rnw) | (rmx & (rna | rnw));
   wire [PW-1:0] s3 = rna ^ rnw ^ rmx;
   wire [PW-1:0] c3 = {g3[PW-2:0], 1'b0};

   wire [PW-1:0] g4 = (q6 & rps) | (rpc & (q6 | rps));
   wire [PW-1:0] s4 = q6 ^ rps ^ rpc;
   wire [PW-1:0] c4 = {g4[PW-2:0], 1'b0};

   wire [PW-1:0] g5 = (s1 & s2) | (s3 & (s1 | s2));
   wire [PW-1:0] s5 = s1 ^ s2 ^ s3;
   wire [PW-1:0] c5 = {g5[PW-2:0], 1'b0};

   wire [PW-1:0] g6 = (c1 & c2) | (c3 & (c1 | c2));
   wire [PW-1:0] s6 = c1 ^ c2 ^ c3;
   wire [PW-1:0] c6 = {g6[PW-2:0], 1'b0};

   wire [PW-1:0] g7 = (s5 & s6) | (s4 & (s5 | s6));
   wire [PW-1:0] s7 = s5 ^ s6 ^ s4;
   wire [PW-1:0] c7 = {g7[PW-2:0], 1'b0};

   wire [PW-1:0] g8 = (c5 & c6) | (c4 & (c5 | c6));
   wire [PW-1:0] s8 = c5 ^ c6 ^ c4;
   wire [PW-1:0] c8 = {g8[PW-2:0], 1'b0};

   wire [PW-1:0] g9 = (s7 & s8) | (c7 & (s7 | s8));
   wire [PW-1:0] s9 = s7 ^ s8 ^ c7;
   wire [PW-1:0] c9 = {g9[PW-2:0], 1'b0};

   wire [PW-1:0] gA = (s9 & c9) | (c8 & (s9 | c9));
   wire [PW-1:0] sA = s9 ^ c9 ^ c8;
   wire [PW-1:0] cA = {gA[PW-2:0], 1'b0};

   always @(posedge clk) begin
      if (rst) begin
	 w	   <= {DW{1'b0}};
	 ps_south_s <= {PW{1'b0}};
	 ps_south_c <= {(PW-1){1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    ps_south_s <= sA;
	    ps_south_c <= cA[PW-1:1];
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
   // The systolic array: carry-save partial sums flow south, raw weight bytes
   // shift south during load.
   //##########################################################################
   wire [PW-1:0]  ps_s [0:N][0:N-1];
   wire [PW-2:0]  ps_c [0:N][0:N-1];
   wire [DW-1:0]  w_v  [0:N][0:N-1];
   wire [PW-1:0]  res  [0:N-1];

   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign ps_s[0][j] = {PW{1'b0}};
	 assign ps_c[0][j] = {(PW-1){1'b0}};
	 assign w_v[0][j]  = w_data[j*DW +: DW];
	 // the only carry-propagate add per column: on the output side
	 assign res[j] = ps_s[N][j] + {ps_c[N][j], 1'b0};
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe_nb #(.DW(DW), .PW(PW)) u_pe
	      (.clk(clk), .rst(rst), .ld(w_valid), .en(en),
	       .w_north(w_v[i][j]), .a_west(atap[i][j]),
	       .ps_north_s(ps_s[i][j]), .ps_north_c(ps_c[i][j]),
	       .w_south(w_v[i+1][j]),
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

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
// Datapath (v9) - carry-save accumulation, no carry-propagate adder in the
// loop, no Booth select chain:
//   * Partial sums are exact at PW = 2*DW + clog2(N) = 19 signed bits
//     (extreme vectors reach |C| = 8*128^2 = 131072); results are
//     sign-extended to 32 bits at the outputs only.
//   * Signed product by sign-bit correction, verified on corner cases:
//        a*w = Al*Wl - 2^7 * corr,   corr = a7*W_ext + w7*A_ext  (< 2^9)
//     with Al/A_ext the low-7 / full unsigned patterns.  Rows: 7 dense AND
//     rows (a_i & w[6:0]) << i, one dense correction row {(2^12-corr, 7'0)}
//     (exactly PW bits), plus the incoming carry-save pair = 10 rows ->
//     4 levels of 3:2 compressors.  The per-cycle front end is a single
//     AND level (the correction row's negate runs in parallel, 2 levels).
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

   // unsigned patterns
   wire [8:0] a_ext = {a7, a_west[DW-2:0]};
   wire [8:0] w_ext = {w7, w[DW-2:0]};

   // 7 unsigned product rows: a_i * w[6:0] placed at bit i
   wire [PW-1:0] q0 = {{(PW-7){1'b0}},  {7{a_west[0]}} & w[DW-2:0]};
   wire [PW-1:0] q1 = {{(PW-8){1'b0}},  {7{a_west[1]}} & w[DW-2:0], 1'b0};
   wire [PW-1:0] q2 = {{(PW-9){1'b0}},  {7{a_west[2]}} & w[DW-2:0], 2'b0};
   wire [PW-1:0] q3 = {{(PW-10){1'b0}}, {7{a_west[3]}} & w[DW-2:0], 3'b0};
   wire [PW-1:0] q4 = {{(PW-11){1'b0}}, {7{a_west[4]}} & w[DW-2:0], 4'b0};
   wire [PW-1:0] q5 = {{(PW-12){1'b0}}, {7{a_west[5]}} & w[DW-2:0], 5'b0};
   wire [PW-1:0] q6 = {{(PW-13){1'b0}}, {7{a_west[6]}} & w[DW-2:0], 6'b0};

   // correction row: -2^7 * corr (mod 2^PW), corr = a7*W_ext + w7*A_ext
   wire [8:0]  corr  = ({9{a7}} & w_ext) + ({9{w7}} & a_ext);
   wire [11:0] ncorr = 12'd0 - corr;
   wire [PW-1:0] rcorr = {ncorr, 7'b0};	 // exactly PW bits

   // incoming partial sum, carry-save form
   wire [PW-1:0] rps = ps_north_s;
   wire [PW-1:0] rpc = {ps_north_c, 1'b0};

   //#########################################################################
   // 4-level 3:2 compressor tree: 10 rows -> 2 rows (mod 2^PW)
   //#########################################################################
   wire [PW-1:0] g1 = (q0 & q1) | (q2 & (q0 | q1));
   wire [PW-1:0] s1 = q0 ^ q1 ^ q2;
   wire [PW-1:0] c1 = {g1[PW-2:0], 1'b0};

   wire [PW-1:0] g2 = (q3 & q4) | (q5 & (q3 | q4));
   wire [PW-1:0] s2 = q3 ^ q4 ^ q5;
   wire [PW-1:0] c2 = {g2[PW-2:0], 1'b0};

   wire [PW-1:0] g3 = (q6 & rcorr) | (rps & (q6 | rcorr));
   wire [PW-1:0] s3 = q6 ^ rcorr ^ rps;
   wire [PW-1:0] c3 = {g3[PW-2:0], 1'b0};

   wire [PW-1:0] g4 = (s1 & s2) | (s3 & (s1 | s2));
   wire [PW-1:0] s4 = s1 ^ s2 ^ s3;
   wire [PW-1:0] c4 = {g4[PW-2:0], 1'b0};

   wire [PW-1:0] g5 = (c1 & c2) | (c3 & (c1 | c2));
   wire [PW-1:0] s5 = c1 ^ c2 ^ c3;
   wire [PW-1:0] c5 = {g5[PW-2:0], 1'b0};

   wire [PW-1:0] g6 = (s4 & s5) | (rpc & (s4 | s5));
   wire [PW-1:0] s6 = s4 ^ s5 ^ rpc;
   wire [PW-1:0] c6 = {g6[PW-2:0], 1'b0};

   wire [PW-1:0] g7 = (c4 & c5) | (c6 & (c4 | c5));
   wire [PW-1:0] s7 = c4 ^ c5 ^ c6;
   wire [PW-1:0] c7 = {g7[PW-2:0], 1'b0};

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
	    ps_south_s <= s7;
	    ps_south_c <= c7[PW-1:1];
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

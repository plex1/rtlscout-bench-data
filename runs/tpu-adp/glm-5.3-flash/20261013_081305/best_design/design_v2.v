//#############################################################################
// Optimized weight-stationary systolic matrix-multiply tile (TPU MXU style).
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Architecture notes (vs. original):
//  * Internal accumulator width PW=19: the exact sum of eight 8x8 signed
//    products is bounded by 8*2^14 = 2^17, which fits in 19-bit two's
//    complement. Results are sign-extended to 32 bits at the output, so the
//    32-bit contract is preserved while every adder shrinks.
//  * Each row's activation is BROADCAST to all columns; row i is delayed by
//    i+K (K=7) compute cycles. PE(i,j) at compute cycle t multiplies
//    A[t-i-K][i] * B[i][j]; the psum chain then delivers dot(A[t-8-K],B) at
//    the bottom of every column simultaneously at t = m+15 = 2N-1 after row m
//    was accepted. No eastward activation registers, no output de-skew.
//  * Weights load exactly as the original: a per-column southward shift
//    chain fed by w_data; cycle k carries B[N-1-k][*].
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // internal accumulator width (signed)
   (
    input		  clk,
    input		  rst,	// synchronous, active high
    input		  ld,	// weight-load shift enable
    input		  en,	// compute enable
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_row,      // broadcast activation for this row
    input [PW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [PW-1:0]	  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_row) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en)
	   psum_south <= psum_north + {{(PW-2*DW){prod[2*DW-1]}}, prod};
      end
   end

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;   // 15 compute cycles injection-to-output
   localparam PW  = 19;	 // exact accumulator width
   localparam K   = N;	 // uniform input delay aligning all columns
			 // D_i = i + N: row i emerges from the psum chain
			 // (7-i) cycles after its product is formed, so the
			 // bottom psum at compute time t is exactly
			 // dot(A[t-2N+1], B) = C[t-15][*].

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i delayed by (i + K) compute cycles, zeros injected while
   // a_valid is low. The delayed row is broadcast to all N columns.
   //##########################################################################
   wire [DW-1:0] a_row [0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 localparam DEPTH = i + K;
	 reg [DW-1:0] sh [0:DEPTH-1];
	 integer k;
	 wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
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
	 assign a_row[i] = sh[DEPTH-1];
      end
   endgenerate

   //##########################################################################
   // Weight-stationary array: psum flows south; activations are broadcast.
   //##########################################################################
   wire [PW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign psum_v[0][j] = {PW{1'b0}};
	 assign w_v[0][j]    = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe #(.DW(DW), .PW(PW)) u_pe
	      (.clk	     (clk),
	       .rst	     (rst),
	       .ld	     (ld),
	       .en	     (en),
	       .w_north	     (w_v[i][j]),
	       .a_row	     (a_row[i]),
	       .psum_north   (psum_v[i][j]),
	       .w_south	     (w_v[i+1][j]),
	       .psum_south   (psum_v[i+1][j]));
	 end
      end
   endgenerate

   //##########################################################################
   // Output: all columns are time-aligned; sign-extend PW -> ACCW.
   //##########################################################################
   genvar jj;
   generate
      for (jj = 0; jj < N; jj = jj + 1) begin : g_out
	 assign c_data[jj*ACCW +: ACCW] =
	   {{(ACCW-PW){psum_v[N][jj][PW-1]}}, psum_v[N][jj]};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
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

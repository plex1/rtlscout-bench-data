//#############################################################################
// Weight-stationary 8x8 matmul tile (TPU-style MXU), cycle-accurate rewrite.
//
// Architecture: instead of a systolic PE mesh, each output column j gets a
// dot-product unit: C[m][j] = sum_i A[m][i] * B[i][j]. All 64 8x8 multipliers
// are kept (mandatory: one full result row per compute cycle). Narrow adder
// trees (sum of 8 signed 16-bit products fits in 19 bits) replace the 32-bit
// psum chain, and a 64-bit-wide activation delay line provides the exact
// 2N-1 = 15 compute-cycle latency. Rows are processed independently, so no
// zero injection is needed; rows whose a_valid was low simply produce
// unchecked (c_valid low) results.
//
// Timing: a_data of row m is visible during cycle m. The delay chain (14
// stages, advanced only on en = ~w_valid) presents it during cycle m+14, the
// product/sum result register lands it during cycle m+15, matching c_valid
// (a_valid delayed 15 en cycles). Weight load replicates the original
// per-column shift chain exactly (bottom matrix row first).
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (input	    clk,
    input	    rst,
    input	    w_valid,
    input [N*DW-1:0] w_data,
    input	    a_valid,
    input [N*DW-1:0] a_data,
    output	    c_valid,
    output [N*ACCW-1:0] c_data);

   localparam LAT = 2*N - 1;   // 15
   localparam D   = LAT - 1;   // 14 delay stages before the result register

   wire en = ~w_valid;

   //##########################################################################
   // Activation row delay chain (64 bits/stage, advances only on en).
   //##########################################################################
   reg [N*DW-1:0] ad [0:D-1];
   integer s;
   always @(posedge clk) begin
      if (rst) begin
	 for (s = 0; s < D; s = s + 1)
	   ad[s] <= {N*DW{1'b0}};
      end
      else if (en) begin
	 ad[0] <= a_data;
	 for (s = 1; s < D; s = s + 1)
	   ad[s] <= ad[s-1];
      end
   end

   //##########################################################################
   // Stationary weights: per-column shift chain, exactly like the original
   // array (w_top enters stage 0, shifts south on every load cycle).
   // After 8 load cycles stage i of column j holds B[i][j].
   //##########################################################################
   reg [DW-1:0] b [0:N-1][0:N-1]; // b[j][i]
   genvar gj, gi;
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_wcol
	 always @(posedge clk) begin
	    if (rst)
	      b[gj][0] <= {DW{1'b0}};
	    else if (w_valid)
	      b[gj][0] <= w_data[gj*DW +: DW];
	 end
	 for (gi = 1; gi < N; gi = gi + 1) begin : g_wstage
	    always @(posedge clk) begin
	       if (rst)
		 b[gj][gi] <= {DW{1'b0}};
	       else if (w_valid)
		 b[gj][gi] <= b[gj][gi-1];
	    end
	 end
      end
   endgenerate

   //##########################################################################
   // Dot product per column: balanced adder tree over 8 signed products,
   // registered once. Internal widths kept minimal (19-bit sums).
   //##########################################################################
   genvar j;
   wire [ACCW-1:0] sum [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
	 wire signed [2*DW-1:0] p0 = $signed(ad[D-1][0*DW +: DW]) * $signed(b[j][0]);
	 wire signed [2*DW-1:0] p1 = $signed(ad[D-1][1*DW +: DW]) * $signed(b[j][1]);
	 wire signed [2*DW-1:0] p2 = $signed(ad[D-1][2*DW +: DW]) * $signed(b[j][2]);
	 wire signed [2*DW-1:0] p3 = $signed(ad[D-1][3*DW +: DW]) * $signed(b[j][3]);
	 wire signed [2*DW-1:0] p4 = $signed(ad[D-1][4*DW +: DW]) * $signed(b[j][4]);
	 wire signed [2*DW-1:0] p5 = $signed(ad[D-1][5*DW +: DW]) * $signed(b[j][5]);
	 wire signed [2*DW-1:0] p6 = $signed(ad[D-1][6*DW +: DW]) * $signed(b[j][6]);
	 wire signed [2*DW-1:0] p7 = $signed(ad[D-1][7*DW +: DW]) * $signed(b[j][7]);

	 wire signed [2*DW:0]   q0 = p0 + p1;
	 wire signed [2*DW:0]   q1 = p2 + p3;
	 wire signed [2*DW:0]   q2 = p4 + p5;
	 wire signed [2*DW:0]   q3 = p6 + p7;

	 wire signed [2*DW+1:0] r0 = q0 + q1;
	 wire signed [2*DW+1:0] r1 = q2 + q3;

	 wire signed [2*DW+2:0] t  = r0 + r1;   // 19 bits: enough for 8 products

	 reg [2*DW+2:0] res;
	 always @(posedge clk) begin
	    if (rst)
	      res <= {(2*DW+3){1'b0}};
	    else if (en)
	      res <= t;
	 end
	 assign sum[j] = {{(ACCW-(2*DW+3)){res[2*DW+2]}}, res};
      end
   endgenerate

   //##########################################################################
   // Output packing and valid pipeline (a_valid delayed LAT en cycles).
   //##########################################################################
   genvar k;
   generate
      for (k = 0; k < N; k = k + 1) begin : g_out
	 assign c_data[k*ACCW +: ACCW] = sum[k];
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   integer v;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

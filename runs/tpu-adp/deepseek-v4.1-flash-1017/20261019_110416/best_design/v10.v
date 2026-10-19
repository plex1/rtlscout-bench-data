// v10: TPU tile with multiply and accumulate in SEPARATE pipeline stages.
//
// Golden per-PE behaviour (a reload freezes every partial sum):
//   P(i,j) at end of cycle t+1 = sum_{k<=i} A(t-i-j)[k] * B[k][j]
// so for an A row accepted in cycle T the completed column value must be
// visible in cycle T+8+j, then de-skewed by (7-j) -> output in cycle T+15.
//
// Instead of a systolic mesh the activation is fanned out from per-row 8-bit
// shift registers and the weights sit in plain registers, but the exact
// staging is kept:
//   tap(i,j)  = A row byte i delayed by (i+j-1) compute cycles
//   pr(i,j)   = registered product tap(i,j)*B[i][j]     (1 stage)
//   acc(i,j)  = acc(i-1,j) + pr(i,j)                    (8 stage cascade)
//   de-skew   = (7-j) register stages
// Total = (i+j-1) + 1 + (8-i) + (7-j) = 15 stages for every (i,j).
//
// The product registers are NOT gated by the compute enable: during a weight
// reload the activation taps are frozen while the new weights shift in, so the
// last product formed before compute resumes is exactly (frozen A)*(new B),
// which is what the golden array does (a PE keeps its a_west across the stall).
// The (0,0) element has tap delay -1: its product is combinational, matching
// the original where PE(0,0) sees a_data directly and therefore picks up the
// new A row on the first cycle after a reload.
//
// Every register-to-register path is a single 8x8 signed multiply or a single
// narrow add.

module tpu
  #(parameter N = 8, DW = 8, ACCW = 32)
   (
    input		 clk,
    input		 rst,
    input		 w_valid,
    input [N*DW-1:0]	 w_data,
    input		 a_valid,
    input [N*DW-1:0]	 a_data,
    output		 c_valid,
    output [N*ACCW-1:0]	 c_data
    );

   localparam LAT = 2*N - 1;	// 15

   wire		      ld = w_valid;
   wire		      en = ~w_valid;

   genvar	      i, j;
   integer	      kk;

   //--------------------------------------------------------------------------
   // Activation delay lines: sr[i][d] = A row byte i delayed by d+1 cycles.
   //--------------------------------------------------------------------------
   wire [N*DW-1:0] din = a_valid ? a_data : {N*DW{1'b0}};
   wire [DW-1:0]   sr [0:N-1][0:N+4];

   generate
      for (i = 0; i < N; i = i + 1) begin : g_dly
         reg [DW-1:0] q [0:i+N-2];
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < i + N - 1; kk = kk + 1)
                 q[kk] <= {DW{1'b0}};
            end
            else if (en) begin
               q[0] <= din[i*DW +: DW];
               for (kk = 1; kk < i + N - 1; kk = kk + 1)
                 q[kk] <= q[kk-1];
            end
         end
         for (j = 0; j < N; j = j + 1) begin : g_t
            if (i + j == 0) begin : g_comb
               assign sr[i][j] = din[i*DW +: DW];
            end
            else if (i + j == 1) begin : g_zero
               assign sr[i][j] = din[i*DW +: DW];
            end
            else begin : g_reg
               assign sr[i][j] = q[i+j-2];
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Stationary weights (shift south, bottom matrix row first).
   //--------------------------------------------------------------------------
   reg [DW-1:0] wreg [0:N-1][0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_wload
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1)
                 wreg[kk][j] <= {DW{1'b0}};
            end
            else if (ld) begin
               wreg[0][j] <= w_data[j*DW +: DW];
               for (kk = 1; kk < N; kk = kk + 1)
                 wreg[kk][j] <= wreg[kk-1][j];
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Per column: product registers, growing-width accumulate cascade, de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [15:0] pc0 = $signed(din[0 +: DW]) * $signed(wreg[0][0]);

         reg signed [15:0] pr [0:N-1];
         integer	      m;
         always @(posedge clk) begin
            if (rst) begin
               for (m = 0; m < N; m = m + 1)
                 pr[m] <= 16'sd0;
            end
            else begin
               if (!(j == 0)) begin : dummy_never
               end
               for (m = 1; m < N; m = m + 1)
                 pr[m] <= $signed(sr[m][j]) * $signed(wreg[m][j]);
               if (j != 0)
                 pr[0] <= $signed(sr[0][j]) * $signed(wreg[0][j]);
            end
         end

         // accumulate cascade (gated by the compute enable)
         reg signed [15:0] a0;
         reg signed [16:0] a1, a2;
         reg signed [17:0] a3, a4, a5, a6;
         reg signed [18:0] a7;

         always @(posedge clk) begin
            if (rst) begin
               a0 <= 16'sd0;
               a1 <= 17'sd0; a2 <= 17'sd0;
               a3 <= 18'sd0; a4 <= 18'sd0; a5 <= 18'sd0; a6 <= 18'sd0;
               a7 <= 19'sd0;
            end
            else if (en) begin
               if (j == 0) a0 <= pc0;
               else        a0 <= pr[0];
               a1 <= a0           + pr[1];
               a2 <= $signed(a1)  + pr[2];
               a3 <= $signed(a2)  + pr[3];
               a4 <= $signed(a3)  + pr[4];
               a5 <= $signed(a4)  + pr[5];
               a6 <= $signed(a5)  + pr[6];
               a7 <= $signed(a6)  + pr[7];
            end
         end

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-19){a7[18]}}, a7};
         end
         else begin : g_dsn
            reg signed [18:0] ds [0:D-1];
            always @(posedge clk) begin
               if (rst) begin
                  for (m = 0; m < D; m = m + 1)
                    ds[m] <= 19'sd0;
               end
               else if (en) begin
                  ds[0] <= a7;
                  for (m = 1; m < D; m = m + 1)
                    ds[m] <= ds[m-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-19){ds[D-1][18]}}, ds[D-1]};
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Valid pipeline.
   //--------------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// v8: multiply and accumulate in separate pipeline stages.
//
// Golden per-PE staging:  P(i,j)(T) = P(i-1,j)(T-1) + A(T-i-j)[i]*B[i][j].
// Here the products live in their own registers (taps delayed by i+j-1) so the
// multiply and the add never share a cycle.  The product registers are NOT
// gated by the compute enable: during a weight reload the activation taps are
// frozen while the new weights shift in, so the last product formed before
// compute resumes is exactly (frozen A) * (new weight) -- matching the golden
// array where a stalled PE keeps its a_west and picks up the reloaded weight.

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
   localparam AW  = 19;
   localparam MD  = 2*N - 3;	// deepest activation tap (13)

   wire		      ld = w_valid;
   wire		      en = ~w_valid;

   genvar	      i, j;
   integer	      kk;

   //--------------------------------------------------------------------------
   // Activation delay line: dd[d] = A delayed by d compute cycles.
   //--------------------------------------------------------------------------
   wire [N*DW-1:0] din = a_valid ? a_data : {N*DW{1'b0}};
   wire [N*DW-1:0] dd [0:MD];

   assign dd[0] = din;
   generate
      for (i = 1; i <= MD; i = i + 1) begin : g_dly
         reg [N*DW-1:0] q;
         always @(posedge clk) begin
            if (rst)     q <= {N*DW{1'b0}};
            else if (en) q <= dd[i-1];
         end
         assign dd[i] = q;
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
   // Per column: product registers + prefix-sum cascade + de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         reg signed [AW-1:0] acc [0:N-1];
         reg signed [AW-1:0] pr  [0:N-1];

         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1) begin
                  acc[kk] <= {AW{1'b0}};
                  pr[kk]  <= {AW{1'b0}};
               end
            end
            else begin
               // ungated product registers (row 0 handled in the cascade)
               for (kk = 1; kk < N; kk = kk + 1)
                 pr[kk] <= $signed(dd[kk+j-1][kk*DW +: DW]) * $signed(wreg[kk][j]);
               if (en) begin
                  // acc[0] is itself a product register (row 0 of the cascade)
                  acc[0] <= $signed(dd[j][0 +: DW]) * $signed(wreg[0][j]);
                  for (kk = 1; kk < N; kk = kk + 1)
                    acc[kk] <= acc[kk-1] + pr[kk];
               end
            end
         end

         // acc[0] must be built from the tap delayed by (j-1); for j == 0 the
         // tap is the combinational input, so dd[0] is used above (j=0 -> dd[0]).
         localparam D = N - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){acc[N-1][AW-1]}}, acc[N-1]};
         end
         else begin : g_dsn
            reg signed [AW-1:0] ds [0:D-1];
            integer	      m;
            always @(posedge clk) begin
               if (rst) begin
                  for (m = 0; m < D; m = m + 1)
                    ds[m] <= {AW{1'b0}};
               end
               else if (en) begin
                  ds[0] <= acc[N-1];
                  for (m = 1; m < D; m = m + 1)
                    ds[m] <= ds[m-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){ds[D-1][AW-1]}}, ds[D-1]};
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

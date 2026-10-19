// v5: broadcast-weight TPU tile with exact reload semantics.
//
// Reference behavior: P(i+1,j,t+1) = P(i,j,t) + a(t-i-j)[i]*w(i,j).  A weight
// reload freezes every partial sum, so the 8-stage accumulation chain must be
// preserved bit-exactly.  Here the weights are held in plain registers (no
// south-shifting mesh) and the activation skew is built with cheap 8-bit shift
// registers instead of 8-bit registers travelling through 64 PEs.

module tpu
  #(parameter N = 8, DW = 8, AW = 19, ACCW = 32)
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
   // Activation skew: a_tap[i][j] = A row byte i delayed by (i+j) compute cycles.
   //--------------------------------------------------------------------------
   wire [DW-1:0] a_tap [0:N-1][0:N-1];
   generate
      for (i = 0; i < N; i = i + 1) begin : g_adly
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i + N - 1 == 0) begin : g_none
            assign a_tap[i][0] = din;
         end
         else begin : g_sr
            reg [DW-1:0] sr [0:i+N-2];
            always @(posedge clk) begin
               if (rst) begin
                  for (kk = 0; kk < i + N - 1; kk = kk + 1)
                    sr[kk] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (kk = 1; kk < i + N - 1; kk = kk + 1)
                    sr[kk] <= sr[kk-1];
               end
            end
            for (j = 0; j < N; j = j + 1) begin : g_tap
               if (i + j == 0) begin : g_t0
                  assign a_tap[i][j] = din;
               end
               else begin : g_tn
                  assign a_tap[i][j] = sr[i+j-1];
               end
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Stationary weights: B[i][j] loaded bottom row first, shifted south.
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
   // 8-stage accumulation chain per column + output de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [AW-1:0] prod [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_prod
            assign prod[i] = $signed(a_tap[i][j]) * $signed(wreg[i][j]);
         end

         reg signed [AW-1:0] acc [0:N-1];
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1)
                 acc[kk] <= {AW{1'b0}};
            end
            else if (en) begin
               acc[0] <= prod[0];
               for (kk = 1; kk < N; kk = kk + 1)
                 acc[kk] <= acc[kk-1] + prod[kk];
            end
         end

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-AW){acc[N-1][AW-1]}}, acc[N-1]};
         end
         else begin : g_dsn
            reg signed [AW-1:0] ds [0:D-1];
            always @(posedge clk) begin
               if (rst) begin
                  for (kk = 0; kk < D; kk = kk + 1)
                    ds[kk] <= {AW{1'b0}};
               end
               else if (en) begin
                  ds[0] <= acc[N-1];
                  for (kk = 1; kk < D; kk = kk + 1)
                    ds[kk] <= ds[kk-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-AW){ds[D-1][AW-1]}}, ds[D-1]};
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

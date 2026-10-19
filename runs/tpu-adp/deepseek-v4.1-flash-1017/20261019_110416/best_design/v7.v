// v7: exact golden staging, but each accumulate stage is a fused multiply-add.
//
// A signed 8x8 product is expanded into 8 masked partial products
//   A*W = sum_{b=0..6} A[b]*(W<<b) - A[7]*(W<<7)
// and the previous partial sum is added as a ninth addend, so the whole stage
// becomes one 9-input carry-save reduction with a single final carry-propagate
// adder instead of (multiplier CPA) + (accumulator CPA).

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
   // Per-column cascade of fused multiply-adds + output de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [AW-1:0] pp [0:N-1][0:DW-1];
         for (i = 0; i < N; i = i + 1) begin : g_pp
            wire signed [AW-1:0] wse = $signed(wreg[i][j]);
            for (genvar b = 0; b < DW; b = b + 1) begin : g_b
               if (b == DW - 1)
                 assign pp[i][b] = a_tap[i][j][b] ? -(wse <<< b) : {AW{1'b0}};
               else
                 assign pp[i][b] = a_tap[i][j][b] ?  (wse <<< b) : {AW{1'b0}};
            end
         end

         reg signed [AW-1:0] acc [0:N-1];
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1)
                 acc[kk] <= {AW{1'b0}};
            end
            else if (en) begin
               acc[0] <= pp[0][0]+pp[0][1]+pp[0][2]+pp[0][3]+
                         pp[0][4]+pp[0][5]+pp[0][6]+pp[0][7];
               for (kk = 1; kk < N; kk = kk + 1)
                 acc[kk] <= acc[kk-1] + pp[kk][0]+pp[kk][1]+pp[kk][2]+pp[kk][3]+
                            pp[kk][4]+pp[kk][5]+pp[kk][6]+pp[kk][7];
            end
         end

         localparam D = N - 1 - j;
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

// v9: exact golden staging with a REDUNDANT (carry-save) accumulation cascade.
//
// Golden per-PE staging:  P(i,j)(T) = P(i-1,j)(T-1) + A(T-i-j)[i]*B[i][j],
// with a weight reload freezing every partial sum.  The numeric value of the
// intermediate partial sums is internal, so they can be carried in redundant
// (sum,carry) form: each stage becomes a single 3:2 compressor and only one
// carry-propagate adder is needed per column at the end of the cascade.

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
   // Per-column redundant cascade + single CPA + output de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [AW-1:0] p [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_prod
            assign p[i] = $signed(a_tap[i][j]) * $signed(wreg[i][j]);
         end

         reg signed [AW-1:0] s [0:N-1];
         reg signed [AW-1:0] c [1:N-1];

         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1)
                 s[kk] <= {AW{1'b0}};
               for (kk = 1; kk < N; kk = kk + 1)
                 c[kk] <= {AW{1'b0}};
            end
            else if (en) begin
               s[0] <= p[0];
               s[1] <= s[0] ^ p[1];
               c[1] <= {s[0] & p[1], 1'b0};
               for (kk = 2; kk < N; kk = kk + 1) begin
                  s[kk] <= s[kk-1] ^ c[kk-1] ^ p[kk];
                  c[kk] <= { (s[kk-1] & c[kk-1]) |
                             (s[kk-1] & p[kk])    |
                             (c[kk-1] & p[kk]), 1'b0 };
               end
            end
         end

         wire signed [AW-1:0] sum = s[N-1] + c[N-1];

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){sum[AW-1]}}, sum};
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
                  ds[0] <= sum;
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

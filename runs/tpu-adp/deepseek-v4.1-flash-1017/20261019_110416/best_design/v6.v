// v6: TPU tile, mult/add pipelined apart.
//
// Golden per-PE behavior:  P(i,j)(t+1) = P(i-1,j)(t) + A(t-i-j)[i]*B[i][j]
// (A word arriving at time t is the row accepted at time t; a reload freezes
// everything).  Instead of a systolic mesh the activation is fanned out from a
// single 8-bit delay line and every weight is held in a plain register, but the
// 8-deep accumulation cascade is preserved exactly.
//
// The multiply and the accumulate are placed in DIFFERENT cycles:  products for
// i>=1 are registered one cycle ahead and read from delay taps shifted by one,
// so every register-to-register path is either a single 8x8 multiply or a
// single narrow add.

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
   localparam ND  = 2*N - 3;	// 13 delay-line registers

   wire		      ld = w_valid;
   wire		      en = ~w_valid;

   genvar	      i, j;
   integer	      k;

   //--------------------------------------------------------------------------
   // Activation delay line:  dd[d] holds the A word delayed by d cycles.
   //--------------------------------------------------------------------------
   wire [N*DW-1:0] din = a_valid ? a_data : {N*DW{1'b0}};
   wire [N*DW-1:0] dd [0:ND];

   assign dd[0] = din;
   generate
      for (i = 1; i <= ND; i = i + 1) begin : g_dly
         reg [N*DW-1:0] q;
         always @(posedge clk) begin
            if (rst)   q <= {N*DW{1'b0}};
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
               for (k = 0; k < N; k = k + 1)
                 wreg[k][j] <= {DW{1'b0}};
            end
            else if (ld) begin
               wreg[0][j] <= w_data[j*DW +: DW];
               for (k = 1; k < N; k = k + 1)
                 wreg[k][j] <= wreg[k-1][j];
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Per-column cascade: 8 registered products + 7 narrow adds.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [DW-1:0] a0 = dd[j];
         wire signed [DW-1:0] a1 = dd[j+0];
         wire signed [DW-1:0] a2 = dd[j+1];
         wire signed [DW-1:0] a3 = dd[j+2];
         wire signed [DW-1:0] a4 = dd[j+3];
         wire signed [DW-1:0] a5 = dd[j+4];
         wire signed [DW-1:0] a6 = dd[j+5];
         wire signed [DW-1:0] a7 = dd[j+6];

         reg signed [15:0] r0, r1, r2, r3, r4, r5, r6, r7;
         reg signed [16:0] s1, s2;
         reg signed [17:0] s3, s4, s5, s6;
         reg signed [18:0] s7;

         always @(posedge clk) begin
            if (rst) begin
               r0 <= 16'sd0; r1 <= 16'sd0; r2 <= 16'sd0; r3 <= 16'sd0;
               r4 <= 16'sd0; r5 <= 16'sd0; r6 <= 16'sd0; r7 <= 16'sd0;
               s1 <= 17'sd0; s2 <= 17'sd0;
               s3 <= 18'sd0; s4 <= 18'sd0; s5 <= 18'sd0; s6 <= 18'sd0;
               s7 <= 19'sd0;
            end
            else if (en) begin
               r0 <= a0 * $signed(wreg[0][j]);
               r1 <= a1 * $signed(wreg[1][j]);
               r2 <= a2 * $signed(wreg[2][j]);
               r3 <= a3 * $signed(wreg[3][j]);
               r4 <= a4 * $signed(wreg[4][j]);
               r5 <= a5 * $signed(wreg[5][j]);
               r6 <= a6 * $signed(wreg[6][j]);
               r7 <= a7 * $signed(wreg[7][j]);
               s1 <= r0 + r1;
               s2 <= s1 + r2;
               s3 <= s2 + r3;
               s4 <= s3 + r4;
               s5 <= s4 + r5;
               s6 <= s5 + r6;
               s7 <= s6 + r7;
            end
         end

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-19){s7[18]}}, s7};
         end
         else begin : g_dsn
            reg signed [18:0] ds [0:D-1];
            integer	      m;
            always @(posedge clk) begin
               if (rst) begin
                  for (m = 0; m < D; m = m + 1)
                    ds[m] <= 19'sd0;
               end
               else if (en) begin
                  ds[0] <= s7;
                  for (m = 1; m < D; m = m + 1)
                    ds[m] <= ds[m-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-19){ds[D-1][18]}}, ds[D-1]};
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

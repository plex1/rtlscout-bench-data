// v3: broadcast-A TPU tile with pipelined adder tree.
//
// All 64 products for one A row are formed in parallel (one 8x8 signed
// multiply per PE site), registered, then reduced by a balanced 3-level adder
// tree. Pipeline budget 15 = 11 (8-bit A delay line) + 1 (product reg)
// + 3 (tree levels). Critical path per stage is a single multiply or a single
// narrow add.

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
   localparam ND  = 11;		// activation delay depth

   wire ld = w_valid;
   wire en = ~w_valid;

   //--------------------------------------------------------------------------
   // Activation delay line (zeros injected while a_valid is low).
   //--------------------------------------------------------------------------
   reg [N*DW-1:0] adel [0:ND-1];
   integer	  k;
   always @(posedge clk) begin
      if (rst) begin
         for (k = 0; k < ND; k = k + 1)
           adel[k] <= {N*DW{1'b0}};
      end
      else if (en) begin
         adel[0] <= a_valid ? a_data : {N*DW{1'b0}};
         for (k = 1; k < ND; k = k + 1)
           adel[k] <= adel[k-1];
      end
   end

   wire [N*DW-1:0] a_cur = adel[ND-1];

   //--------------------------------------------------------------------------
   // Stationary weights: vertical shift register loaded from the top.
   //--------------------------------------------------------------------------
   reg [DW-1:0] wreg [0:N-1][0:N-1];
   integer	i, j;
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             wreg[i][j] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (j = 0; j < N; j = j + 1)
           wreg[0][j] <= w_data[j*DW +: DW];
         for (i = 1; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             wreg[i][j] <= wreg[i-1][j];
      end
   end

   //--------------------------------------------------------------------------
   // Product registers (16-bit signed).
   //--------------------------------------------------------------------------
   reg signed [2*DW-1:0] preg [0:N-1][0:N-1];
   integer		 i2, j2;
   always @(posedge clk) begin
      if (rst) begin
         for (i2 = 0; i2 < N; i2 = i2 + 1)
           for (j2 = 0; j2 < N; j2 = j2 + 1)
             preg[i2][j2] <= {(2*DW){1'b0}};
      end
      else if (en) begin
         for (i2 = 0; i2 < N; i2 = i2 + 1)
           for (j2 = 0; j2 < N; j2 = j2 + 1)
             preg[i2][j2] <= $signed(a_cur[i2*DW +: DW]) * $signed(wreg[i2][j2]);
      end
   end

   //--------------------------------------------------------------------------
   // Adder tree, level 1: 8 -> 4 sums of two 16-bit products (17 bits).
   //--------------------------------------------------------------------------
   reg signed [2*DW:0] s1 [0:N/2-1][0:N-1];
   integer	      a1, j3;
   always @(posedge clk) begin
      if (rst) begin
         for (a1 = 0; a1 < N/2; a1 = a1 + 1)
           for (j3 = 0; j3 < N; j3 = j3 + 1)
             s1[a1][j3] <= {2*DW+1{1'b0}};
      end
      else if (en) begin
         for (a1 = 0; a1 < N/2; a1 = a1 + 1)
           for (j3 = 0; j3 < N; j3 = j3 + 1)
             s1[a1][j3] <= preg[2*a1][j3] + preg[2*a1+1][j3];
      end
   end

   //--------------------------------------------------------------------------
   // Level 2: 4 -> 2 sums (18 bits).
   //--------------------------------------------------------------------------
   reg signed [2*DW+1:0] s2 [0:N/4-1][0:N-1];
   integer		 a2, j4;
   always @(posedge clk) begin
      if (rst) begin
         for (a2 = 0; a2 < N/4; a2 = a2 + 1)
           for (j4 = 0; j4 < N; j4 = j4 + 1)
             s2[a2][j4] <= {2*DW+2{1'b0}};
      end
      else if (en) begin
         for (a2 = 0; a2 < N/4; a2 = a2 + 1)
           for (j4 = 0; j4 < N; j4 = j4 + 1)
             s2[a2][j4] <= s1[2*a2][j4] + s1[2*a2+1][j4];
      end
   end

   //--------------------------------------------------------------------------
   // Level 3: 2 -> 1 (19 bits) and output sign extension.
   //--------------------------------------------------------------------------
   reg signed [AW-1:0] s3 [0:N-1];
   integer	      j5;
   always @(posedge clk) begin
      if (rst) begin
         for (j5 = 0; j5 < N; j5 = j5 + 1)
           s3[j5] <= {AW{1'b0}};
      end
      else if (en) begin
         for (j5 = 0; j5 < N; j5 = j5 + 1)
           s3[j5] <= s2[0][j5] + s2[1][j5];
      end
   end

   genvar jj;
   generate
      for (jj = 0; jj < N; jj = jj + 1) begin : g_out
         assign c_data[jj*ACCW +: ACCW] = {{(ACCW-AW){s3[jj][AW-1]}}, s3[jj]};
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

// v2: broadcast-A TPU tile.
// Instead of a systolic mesh, A is delayed and broadcast; all 64 products are
// formed in parallel and summed by a per-column pipelined adder chain.
// Latency budget 15 = 7 (A delay) + 1 (product reg) + 7 (adder chain).
// Critical path per stage is either one 8x8 multiply or one 19-bit add.

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
   localparam ND  = 7;		// activation delay stages
   localparam NC  = 7;		// adder chain stages

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
   // Product registers (one 8x8 signed multiply each, then registered).
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
   // Per-column pipelined adder chain: acc[0] = p0+p1, acc[k] = acc[k-1]+p(k+1).
   //--------------------------------------------------------------------------
   reg signed [AW-1:0] acc [0:NC-1][0:N-1];
   integer	       s, j3;
   always @(posedge clk) begin
      if (rst) begin
         for (s = 0; s < NC; s = s + 1)
           for (j3 = 0; j3 < N; j3 = j3 + 1)
             acc[s][j3] <= {AW{1'b0}};
      end
      else if (en) begin
         for (j3 = 0; j3 < N; j3 = j3 + 1) begin
            acc[0][j3] <= $signed(preg[0][j3]) + $signed(preg[1][j3]);
            for (s = 1; s < NC; s = s + 1)
              acc[s][j3] <= acc[s-1][j3] + $signed(preg[s+1][j3]);
         end
      end
   end

   genvar jj;
   generate
      for (jj = 0; jj < N; jj = jj + 1) begin : g_out
         wire signed [AW-1:0] r = acc[NC-1][jj];
         assign c_data[jj*ACCW +: ACCW] = {{(ACCW-AW){r[AW-1]}}, r};
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

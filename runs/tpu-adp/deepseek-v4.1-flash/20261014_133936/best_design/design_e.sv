// Variant E: smallest datapath.  All activations of a row are aligned to
// arrive at the multiply site simultaneously (cycle n+7+j at column j) via a
// uniform 7-cycle delay plus the systolic east chain.  The whole column sum is
// then one combinational multiply + adder tree feeding a single register per
// column.  Same cycle-exact result timing as the original systolic array.
module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
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

   localparam LAT = 2*N - 1;
   localparam IW  = 19;
   localparam PW  = 2*DW;

   wire ld = w_valid;
   wire en = ~w_valid;

   integer i, j, k;
   genvar gi, gj;

   // ---- stationary weights (shift south while loading) --------------------
   reg [DW-1:0] wq [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) wq[i][j] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (j = 0; j < N; j = j + 1) wq[0][j] <= w_data[j*DW +: DW];
         for (i = 1; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) wq[i][j] <= wq[i-1][j];
      end
   end

   // ---- activation injection (zeros while a_valid low) --------------------
   wire [DW-1:0] din [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   // ---- uniform (N-1)-cycle delay: aligned[i] = din[i] delayed N-1 --------
   reg  [DW-1:0] de [0:N-1][0:N-2];
   wire [DW-1:0] aligned [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_align
         assign aligned[gi] = de[gi][N-2];
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (k = 0; k < N-1; k = k + 1) de[i][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1) begin
            de[i][0] <= din[i];
            for (k = 1; k < N-1; k = k + 1) de[i][k] <= de[i][k-1];
         end
      end
   end

   // ---- east chain (one register per PE, N-1 per row) ---------------------
   wire [DW-1:0] a_h [0:N-1][0:N];
   reg  [DW-1:0] areg [0:N-1][0:N-2];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ah
         assign a_h[gi][0] = aligned[gi];
         for (gj = 0; gj < N-1; gj = gj + 1) begin : g_ahj
            assign a_h[gi][gj+1] = areg[gi][gj];
         end
         assign a_h[gi][N] = {DW{1'b0}};
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N-1; j = j + 1) areg[i][j] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N-1; j = j + 1) areg[i][j] <= a_h[i][j];
      end
   end

   // ---- combinational multiply + column sum tree, one register per column -
   reg [IW-1:0] sumreg [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            sumreg[j] <= $signed(a_h[0][j]) * $signed(wq[0][j])
                       + $signed(a_h[1][j]) * $signed(wq[1][j])
                       + $signed(a_h[2][j]) * $signed(wq[2][j])
                       + $signed(a_h[3][j]) * $signed(wq[3][j])
                       + $signed(a_h[4][j]) * $signed(wq[4][j])
                       + $signed(a_h[5][j]) * $signed(wq[5][j])
                       + $signed(a_h[6][j]) * $signed(wq[6][j])
                       + $signed(a_h[7][j]) * $signed(wq[7][j]);
         end
      end
   end

   // ---- output de-skew: column j delayed N-1-j ----------------------------
   reg [IW-1:0] ds [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (k = 0; k < N-1-j; k = k + 1) ds[j][k] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            if (N-1-j > 0) begin
               ds[j][0] <= sumreg[j];
               for (k = 1; k < N-1-j; k = k + 1) ds[j][k] <= ds[j][k-1];
            end
         end
      end
   end

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_out
         if (N-1-gi == 0)
           assign c_data[gi*ACCW +: ACCW] = {{(ACCW-IW){sumreg[gi][IW-1]}}, sumreg[gi]};
         else
           assign c_data[gi*ACCW +: ACCW] =
                    {{(ACCW-IW){ds[gi][N-2-gi][IW-1]}}, ds[gi][N-2-gi]};
      end
   endgenerate

   // ---- valid pipeline ----------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// Variant D: products are registered (one pipeline stage), then the 8 column
// products are summed with a COMBINATIONAL adder tree (no accumulator
// registers at all).  Same cycle-exact product/weight timing as the original
// systolic array; column j settles at cycle n+7+j, de-skewed by 7-j cycles.
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

   // ---- row skew: arow[i] = din[i] delayed i cycles -----------------------
   reg  [DW-1:0] sk [0:N-1][0:N-1];
   wire [DW-1:0] arow [0:N-1];
   assign arow[0] = din[0];
   generate
      for (gi = 1; gi < N; gi = gi + 1) begin : g_arow
         assign arow[gi] = sk[gi][gi-1];
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 1; i < N; i = i + 1)
           for (k = 0; k < i; k = k + 1) sk[i][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 1; i < N; i = i + 1) begin
            sk[i][0] <= din[i];
            for (k = 1; k < i; k = k + 1) sk[i][k] <= sk[i][k-1];
         end
      end
   end

   // ---- east activation chain (one register per PE) ----------------------
   wire [DW-1:0] a_h [0:N-1][0:N];
   reg  [DW-1:0] areg [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ah
         assign a_h[gi][0] = arow[gi];
         for (gj = 0; gj < N; gj = gj + 1) begin : g_ahj
            assign a_h[gi][gj+1] = (gj == N-1) ? {DW{1'b0}} : areg[gi][gj];
         end
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

   // ---- products, registered at the PE stage ------------------------------
   reg [PW-1:0] P [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) P[i][j] <= {PW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             P[i][j] <= $signed(a_h[i][j]) * $signed(wq[i][j]);
      end
   end

   // ---- combinational column sum tree -------------------------------------
   wire [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_col
         wire signed [IW-1:0] t0 = $signed(P[0][gi]);
         wire signed [IW-1:0] t1 = $signed(P[1][gi]);
         wire signed [IW-1:0] t2 = $signed(P[2][gi]);
         wire signed [IW-1:0] t3 = $signed(P[3][gi]);
         wire signed [IW-1:0] t4 = $signed(P[4][gi]);
         wire signed [IW-1:0] t5 = $signed(P[5][gi]);
         wire signed [IW-1:0] t6 = $signed(P[6][gi]);
         wire signed [IW-1:0] t7 = $signed(P[7][gi]);
         assign csum[gi] = t0 + t1 + t2 + t3 + t4 + t5 + t6 + t7;
      end
   endgenerate

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
               ds[j][0] <= csum[j];
               for (k = 1; k < N-1-j; k = k + 1) ds[j][k] <= ds[j][k-1];
            end
         end
      end
   end

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_out
         if (N-1-gi == 0)
           assign c_data[gi*ACCW +: ACCW] = {{(ACCW-IW){csum[gi][IW-1]}}, csum[gi]};
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

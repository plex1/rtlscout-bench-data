// Variant G: aligned-activation datapath with TWO clean pipeline stages:
//   stage 1: one multiply per PE -> product registers
//   stage 2: explicit 3-level adder tree -> column-sum register
// then a (N-1-j)-cycle output de-skew.  Cycle-exact vs. the original array.
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

   // ---- per-row weight delay of (N-1-i) cycles ----------------------------
   reg  [DW-1:0] wd [0:N-1][0:N-1][0:N-2];
   wire [DW-1:0] wsel [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_wd
         for (gj = 0; gj < N; gj = gj + 1) begin : g_wdj
            if (N-1-gi == 0)
              assign wsel[gi][gj] = wq[gi][gj];
            else
              assign wsel[gi][gj] = wd[gi][gj][N-2-gi];
         end
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N-1; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             for (k = 0; k < N-1-i; k = k + 1) wd[i][j][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N-1; i = i + 1)
           for (j = 0; j < N; j = j + 1) begin
              wd[i][j][0] <= wq[i][j];
              for (k = 1; k < N-1-i; k = k + 1)
                wd[i][j][k] <= wd[i][j][k-1];
           end
      end
   end

   // ---- activation injection -------------------------------------------------
   wire [DW-1:0] din [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   // ---- uniform (N-1)-cycle delay -----------------------------------------
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

   // ---- east chain ----------------------------------------------------------
   wire [DW-1:0] a_h [0:N-1][0:N];
   reg  [DW-1:0] areg [0:N-1][0:N-2];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ah
         assign a_h[gi][0] = aligned[gi];
         for (gj = 0; gj < N-1; gj = gj + 1) begin : g_ahj
            assign a_h[gi][gj+1] = areg[gi][gj];
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

   // ---- stage 1: product registers -----------------------------------------
   reg [PW-1:0] P [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) P[i][j] <= {PW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             P[i][j] <= $signed(a_h[i][j]) * $signed(wsel[i][j]);
      end
   end

   // ---- stage 2: explicit 3-level signed adder tree -------------------------
   wire signed [IW-1:0] s1 [0:N-1];
   reg [IW-1:0] sumreg [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= s1[j];
      end
   end

   // tree wiring
   wire signed [IW-1:0] t01 [0:N-1], t23 [0:N-1], t45 [0:N-1], t67 [0:N-1];
   wire signed [IW-1:0] u0 [0:N-1], u1 [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_tr
         assign t01[gi] = $signed(P[0][gi]) + $signed(P[1][gi]);
         assign t23[gi] = $signed(P[2][gi]) + $signed(P[3][gi]);
         assign t45[gi] = $signed(P[4][gi]) + $signed(P[5][gi]);
         assign t67[gi] = $signed(P[6][gi]) + $signed(P[7][gi]);
         assign u0[gi]  = t01[gi] + t23[gi];
         assign u1[gi]  = t45[gi] + t67[gi];
         assign s1[gi]  = u0[gi] + u1[gi];
      end
   endgenerate

   // ---- output de-skew -----------------------------------------------------
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

   // ---- valid pipeline -----------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// Variant I: design_b alignment (products registered, all 8 products of a
// column in the same cycle) + a combinational 3-level adder tree feeding ONE
// sum register per column (replaces the 7-deep accumulator chain). Then the
// (N-1-j)-cycle output de-skew.  Cycle-exact vs. the original array.
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

   // ---- stationary weights -------------------------------------------------
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

   // ---- activation injection ----------------------------------------------
   wire [DW-1:0] din [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   // ---- row skew (i cycles) + per-column delay (j cycles) ------------------
   // arow[i] = din[i] delayed i cycles
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

   // cdin[j][i] = arow[i] delayed j cycles  (activates multiply at n+i+j)
   reg  [DW-1:0] cd [0:N-1][0:N-1][0:N-1];
   wire [DW-1:0] cdin [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) assign cdin[0][gi] = arow[gi];
      for (gi = 1; gi < N; gi = gi + 1) begin : g_cdj
         for (gj = 0; gj < N; gj = gj + 1) assign cdin[gi][gj] = cd[gi][gj][gi-1];
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (j = 1; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1)
             for (k = 0; k < j; k = k + 1) cd[j][i][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (j = 1; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) begin
              cd[j][i][0] <= arow[i];
              for (k = 1; k < j; k = k + 1) cd[j][i][k] <= cd[j][i][k-1];
           end
      end
   end

   // ---- products, registered (all column products aligned) ----------------
   reg [PW-1:0] P [0:N-1][0:N-1];	// P[j][i] for column j, row i
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) P[j][i] <= {PW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1)
             P[j][i] <= $signed(cdin[j][i]) * $signed(wq[i][j]);
      end
   end

   // ---- combinational 3-level signed adder tree ---------------------------
   wire signed [IW-1:0] t0 [0:N-1], t1 [0:N-1], t2 [0:N-1], t3 [0:N-1];
   wire signed [IW-1:0] u0 [0:N-1], u1 [0:N-1];
   wire signed [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_tr
         assign t0[gi] = $signed(P[gi][0]) + $signed(P[gi][1]);
         assign t1[gi] = $signed(P[gi][2]) + $signed(P[gi][3]);
         assign t2[gi] = $signed(P[gi][4]) + $signed(P[gi][5]);
         assign t3[gi] = $signed(P[gi][6]) + $signed(P[gi][7]);
         assign u0[gi] = t0[gi] + t1[gi];
         assign u1[gi] = t2[gi] + t3[gi];
         assign csum[gi] = u0[gi] + u1[gi];
      end
   endgenerate

   // ---- single sum register per column ------------------------------------
   reg [IW-1:0] sumreg [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= csum[j];
      end
   end

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

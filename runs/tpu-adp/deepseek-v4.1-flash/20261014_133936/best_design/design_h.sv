// Variant H: aligned-activation datapath, explicit subtract-free structure.
//   activations aligned to a common time, per-row weight delay (N-1-i),
//   products as explicit wires, explicit 3-level adder tree, ONE sum register
//   per column, then (N-1-j)-cycle output de-skew.  Cycle-exact.
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

   // ---- explicit products --------------------------------------------------
   wire signed [PW-1:0] prod [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_pi
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pj
            assign prod[gi][gj] = $signed(a_h[gi][gj]) * $signed(wsel[gi][gj]);
         end
      end
   endgenerate

   // ---- explicit signed adder tree (3 levels) ------------------------------
   wire signed [IW-1:0] t0 [0:N-1], t1 [0:N-1], t2 [0:N-1], t3 [0:N-1];
   wire signed [IW-1:0] u0 [0:N-1], u1 [0:N-1];
   wire signed [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_tr
         assign t0[gi] = $signed(prod[0][gi]) + $signed(prod[1][gi]);
         assign t1[gi] = $signed(prod[2][gi]) + $signed(prod[3][gi]);
         assign t2[gi] = $signed(prod[4][gi]) + $signed(prod[5][gi]);
         assign t3[gi] = $signed(prod[6][gi]) + $signed(prod[7][gi]);
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

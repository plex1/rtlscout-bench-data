// Variant L: design_b semantics (products formed at the exact systolic cycle
// n+i+j so mid-compute weight reloads behave identically) but the per-column
// activation delay chains are replaced by a shared systolic east chain (same
// values, 4x fewer flops), and the accumulator stages use just-wide-enough
// widths.
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

   // ---- row skew: arow[i] = din[i] delayed i ------------------------------
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

   // ---- systolic east chain: a_h[i][j] = din[i] delayed i+j ---------------
   wire [DW-1:0] a_h [0:N-1][0:N];
   reg  [DW-1:0] areg [0:N-1][0:N-2];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ah
         assign a_h[gi][0] = arow[gi];
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

   // ---- products registered (identical cycle to the original PE multiply) --
   reg [PW-1:0] P [0:N-1][0:N-1];	// P[j][i]
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) P[j][i] <= {PW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1)
             P[j][i] <= $signed(a_h[i][j]) * $signed(wq[i][j]);
      end
   end

   // ---- accumulator chain, progressively sized ----------------------------
   // acc[j][i] = sum of products 0..i  -> needs 16+ceil(log2(i+1)) signed bits
   reg [PW-1:0]   acc0 [0:N-1];			// 16b, k=1
   reg [PW:0]     acc1 [0:N-1];			// 17b, k=2
   reg [PW:0]     acc2 [0:N-1];			// 17b, k=3
   reg [PW+1:0]   acc3 [0:N-1];			// 18b, k=4
   reg [PW+1:0]   acc4 [0:N-1];			// 18b, k=5
   reg [PW+1:0]   acc5 [0:N-1];			// 18b, k=6
   reg [PW+1:0]   acc6 [0:N-1];			// 18b, k=7
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1) begin
            acc0[j] <= 0; acc1[j] <= 0; acc2[j] <= 0; acc3[j] <= 0;
            acc4[j] <= 0; acc5[j] <= 0; acc6[j] <= 0;
         end
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            acc0[j] <= P[j][0];
            acc1[j] <= $signed(acc0[j]) + $signed(P[j][1]);
            acc2[j] <= $signed(acc1[j]) + $signed(P[j][2]);
            acc3[j] <= $signed(acc2[j]) + $signed(P[j][3]);
            acc4[j] <= $signed(acc3[j]) + $signed(P[j][4]);
            acc5[j] <= $signed(acc4[j]) + $signed(P[j][5]);
            acc6[j] <= $signed(acc5[j]) + $signed(P[j][6]);
         end
      end
   end

   wire signed [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_csum
         assign csum[gi] = $signed(acc6[gi]) + $signed(P[gi][7]);
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

   // ---- valid pipeline -----------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

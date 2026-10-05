//#############################################################################
// Column-parallel dot-product MXU. Weights stationary. Each column registers
// the skewed activations as they stream west->east so all 8 products for that
// column are time-aligned, sums them with a balanced adder tree, then delays
// column j by j cycles to reproduce the systolic bottom-row + de-skew timing.
//#############################################################################

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
   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   genvar i, j;

   // ------------------------------------------------------------------
   // Weight registers: per-column vertical shift chain (load from top).
   // ------------------------------------------------------------------
   wire [DW-1:0] w [0:N][0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_wedge
         assign w[0][j] = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_wrow
         for (j = 0; j < N; j = j + 1) begin : g_wcol
            reg [DW-1:0] wr;
            always @(posedge clk) begin
               if (rst) wr <= {DW{1'b0}};
               else if (ld) wr <= w[i][j];
            end
            assign w[i+1][j] = wr;
         end
      end
   endgenerate

   // ------------------------------------------------------------------
   // Input skew: row i delayed i cycles (zeros when a_valid low).
   // ------------------------------------------------------------------
   wire [N*DW-1:0] a_left;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   // ------------------------------------------------------------------
   // Per-column registered activations (west->east), so column j sees all 8
   // products time-aligned. ac[i][j] is row i's activation as it reaches
   // column j (the value the original PE(i,j) would multiply).
   // ------------------------------------------------------------------
   wire [DW-1:0] ac [0:N-1][0:N-1];
   generate
      for (i = 0; i < N; i = i + 1) begin : g_arow
         // column 0 registers the west-edge activation
         reg [DW-1:0] ar0;
         always @(posedge clk) begin
            if (rst) ar0 <= {DW{1'b0}};
            else if (en) ar0 <= a_left[i*DW +: DW];
         end
         assign ac[i][0] = ar0;
         for (j = 1; j < N; j = j + 1) begin : g_acol
            reg [DW-1:0] arj;
            always @(posedge clk) begin
               if (rst) arj <= {DW{1'b0}};
               else if (en) arj <= ac[i][j-1];
            end
            assign ac[i][j] = arj;
         end
      end
   endgenerate

   // ------------------------------------------------------------------
   // Column dot products + per-column alignment.
   // Column j: products use ac[i][j] (all aligned at column j), tree-sum,
   // register once (colr), then delay j cycles to hit the bottom-row time.
   // ------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire [2*DW-1:0] pp [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_pp
            assign pp[i] = $signed(ac[i][j]) * $signed(w[i+1][j]);
         end
         wire [ACCW-1:0] s [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_sext
            assign s[i] = {{(ACCW-2*DW){pp[i][2*DW-1]}}, pp[i]};
         end
         wire [ACCW-1:0] l1 [0:3];
         for (i = 0; i < 4; i = i + 1) begin : g_l1
            assign l1[i] = s[2*i] + s[2*i+1];
         end
         wire [ACCW-1:0] l2 [0:1];
         assign l2[0] = l1[0] + l1[1];
         assign l2[1] = l1[2] + l1[3];
         wire [ACCW-1:0] dot = l2[0] + l2[1];

         reg [ACCW-1:0] colr;
         always @(posedge clk) begin
            if (rst) colr <= {ACCW{1'b0}};
            else if (en) colr <= dot;
         end

         // Column j: delay the registered tree result by j cycles to align.
         if (j == 0) begin : g_c0
            assign c_data[j*ACCW +: ACCW] = colr;
         end
         else begin : g_cn
            reg [ACCW-1:0] rr [0:j-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < j; k = k + 1) rr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  rr[0] <= colr;
                  for (k = 1; k < j; k = k + 1) rr[k] <= rr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = rr[j-1];
         end
      end
   endgenerate

   // ------------------------------------------------------------------
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   // ------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

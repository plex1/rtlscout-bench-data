// Variant B: restructured datapath. Products and accumulation are placed in
// SEPARATE pipeline stages (short critical path), while the product for each
// (row i, col j) is still computed at exactly the same compute-cycle as the
// original systolic array (so the weight value sampled is identical, keeping
// reload-during-drain semantics exact).
//
// Original timing (compute cycle n = row presented):
//   product P[i][j] for row n is formed at cycle n+i+j and registered.
//   column sums flow south; result row n appears at cycle n+15.
//
// Here: A rows are skewed by i (row skew) and then delayed by j per column, so
// A[n][i] reaches column j at cycle n+i+j -> product register identical to the
// original. The 8 products of a column arrive one per cycle; a single-adder
// per stage chain (acc_i = acc_{i-1} + P_i, registered) sums them, and the
// (7-j) cycle output de-skew aligns the columns at cycle n+15 exactly.
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
   localparam IW  = 19;	// signed accumulator width (all true sums fit)
   localparam PW  = 2*DW;	// stored product width

   wire ld = w_valid;
   wire en = ~w_valid;

   integer i, j, k;

   //##########################################################################
   // Stationary weights: wq[i][j], shifted south while loading (exact copy of
   // the original PE weight shift register).
   //##########################################################################
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

   //##########################################################################
   // Activation injection: zeros while a_valid is low.
   //##########################################################################
   wire [DW-1:0] din [0:N-1];
   genvar gi;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   // row skew: arow[i] = din[i] delayed i compute cycles
   reg  [DW-1:0] sk [0:N-1][0:N-1];
   wire [DW-1:0] arow [0:N-1];
   generate
      assign arow[0] = din[0];
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

   // per-column activation delay: cdin[j][i] = arow[i] delayed j compute cycles
   reg  [DW-1:0] cd [0:N-1][0:N-1][0:N-1];
   wire [DW-1:0] cdin [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) assign cdin[0][gi] = arow[gi];
      for (gi = 1; gi < N; gi = gi + 1) begin : g_cdj
         genvar gj;
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

   //##########################################################################
   // Products (registered) - same cycle as the original PE multiply.
   //##########################################################################
   wire [PW-1:0] prod [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_prod_i
         genvar gj;
         for (gj = 0; gj < N; gj = gj + 1) begin : g_prod_j
            assign prod[gj][gi] = $signed(cdin[gj][gi]) * $signed(wq[gi][gj]);
         end
      end
   endgenerate

   reg [PW-1:0] P [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) P[j][i] <= {PW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) P[j][i] <= prod[j][i];
      end
   end

   //##########################################################################
   // Per-column accumulation chain: one adder per registered stage.
   // acc[j][i] = sum of P[0..i]; acc[j][N-2] registered, last add combinational.
   //##########################################################################
   reg [IW-1:0] acc [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N-1; i = i + 1) acc[j][i] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            acc[j][0] <= {{(IW-PW){P[j][0][PW-1]}}, P[j][0]};
            for (i = 1; i < N-1; i = i + 1)
              acc[j][i] <= acc[j][i-1] + {{(IW-PW){P[j][i][PW-1]}}, P[j][i]};
         end
      end
   end

   wire [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_csum
         assign csum[gi] = acc[gi][N-2] + {{(IW-PW){P[gi][N-1][PW-1]}}, P[gi][N-1]};
      end
   endgenerate

   //##########################################################################
   // Output de-skew: column j delayed N-1-j cycles -> all align at cycle n+15
   //##########################################################################
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

   //##########################################################################
   // Valid pipeline (advances only while the array computes).
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

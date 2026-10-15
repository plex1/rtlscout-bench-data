// TPU MXU tile, cycle-accurate to the logikbench 'tpu' benchmark, restructured.
// Weights held in a south-shifting register array (identical load semantics);
// activation rows delayed in a line; per-column dot product in two registered
// stages: registered products, then a balanced sum tree. Small cones per stage.
//
// Timing: an A row accepted on compute cycle s produces its result row on
// cycle s+14 (c_valid from a 15-deep valid pipe, inclusive count).
//   ap[0] .. ap[12]  : 13-stage activation delay line
//   pr               : registered products (cycle s+13)
//   r  -> c_data     : registered column sums (cycle s+14)
module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter AW = 19)	 // accumulator width (signed; |sum| <= 131072)
   (
    input		    clk,
    input		    rst,	// synchronous, active high
    input		    w_valid,
    input [N*DW-1:0]	    w_data,
    input		    a_valid,
    input [N*DW-1:0]	    a_data,
    output		    c_valid,
    output [N*32-1:0]	    c_data
    );

   localparam LAT = 2*N - 1;	// 15: valid-pipe depth (inclusive latency)
   localparam D   = LAT - 2;	// 13: activation delay-line depth
   localparam PW  = 18;		// partial-sum width (signed)

   // ------------------------------------------------------------------
   // Weight array: column-wise south shift on w_valid, held otherwise.
   // ------------------------------------------------------------------
   reg [DW-1:0] wv [0:N*N-1];
   integer      i, j, k, m;
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N*N; i = i + 1) wv[i] <= {DW{1'b0}};
      end
      else if (w_valid) begin
         for (j = 0; j < N; j = j + 1) wv[j] <= w_data[j*DW +: DW];
         for (i = 1; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) wv[i*N+j] <= wv[(i-1)*N+j];
      end
   end

   // ------------------------------------------------------------------
   // Activation delay line (compute cycles only; zeros when a_valid low;
   // rows arriving during w_valid are dropped, as in the original).
   // ------------------------------------------------------------------
   reg [N*DW-1:0] ap [0:D-1];
   always @(posedge clk) begin
      if (rst) begin
         for (k = 0; k < D; k = k + 1) ap[k] <= {N*DW{1'b0}};
      end
      else if (~w_valid) begin
         ap[0] <= a_valid ? a_data : {N*DW{1'b0}};
         for (k = 1; k < D; k = k + 1) ap[k] <= ap[k-1];
      end
   end

   // ------------------------------------------------------------------
   // Stage 1: registered products a[i] * w[i][j].
   // ------------------------------------------------------------------
   wire signed [2*DW-1:0] prod [0:N*N-1];
   reg  signed [2*DW-1:0] pr   [0:N*N-1];
   genvar gi, gj, gk;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_mul
         for (gj = 0; gj < N; gj = gj + 1) begin : g_mulj
            assign prod[gi*N+gj] = $signed(ap[D-1][gi*DW +: DW]) *
                                   $signed(wv[gi*N+gj]);
         end
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) for (m = 0; m < N*N; m = m + 1) pr[m] <= {(2*DW){1'b0}};
      else if (~w_valid) for (m = 0; m < N*N; m = m + 1) pr[m] <= prod[m];
   end

   // ------------------------------------------------------------------
   // Stage 2: balanced sum-of-products tree per column.
   // ------------------------------------------------------------------
   wire signed [PW-1:0] t0 [0:4*N-1];   // sums of 2 products
   wire signed [PW-1:0] t1 [0:2*N-1];   // sums of 4 products
   wire signed [AW-1:0] t2 [0:N-1];     // sums of 8 products
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_tree
         for (gk = 0; gk < 4; gk = gk + 1) begin : g_t0
            assign t0[gj*4+gk] = pr[(2*gk)*N+gj] + pr[(2*gk+1)*N+gj];
         end
         for (gk = 0; gk < 2; gk = gk + 1) begin : g_t1
            assign t1[gj*2+gk] = t0[gj*4+2*gk] + t0[gj*4+2*gk+1];
         end
         assign t2[gj] = t1[gj*2] + t1[gj*2+1];
      end
   endgenerate

   reg signed [AW-1:0] r [0:N-1];
   always @(posedge clk) begin
      if (rst) for (m = 0; m < N; m = m + 1) r[m] <= {AW{1'b0}};
      else if (~w_valid) for (m = 0; m < N; m = m + 1) r[m] <= t2[m];
   end

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_out
         assign c_data[gj*32 +: 32] = {{(32-AW){r[gj][AW-1]}}, r[gj]};
      end
   endgenerate

   // ------------------------------------------------------------------
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   // ------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (~w_valid) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

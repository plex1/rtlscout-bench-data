//#############################################################################
// Weight-stationary systolic matmul tile, 8x8, 8-bit signed operands.
//
// Architecture: the 8 matrix rows are FOLDED into a single row of 8 PEs
// (PE(j) owns column j of B and all 8 of its weights).  The activation row
// streams east through a 64-bit shift chain, one PE per column.  Each PE
// forms the 8 products of the current activation row with its stationary
// weights in a registered stage, then a balanced adder tree reduces them in
// the next stage.  The 15-cycle contract latency is filled by an input
// pipeline delay (5 cycles) plus the column de-skew (7-j cycles).
//
//   input row m at cycle m
//   ac[j] holds row m at cycle m+6+j          (5 skew + 1 + j chain hops)
//   product regs valid at cycle m+7+j
//   tree result (19-bit) valid at cycle m+8+j
//   de-skew 7-j  ->  c_data row m valid at cycle m+15  == c_valid
//#############################################################################
module tpu
  (input		 clk,
   input		 rst,
   input		 w_valid,
   input [63:0]	 w_data,
   input		 a_valid,
   input [63:0]	 a_data,
   output		 c_valid,
   output [255:0] c_data);

   localparam N    = 8;
   localparam DW   = 8;
   localparam AW   = 19;	  // exact accumulator width: |sum| <= 8*16384 = 2^17
   localparam ACCW = 32;
   localparam LAT  = 2*N - 1;
   localparam SK   = 5;		  // input pipeline delay in compute cycles

   wire	      en = ~w_valid;
   wire	      ld = w_valid;

   integer    i, j;

   //------------------------------------------------------------------
   // stationary weights: per column an 8-stage byte shift register
   // (cycle k carries B[7-k][*], so depth r ends up holding B[r][j])
   //------------------------------------------------------------------
   reg [DW-1:0] wq [0:N-1][0:N-1];	// wq[col][row]
   always @(posedge clk) begin
      if (rst) begin
	 for (j = 0; j < N; j = j + 1)
	   for (i = 0; i < N; i = i + 1)
	     wq[j][i] <= {DW{1'b0}};
      end
      else if (ld) begin
	 for (j = 0; j < N; j = j + 1) begin
	    wq[j][0] <= w_data[j*DW +: DW];
	    for (i = 1; i < N; i = i + 1)
	      wq[j][i] <= wq[j][i-1];
	 end
      end
   end

   //------------------------------------------------------------------
   // input pipeline delay (fills the latency budget)
   //------------------------------------------------------------------
   reg [N*DW-1:0] sk [0:SK-1];
   always @(posedge clk) begin
      if (rst) begin
	 for (i = 0; i < SK; i = i + 1) sk[i] <= {N*DW{1'b0}};
      end
      else if (en) begin
	 sk[0] <= a_valid ? a_data : {N*DW{1'b0}};
	 for (i = 1; i < SK; i = i + 1) sk[i] <= sk[i-1];
      end
   end

   //------------------------------------------------------------------
   // east propagation chain: ac[j] carries the activation row that PE(j)
   // consumes this cycle.
   //------------------------------------------------------------------
   reg [N*DW-1:0] ac [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
	 for (i = 0; i < N; i = i + 1) ac[i] <= {N*DW{1'b0}};
      end
      else if (en) begin
	 ac[0] <= sk[SK-1];
	 for (i = 1; i < N; i = i + 1) ac[i] <= ac[i-1];
      end
   end

   //------------------------------------------------------------------
   // products (registered) and per-column reduction tree
   //------------------------------------------------------------------
   reg signed [2*DW-1:0] pr [0:N-1][0:N-1];
   reg signed [AW-1:0]	 acc [0:N-1];

   genvar gj;
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
	 integer k;
	 always @(posedge clk) begin
	    if (rst) begin
	       for (k = 0; k < N; k = k + 1) pr[gj][k] <= {2*DW{1'b0}};
	       acc[gj] <= {AW{1'b0}};
	    end
	    else if (en) begin
	       for (k = 0; k < N; k = k + 1)
		 pr[gj][k] <= $signed(ac[gj][k*DW +: DW]) * $signed(wq[gj][k]);
	       acc[gj] <= (  {{3{pr[gj][0][2*DW-1]}}, pr[gj][0]}
			   + {{3{pr[gj][1][2*DW-1]}}, pr[gj][1]}
			   + {{3{pr[gj][2][2*DW-1]}}, pr[gj][2]}
			   + {{3{pr[gj][3][2*DW-1]}}, pr[gj][3]})
			+ (  {{3{pr[gj][4][2*DW-1]}}, pr[gj][4]}
			   + {{3{pr[gj][5][2*DW-1]}}, pr[gj][5]}
			   + {{3{pr[gj][6][2*DW-1]}}, pr[gj][6]}
			   + {{3{pr[gj][7][2*DW-1]}}, pr[gj][7]});
	    end
	 end
      end
   endgenerate

   //------------------------------------------------------------------
   // output de-skew (column j delayed N-1-j) + sign extension to 32 bits
   //------------------------------------------------------------------
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_desk
	 localparam D = N - 1 - gj;
	 if (D == 0) begin : g_d0
	    assign c_data[gj*ACCW +: ACCW] = {{(ACCW-AW){acc[gj][AW-1]}}, acc[gj]};
	 end
	 else begin : g_dn
	    reg signed [AW-1:0] dr [0:D-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1) dr[k] <= {AW{1'b0}};
	       end
	       else if (en) begin
		  dr[0] <= acc[gj];
		  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[gj*ACCW +: ACCW] = {{(ACCW-AW){dr[D-1][AW-1]}}, dr[D-1]};
	 end
      end
   endgenerate

   //------------------------------------------------------------------
   // valid pipeline
   //------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)		vpipe <= {LAT{1'b0}};
      else if (en)	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// v8: exact-value two-stage PE.
//
// The reference PE computes  psum_south(t+1) = psum_north(t) + a(t)*w(t),
// so the FF->FF path is "multiply + 19-bit add".  We keep the *exact same
// values* but split the datapath into two registers:
//
//   pr(t+1)         = a(t+1) * w(t+1)          <- product register
//   psum_south(t+1) = psum_north(t) + pr(t)
//
// a(t+1) and w(t+1) are known one cycle early because
//   a(t+1) = en ? a_prev_stage(t) : a(t)      (activation chain, frozen on load)
//   w(t+1) = ld ? w_north(t)     : w(t)       (weight shift chain)
// Both are pure FF outputs / registered shifts, so the prediction is exact.
//
// Row 0 needs no prediction: psum_v[1][j](t+1) = a_h[0][j](t) * w[0][j](t).
//
// So the FF->FF critical path becomes max(multiply, 19-bit add) instead of
// multiply+add, with no change to the 15-cycle contract or to the transient
// (garbage) values that the golden trace checks.
//#############################################################################
module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19,
    parameter ACCW = 32)
   (
    input		  clk,
    input		  rst,
    input		  w_valid,
    input [N*DW-1:0]	  w_data,
    input		  a_valid,
    input [N*DW-1:0]	  a_data,
    output		  c_valid,
    output [N*ACCW-1:0]  c_data
    );

   localparam LAT = 2*N - 1;
   localparam PW  = 2*DW;	  // product width

   wire ld = w_valid;
   wire en = ~w_valid;

   integer i, j, d;

   //------------------------------------------------------------------
   // stationary weights (south-shifting chain, one byte per row)
   //------------------------------------------------------------------
   reg [DW-1:0] wq [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
	 for (i = 0; i < N; i = i + 1)
	   for (j = 0; j < N; j = j + 1) wq[i][j] <= {DW{1'b0}};
      end
      else if (ld) begin
	 for (j = 0; j < N; j = j + 1) begin
	    wq[0][j] <= w_data[j*DW +: DW];
	    for (i = 1; i < N; i = i + 1) wq[i][j] <= wq[i-1][j];
	 end
      end
   end

   //------------------------------------------------------------------
   // activation delay chains: row i needs a_data[i] delayed by i+j.
   // dly[i][d] = din[i] delayed d en-steps, d = 1 .. i+7
   //------------------------------------------------------------------
   wire [DW-1:0] din [0:N-1];
   reg  [DW-1:0] dly [0:N-1][0:N+6];

   genvar gi, gj;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
	 assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   always @(posedge clk) begin
      if (rst) begin
	 for (i = 0; i < N; i = i + 1)
	   for (d = 0; d < N+7; d = d + 1) dly[i][d] <= {DW{1'b0}};
      end
      else if (en) begin
	 for (i = 0; i < N; i = i + 1) begin
	    dly[i][1] <= din[i];
	    for (d = 2; d <= i+7; d = d + 1) dly[i][d] <= dly[i][d-1];
	 end
      end
   end

   //------------------------------------------------------------------
   // per-PE combinational operands
   //   ah  = activation presented to PE(i,j) this cycle
   //   ahp = activation PE(i,j) will see next cycle
   //   wn  = weight PE(i,j) will hold next cycle
   //------------------------------------------------------------------
   wire [DW-1:0] ah  [0:N-1][0:N-1];
   wire [DW-1:0] ahp [0:N-1][0:N-1];
   wire [DW-1:0] wn  [0:N-1][0:N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_r
	 for (gj = 0; gj < N; gj = gj + 1) begin : g_c
	    localparam integer D = gi + gj;
	    if (D == 0) begin : g_d0
	       assign ah[gi][gj]  = din[0];
	       assign ahp[gi][gj] = din[0];
	    end
	    else begin : g_dn
	       assign ah[gi][gj]  = dly[gi][D];
	       assign ahp[gi][gj] = en ? ((D == 1) ? din[gi] : dly[gi][D-1])
				       : ah[gi][gj];
	    end
	    if (gi == 0) begin : g_w0
	       assign wn[gi][gj] = ld ? w_data[gj*DW +: DW] : wq[gi][gj];
	    end
	    else begin : g_wn
	       assign wn[gi][gj] = ld ? wq[gi-1][gj] : wq[gi][gj];
	    end
	 end
      end
   endgenerate

   //------------------------------------------------------------------
   // product registers (rows 1..N-1) and column accumulators
   //------------------------------------------------------------------
   reg [PW-1:0]  pr [1:N-1][0:N-1];
   reg [AW-1:0]  psum [1:N][0:N-1];

   always @(posedge clk) begin
      if (rst) begin
	 for (i = 1; i < N; i = i + 1)
	   for (j = 0; j < N; j = j + 1) pr[i][j] <= {PW{1'b0}};
	 for (i = 1; i <= N; i = i + 1)
	   for (j = 0; j < N; j = j + 1) psum[i][j] <= {AW{1'b0}};
      end
      else begin
	 // row 0: product lands straight in the first accumulator
	 for (j = 0; j < N; j = j + 1)
	   psum[1][j] <= en ? $signed(ah[0][j]) * $signed(wq[0][j])
			    : psum[1][j];
	 // rows 1..N-1: product register, then add
	 for (i = 1; i < N; i = i + 1)
	   for (j = 0; j < N; j = j + 1)
	     pr[i][j] <= $signed(ahp[i][j]) * $signed(wn[i][j]);
	 for (i = 1; i < N; i = i + 1)
	   for (j = 0; j < N; j = j + 1)
	     psum[i+1][j] <= en
			     ? (psum[i][j] + {{(AW-PW){pr[i][j][PW-1]}}, pr[i][j]})
			     : psum[i+1][j];
      end
   end

   //------------------------------------------------------------------
   // output de-skew + sign extension
   //------------------------------------------------------------------
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_deskew
	 localparam integer DD = N - 1 - gj;
	 if (DD == 0) begin : g_d0
	    assign c_data[gj*ACCW +: ACCW] =
		   {{(ACCW-AW){psum[N][gj][AW-1]}}, psum[N][gj]};
	 end
	 else begin : g_dn
	    reg [AW-1:0] dr [0:DD-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < DD; k = k + 1) dr[k] <= {AW{1'b0}};
	       end
	       else if (en) begin
		  dr[0] <= psum[N][gj];
		  for (k = 1; k < DD; k = k + 1) dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[gj*ACCW +: ACCW] =
		   {{(ACCW-AW){dr[DD-1][AW-1]}}, dr[DD-1]};
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

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
//
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Architecture notes (vs. the original per-PE systolic mesh):
//  * Partial sums are held at exactly PW = 2*DW + clog2(N) = 19 bits: the sum
//    of N signed DWxDW products lies in [-N*2^(2DW-2), N*2^(2DW-2)] which
//    always fits (N*2^14 = 2^17 <= 2^18-1 ... the extreme +131072 needs the
//    19th bit).  Results are sign-extended to ACCW at the output only.
//  * Column j of the original array finishes (8 + j) cycles after its A row
//    enters, so the original needs a per-column de-skew of (N-1-j) cycles.
//    Here instead every column's activation is delayed by exactly N-1 = 7
//    cycles, so ALL columns finish simultaneously at cycle 2N-1 = 15 and the
//    de-skew network disappears completely.
//  * Row i's activation byte is therefore delayed by (i + N-1) compute
//    cycles: a single shift line per row (products stay staggered by one
//    cycle per row, which is what the south-flowing partial-sum chain needs).
//  * The valid pipe delays a_valid by LAT = 2N-1 compute cycles; everything
//    (delay lines, accumulators, valid pipe) stalls while w_valid is high.
//#############################################################################

module tpu
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		      clk,
    input		      rst,     // synchronous, active high
    input		      w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]	      w_data,  // weight row, packed { col N-1 ... col 0 }
    input		      a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]	      a_data,  // activation row, packed { col N-1 ... col 0 }
    output		      c_valid, // result-row valid
    output [N*ACCW-1:0]	      c_data   // result row, packed { col N-1 ... col 0 }
    );

   // exact partial-sum width: N signed 2*DW-bit products
   localparam PW = 2*DW + $clog2(N);
   localparam LAT = 2*N - 1;

   // the array advances (and the weight registers shift) on disjoint phases
   wire en = ~w_valid;

   //##########################################################################
   // Activation distribution: row i carries byte i of a_data delayed by
   // (i + N-1) compute cycles; zeros are injected while a_valid is low so
   // idle/drain cycles contribute nothing to the accumulators.
   //##########################################################################
   wire [DW-1:0] arow [0:N-1];

   genvar	   gi, gj;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_aline
	 localparam DEPTH = gi + N - 1;
	 reg [DW-1:0] sh [0:DEPTH-1];
	 integer	    k;
	 wire [DW-1:0]	    din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
	 always @(posedge clk) begin
	    if (rst) begin
	       for (k = 0; k < DEPTH; k = k + 1)
		 sh[k] <= {DW{1'b0}};
	    end
	    else if (en) begin
	       sh[0] <= din;
	       for (k = 1; k < DEPTH; k = k + 1)
		 sh[k] <= sh[k-1];
	    end
	 end
	 assign arow[gi] = sh[DEPTH-1];
      end
   endgenerate

   //##########################################################################
   // The PE mesh: stationary weight per cell, product accumulated into the
   // south-flowing partial-sum chain.  Row 0's chain input is hardwired zero.
   //##########################################################################
   reg [DW-1:0]	  wreg [0:N-1][0:N-1];
   reg [PW-1:0]	  psum [0:N-1][0:N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_row
	 for (gj = 0; gj < N; gj = gj + 1) begin : g_col
	    wire [DW-1:0] w_north =
			  (gi == 0) ? w_data[gj*DW +: DW] : wreg[gi-1][gj];
	    wire [PW-1:0] p_north =
			  (gi == 0) ? {PW{1'b0}} : psum[gi-1][gj];
	    wire [2*DW-1:0] prod = $signed(arow[gi]) * $signed(wreg[gi][gj]);

	    always @(posedge clk) begin
	       if (rst)
		 wreg[gi][gj] <= {DW{1'b0}};
	       else if (w_valid)
		 wreg[gi][gj] <= w_north;
	    end

	    always @(posedge clk) begin
	       if (rst)
		 psum[gi][gj] <= {PW{1'b0}};
	       else if (en)
		 psum[gi][gj] <= p_north +
				 {{(PW-2*DW){prod[2*DW-1]}}, prod};
	    end
	 end
      end
   endgenerate

   //##########################################################################
   // Outputs: bottom-row partial sums are already complete results.
   //##########################################################################
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_out
	 assign c_data[gj*ACCW +: ACCW] =
	   {{(ACCW-PW){psum[N-1][gj][PW-1]}}, psum[N-1][gj]};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

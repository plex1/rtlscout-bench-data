//#############################################################################
// Cycle-exact rewrite of the logikbench TPU (weight-stationary MXU tile).
//
// Structure (equivalent to the 64-PE mesh, but much smaller):
//   * weights: 8 rows x 8 bytes, shifting south on every load cycle.
//   * activation skew: element i of an accepted row is delayed i+j cycles
//     before reaching column j's multiplier (one shared chain of length i+7
//     per element, tapped at each column) -> product (i,j) of the row
//     accepted at compute step t is computed at step t+i+j, with the weight
//     value in force at that step, exactly like the original array.
//   * per-column psum chain: one narrow add per row (a partial sum of at
//     most i+1 8x8 products needs only pw(i) = 16..19 bits, not 32).
//   * per-column de-skew of 7-j cycles and a 15-deep valid pipeline.
//
// Row accepted at compute step t => C row on c_data during step t+15 with
// c_valid = a_valid delayed 15 compute cycles. w_valid freezes the whole
// compute pipeline (outputs hold) and reloads the weights.
//#############################################################################

module tpu (
   input		 clk,
   input		 rst,
   input		 w_valid,
   input [63:0]		 w_data,
   input		 a_valid,
   input [63:0]		 a_data,
   output		 c_valid,
   output [255:0]	 c_data
   );

   wire en = ~w_valid;

   //-------------------------------------------------------------------
   // Weight storage: row i at W[i*64 +: 64], byte j at + j*8. Rows shift
   // away from row 0 (which captures w_data) on every load cycle.
   //-------------------------------------------------------------------
   reg [511:0] W;
   always @(posedge clk) begin
      if (rst)	       W <= 512'd0;
      else if (w_valid) W <= {W[447:0], w_data};
   end

   //-------------------------------------------------------------------
   // Activation skew chains: tap for (i,j) = element i delayed (i+j).
   //-------------------------------------------------------------------
   wire [511:0] tap;

   genvar i, j;
   generate
      for (i = 0; i < 8; i = i + 1) begin : g_arow
	 localparam L = i + 7;
	 reg [8*L-1:0] sr;
	 wire [7:0]    din = a_valid ? a_data[i*8 +: 8] : 8'd0;
	 always @(posedge clk) begin
	    if (rst)		 sr <= {8*L{1'b0}};
	    else if (en)	 sr <= {sr[8*L-9:0], din};
	 end
	 for (j = 0; j < 8; j = j + 1) begin : g_atap
	    if (i + j == 0) begin : g_z
	       assign tap[7:0] = din;
	    end
	    else begin : g_n
	       assign tap[(i*8+j)*8 +: 8] = sr[8*(i+j-1) +: 8];
	    end
	 end
      end
   endgenerate

   //-------------------------------------------------------------------
   // Per-column psum chains.
   //-------------------------------------------------------------------
   function integer pw(input integer r);
      pw = 16 + (r > 0) + (r > 2) + (r > 6);
   endfunction

   generate
      for (j = 0; j < 8; j = j + 1) begin : g_col
	 wire [151:0] pv;	// pv[i*19 +: 19] = row i partial sum
	 for (i = 0; i < 8; i = i + 1) begin : g_pe
	    localparam Wd = pw(i);
	    wire [15:0] pr = $signed(tap[(i*8+j)*8 +: 8]) *
			     $signed(W[i*64 + j*8 +: 8]);
	    reg  signed [Wd-1:0] p;
	    wire signed [Wd-1:0] pn;
	    if (i == 0) begin : g_p0
	       assign pn = {Wd{1'b0}};
	    end
	    else begin : g_pn
	       assign pn = pv[(i-1)*19 +: Wd];
	    end
	    always @(posedge clk) begin
	       if (rst)	    p <= {Wd{1'b0}};
	       else if (en) p <= pn + $signed(pr);
	    end
	    assign pv[i*19 +: 19] = p;
	 end

	 // de-skew column j by 7-j cycles
	 localparam D = 7 - j;
	 if (D == 0) begin : g_d0
	    assign c_data[j*32 +: 32] = $signed(pv[7*19 +: 19]);
	 end
	 else begin : g_dn
	    reg signed [18:0] d [0:D-1];
	    integer t;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (t = 0; t < D; t = t + 1)
		    d[t] <= 19'sd0;
	       end
	       else if (en) begin
		  d[0] <= pv[7*19 +: 19];
		  for (t = 1; t < D; t = t + 1)
		    d[t] <= d[t-1];
	       end
	    end
	    assign c_data[j*32 +: 32] = $signed(d[D-1]);
	 end
      end
   endgenerate

   //-------------------------------------------------------------------
   // c_valid = a_valid delayed 15 compute cycles.
   //-------------------------------------------------------------------
   reg [14:0] vp;
   always @(posedge clk) begin
      if (rst)	  vp <= 15'd0;
      else if (en) vp <= {vp[13:0], a_valid};
   end
   assign c_valid = vp[14];

endmodule

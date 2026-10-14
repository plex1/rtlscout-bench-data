//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- narrow datapath
// with sign-magnitude operand encoding.
//
// Cycle-accurate vs the logikbench `tpu` golden trace (anti-diagonal term
// timing u = i+j, exact under mid-drain weight reloads).
//
// Key ideas:
//  * partial sums / de-skew store only the EXACT signed width (16..19 bits).
//  * synchronous reset only on the valid pipeline.
//  * Both operands are carried through their pipelines in sign-magnitude form
//    {sign, mag[7:0]} (9 bits), encoded ONCE at the injection points (a_data
//    into the input skew, w_data into the weight chain) -- i.e. outside the
//    critical cone.  The multiply becomes a fully UNSIGNED 8x8 of the
//    magnitudes plus a single conditional negate (neg = sign_a ^ sign_w),
//    whose invert/carry-in folds into the final addition:
//        a*w = (-1)^(sa^sw) * magA * magW      (exact integer identity)
//    Values are unchanged, so all accumulators keep exact-width semantics.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter EW = 9,	 // encoded operand width {sign, mag[DW-1:0]}
    parameter PIN = 16,	 // psum_north width (signed)
    parameter POUT = 16) // psum_south width (signed)
   (
    input		      clk,
    input		      rst,	      // (used only by the valid pipe)
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [EW-1:0]	      w_north,    // encoded weight in from PE above
    input [EW-1:0]	      a_west,     // encoded activation in from west
    input [PIN-1:0]	      psum_north, // partial sum in from above
    output [EW-1:0]	      w_south,    // encoded weight out to PE below
    output reg [EW-1:0]	      a_east,     // encoded activation out to east
    output reg [POUT-1:0]     psum_south  // partial sum out to below
    );

   reg [EW-1:0] w;
   assign w_south = w;

   // unsigned magnitude product (top bits are constant 0 -> reduced by opt)
   wire [2*DW-1:0] magprod = {1'b0, a_west[DW-1:0]} * {1'b0, w[DW-1:0]};
   // conditional negate: a*w = (-1)^(sa^sw) * magA * magW  (exact mod 2^16)
   wire	     neg = a_west[DW] ^ w[DW];
   wire [2*DW-1:0] prod = neg ? (~magprod + {{(2*DW-1){1'b0}}, 1'b1}) : magprod;

   always @(posedge clk) begin
      if (ld)
	w <= w_north;
      if (en) begin
	 a_east     <= a_west;
	 // context width = POUT (>= PIN and >= 2*DW); both operands signed,
	 // so they sign-extend into the POUT-bit exact-width accumulator.
	 psum_south <= $signed(psum_north) + $signed(prod);
      end
   end

endmodule

//#############################################################################

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		clk,
    input		rst,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*(DW+1)-1:0] a_left,	// encoded activation inputs, one per row
    output [N*ACCW-1:0]	c_bot	// column results, one per column
    );

   localparam EW = DW + 1; // {sign, magnitude} encoding

   // exact signed width of a sum of k products of DW-bit signed numbers
   function integer pwidth;
      input integer k;
      integer	       w, pmax, pmin;
      begin
	 pmax = (1 << (2*DW - 2));			// max product
	 pmin = -((1 << (DW-1)) * ((1 << (DW-1)) - 1)); // min product
	 w = 2*DW;
	 while ((k*pmax > ((1 << (w-1)) - 1)) || (k*pmin < -((1 << (w-1)))))
	   w = w + 1;
	 pwidth = w;
      end
   endfunction

   // sign-magnitude encoding: {w[7], |w|} with |w| as UNSIGNED DW bits
   // (|−128| = 128 fits in 8 unsigned bits)
   function [EW-1:0] sm_enc;
      input [DW-1:0] w;
      begin
	 sm_enc = w[DW-1] ? {1'b1, (~w + {{(DW-1){1'b0}}, 1'b1})}
			  : {1'b0, w};
      end
   endfunction

   // partial-sum nets: (N+1) rows x N columns, ACCW-wide slots (wires only;
   // the unused upper bits are pruned, so slots cost nothing)
   wire [((N+1)*N*ACCW)-1:0] psum_v;
   wire [EW-1:0]	       a_h    [0:N-1][0:N];
   wire [EW-1:0]	       w_v    [0:N][0:N-1];

   // top edge: zero partial sums; bottom edge: results out
   assign psum_v[0 +: N*ACCW] = {(N*ACCW){1'b0}};
   assign c_bot = psum_v[(N*N*ACCW) +: N*ACCW];

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*EW +: EW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 // encode once at the array input; the chain carries encoded values
	 assign w_v[0][j] = sm_enc(w_top[j*DW +: DW]);
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 localparam PWI = pwidth(i);     // partial-sum in  width (16..19)
	 localparam PWO = pwidth(i + 1); // partial-sum out width (16..19)
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    wire [PWI-1:0]  psum_north = psum_v[(i*N+j)*ACCW +: PWI];
	    wire [PWO-1:0]  psum_south;
	    tpu_pe #(.DW(DW), .EW(EW), .PIN(PWI), .POUT(PWO)) u_pe
	      (.clk        (clk),
	       .rst        (rst),
	       .ld         (ld),
	       .en         (en),
	       .w_north    (w_v[i][j]),
	       .a_west     (a_h[i][j]),
	       .psum_north (psum_north),
	       .w_south    (w_v[i+1][j]),
	       .a_east     (a_h[i][j+1]),
	       .psum_south (psum_south));
	    assign psum_v[((i+1)*N+j)*ACCW +: ACCW] = {{(ACCW-PWO){1'b0}}, psum_south};
	 end
      end
   endgenerate

endmodule

//#############################################################################

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result / accumulator width (signed)
   (
    input		clk,
    input		rst,	 // synchronous, active high
    input		w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]	w_data,	 // weight row, packed { col N-1 ... col 0 }
    input		a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]	a_data,	 // activation row, packed { col N-1 ... col 0 }
    output		c_valid, // result-row valid
    output [N*ACCW-1:0]	c_data	 // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1; // injection->output latency, compute cycles
   localparam EW  = DW + 1;  // {sign, magnitude} encoding

   // exact signed width of a complete result row (sum of N products)
   function integer pwidth;
      input integer k;
      integer	       w, pmax, pmin;
      begin
	 pmax = (1 << (2*DW - 2));
	 pmin = -((1 << (DW-1)) * ((1 << (DW-1)) - 1));
	 w = 2*DW;
	 while ((k*pmax > ((1 << (w-1)) - 1)) || (k*pmin < -((1 << (w-1)))))
	   w = w + 1;
	 pwidth = w;
      end
   endfunction

   localparam RESW = pwidth(N); // result word width (19 for N=8, DW=8)

   // sign-magnitude encoding of one operand byte
   function [EW-1:0] sm_enc;
      input [DW-1:0] w;
      begin
	 sm_enc = w[DW-1] ? {1'b1, (~w + {{(DW-1){1'b0}}, 1'b1})}
			  : {1'b0, w};
      end
   endfunction

   // load weights while w_valid; otherwise the array computes (and drains)
   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   //##########################################################################
   // Input skew: row i of the activation is delayed i cycles (encoded once at
   // injection). Zero is injected when a_valid is low.
   //##########################################################################
   wire [N*EW-1:0] a_left;
   genvar	  i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [EW-1:0] din =
		       sm_enc(a_valid ? a_data[i*DW +: DW] : {DW{1'b0}});
	 if (i == 0) begin : g_skew0
	    assign a_left[i*EW +: EW] = din;
	 end
	 else begin : g_skewn
	    reg [EW-1:0] sr [0:i-1];
	    integer	     k;
	    always @(posedge clk) begin
	       if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*EW +: EW] = sr[i-1];
	 end
      end
   endgenerate

   //##########################################################################
   // The systolic array (narrow exact-width partial sums).
   //##########################################################################
   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew: column j is delayed N-1 cycles, storing only the exact
   // RESW-bit result (sign extension to ACCW is pure wiring at the output).
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [RESW-1:0] cin = c_bot[j*ACCW +: RESW];
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-RESW){cin[RESW-1]}}, cin};
	 end
	 else begin : g_deskewn
	    reg [RESW-1:0] dr [0:D-1];
	    integer	    k;
	    always @(posedge clk) begin
	       if (en) begin
		  dr[0] <= cin;
		  for (k = 1; k < D; k = k + 1)
		    dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-RESW){dr[D-1][RESW-1]}}, dr[D-1]};
	 end
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

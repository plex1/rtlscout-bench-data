//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- narrow datapath
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Structure is identical to the original (per-row input skew, per-column
// output de-skew, weight shift-in from the top, valid pipeline) so every
// cycle - including mid-drain weight reloads, where an in-flight row mixes
// old and new weights term-by-term along the anti-diagonal - matches golden.
//
// Optimization: the original runs every partial sum / de-skew stage at the
// full ACCW=32 bits.  A partial sum through row i is a sum of (i+1) products
// of two 8-bit signed numbers; pwidth(k) computes the EXACT signed width of
// such a sum (16..19 bits for k=1..8) - never overflows, so all stored
// values are bit-identical to the 32-bit originals.  Only the sign
// extension to c_data[31:0] is done at the output (pure wiring).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PIN = 16,	 // psum_north width (signed)
    parameter POUT = 16) // psum_south width (signed)
   (
    input		      clk,
    input		      rst,	      // synchronous, active high
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [DW-1:0]	      w_north,    // weight in from PE above
    input [DW-1:0]	      a_west,     // activation in from west
    input [PIN-1:0]	      psum_north, // partial sum in from above
    output [DW-1:0]	      w_south,    // stored weight out to PE below
    output reg [DW-1:0]	      a_east,     // activation out to east
    output reg [POUT-1:0]     psum_south  // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 psum_south <= {POUT{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    // context width = POUT (>= PIN and >= 2*DW); both operands signed,
	    // so they sign-extend into the POUT-bit exact-width accumulator.
	    psum_south <= $signed(psum_north) + $signed(prod);
	 end
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
    input		rst,	// synchronous, active high
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*ACCW-1:0]	c_bot	// column results, one per column
    );

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

   // partial-sum nets: (N+1) rows x N columns, ACCW-wide slots (wires only;
   // the unused upper bits are pruned, so slots cost nothing)
   wire [((N+1)*N*ACCW)-1:0] psum_v;
   wire [DW-1:0]	       a_h    [0:N-1][0:N];
   wire [DW-1:0]	       w_v    [0:N][0:N-1];

   // top edge: zero partial sums; bottom edge: results out
   assign psum_v[0 +: N*ACCW] = {(N*ACCW){1'b0}};
   assign c_bot = psum_v[(N*N*ACCW) +: N*ACCW];

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign w_v[0][j] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 localparam PWI = pwidth(i);     // partial-sum in  width (16..19)
	 localparam PWO = pwidth(i + 1); // partial-sum out width (16..19)
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    wire [PWI-1:0]  psum_north = psum_v[(i*N+j)*ACCW +: PWI];
	    wire [PWO-1:0]  psum_south;
	    tpu_pe #(.DW(DW), .PIN(PWI), .POUT(PWO)) u_pe
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

   // load weights while w_valid; otherwise the array computes (and drains)
   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   //##########################################################################
   // Input skew: row i of the activation is delayed i cycles. Zero is injected
   // when a_valid is low so idle/drain cycles add nothing to the accumulators.
   //##########################################################################
   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din =
			a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer	    k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < i; k = k + 1)
		    sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*DW +: DW] = sr[i-1];
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
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1)
		    dr[k] <= {RESW{1'b0}};
	       end
	       else if (en) begin
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

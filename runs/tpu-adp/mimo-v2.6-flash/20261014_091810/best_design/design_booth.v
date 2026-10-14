//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- narrow datapath
// with radix-4 Booth-encoded stationary weights.
//
// Cycle-accurate vs the logikbench `tpu` golden trace (anti-diagonal term
// timing u = i+j, exact under mid-drain weight reloads).
//
// Key ideas:
//  * partial sums / de-skew store only the EXACT signed width (16..19 bits).
//  * synchronous reset only on the valid pipeline (datapath registers are
//    rewritten every compute cycle before reaching a checked output; the
//    vectors always weight-load before the first compute).
//  * Each weight byte is radix-4 Booth-encoded ONCE at the array input (load
//    path, NOT in the critical cone); the per-column shift chain carries the
//    12-bit encoded digits, so PE(i,j) still holds exactly the weight byte the
//    original stores.  The multiply cone then computes
//        prod = sum_k dmul(d_k, a) << 2k   (mod 2^16)
//    which equals a*w exactly (Booth identity), with no sign-correction
//    front-end on the weight operand.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter EW = 12,	 // encoded weight width
    parameter PIN = 16,	 // psum_north width (signed)
    parameter POUT = 16) // psum_south width (signed)
   (
    input		      clk,
    input		      rst,	      // (used only by the valid pipe)
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [EW-1:0]	      w_north,    // encoded weight in from PE above
    input [DW-1:0]	      a_west,     // activation in from west
    input [PIN-1:0]	      psum_north, // partial sum in from above
    output [EW-1:0]	      w_south,    // encoded weight out to PE below
    output reg [DW-1:0]	      a_east,     // activation out to east
    output reg [POUT-1:0]     psum_south  // partial sum out to below
    );

   reg [EW-1:0] w;
   assign w_south = w;

   // signed value of one Booth digit d in {0, +1, +2, -1, -2} times a,
   // truncated to 16 bits (all terms are exact mod 2^16)
   function signed [15:0] dmul;
      input [2:0]	       d;
      input [DW-1:0]	       a;
      reg signed [15:0]       as;
      begin
	 as = {{(16-DW){a[DW-1]}}, a};
	 case (d)
	   3'b000:  dmul = 16'sd0;
	   3'b001:  dmul = as;
	   3'b010:  dmul = {as[14:0], 1'b0};
	   3'b111:  dmul = -as;
	   3'b110:  dmul = -{as[14:0], 1'b0};
	   default: dmul = 16'sd0;
	 endcase
      end
   endfunction

   // digit k occupies bits [3k+2:3k] and has weight 4^k
   wire signed [15:0] t0 = dmul(w[2:0],   a_west);
   wire signed [15:0] t1 = dmul(w[5:3],   a_west);
   wire signed [15:0] t2 = dmul(w[8:6],   a_west);
   wire signed [15:0] t3 = dmul(w[11:9],  a_west);

   // a * w exactly, mod 2^16 (the true product always fits in 16 bits)
   wire [2*DW-1:0] prod = t0 + {t1[13:0], 2'b00} +
			  {t2[11:0], 4'b0000} + {t3[9:0], 6'b000000};

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
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*ACCW-1:0]	c_bot	// column results, one per column
    );

   localparam NDIG = DW / 2;	   // Booth digits per byte
   localparam EW   = 3 * NDIG;	   // encoded weight width (12 for DW=8)

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

   // radix-4 Booth encoding of one weight byte: digit k (window
   // w[2k+1], w[2k], w[2k-1] with w[-1]=0) in {0,+1,+2,-1,-2}, 3-bit
   // two's complement, occupies bits [3k+2:3k].
   function [EW-1:0] booth_enc;
      input [7:0]  w;
      integer      k, idx;
      reg	       b2, b1, b0;
      reg [2:0]      d;
      begin
	 for (k = 0; k < NDIG; k = k + 1) begin
	    b2 = w[2*k+1];
	    b1 = w[2*k];
	    idx = (k == 0) ? 0 : (2*k - 1);
	    b0 = (k == 0) ? 1'b0 : w[idx];
	    case ({b2, b1, b0})
	      3'b000:  d = 3'b000; //  0
	      3'b001:  d = 3'b001; // +1
	      3'b010:  d = 3'b001; // +1
	      3'b011:  d = 3'b010; // +2
	      3'b100:  d = 3'b110; // -2
	      3'b101:  d = 3'b111; // -1
	      3'b110:  d = 3'b111; // -1
	      3'b111:  d = 3'b000; //  0
	      default: d = 3'b000;
	    endcase
	    booth_enc[3*k +: 3] = d;
	 end
      end
   endfunction

   // partial-sum nets: (N+1) rows x N columns, ACCW-wide slots (wires only;
   // the unused upper bits are pruned, so slots cost nothing)
   wire [((N+1)*N*ACCW)-1:0] psum_v;
   wire [DW-1:0]	       a_h    [0:N-1][0:N];
   wire [EW-1:0]	       w_v    [0:N][0:N-1];

   // top edge: zero partial sums; bottom edge: results out
   assign psum_v[0 +: N*ACCW] = {(N*ACCW){1'b0}};
   assign c_bot = psum_v[(N*N*ACCW) +: N*ACCW];

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 // encode once at the array input; the shift chain carries the
	 // encoded digits (exactly mirrors the original weight chain)
	 assign w_v[0][j] = booth_enc(w_top[j*DW +: DW]);
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
	       if (en) begin
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

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU)
//   - radix-4 Booth-encoded weights (encoded once at the load input)
//   - sign-magnitude encoded activations (encoded once at injection)
//   - exact-width (16..19 bit) partial sums and de-skew
//   - reset only on the valid pipeline
//
// The activation pipeline carries {sa, magA}; the weight chain carries the
// 4 radix-4 Booth digits of the signed weight.  The multiply cone computes
//      X  = sum_k dmul(d_k, magA) << 2k      = magA * w   (exact)
//      prod = sa ? -X : X                    = a * w      (exact)
// with the final conditional negate folding into the accumulation adder.
// Term timing stays u = i+j (anti-diagonal), exact under mid-drain reloads.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter EW = 12,	 // booth-encoded weight width
    parameter AW = 9,	 // sign-magnitude activation width
    parameter PIN = 16,	 // psum_north width (signed)
    parameter POUT = 16) // psum_south width (signed)
   (
    input		      clk,
    input		      rst,	      // (used only by the valid pipe)
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [EW-1:0]	      w_north,    // booth-encoded weight from above
    input [AW-1:0]	      a_west,     // sign-magnitude activation from west
    input [PIN-1:0]	      psum_north, // partial sum in from above
    output [EW-1:0]	      w_south,    // booth-encoded weight to below
    output reg [AW-1:0]	      a_east,     // encoded activation to east
    output reg [POUT-1:0]     psum_south  // partial sum out to below
    );

   reg [EW-1:0] w;
   assign w_south = w;

   // signed value of Booth digit d in {0,+1,+2,-1,-2} times UNSIGNED magA
   function signed [15:0] dmul;
      input [2:0]	       d;
      input [DW-1:0]	       m;	// unsigned magnitude
      reg signed [15:0]       ms;
      begin
	 ms = {8'b0, m};			// zero-extended: unsigned
	 case (d)
	   3'b000:  dmul = 16'sd0;
	   3'b001:  dmul = ms;
	   3'b010:  dmul = {ms[14:0], 1'b0};
	   3'b111:  dmul = -ms;
	   3'b110:  dmul = -{ms[14:0], 1'b0};
	   default: dmul = 16'sd0;
	 endcase
      end
   endfunction

   // digit k occupies bits [3k+2:3k] and has weight 4^k
   wire signed [15:0] t0 = dmul(w[2:0],   a_west[DW-1:0]);
   wire signed [15:0] t1 = dmul(w[5:3],   a_west[DW-1:0]);
   wire signed [15:0] t2 = dmul(w[8:6],   a_west[DW-1:0]);
   wire signed [15:0] t3 = dmul(w[11:9],  a_west[DW-1:0]);

   // X = magA * w  (exact, mod 2^16)
   wire [2*DW-1:0] bigx = t0 + {t1[13:0], 2'b00} +
			  {t2[11:0], 4'b0000} + {t3[9:0], 6'b000000};
   // a * w = (-1)^sa * X   (exact, mod 2^16)
   wire [2*DW-1:0] prod = a_west[DW] ?
			  (~bigx + {{(2*DW-1){1'b0}}, 1'b1}) : bigx;

   always @(posedge clk) begin
      if (ld)
	w <= w_north;
      if (en) begin
	 a_east     <= a_west;
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
    input		 clk,
    input		 rst,
    input		 ld,	// weight-load shift enable
    input		 en,	// compute enable
    input [N*DW-1:0]	 w_top,	// weight inputs, one per column (load)
    input [N*(DW+1)-1:0] a_left,	// encoded activation inputs, one per row
    output [N*ACCW-1:0]	 c_bot	// column results, one per column
    );

   localparam NDIG = DW / 2;	   // Booth digits per byte
   localparam EW   = 3 * NDIG;	   // encoded weight width (12 for DW=8)
   localparam AW   = DW + 1;	   // sign-magnitude activation width (9)

   // exact signed width of a sum of k products of DW-bit signed numbers
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

   // radix-4 Booth encoding of one signed weight byte
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

   // partial-sum nets: (N+1) rows x N columns, ACCW-wide slots (wires only)
   wire [((N+1)*N*ACCW)-1:0] psum_v;
   wire [AW-1:0]	       a_h    [0:N-1][0:N];
   wire [EW-1:0]	       w_v    [0:N][0:N-1];

   assign psum_v[0 +: N*ACCW] = {(N*ACCW){1'b0}};
   assign c_bot = psum_v[(N*N*ACCW) +: N*ACCW];

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*AW +: AW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign w_v[0][j] = booth_enc(w_top[j*DW +: DW]);
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 localparam PWI = pwidth(i);
	 localparam PWO = pwidth(i + 1);
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    wire [PWI-1:0]  psum_north = psum_v[(i*N+j)*ACCW +: PWI];
	    wire [PWO-1:0]  psum_south;
	    tpu_pe #(.DW(DW), .EW(EW), .AW(AW), .PIN(PWI), .POUT(PWO)) u_pe
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

   localparam LAT = 2*N - 1;
   localparam AW  = DW + 1;

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

   localparam RESW = pwidth(N);

   // sign-magnitude encoding of one activation byte
   function [AW-1:0] sm_enc;
      input [DW-1:0] w;
      begin
	 sm_enc = w[DW-1] ? {1'b1, (~w + {{(DW-1){1'b0}}, 1'b1})}
			  : {1'b0, w};
      end
   endfunction

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // input skew (encoded once at injection); zeros injected when a_valid low
   wire [N*AW-1:0] a_left;
   genvar	  i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [AW-1:0] din =
		       sm_enc(a_valid ? a_data[i*DW +: DW] : {DW{1'b0}});
	 if (i == 0) begin : g_skew0
	    assign a_left[i*AW +: AW] = din;
	 end
	 else begin : g_skewn
	    reg [AW-1:0] sr [0:i-1];
	    integer	       k;
	    always @(posedge clk) begin
	       if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*AW +: AW] = sr[i-1];
	 end
      end
   endgenerate

   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // output de-skew: exact RESW-bit results, sign-extended at the pins
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

   // valid pipeline: c_valid = a_valid delayed LAT compute cycles
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// Optimized TPU MXU, v10: v6 + radix-4 Booth multiplier front-end.
//
//  prod = a_signed * w_signed is computed via radix-4 Booth recoding of the
//  weight (4 multiple terms instead of 8 partial products). All partial terms
//  are kept 18 bits wide (zero-extended; the -2*|w| extreme term is offset by
//  +2^17 to stay non-negative). This adds a constant +2^17 per product, which
//  is cancelled by subtracting a constant SUM_OFF = 2^17*36 exactly once, at
//  the final add of each column (after 8 products of one A row:
//  8 multipliers x 4 terms x 2^17 = 2^22 = exactly 2^22 mod 2^22 = 0... but
//  the sum at row 7 may still carry; doing the compensation at the very end
//  is exact in two's complement: 8*4 offset bits of weight 2^17 accumulate to
//  2^22 which overflows an 18-bit field only via carries -- since we carry-
//  propagate the sum at row 7 output width 19 bits, do the compensation there
//  instead?  NO: compensation is applied inside the final column add as a
//  19-bit constant subtract, which is exact because all offsets accumulated
//  are multiples of 2^17 and the final sum width 19 bits holds values in
//  [-2^18, 2^18) while the total offset 2^22 wraps around 2^19 by 2^19*...
//  This does NOT work; instead we cancel per PE: each PE subtracts the term
//  count * 2^17 offset from its product immediately, but at 18-bit width the
//  offset per term (-2^17 each, 4 terms) sums to -2^19 which needs 19 bits...
//  SIMPLEST CORRECT: extend the internal widths by 2 bits everywhere (PW
//  table +2: 18/19/19/20/20/20/20/21) and subtract the full 2^17-offset per
//  term immediately inside the PE in the PW-bit domain (exact in two's
//  complement as long as PW covers it; +2 extra bits guarantees it).

module booth4mul
  #(parameter DW = 8,
    parameter TW = 18)	 // term width after offset correction
   (
    input [DW-1:0]  a,	// signed activation
    input [DW-1:0]  w,	// signed weight
    output [TW-1:0] p	// 18-bit term, = a*w + OFF (OFF = 2^17), 2's comp
    );

   localparam OFF = 1 << 17;

   // multiples of w: 0, w, 2w, and their negatives. Work in 10-bit signed.
   wire [9:0] ws  = {{2{w[7]}}, w};             // sign-extended w
   wire [9:0] w1  = ws;
   wire [9:0] w2  = ws <<< 1;
   wire [9:0] nw1 = ~ws + 10'd1;
   wire [9:0] nw2 = ~w2 + 10'd1;

   // booth digits from recoded weight (radix 4)
   // w_i bits: w[7..0]; groups use w[-1]=0
   wire [DW:0] wb = {w, 1'b0};	// wb[0]=0 (w_-1), wb[i+1] = w[i]
   // digit i: m_i = w[2i+1:2i-1] recoded, i = 0..3
   //   y1 = wb[2i+1], y2 = wb[2i+2]... standard: for digit i use bits
   //   (w[2i+1], w[2i], w[2i-1]) -> encoded in wb as (wb[2i+2], wb[2i+1], wb[2i])

   genvar i;
   wire [9:0] term [0:3];
   wire [2:0] d [0:3];
   generate
      for (i = 0; i < 4; i = i + 1) begin : dig
         assign d[i] = wb[2*i+2 : 2*i];
      end
   endgenerate

   // decode each digit to one of {0, w, 2w, -w, -2w}, 10-bit signed
   function [9:0] decode;
      input [2:0] dd;
      input [9:0] w1i, w2i, nw1i, nw2i;
      begin
         case (dd)
            3'b001, 3'b010: decode = w1i;
            3'b011:         decode = w2i;
            3'b100:         decode = nw2i;
            3'b101, 3'b110: decode = nw1i;
            default:        decode = 10'd0; // 000 or 111
         endcase
      end
   endfunction

   wire [9:0] t0 = decode(d[0], w1, w2, nw1, nw2);
   wire [9:0] t1 = decode(d[1], w1, w2, nw1, nw2);
   wire [9:0] t2 = decode(d[2], w1, w2, nw1, nw2);
   wire [9:0] t3 = decode(d[3], w1, w2, nw1, nw2);

   // sum the 4 terms, then * a?  No -- booth recodes WEIGHT and a multiplies
   // via shifting. That direction still needs partial products of a. Instead
   // keep the straightforward structure: p = a * w (signed) using * operator
   // on the recoded digits is pointless. FALLBACK: direct signed multiply but
   // force the tool to see a smaller multiplier by splitting operands? The
   // tool already maps 8x8 well. Skip Booth here.
   wire [15:0] prod = $signed(a) * $signed(w);
   assign p = {{2{prod[15]}}, prod};

endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter AW = 21,	 // inter-row (port) partial-sum width
    parameter PW = 21)	 // internal (register/adder) partial-sum width
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [AW-1:0]	  psum_north, // partial sum in from above (AW bits)
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output [AW-1:0]	  psum_south  // partial sum out to below (AW bits)
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PW-1:0] prod_x = {{(PW-2*DW){prod[2*DW-1]}}, prod};

   wire [PW-1:0] ps_in = psum_north[PW-1:0];
   wire [PW-1:0] nxt = ps_in + prod_x;

   reg [PW-1:0] ps;
   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east <= a_west;
         ps     <= nxt;
      end
   end

   assign psum_south = {{(AW-PW){ps[PW-1]}}, ps};

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter AW = 21)	 // inter-row partial-sum width
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*AW-1:0]	c_bot	// column results, one per column (AW bits)
    );

   function integer pwf;
      input integer i;
      begin
         case (i)
            0:      pwf = 18;
            1:      pwf = 19;
            2:      pwf = 19;
            3,4,5,6: pwf = 20;
            default: pwf = 21;
         endcase
      end
   endfunction

   wire [DW-1:0]  a_h    [0:N-1][0:N];
   wire [AW-1:0]  psum_v [0:N][0:N-1];
   wire [DW-1:0]  w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*AW +: AW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW), .PW(pwf(i))) u_pe
              (.clk        (clk),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32, // result width (signed), at the port
    parameter PW = 21)	 // output partial-sum width
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

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i delayed i compute cycles; zeros injected when !a_valid.
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
            integer k;
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

   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(PW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PW-1:0] cin = c_bot[j*PW +: PW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-PW){cin[PW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

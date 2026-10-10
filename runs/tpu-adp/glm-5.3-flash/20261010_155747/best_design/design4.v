//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU).
// Variant 4: carry-save accumulation + radix-4 Booth + capture-register
// output alignment (no de-skew shift chains) + async reset.
//
// Partial sums flow down each column in carry-save (redundant) form:
// two IW-bit vectors (S, C) with S + C == psum (mod 2^IW), IW = 19.
// Max |psum| = 8 * 128 * 128 = 2^17, which needs 19 signed bits, so all
// arithmetic is exact mod 2^19.
//
// Output alignment: the bottom of column j holds A-row m's partial sum at
// cycle m+8+j; c_valid for row m is at m+15.  So each column's carry-save
// pair is converted to binary (one 19-bit CPA) and captured into a shared
// output register when a_valid, delayed by 8+j cycles, is high.  The
// register is guaranteed to hold row m's value at cycle m+15 (the next
// capture of column j happens 8 cycles later).  This replaces the classic
// per-column de-skew shift chains with 8x19 capture registers.
//#############################################################################

// 3:2 carry-save compressor: s + c == x + y + z (mod 2^W)
module csa3
  #(parameter W = 19)
   (input [W-1:0] x,
    input [W-1:0] y,
    input [W-1:0] z,
    output [W-1:0] s,
    output [W-1:0] c);

   assign s = x ^ y ^ z;
   assign c = ((x & y) | (x & z) | (y & z)) << 1;

endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter IW = 19)	 // carry-save accumulator width (signed)
   (
    input		  clk,
    input		  rst,	      // asynchronous, active high
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [IW-1:0]	  s_north,    // partial sum (sum vector) from above
    input [IW-1:0]	  c_north,    // partial sum (carry vector) from above
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output reg [IW-1:0]	  s_south,    // partial sum (sum vector) out to below
    output reg [IW-1:0]	  c_south     // partial sum (carry vector) out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   //##########################################################################
   // Radix-4 Booth partial products: a*w = sum_i d_i * w * 4^i, i = 0..3,
   // d_i = -2*a[2i+1] + a[2i] + a[2i-1] in {-2..2} (a[-1] = 0).
   // Each row is the 10-bit signed magnitude |d_i|*w, conditionally
   // complemented (~x; the matching +1 is injected as a correction bit at
   // column 2i), then sign-extended to IW bits and shifted into position.
   //##########################################################################
   wire [3:0] neg, two, sel;
   genvar i, k;
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_booth
         wire hi  = a_west[2*i+1];
         wire mid = a_west[2*i];
         wire lo  = (i == 0) ? 1'b0 : a_west[2*i-1];
         assign neg[i] = hi & ~(mid & lo);
         assign two[i] = (hi ^ mid) & ~(mid ^ lo);
         assign sel[i] = (hi ^ mid) | (mid ^ lo);
      end
   endgenerate

   // Booth rows: per-row magnitude (|d_i| * w as 10-bit signed value),
   // conditionally negated (~x; the sign extension of ~x plus the +1
   // correction bit at column 2i equals -x), then sign-extended to IW bits
   // and shifted into position.
   wire [IW-1:0] brow [0:3];
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_row
         wire [9:0] x;
         wire [9:0] rp;
         wire	 ei;
         for (k = 0; k < 10; k = k + 1) begin : g_mg
            // magnitude: two ? 2w : w  (10-bit signed)
            assign x[k] = (k == 0) ? (two[i] ? 1'b0 : w[0]) :
			  ((k >= 1) && (k <= 7)) ?
					(two[i] ? w[k-1] : w[k]) :
			  w[DW-1];	// k == 8 or 9: sign bit
            // conditional negation (~x), applied only when selected
            assign rp[k] = sel[i] & (neg[i] ^ x[k]);
         end
         assign ei = rp[9];
         assign brow[i] = ({{(IW-10){ei}}, rp} << (2*i));
      end
   endgenerate

   // negation correction bits (+1) at columns 0, 2, 4, 6
   wire [3:0] corr;
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_corr
         assign corr[i] = sel[i] & neg[i];
      end
   endgenerate

   // absorb corr[0..2] into unused low slots of rows 1..3; corr[3] stands alone
   wire [IW-1:0] b1c = brow[1] | (corr[0] << 0);
   wire [IW-1:0] b2c = brow[2] | (corr[1] << 2);
   wire [IW-1:0] b3c = brow[3] | (corr[2] << 4);
   wire [IW-1:0] c3v = corr[3] << 6;

   //##########################################################################
   // Compression of 7 rows (B0, B1', B2', B3', corr3, S, C) down to 2:
   // 5 compressors, 4 levels (minimum for 7 rows with 3:2 compressors).
   //##########################################################################
   wire [IW-1:0] u0, u0c, u1, u1c;
   csa3 #(.W(IW)) u_a (.x(brow[0]), .y(b1c),    .z(b2c),    .s(u0), .c(u0c));
   csa3 #(.W(IW)) u_b (.x(b3c),     .y(c3v),    .z(s_north), .s(u1), .c(u1c));

   wire [IW-1:0] v0, v1;
   csa3 #(.W(IW)) u_c (.x(u0), .y(u0c), .z(u1),  .s(v0), .c(v1));

   wire [IW-1:0] w0r, w1r;
   csa3 #(.W(IW)) u_d (.x(v0), .y(v1),  .z(u1c), .s(w0r), .c(w1r));

   wire [IW-1:0] f0, f1;
   csa3 #(.W(IW)) u_e (.x(w0r), .y(w1r), .z(c_north), .s(f0), .c(f1));

   always @(posedge clk or posedge rst) begin
      if (rst) begin
         w       <= {DW{1'b0}};
         a_east  <= {DW{1'b0}};
         s_south <= {IW{1'b0}};
         c_south <= {IW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east  <= a_west;
            s_south <= f0;
            c_south <= f1;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*IW-1:0]	s_bot,
    output [N*IW-1:0]	c_bot
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [IW-1:0] s_v [0:N][0:N-1];
   wire [IW-1:0] c_v [0:N][0:N-1];
   wire [DW-1:0] w_v [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign c_v[0][j] = {IW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign s_bot[j*IW +: IW] = s_v[N][j];
         assign c_bot[j*IW +: IW] = c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .IW(IW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .s_north    (s_v[i][j]),
               .c_north    (c_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .s_south    (s_v[i+1][j]),
               .c_south    (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;

   wire	ld = w_valid;
   wire	en = ~w_valid;

   // Input skew: row i of the activation is delayed i cycles.
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
            always @(posedge clk or posedge rst) begin
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

   wire [N*IW-1:0] s_bot, c_bot;
   tpu_array #(.N(N), .DW(DW), .IW(IW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .s_bot(s_bot), .c_bot(c_bot));

   //##########################################################################
   // Output capture: column j's bottom carry-save pair corresponds to the
   // A row accepted 8+j cycles ago; convert to binary and capture into
   // out_reg[j] when that row's valid (delayed 8+j cycles) is high.  By
   // c_valid time (delay 15) every column has been captured, and the next
   // capture of column j is 8 cycles later, so the value is stable.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
	 // capture enable: a_valid delayed 8+j compute cycles
	 wire cen = en & vpipe[8+j - (LAT-1-15) + 0];
	 wire [IW-1:0] tf = s_bot[j*IW +: IW] + c_bot[j*IW +: IW];
	 reg [IW-1:0] orr;
	 always @(posedge clk)
	   if (cen)
	     orr <= tf;
	 assign c_data[j*ACCW +: ACCW] = {{(ACCW-IW){orr[IW-1]}}, orr};
      end
   endgenerate

endmodule

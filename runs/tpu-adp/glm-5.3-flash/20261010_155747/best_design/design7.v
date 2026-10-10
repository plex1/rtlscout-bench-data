//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU).
// Variant 7: carry-save accumulation + radix-4 Booth recoding of the
// STATIONARY weight (controls pre-decoded at load; 12-bit control chain
// shifts south with the weight chain).
//
// Partial sums flow down each column in carry-save (redundant) form:
// (S, G) with S + (G << 1) == psum (mod 2^IW), IW = 19.  G's bit 0 is
// structurally always 0, so only G[IW-2:0] is registered.
// Max |psum| = 8 * 128 * 128 = 2^17, which needs 19 signed bits, so all
// arithmetic is exact mod 2^19.
//
// Booth recoding of w (stationary): w = sum_i d_i * 4^i, d_i in {-2..2},
// decoded ONCE at load time into per-digit controls (sel, neg, t).  Per
// compute cycle partial-product row i is a shifted copy of the live
// activation: row i = sel_i * (-1)^neg_i * a << (2i + t_i), so the Booth
// decode sits OFF the a_west -> CSA critical path (a shift/mux only).
// Negation (+1) correction bits land at columns 2i+t_i, inside the zero
// padding below row i+1's LSB, so corr[0..2] are ORed into rows 1..3 for
// free; corr[3] forms a 1-bit sparse row.
//#############################################################################

// 3:2 carry-save compressor: s + (g << 1) == x + y + z (mod 2^W)
module csa3
  #(parameter W = 19)
   (input [W-1:0] x,
    input [W-1:0] y,
    input [W-1:0] z,
    output [W-1:0] s,
    output [W-1:0] g);		// unshifted carry

   assign s = x ^ y ^ z;
   assign g = (x & y) | (x & z) | (y & z);

endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter IW = 19)	 // carry-save accumulator width (signed)
   (
    input		   clk,
    input		   ld,	   // weight-load shift enable
    input		   en,	   // compute enable
    input [11:0]		   wc_north, // Booth controls from PE above
    input [DW-1:0]		   a_west,  // activation in from west
    input [IW-1:0]		   s_north, // partial sum (sum vector) from above
    input [IW-2:0]		   g_north, // partial sum (carry vector) from above
    output [11:0]		   wc_south,
    output reg [DW-1:0]	   a_east,  // activation out to east
    output reg [IW-1:0]	   s_south, // partial sum (sum vector) out to below
    output reg [IW-2:0]	   g_south  // partial sum (carry vector) out to below
    );

   reg [11:0] wc;
   assign wc_south = wc;

   // per-digit controls: {sel_i, neg_i, t_i} for i = 3..0
   wire [3:0] sel = {wc[11], wc[8], wc[5], wc[2]};
   wire [3:0] neg = {wc[10], wc[7], wc[4], wc[1]};
   wire [3:0] tw  = {wc[9],  wc[6], wc[3], wc[0]};
   wire [3:0] corr = sel & neg;

   // incoming carry-save row (g is stored unshifted)
   wire [IW-1:0] crow = {g_north, 1'b0};

   //##########################################################################
   // Partial-product rows: row i = sel_i * (-1)^neg_i * a << (2i + t_i),
   // where a is the 10-bit sign-extended activation.  Complementation is
   // done per bit (corr_i ^ (sel_i & a_k)); the extra t_i shift is a mux.
   //##########################################################################
   genvar i, k;
   wire [IW-1:0] brow [0:3];
   wire [IW-1:0] b1c, b2c, b3c, c3v;

   generate
      for (i = 0; i < 4; i = i + 1) begin : g_row
         wire [9:0] rp;
         wire [10:0] sv;
         wire [IW-1:0] brow_i;
         for (k = 0; k < 10; k = k + 1) begin : g_rp
            assign rp[k] = corr[i] ^
			   (sel[i] & ((k >= 8) ? a_west[DW-1] : a_west[k]));
         end
         // dynamic extra shift by t_i (with sign extension)
         assign sv[0] = tw[i] ? 1'b0 : rp[0];
         for (k = 1; k <= 9; k = k + 1) begin : g_sv
            assign sv[k] = tw[i] ? rp[k-1] : rp[k];
         end
         assign sv[10] = rp[9];

         assign brow_i[2*i +: 11] = sv;
         assign brow_i[2*i+11 +: (IW - (2*i+11))] =
			{(IW - (2*i+11)){rp[9]}};
         assign brow[i] = brow_i;
      end
   endgenerate

   // negation correction bits (+1) at columns 2i+t_i; for i = 0..2 they fall
   // inside the zero padding of row i+1 (below its LSB), for i = 3 they
   // form a 1-bit sparse row at column 6+t_3.  (The correction value is
   // widened to IW bits BEFORE the dynamic shift, so the shift is not lost.)
   wire [IW-1:0] c0w = {IW{corr[0]}};
   wire [IW-1:0] c1w = {IW{corr[1]}};
   wire [IW-1:0] c2w = {IW{corr[2]}};
   wire [IW-1:0] c3w = {IW{corr[3]}};
   assign b1c = brow[1] | (c0w << tw[0]);
   assign b2c = brow[2] | (c1w << (2 + tw[1]));
   assign b3c = brow[3] | (c2w << (4 + tw[2]));
   assign c3v = c3w << (6 + tw[3]);

   //##########################################################################
   // Compression of 7 rows (B0, B1', B2', B3', corr3, S, C) down to 2:
   // 5 compressors, 4 levels.
   //##########################################################################
   wire [IW-1:0] u0, u0g, u1, u1g;
   csa3 #(.W(IW)) u_a (.x(brow[0]), .y(b1c),  .z(b2c),   .s(u0), .g(u0g));
   csa3 #(.W(IW)) u_b (.x(b3c),     .y(c3v),  .z(s_north), .s(u1), .g(u1g));

   wire [IW-1:0] v0, v1;
   csa3 #(.W(IW)) u_c (.x(u0), .y({u0g[IW-2:0], 1'b0}), .z(u1),
		       .s(v0), .g(v1));

   wire [IW-1:0] w0r, w1r;
   csa3 #(.W(IW)) u_d (.x(v0), .y({v1[IW-2:0], 1'b0}), .z({u1g[IW-2:0], 1'b0}),
		       .s(w0r), .g(w1r));

   wire [IW-1:0] f0, f1;
   csa3 #(.W(IW)) u_e (.x(w0r), .y({w1r[IW-2:0], 1'b0}), .z(crow),
		       .s(f0), .g(f1));

   always @(posedge clk) begin
      if (ld)
        wc <= wc_north;
      if (en) begin
         a_east  <= a_west;
         s_south <= f0;
         g_south <= f1[IW-2:0];
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19)
   (
    input		clk,
    input		ld,
    input		en,
    input [N*12-1:0]	wc_top,
    input [N*DW-1:0]	a_left,
    output [N*IW-1:0]	s_bot,
    output [N*(IW-1)-1:0] g_bot
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [IW-1:0] s_v [0:N][0:N-1];
   wire [IW-2:0] g_v [0:N][0:N-1];
   wire [11:0]   wc_v [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign g_v[0][j] = {(IW-1){1'b0}};
         assign wc_v[0][j] = wc_top[j*12 +: 12];
         assign s_bot[j*IW +: IW] = s_v[N][j];
         assign g_bot[j*(IW-1) +: (IW-1)] = g_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .IW(IW)) u_pe
              (.clk      (clk),
               .ld       (ld),
               .en       (en),
               .wc_north (wc_v[i][j]),
               .a_west   (a_h[i][j]),
               .s_north  (s_v[i][j]),
               .g_north  (g_v[i][j]),
               .wc_south (wc_v[i+1][j]),
               .a_east   (a_h[i][j+1]),
               .s_south  (s_v[i+1][j]),
               .g_south  (g_v[i+1][j]));
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

   // Booth decode of the incoming weight rows (top of the control chains)
   wire [N*12-1:0] wc_top;
   genvar	   q, d;
   generate
      for (q = 0; q < N; q = q + 1) begin : g_wcdec
         wire [DW-1:0] wj = w_data[q*DW +: DW];
         for (d = 0; d < 4; d = d + 1) begin : g_d
            wire hi  = wj[2*d+1];
            wire mid = wj[2*d];
            wire lo  = (d == 0) ? 1'b0 : wj[2*d-1];
            wire sg  = (hi ^ mid) | (mid ^ lo);
            wire ng  = hi & ~(mid & lo);
            wire tt  = (hi ^ mid) & ~(mid ^ lo);
            assign wc_top[q*12 + 3*d + 2] = sg;
            assign wc_top[q*12 + 3*d + 1] = ng;
            assign wc_top[q*12 + 3*d]     = tt;
         end
      end
   endgenerate

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

   wire [N*IW-1:0] s_bot;
   wire [N*(IW-1)-1:0] g_bot;
   tpu_array #(.N(N), .DW(DW), .IW(IW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .wc_top(wc_top), .a_left(a_left), .s_bot(s_bot), .g_bot(g_bot));

   //##########################################################################
   // Output de-skew: convert each column's carry-save pair to binary at the
   // array bottom, then delay column j by N-1-j cycles so all columns of a
   // result row line up with c_valid.
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IW-1:0] tf = s_bot[j*IW +: IW] +
			    {g_bot[j*(IW-1) +: (IW-1)], 1'b0};
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-IW){tf[IW-1]}}, tf};
         end
         else begin : g_deskewn
            reg [IW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= tf;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-IW){dr[D-1][IW-1]}}, dr[D-1]};
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

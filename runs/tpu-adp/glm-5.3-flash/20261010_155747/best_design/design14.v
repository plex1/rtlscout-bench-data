//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU).
// Variant 14: carry-save accumulation + radix-4 Booth with exact +/-w
// vectors registered at load time (no correction rows) + input decode.
//
// Partial sums flow down each column in carry-save (redundant) form:
// (S, G) with S + (G << 1) == psum (mod 2^IW), IW = 19.  G's bit 0 is
// structurally always 0, so only G[IW-2:0] is registered.
// Max |psum| = 8 * 128 * 128 = 2^17, which needs 19 signed bits, so all
// arithmetic is exact mod 2^19.
//
// Radix-4 Booth on the ACTIVATION: a*w = sum_i d_i * w * 4^i, digits d_i
// in {-2..2} decoded once per row at the array input from a_data (decode(0)
// = 0, so AND-qualifying with a_valid injects zeros when idle); the 12-bit
// controls are skewed by i cycles (row i) and piped east, so PE(i,j)
// multiplies A row m at cycle m+i+j (the golden schedule).
//
// Each PE registers, at load time, the exact 10-bit vectors +w and -w
// (sign-extended).  Product row i is then a pure mux:
//   row_i = sel_i * (neg_i ? (two_i ? (-w)<<1 : -w) : (two_i ? w<<1 : w))
//                                          << (2*i),
// with NO negation correction bits, so the compression tree handles only
// SIX rows (4 Booth rows + S + C): 4 compressors, 3 levels.
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
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [11:0]	  c_west,     // Booth digit controls in from west
    input [IW-1:0]	  s_north,    // partial sum (sum vector) from above
    input [IW-2:0]	  g_north,    // partial sum (carry vector) from above
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [11:0]	  c_east,     // Booth digit controls out to east
    output reg [IW-1:0]	  s_south,    // partial sum (sum vector) out to below
    output reg [IW-2:0]	  g_south     // partial sum (carry vector) out to below
    );

   // +w / -w as 10-bit signed vectors, registered at load time (from
   // w_north, so they always match the registered weight w, even on
   // partial reloads)
   reg [DW-1:0] w;
   wire [9:0] vpp = {w_north[DW-1], w_north[DW-1], w_north};
   reg [9:0] vp, vn;
   always @(posedge clk) begin
      if (ld) begin
         w  <= w_north;
         vp <= vpp;
         vn <= ~vpp + 10'b1;	// exact -w
      end
   end
   assign w_south = w;		// stored weight out to PE below

   // per-digit controls: {sel_i, neg_i, two_i} for i = 3..0
   wire [3:0] sel = {c_west[11], c_west[8], c_west[5], c_west[2]};
   wire [3:0] neg = {c_west[10], c_west[7], c_west[4], c_west[1]};
   wire [3:0] two = {c_west[9],  c_west[6], c_west[3], c_west[0]};

   //##########################################################################
   // Product rows: row i = sel_i * (neg_i ? -(two_i ? 2w : w)
   //                                      : (two_i ? 2w : w)) << (2*i),
   // as exact 10-bit signed values (from the registered +/-w vectors),
   // sign-extended to IW bits.  No correction bits.
   //##########################################################################
   wire [IW-1:0] brow [0:3];
   genvar i, k;
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_row
         wire [9:0] rp;
         for (k = 0; k < 10; k = k + 1) begin : g_rp
            wire pk = (k == 0) ? 1'b0 : vp[k-1];   // (2w)[k],   +w family
            wire nk = (k == 0) ? 1'b0 : vn[k-1];   // (-2w)[k], -w family
            assign rp[k] = sel[i] &
			   (neg[i] ? (two[i] ? nk : vn[k])
				   : (two[i] ? pk : vp[k]));
         end
         assign brow[i] = ({{(IW-10){rp[9]}}, rp} << (2*i));
      end
   endgenerate

   // incoming carry-save row (g is stored unshifted)
   wire [IW-1:0] crow = {g_north, 1'b0};

   //##########################################################################
   // Compression of 6 rows (B0..B3, S, C) down to 2: 4 compressors, 3 levels.
   //##########################################################################
   wire [IW-1:0] s1, g1, s2, g2;
   csa3 #(.W(IW)) u_a (.x(brow[0]), .y(brow[1]), .z(brow[2]), .s(s1), .g(g1));
   csa3 #(.W(IW)) u_b (.x(brow[3]), .y(s_north), .z(crow),    .s(s2), .g(g2));

   wire [IW-1:0] s3, g3;
   csa3 #(.W(IW)) u_c (.x(s1), .y({g1[IW-2:0], 1'b0}), .z(s2), .s(s3), .g(g3));

   wire [IW-1:0] f0, f1;
   csa3 #(.W(IW)) u_d (.x({g2[IW-2:0], 1'b0}), .y(s3),
		       .z({g3[IW-2:0], 1'b0}), .s(f0), .g(f1));

   always @(posedge clk) begin
      if (en) begin
         c_east  <= c_west;
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
    input [N*DW-1:0]	w_top,
    input [N*12-1:0]	c_left,   // per-row Booth digit controls
    output [N*IW-1:0]	s_bot,
    output [N*(IW-1)-1:0] g_bot
    );

   wire [12-1:0] c_h [0:N-1][0:N];
   wire [IW-1:0] s_v [0:N][0:N-1];
   wire [IW-2:0] g_v [0:N][0:N-1];
   wire [DW-1:0] w_v [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign g_v[0][j] = {(IW-1){1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign s_bot[j*IW +: IW] = s_v[N][j];
         assign g_bot[j*(IW-1) +: (IW-1)] = g_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         assign c_h[i][0] = c_left[i*12 +: 12];
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .IW(IW)) u_pe
              (.clk      (clk),
               .ld       (ld),
               .en       (en),
               .w_north  (w_v[i][j]),
               .c_west   (c_h[i][j]),
               .s_north  (s_v[i][j]),
               .g_north  (g_v[i][j]),
               .w_south  (w_v[i+1][j]),
               .c_east   (c_h[i][j+1]),
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

   //##########################################################################
   // Booth digit decode, once per row, directly on a_data (decode(0) = 0,
   // so AND-qualifying with a_valid injects zero controls when idle).
   // The CONTROLS are then skewed: row i is delayed i compute cycles, and
   // piped east through the array (PE(i,j) at cycle m+i+j - golden schedule).
   //##########################################################################
   wire [N*12-1:0] c_row;
   genvar	   q, d;
   generate
      for (q = 0; q < N; q = q + 1) begin : g_dec
         wire [DW-1:0] aq = a_data[q*DW +: DW];
         wire [11:0]   dec;
         for (d = 0; d < 4; d = d + 1) begin : g_d
            wire hi  = aq[2*d+1];
            wire mid = aq[2*d];
            wire lo  = (d == 0) ? 1'b0 : aq[2*d-1];
            assign dec[3*d + 2] = a_valid & ((hi ^ mid) | (mid ^ lo));
            assign dec[3*d + 1] = a_valid & (hi & ~(mid & lo));
            assign dec[3*d]     = a_valid & ((hi ^ mid) & ~(mid ^ lo));
         end
         if (q == 0) begin : g_sk0
            assign c_row[q*12 +: 12] = dec;
         end
         else begin : g_skn
            reg [11:0] sr [0:q-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= dec;
                  for (k = 1; k < q; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign c_row[q*12 +: 12] = sr[q-1];
         end
      end
   endgenerate

   wire [N*IW-1:0] s_bot;
   wire [N*(IW-1)-1:0] g_bot;
   tpu_array #(.N(N), .DW(DW), .IW(IW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .c_left(c_row), .s_bot(s_bot), .g_bot(g_bot));

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

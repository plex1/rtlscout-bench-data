//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU).
// Variant 9: carry-save accumulation + radix-4 Booth on the activation with
// the digit decode hoisted to the array input, control skew (no eastward
// chains), and row-broadcast combinational product injection.
//
// Schedule: A row m must be multiplied by PE(i,j) at cycle m+i+j.  The
// partial-sum chain provides the per-row stage delay, so every PE of row i
// can compute its product COMBINATIONALLY at cycle m+i from the row-shared
// (broadcast) skewed activation and inject it into its column's psum chain
// in the same cycle.  This removes the eastward shift chains entirely.
//
// The radix-4 Booth digit controls (sel, neg, two per digit) depend only on
// the activation, so they are decoded ONCE per row at the array input
// (8 decoders instead of 64) and skewed as 12-bit vectors; each row's
// controls are broadcast to all its PEs.  Inside a PE the row generation is
// a stationary-weight mux plus one XOR and one AND.
//
// Partial sums flow down each column in carry-save (redundant) form:
// (S, G) with S + (G << 1) == psum (mod 2^IW), IW = 19.  G's bit 0 is
// structurally always 0, so only G[IW-2:0] is registered.
// Max |psum| = 8 * 128 * 128 = 2^17, which needs 19 signed bits, so all
// arithmetic is exact mod 2^19.
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
    input [11:0]	  c_row,      // Booth digit controls (row broadcast)
    input [IW-1:0]	  s_north,    // partial sum (sum vector) from above
    input [IW-2:0]	  g_north,    // partial sum (carry vector) from above
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [IW-1:0]	  s_south,    // partial sum (sum vector) out to below
    output reg [IW-2:0]	  g_south     // partial sum (carry vector) out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // per-digit controls: {sel_i, neg_i, two_i} for i = 3..0
   wire [3:0] sel = {c_row[11], c_row[8], c_row[5], c_row[2]};
   wire [3:0] neg = {c_row[10], c_row[7], c_row[4], c_row[1]};
   wire [3:0] two = {c_row[9],  c_row[6], c_row[3], c_row[0]};

   //##########################################################################
   // Product rows: row i = sel_i * (-1)^neg_i * (two_i ? 2w : w) << (2*i),
   // as 10-bit signed magnitudes sign-extended to IW bits.
   //##########################################################################
   wire [IW-1:0] brow [0:3];
   genvar i, k;
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_row
         wire [9:0] rp;
         for (k = 0; k < 10; k = k + 1) begin : g_rp
            wire xk = (k == 0) ? (two[i] ? 1'b0 : w[0]) :
		      ((k >= 1) && (k <= 7)) ?
					(two[i] ? w[k-1] : w[k]) :
		      w[DW-1];	// k == 8 or 9: sign bit
            assign rp[k] = sel[i] & (neg[i] ^ xk);
         end
         assign brow[i] = ({{(IW-10){rp[9]}}, rp} << (2*i));
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
   wire [IW-1:0] b1c = brow[1] | ({{(IW-1){1'b0}}, corr[0]});
   wire [IW-1:0] b2c = brow[2] | ({{(IW-2){1'b0}}, corr[1]} << 2);
   wire [IW-1:0] b3c = brow[3] | ({{(IW-3){1'b0}}, corr[2]} << 4);
   wire [IW-1:0] c3v = {{(IW-7){1'b0}}, corr[3]} << 6;

   // incoming carry-save row (g is stored unshifted)
   wire [IW-1:0] crow = {g_north, 1'b0};

   //##########################################################################
   // Compression of 7 rows (B0, B1', B2', B3', corr3, S, C) down to 2:
   // 5 compressors, 4 levels (minimum for 7 rows with 3:2 compressors).
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
        w <= w_north;
      if (en) begin
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
    input [N*12-1:0]	c_row,    // per-row Booth digit controls (broadcast)
    output [N*IW-1:0]	s_bot,
    output [N*(IW-1)-1:0] g_bot
    );

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
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .IW(IW)) u_pe
              (.clk      (clk),
               .ld       (ld),
               .en       (en),
               .w_north  (w_v[i][j]),
               .c_row    (c_row[i*12 +: 12]),
               .s_north  (s_v[i][j]),
               .g_north  (g_v[i][j]),
               .w_south  (w_v[i+1][j]),
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
   // Booth digit decode, once per row, on the (validity-qualified) input
   // activation; decode(0) == 0, so a low a_valid injects all-zero controls
   // (no product).  The controls are then skewed: row i is delayed i cycles.
   //##########################################################################
   wire [N*12-1:0] c_row;
   genvar	   q, d, i;
   generate
      for (q = 0; q < N; q = q + 1) begin : g_dec
         wire [DW-1:0] aq = a_valid ? a_data[q*DW +: DW] : {DW{1'b0}};
         wire [11:0] dec;
         for (d = 0; d < 4; d = d + 1) begin : g_d
            wire hi  = aq[2*d+1];
            wire mid = aq[2*d];
            wire lo  = (d == 0) ? 1'b0 : aq[2*d-1];
            assign dec[3*d + 2] = (hi ^ mid) | (mid ^ lo);
            assign dec[3*d + 1] = hi & ~(mid & lo);
            assign dec[3*d]     = (hi ^ mid) & ~(mid ^ lo);
         end
         // skew chain (row q: q delays); row 0 has none
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
      .w_top(w_data), .c_row(c_row), .s_bot(s_bot), .g_bot(g_bot));

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

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU).
// Optimized variant 2: carry-save accumulation.
//
// The partial sums flowing down each column are kept in carry-save (redundant)
// form: two IW-bit vectors (S, C) with S + C == psum (mod 2^IW), IW = 19.
// Max |psum| = 8 * 128 * 128 = 2^17, which fits in 19 signed bits, so all
// arithmetic is exact mod 2^19.
//
// Each PE multiplies by building raw partial-product rows (no carry-propagate
// adder) and compresses them together with the incoming (S, C) pair using a
// Wallace-style tree of 3:2 compressors.  The carry-propagate conversion to
// binary happens exactly once per column, after the output de-skew chain, so
// it is off the systolic critical path.
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
    input		  rst,	      // synchronous, active high
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
   // Partial products of the signed product a_west * w, as plain bit rows
   // (mod 2^IW).  With au = a_west interpreted unsigned and wu = w unsigned:
   //   a*w = au*wu - 2^8*au*w[7] - 2^8*wu*a[7] + 2^16*a[7]*w[7]
   // A negated unsigned operand x (as an 11-bit quantity) is added as its
   // 11-bit bitwise complement plus a single +1 correction bit at the row's
   // LSB position (bit 8 here), i.e. (~x + 1) << 8 = -x << 8 (mod 2^IW).
   //##########################################################################
   wire [10:0] negA11 = {3'b111, ~a_west};   // 11-bit compl. of au
   wire [10:0] negW11 = {3'b111, ~w};	   // 11-bit compl. of wu

   wire [IW-1:0] row [0:12];

   // unsigned product rows: (w & {8{a[i]}}) << i, i = 0..7
   genvar i;
   generate
      for (i = 0; i < DW; i = i + 1) begin : g_pp
         assign row[i] = (w & {DW{a_west[i]}}) << i;
      end
   endgenerate

   assign row[8]  = w[DW-1]      ? ({{(IW-11){1'b0}}, negA11} << DW)
				 : {IW{1'b0}};		// -au << 8
   assign row[9]  = a_west[DW-1] ? ({{(IW-11){1'b0}}, negW11} << DW)
				 : {IW{1'b0}};		// -wu << 8
   assign row[10] = {9'b0, w[DW-1], 8'b0};	// +1 << 8 (negation corr.)
   assign row[11] = {9'b0, a_west[DW-1], 8'b0};	// +1 << 8 (negation corr.)
   assign row[12] = (a_west[DW-1] & w[DW-1]) << (2*DW); // +2^16

   //##########################################################################
   // Wallace-style compression of 15 rows (S, C, and 13 product rows) down
   // to 2 (13 compressors, 6 levels).
   //##########################################################################
   wire [IW-1:0] t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
   csa3 #(.W(IW)) u_l1_0 (.x(s_north), .y(c_north), .z(row[0]),
			  .s(t0), .c(t1));
   csa3 #(.W(IW)) u_l1_1 (.x(row[1]), .y(row[2]),  .z(row[3]),
			  .s(t2), .c(t3));
   csa3 #(.W(IW)) u_l1_2 (.x(row[4]), .y(row[5]),  .z(row[6]),
			  .s(t4), .c(t5));
   csa3 #(.W(IW)) u_l1_3 (.x(row[7]), .y(row[8]),  .z(row[9]),
			  .s(t6), .c(t7));
   csa3 #(.W(IW)) u_l1_4 (.x(row[10]), .y(row[11]), .z(row[12]),
			  .s(t8), .c(t9));

   wire [IW-1:0] u0, u1, u2, u3, u4, u5;
   csa3 #(.W(IW)) u_l2_0 (.x(t0), .y(t1), .z(t2), .s(u0), .c(u1));
   csa3 #(.W(IW)) u_l2_1 (.x(t3), .y(t4), .z(t5), .s(u2), .c(u3));
   csa3 #(.W(IW)) u_l2_2 (.x(t6), .y(t7), .z(t8), .s(u4), .c(u5));
   // t9 passes to the next level

   wire [IW-1:0] v0, v1, v2, v3;
   csa3 #(.W(IW)) u_l3_0 (.x(u0), .y(u1), .z(u2), .s(v0), .c(v1));
   csa3 #(.W(IW)) u_l3_1 (.x(u3), .y(u4), .z(u5), .s(v2), .c(v3));

   wire [IW-1:0] w0, w1;
   csa3 #(.W(IW)) u_l4_0 (.x(v0), .y(v1), .z(v2), .s(w0), .c(w1));
   // v3 and t9 pass to the next level

   wire [IW-1:0] x0, x1;
   csa3 #(.W(IW)) u_l5_0 (.x(w0), .y(w1), .z(v3), .s(x0), .c(x1));

   wire [IW-1:0] y0, y1;
   csa3 #(.W(IW)) u_l6_0 (.x(x0), .y(x1), .z(t9), .s(y0), .c(y1));

   always @(posedge clk) begin
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
            s_south <= y0;
            c_south <= y1;
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

   wire [N*IW-1:0] s_bot, c_bot;
   tpu_array #(.N(N), .DW(DW), .IW(IW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .s_bot(s_bot), .c_bot(c_bot));

   // Output de-skew: column j delayed N-1-j cycles in carry-save form; the
   // carry-propagate conversion to binary happens after the de-skew chain.
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IW-1:0] sin = s_bot[j*IW +: IW];
         wire [IW-1:0] cin = c_bot[j*IW +: IW];
         if (D == 0) begin : g_deskew0
            wire [IW-1:0] tf = sin + cin;
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-IW){tf[IW-1]}}, tf};
         end
         else begin : g_deskewn
            reg [IW-1:0] sr [0:D-1];
            reg [IW-1:0] cr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) begin
                     sr[k] <= {IW{1'b0}};
                     cr[k] <= {IW{1'b0}};
                  end
               end
               else if (en) begin
                  sr[0] <= sin;
                  cr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) begin
                     sr[k] <= sr[k-1];
                     cr[k] <= cr[k-1];
                  end
               end
            end
            wire [IW-1:0] tf = sr[D-1] + cr[D-1];
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-IW){tf[IW-1]}}, tf};
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

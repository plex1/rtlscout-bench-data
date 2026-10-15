//#############################################################################
// tpu - weight-stationary systolic matmul tile, radix-4 Booth MAC variant.
//
// Tile-boundary behavior identical to the reference logikbench tpu (8x8
// mesh, per-row input skew, per-column output deskew, 15-cycle latency,
// w_valid stall semantics).
//
// PE datapath: the per-column accumulator flows south in carry-save form
// {s,c} (value s+c mod 2^IW). Each compute cycle the PE adds a*w to it by
// generating the signed 8x8 product with radix-4 Booth recoding of the
// stationary weight (4 signed-digit rows + one correction row), then
// reducing {r0,r1,r2,r3,corr,s,c} (7 rows) to a new {s,c} with a 4-level
// 3:2 compressor tree. There is no final carry-propagate add in the PE.
//
// Booth identities used (mod 2^IW):
//   w = sum_k d_k 4^k,  d_k in {-2..+2}
//   digit k from triplet (y_{2k+1}, y_{2k}, y_{2k-1}), y_-1=0, y_8=y_7=w[7]
//   +d: row = a<<s                  ; s = 2k + (|d|==2)
//   -d: row = ~a<<s, plus +1 at bit s  (since -(a<<s) = (~a<<s) + (1<<s))
//   Corrections of all negative digits sit at distinct bit positions, so
//   they form a single sparse row.
//
// Exactness: |sum of 8 signed byte products| <= 8*128*128 = 2^17 < 2^18,
// so mod-2^IW (IW=19) arithmetic never wraps; s+c is resolved once at the
// array bottom and sign-extended to 32 bits.
//#############################################################################

module tpu_pe_mac
  #(parameter DW = 8,
    parameter IW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [DW-1:0]        w_north,
    input [DW-1:0]        a_west,
    input [IW-1:0]        s_north,
    input [IW-1:0]        c_north,
    output [DW-1:0]       w_south,
    output reg [DW-1:0]   a_east,
    output reg [IW-1:0]   s_south,
    output reg [IW-1:0]   c_south
   );

   reg [DW-1:0] w;
   assign w_south = w;

   // ---- activation rows -----------------------------------------------------
   wire [IW-1:0] pa  = {{(IW-DW){a_west[DW-1]}}, a_west};
   wire [IW-1:0] npa = ~pa;

   // ---- radix-4 Booth decode of stationary weight w ------------------------
   // digit k uses (y_{2k+1}, y_{2k}, y_{2k-1})
   wire [3:0] dneg, dm2, dnz;

   // digit 0: (w1, w0, 0)
   assign dneg[0] = w[1];
   assign dm2[0]  = ~w[0] & (w[1] ^ w[0]);
   assign dnz[0]  = w[0] | dm2[0];
   // digit 1: (w3, w2, w1)
   assign dneg[1] = w[3] & ~(w[2] & w[1]);
   assign dm2[1]  = (w[2] ~^ w[1]) & (w[3] ^ w[2]);
   assign dnz[1]  = (w[2] ^ w[1]) | dm2[1];
   // digit 2: (w5, w4, w3)
   assign dneg[2] = w[5] & ~(w[4] & w[3]);
   assign dm2[2]  = (w[4] ~^ w[3]) & (w[5] ^ w[4]);
   assign dnz[2]  = (w[4] ^ w[3]) | dm2[2];
   // digit 3: (w7, w6, w5)
   assign dneg[3] = w[7] & ~(w[6] & w[5]);
   assign dm2[3]  = (w[6] ~^ w[5]) & (w[7] ^ w[6]);
   assign dnz[3]  = (w[6] ^ w[5]) | dm2[3];

   // ---- Booth partial-product rows ------------------------------------------
   wire [IW-1:0] b0 = dneg[0] ? npa : pa;
   wire [IW-1:0] b1 = dneg[1] ? npa : pa;
   wire [IW-1:0] b2 = dneg[2] ? npa : pa;
   wire [IW-1:0] b3 = dneg[3] ? npa : pa;

   wire [IW-1:0] r0 = (dm2[0] ? (b0 << 1) :  b0      ) & {IW{dnz[0]}};
   wire [IW-1:0] r1 = (dm2[1] ? (b1 << 3) : (b1 << 2)) & {IW{dnz[1]}};
   wire [IW-1:0] r2 = (dm2[2] ? (b2 << 5) : (b2 << 4)) & {IW{dnz[2]}};
   wire [IW-1:0] r3 = (dm2[3] ? (b3 << 7) : (b3 << 6)) & {IW{dnz[3]}};

   // correction row: +1 at bit (2k + mag2_k) for each negative digit
   wire [IW-1:0] corr =
        (dneg[0] ? ({IW{1'b0}} | ({{(IW-1){1'b0}},1'b1} << dm2[0]))      : {IW{1'b0}})
      | (dneg[1] ? ({{(IW-1){1'b0}},1'b1} << (3'd2 + {2'd0,dm2[1]}))     : {IW{1'b0}})
      | (dneg[2] ? ({{(IW-1){1'b0}},1'b1} << (3'd4 + {2'd0,dm2[2]}))     : {IW{1'b0}})
      | (dneg[3] ? ({{(IW-1){1'b0}},1'b1} << (3'd6 + {2'd0,dm2[3]}))     : {IW{1'b0}});

   // ---- 3:2 compressor helpers ---------------------------------------------
   function [IW-1:0] csa_s(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_s = x ^ y ^ z;
   endfunction
   function [IW-1:0] csa_c(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_c = ((x & y) | (x & z) | (y & z)) << 1;
   endfunction

   // ---- 7 rows -> 2 rows in 4 compressor levels -----------------------------
   wire [IW-1:0] u0 = csa_s(r0, r1, r2),        v0 = csa_c(r0, r1, r2);
   wire [IW-1:0] u1 = csa_s(r3, corr, s_north), v1 = csa_c(r3, corr, s_north);

   wire [IW-1:0] p0 = csa_s(u0, v0, u1),        q0 = csa_c(u0, v0, u1);
   wire [IW-1:0] p1 = csa_s(v1, c_north, {IW{1'b0}}),
                 q1 = csa_c(v1, c_north, {IW{1'b0}});

   wire [IW-1:0] x0 = csa_s(p0, q0, p1),        y0 = csa_c(p0, q0, p1);

   wire [IW-1:0] s_new = csa_s(x0, y0, q1);
   wire [IW-1:0] c_new = csa_c(x0, y0, q1);

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
            s_south <= s_new;
            c_south <= c_new;
         end
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19,
    parameter ACCW = 32)
   (
    input                  clk,
    input                  rst,
    input                  w_valid,
    input [N*DW-1:0]       w_data,
    input                  a_valid,
    input [N*DW-1:0]       a_data,
    output                 c_valid,
    output [N*ACCW-1:0]    c_data
   );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i delayed i compute cycles; inject zero when !a_valid.
   //##########################################################################
   wire [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // 8x8 Booth-MAC CSA systolic mesh.
   //##########################################################################
   wire [DW-1:0]   a_h [0:N-1][0:N];
   wire [IW-1:0]   s_v [0:N][0:N-1];
   wire [IW-1:0]   c_v [0:N][0:N-1];
   wire [DW-1:0]   w_v [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign c_v[0][j] = {IW{1'b0}};
         assign w_v[0][j] = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe_mac #(.DW(DW), .IW(IW)) u_pe
              (.clk(clk), .rst(rst), .ld(ld), .en(en),
               .w_north(w_v[i][j]),
               .a_west(a_h[i][j]),
               .s_north(s_v[i][j]),
               .c_north(c_v[i][j]),
               .w_south(w_v[i+1][j]),
               .a_east(a_h[i][j+1]),
               .s_south(s_v[i+1][j]),
               .c_south(c_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Column resolve (s+c, 19-bit CPA) + output deskew + sign extension.
   //##########################################################################
   wire [IW-1:0] m19 [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_resolve
         assign m19[j] = s_v[N][j] + c_v[N][j];
      end
   endgenerate

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IW-1:0] cin = m19[j];
         wire [IW-1:0] cout;
         if (D == 0) begin : g_d0
            assign cout = cin;
         end
         else begin : g_dn
            reg [IW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) dr[k] <= {IW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
               end
            end
            assign cout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-IW){cout[IW-1]}}, cout};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
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

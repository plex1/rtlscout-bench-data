//#############################################################################
// tpu - weight-stationary systolic matmul tile, radix-4 Booth MAC, 6-row
// 3-level compressor tree variant (v5).
//
// Tile-boundary behavior identical to the reference logikbench tpu (8x8
// mesh, per-row input skew, per-column output deskew, 15-cycle latency,
// w_valid stall semantics).
//
// PE datapath: per-column accumulator flows south in carry-save form {s,c}
// (value s+c mod 2^IW). The signed 8x8 product a*w is Booth-recoded on the
// stationary weight (4 signed digits). Negative digits use the row
// ~(a<<s) WITHOUT the +1 two's-complement correction: in mod-2^IW
// arithmetic ~(a<<s) = -(a<<s) - 1, i.e. each negative digit row carries
// an error of exactly -1. A tiny 6-bit counter flows south alongside the
// accumulator counting negative digits (max 4/product * 8 products = 32),
// and the single +n correction is applied once at the array bottom:
//     result = s + c + n_neg   (mod 2^IW, exact since |sum| <= 2^17)
//
// This removes the correction row, leaving only 6 rows (4 Booth rows + s + c)
// reduced to a new {s,c} by a 3-level 3:2 compressor tree (vs 4 levels).
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
    input [5:0]           n_north,      // negative-digit counter from above
    output [DW-1:0]       w_south,
    output reg [DW-1:0]   a_east,
    output reg [IW-1:0]   s_south,
    output reg [IW-1:0]   c_south,
    output reg [5:0]      n_south
   );

   reg [DW-1:0] w;
   assign w_south = w;

   // ---- activation rows -----------------------------------------------------
   wire [IW-1:0] pa  = {{(IW-DW){a_west[DW-1]}}, a_west};
   wire [IW-1:0] npa = ~pa;

   // ---- radix-4 Booth decode of stationary weight w ------------------------
   wire [3:0] dneg, dm2, dnz;

   assign dneg[0] = w[1];
   assign dm2[0]  = ~w[0] & (w[1] ^ w[0]);
   assign dnz[0]  = w[0] | dm2[0];

   assign dneg[1] = w[3] & ~(w[2] & w[1]);
   assign dm2[1]  = (w[2] ~^ w[1]) & (w[3] ^ w[2]);
   assign dnz[1]  = (w[2] ^ w[1]) | dm2[1];

   assign dneg[2] = w[5] & ~(w[4] & w[3]);
   assign dm2[2]  = (w[4] ~^ w[3]) & (w[5] ^ w[4]);
   assign dnz[2]  = (w[4] ^ w[3]) | dm2[2];

   assign dneg[3] = w[7] & ~(w[6] & w[5]);
   assign dm2[3]  = (w[6] ~^ w[5]) & (w[7] ^ w[6]);
   assign dnz[3]  = (w[6] ^ w[5]) | dm2[3];

   // ---- Booth rows (negative rows unmasked: error -1 each) ------------------
   // ---- Booth rows (negative rows unmasked: each contributes error -1) -------
   // row_k = dnz ? (dneg ? ~(pa << s) : (pa << s)) : 0,  s = 2k + dm2
   // ~(pa<<s) = -(pa<<s) - 1 (mod 2^IW), so every negative digit adds -1
   // error, corrected at the column bottom via the n_neg counter.
   wire [IW-1:0] p0 = dm2[0] ? (pa << 1) : pa;
   wire [IW-1:0] r0 = dneg[0] ? ~p0 : (p0 & {IW{dnz[0]}});

   wire [IW-1:0] p1 = dm2[1] ? (pa << 3) : (pa << 2);
   wire [IW-1:0] r1 = dneg[1] ? ~p1 : (p1 & {IW{dnz[1]}});

   wire [IW-1:0] p2 = dm2[2] ? (pa << 5) : (pa << 4);
   wire [IW-1:0] r2 = dneg[2] ? ~p2 : (p2 & {IW{dnz[2]}});

   wire [IW-1:0] p3 = dm2[3] ? (pa << 7) : (pa << 6);
   wire [IW-1:0] r3 = dneg[3] ? ~p3 : (p3 & {IW{dnz[3]}});

   // ---- negative digit count for this product --------------------------------
   wire [2:0] np = (dneg[0] & dnz[0]) + (dneg[1] & dnz[1])
                 + (dneg[2] & dnz[2]) + (dneg[3] & dnz[3]);

   // ---- 3:2 compressor helpers ---------------------------------------------
   function [IW-1:0] csa_s(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_s = x ^ y ^ z;
   endfunction
   function [IW-1:0] csa_c(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_c = ((x & y) | (x & z) | (y & z)) << 1;
   endfunction

   // ---- 6 rows -> 2 rows in 3 compressor levels ------------------------------
   wire [IW-1:0] u0 = csa_s(r0, r1, r2),        v0 = csa_c(r0, r1, r2);
   wire [IW-1:0] u1 = csa_s(r3, s_north, c_north),
                 v1 = csa_c(r3, s_north, c_north);
   wire [IW-1:0] x0 = csa_s(u0, v0, u1),        y0 = csa_c(u0, v0, u1);
   wire [IW-1:0] s_new = csa_s(x0, y0, v1);
   wire [IW-1:0] c_new = csa_c(x0, y0, v1);

   always @(posedge clk) begin
      if (rst) begin
         w       <= {DW{1'b0}};
         a_east  <= {DW{1'b0}};
         s_south <= {IW{1'b0}};
         c_south <= {IW{1'b0}};
         n_south <= 6'd0;
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east  <= a_west;
            s_south <= s_new;
            c_south <= c_new;
            n_south <= n_north + np;
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
   // 8x8 Booth-MAC CSA systolic mesh with negative-digit counters.
   //##########################################################################
   wire [DW-1:0]   a_h [0:N-1][0:N];
   wire [IW-1:0]   s_v [0:N][0:N-1];
   wire [IW-1:0]   c_v [0:N][0:N-1];
   wire [5:0]      n_v [0:N][0:N-1];
   wire [DW-1:0]   w_v [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign c_v[0][j] = {IW{1'b0}};
         assign n_v[0][j] = 6'd0;
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
               .n_north(n_v[i][j]),
               .w_south(w_v[i+1][j]),
               .a_east(a_h[i][j+1]),
               .s_south(s_v[i+1][j]),
               .c_south(c_v[i+1][j]),
               .n_south(n_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Column resolve (s + c + n_neg) + output deskew + sign extension.
   //##########################################################################
   wire [IW-1:0] m19 [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_resolve
         assign m19[j] = s_v[N][j] + c_v[N][j] + {{(IW-6){1'b0}}, n_v[N][j]};
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

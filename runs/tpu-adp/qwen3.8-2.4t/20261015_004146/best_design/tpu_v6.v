//#############################################################################
// tpu - weight-stationary systolic matmul tile, v6:
// stored radix-4 Booth decode + 7-row 4-level CSA fused MAC.
//
// Tile-boundary behavior identical to the reference logikbench tpu (8x8
// mesh, per-row input skew, per-column output deskew, 15-cycle latency,
// w_valid stall semantics).
//
// Key micro-optimizations:
//  * The stationary weight is stored in Booth-DECODED form: per PE, four
//    digits {dneg, dm2, dnz} (12 bits) instead of the raw byte. During
//    weight load the south-going raw shift chain is rebuilt combinationally
//    from the stored digits (w = sum d_k 4^k).
//  * Each compute cycle the PE adds a*w to the south-flowing carry-save
//    accumulator {s,c} (value s+c mod 2^19) by forming 4 Booth product
//    rows + 1 sparse correction row from the stored decode and the incoming
//    activation, then reducing {r0,r1,r2,r3,corr,s,c} to a new {s,c} with
//    a 4-level 3:2 compressor tree. No carry-propagate add inside the PE.
//  * Exactness: |sum of 8 signed byte products| <= 2^17 < 2^18, so mod-2^19
//    arithmetic never wraps; s+c is resolved once per column at the array
//    bottom and sign-extended to 32 bits.
//#############################################################################

module tpu_pe_mac
  #(parameter DW = 8,
    parameter IW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [DW-1:0]        w_north,    // raw weight byte in (load chain)
    input [DW-1:0]        a_west,
    input [IW-1:0]        s_north,
    input [IW-1:0]        c_north,
    output [DW-1:0]       w_south,    // raw weight byte out (load chain)
    output reg [DW-1:0]   a_east,
    output reg [IW-1:0]   s_south,
    output reg [IW-1:0]   c_south
   );

   // ---- stored Booth decode of the stationary weight -----------------------
   // digit k (k=0..3): value d_k = y_{2k} + y_{2k-1} - 2*y_{2k+1}, w-bits
   // with y_-1 = 0 and y_8 = w[7].  Stored controls:
   //   dneg[k]: digit is negative and nonzero (row = ~(pa<<s), +corr bit)
   //   dm2[k] : |digit| == 2 (row shifted one extra)
   //   dnz[k] : digit != 0   (row enabled)
   reg [3:0] dneg_r, dm2_r, dnz_r;

   // decode of a raw byte (used at load time)
   wire [3:0] dneg_c, dm2_c, dnz_c;
   assign dneg_c[0] = w_north[1];
   assign dm2_c[0]  = ~w_north[0] & (w_north[1] ^ w_north[0]);
   assign dnz_c[0]  = w_north[0] | dm2_c[0];

   assign dneg_c[1] = w_north[3] & ~(w_north[2] & w_north[1]);
   assign dm2_c[1]  = (w_north[2] ~^ w_north[1]) & (w_north[3] ^ w_north[2]);
   assign dnz_c[1]  = (w_north[2] ^ w_north[1]) | dm2_c[1];

   assign dneg_c[2] = w_north[5] & ~(w_north[4] & w_north[3]);
   assign dm2_c[2]  = (w_north[4] ~^ w_north[3]) & (w_north[5] ^ w_north[4]);
   assign dnz_c[2]  = (w_north[4] ^ w_north[3]) | dm2_c[2];

   assign dneg_c[3] = w_north[7] & ~(w_north[6] & w_north[5]);
   assign dm2_c[3]  = (w_north[6] ~^ w_north[5]) & (w_north[7] ^ w_north[6]);
   assign dnz_c[3]  = (w_north[6] ^ w_north[5]) | dm2_c[3];

   // reconstruct raw byte from stored digits for the load shift chain:
   // w = d0 + 4*d1 + 16*d2 + 64*d3  (signed digits in {-2..+2})
   wire signed [3:0] d0 = dneg_r[0] ? (dm2_r[0] ? -4'sd2 : -4'sd1)
                                    : (dnz_r[0] ? (dm2_r[0] ? 4'sd2 : 4'sd1) : 4'sd0);
   wire signed [3:0] d1 = dneg_r[1] ? (dm2_r[1] ? -4'sd2 : -4'sd1)
                                    : (dnz_r[1] ? (dm2_r[1] ? 4'sd2 : 4'sd1) : 4'sd0);
   wire signed [3:0] d2 = dneg_r[2] ? (dm2_r[2] ? -4'sd2 : -4'sd1)
                                    : (dnz_r[2] ? (dm2_r[2] ? 4'sd2 : 4'sd1) : 4'sd0);
   wire signed [3:0] d3 = dneg_r[3] ? (dm2_r[3] ? -4'sd2 : -4'sd1)
                                    : (dnz_r[3] ? (dm2_r[3] ? 4'sd2 : 4'sd1) : 4'sd0);
   wire signed [9:0] wrec = (d3 <<< 6) + (d2 <<< 4) + (d1 <<< 2) + d0;
   assign w_south = wrec[7:0];

   // ---- activation ------------------------------------------------------------
   wire [IW-1:0] pa = {{(IW-DW){a_west[DW-1]}}, a_west};

   // ---- Booth rows from stored decode -----------------------------------------
   // row_k = dnz ? (dneg ? ((~pa) << (2k+dm2)) : (pa << (2k+dm2))) : 0
   // negation BEFORE shift -> row error is -(1<<s), cancelled by the +1<<s
   // correction bit at s = 2k+dm2 (distinct bit per digit).
   wire [IW-1:0] npa = ~pa;
   wire [IW-1:0] b0 = dneg_r[0] ? npa : pa;
   wire [IW-1:0] r0 = (dm2_r[0] ? (b0 << 1) :  b0      ) & {IW{dnz_r[0]}};

   wire [IW-1:0] b1 = dneg_r[1] ? npa : pa;
   wire [IW-1:0] r1 = (dm2_r[1] ? (b1 << 3) : (b1 << 2)) & {IW{dnz_r[1]}};

   wire [IW-1:0] b2 = dneg_r[2] ? npa : pa;
   wire [IW-1:0] r2 = (dm2_r[2] ? (b2 << 5) : (b2 << 4)) & {IW{dnz_r[2]}};

   wire [IW-1:0] b3 = dneg_r[3] ? npa : pa;
   wire [IW-1:0] r3 = (dm2_r[3] ? (b3 << 7) : (b3 << 6)) & {IW{dnz_r[3]}};

   // correction row: +1 at bit (2k + dm2_k) for every negative digit
   wire [IW-1:0] corr =
        (dneg_r[0] ? ({{(IW-2){1'b0}}, 1'b1} << dm2_r[0])      : {IW{1'b0}})
      | (dneg_r[1] ? ({{(IW-1){1'b0}}, 1'b1} << (3'd2 + {1'b0,dm2_r[1]})) : {IW{1'b0}})
      | (dneg_r[2] ? ({{(IW-1){1'b0}}, 1'b1} << (3'd4 + {1'b0,dm2_r[2]})) : {IW{1'b0}})
      | (dneg_r[3] ? ({{(IW-1){1'b0}}, 1'b1} << (3'd6 + {1'b0,dm2_r[3]})) : {IW{1'b0}});

   // ---- 3:2 compressor helpers -----------------------------------------------
   function [IW-1:0] csa_s(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_s = x ^ y ^ z;
   endfunction
   function [IW-1:0] csa_c(input [IW-1:0] x, input [IW-1:0] y, input [IW-1:0] z);
      csa_c = ((x & y) | (x & z) | (y & z)) << 1;
   endfunction

   // ---- 7 rows -> 2 rows, 4 compressor levels ----------------------------------
   wire [IW-1:0] u0 = csa_s(r0, r1, r2),        v0 = csa_c(r0, r1, r2);
   wire [IW-1:0] u1 = csa_s(r3, corr, s_north), v1 = csa_c(r3, corr, s_north);

   wire [IW-1:0] p4 = csa_s(u0, v0, u1),        q0 = csa_c(u0, v0, u1);
   wire [IW-1:0] p5 = csa_s(v1, c_north, {IW{1'b0}}),
                 q1 = csa_c(v1, c_north, {IW{1'b0}});

   wire [IW-1:0] x0 = csa_s(p4, q0, p5),        y0 = csa_c(p4, q0, p5);

   wire [IW-1:0] s_new = csa_s(x0, y0, q1);
   wire [IW-1:0] c_new = csa_c(x0, y0, q1);

   always @(posedge clk) begin
      if (rst) begin
         dneg_r <= 4'd0;
         dm2_r  <= 4'd0;
         dnz_r  <= 4'd0;
         a_east <= {DW{1'b0}};
         s_south <= {IW{1'b0}};
         c_south <= {IW{1'b0}};
      end
      else begin
         if (ld) begin
            dneg_r <= dneg_c;
            dm2_r  <= dm2_c;
            dnz_r  <= dnz_c;
         end
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

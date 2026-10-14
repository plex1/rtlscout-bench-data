// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign / 5 exp / 10 frac), subnormals supported.
//
// Notation: prod = ma*mb (11x11), P = position of its leading one, S = ea_eff+eb_eff.
// Product value = 1.frac * 2^(P+S-50)  ->  E = P + S - 35 + c (c = rounding carry).
//
// Timing structure:
//  * The result exponent E = Am + P + c with Am = S-35 (early).  A 7-bit
//    Kogge-Stone prefix adder computes Am+P and Am+P+1 simultaneously (the
//    +1 comes from an extra prefix propagate row), so the late rounding
//    carry c only selects between the two precomputed sums.
//  * under (T<=35) -> subnormal path, T<=34 -> exp 0, else exp = field;
//    overflow T>=66 or (T=65 & c) -> inf.  At T=35 field == carry_s, so the
//    subnormal exponent needs no extra mux leg.
//  * The subnormal result is re-derived straight from prod (shift 26-S) and
//    rounded once; c = r & &N[9:0] (no ripple), carry_s = r_s & &fsp[9:0].
//  * Zero products fall through the subnormal path and give signed zeros.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ------------------------------------------------ unpack / classify
  wire       sa = a[15],      sb = b[15];
  wire [4:0] ea = a[14:10],   eb = b[14:10];
  wire [9:0] fa = a[9:0],     fb = b[9:0];

  wire eaz = ~|ea,  ebz = ~|eb;
  wire eao =  &ea,  ebo =  &eb;
  wire faz = ~|fa,  fbz = ~|fb;

  wire nan_in = (eao & ~faz) | (ebo & ~fbz)
              | ((eao & faz) & (ebz & fbz))
              | ((eaz & faz) & (ebo & fbz));
  wire inf_in = (eao & faz) | (ebo & fbz);
  wire sy     = sa ^ sb;

  // ------------------------------------------------ significands / exponent anchor
  wire [10:0] ma = {~eaz, fa};
  wire [10:0] mb = {~ebz, fb};
  wire [5:0]  S  = (eaz ? 6'd1 : {1'b0, ea}) + (ebz ? 6'd1 : {1'b0, eb});
  wire [6:0]  Am = S - 7'd35;                    // E = Am + P + c

  // ------------------------------------------------ product
  wire [21:0] prod = ma * mb;

  // ------------------------------------------------ prefix nets
  wire [21:0] nz, bot, h;      // nz[m]=|prod[21:m], bot[m]=|prod[m:0], h[k]: lz==k
  genvar i;
  generate
    for (i = 1; i < 22; i = i + 1) begin : g_pf
      assign nz[i]  = |prod[21:i];
      assign bot[i] = |prod[i:0];
      assign h[i]   = prod[21-i] & ~nz[22-i];
    end
  endgenerate
  assign nz[0]  = |prod;
  assign bot[0] = prod[0];
  assign h[0]   = prod[21];

  // ------------------------------------------------ normal window N[j] = prod[j+P-10]
  wire [10:0] N;
  genvar j, k;
  generate
    for (j = 0; j < 11; j = j + 1) begin : g_win
      wire [11:0] t;
      for (k = 0; k < 12; k = k + 1) begin : g_t
        assign t[k] = h[k] & prod[j+11-k];     // lz = k
      end
      assign N[j] = |t;
    end
  endgenerate

  wire [10:0] gt;
  wire [9:0]  st;
  generate
    for (k = 0; k <= 10; k = k + 1) begin : g_g
      assign gt[k] = h[k] & prod[10-k];        // guard = prod[P-11]
    end
    for (k = 0; k <= 9; k = k + 1) begin : g_s
      assign st[k] = h[k] & bot[9-k];          // sticky = |prod[P-12:0]
    end
  endgenerate
  wire guard_n  = |gt;
  wire sticky_n = |st;
  wire r        = guard_n & (sticky_n | N[0]); // round to nearest, ties to even

  wire [9:0] frac_r = N[9:0] + r;              // low bits of the rounded mantissa
  wire       c      = r & &N[9:0];             // rounding carry (no ripple needed)

  // ------------------------------------------------ leading-one position
  wire [4:0] P;
  assign P[0] = h[0]|h[2]|h[4]|h[6]|h[8]|h[10]|h[12]|h[14]|h[16]|h[18]|h[20];
  assign P[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign P[2] = h[0]|h[1]|h[6]|h[7]|h[8]|h[9]|h[14]|h[15]|h[16]|h[17];
  assign P[3] = h[6]|h[7]|h[8]|h[9]|h[10]|h[11]|h[12]|h[13];
  assign P[4] = h[0]|h[1]|h[2]|h[3]|h[4]|h[5];

  // ------------------------------------------------ exponent prefix adder (7-bit KS)
  // Va = Am + P, Vb = Am + P + 1  (both computed in parallel)
  wire [6:0] Bp = {2'b00, P};
  wire [6:0] px = Am ^ Bp;
  wire [6:0] gx = Am & Bp;

  wire [6:0] G1, P1, G2, P2, G3, P3;
  assign G1[0] = gx[0];                    assign P1[0] = px[0];
  assign G1[1] = gx[1] | (px[1] & gx[0]);  assign P1[1] = px[1] & px[0];
  assign G1[2] = gx[2] | (px[2] & gx[1]);  assign P1[2] = px[2] & px[1];
  assign G1[3] = gx[3] | (px[3] & gx[2]);  assign P1[3] = px[3] & px[2];
  assign G1[4] = gx[4] | (px[4] & gx[3]);  assign P1[4] = px[4] & px[3];
  assign G1[5] = gx[5] | (px[5] & gx[4]);  assign P1[5] = px[5] & px[4];
  assign G1[6] = gx[6] | (px[6] & gx[5]);  assign P1[6] = px[6] & px[5];

  assign G2[0] = G1[0];                    assign P2[0] = P1[0];
  assign G2[1] = G1[1];                    assign P2[1] = P1[1];
  assign G2[2] = G1[2] | (P1[2] & G1[0]);  assign P2[2] = P1[2] & P1[0];
  assign G2[3] = G1[3] | (P1[3] & G1[1]);  assign P2[3] = P1[3] & P1[1];
  assign G2[4] = G1[4] | (P1[4] & G1[2]);  assign P2[4] = P1[4] & P1[2];
  assign G2[5] = G1[5] | (P1[5] & G1[3]);  assign P2[5] = P1[5] & P1[3];
  assign G2[6] = G1[6] | (P1[6] & G1[4]);  assign P2[6] = P1[6] & P1[4];

  assign G3[0] = G2[0];                    assign P3[0] = P2[0];
  assign G3[1] = G2[1];                    assign P3[1] = P2[1];
  assign G3[2] = G2[2];                    assign P3[2] = P2[2];
  assign G3[3] = G2[3];                    assign P3[3] = P2[3];
  assign G3[4] = G2[4] | (P2[4] & G2[0]);  assign P3[4] = P2[4] & P2[0];
  assign G3[5] = G2[5] | (P2[5] & G2[1]);  assign P3[5] = P2[5] & P2[1];
  assign G3[6] = G2[6] | (P2[6] & G2[2]);  assign P3[6] = P2[6] & P2[2];

  // carry into bit i (Va): G3[i-1]; into Vb: G3[i-1] | P3[i-1] (carry-in 1)
  wire [6:0] Va, Vb;
  assign Va[0] = px[0];                    assign Vb[0] = ~px[0];
  assign Va[1] = px[1] ^ G3[0];            assign Vb[1] = px[1] ^ (G3[0] | P3[0]);
  assign Va[2] = px[2] ^ G3[1];            assign Vb[2] = px[2] ^ (G3[1] | P3[1]);
  assign Va[3] = px[3] ^ G3[2];            assign Vb[3] = px[3] ^ (G3[2] | P3[2]);
  assign Va[4] = px[4] ^ G3[3];            assign Vb[4] = px[4] ^ (G3[3] | P3[3]);
  assign Va[5] = px[5] ^ G3[4];            assign Vb[5] = px[5] ^ (G3[4] | P3[4]);
  assign Va[6] = px[6] ^ G3[5];            assign Vb[6] = px[6] ^ (G3[5] | P3[5]);

  // ------------------------------------------------ exponent decode
  wire under_sel = Va[6] | ~|Va[6:0];            // T <= 35 -> subnormal path
  wire under34   = Va[6];                        // T <= 34 -> exp 0
  wire over_T    = ~Va[6] & (Va[5] | &Va[4:0]);  // T >= 66 -> inf
  wire eq65      = ~Va[6] & Va[5] & &Va[4:1] & ~Va[0];
  wire over      = over_T | (eq65 & c);          // overflow to inf
  wire [4:0] field = c ? Vb[4:0] : Va[4:0];      // E = T + c - 35

  // ------------------------------------------------ subnormal result path
  // value*2^24 = prod * 2^(S-26): right shift by 26-S (>=1 whenever used)
  wire [5:0]  sfull = 6'd26 - S;
  wire [4:0]  s5    = sfull[4:0];
  wire [21:0] shs   = prod >> s5;
  wire [9:0]  fsp   = shs[9:0];

  wire [31:0] pg32 = {9'd0, prod, 1'b0};         // pg32[m] = prod[m-1]
  wire [31:0] bt2  = {8'd0, bot, 2'b00};         // bt2[m]  = bot[m-2]
  wire guard_s  = pg32[s5];
  wire sticky_s = bt2[s5];
  wire r_s      = guard_s & (sticky_s | fsp[0]);

  wire [9:0] frac_sb = fsp + r_s;                // rounded subnormal fraction
  wire       carry_s = r_s & &fsp[9:0];          // round-up to smallest normal

  // ------------------------------------------------ pack
  wire is_nan = nan_in;
  wire is_inf = inf_in | over;
  wire is_sub = under_sel & ~is_nan & ~inf_in;

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_sub            ? {4'd0, carry_s} : field;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      is_inf ? 10'd0 :
                      is_sub ? frac_sb : frac_r;
  wire       sg_y  = is_nan ? 1'b0 : sy;

  assign y = {sg_y, exp_y, frac_y};
endmodule

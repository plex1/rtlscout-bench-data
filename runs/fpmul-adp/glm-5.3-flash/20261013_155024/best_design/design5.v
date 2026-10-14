// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign / 5 exp / 10 frac), subnormals supported.
//
// Notation: prod = ma*mb (exact 11x11), P = position of its leading one,
// S = ea_eff + eb_eff.  Product value = 1.frac * 2^(P+S-50), so the biased
// result exponent is E = P + S - 35 + c (c = rounding carry).
//
// Datapath:
//  * The 11x11 product is reduced with an explicit carry-save tree (5 CSA
//    levels) and the final carry-propagate adder is a Kogge-Stone prefix
//    adder -- much faster than the ripple CPA synthesis infers for "*".
//  * Rounding carry without a ripple: adding r in {0,1} to N[9:0] carries
//    out iff N[9:0] is all ones and r = 1, so c = r & &N[9:0]; the same for
//    the subnormal path (carry_s = r_s & &fsp[9:0]).
//  * Result exponent selects:  under (T<=35) -> subnormal path, T<=34 -> exp 0,
//    T=35 -> exp = c (= carry_s), T>=36 -> exp = field; overflow T>=66 or
//    (T=65 & c) -> inf.  The subnormal result is re-derived straight from
//    prod (shift 26-S) and rounded once.
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

  // ------------------------------------------------ 11x11 carry-save product tree
  wire [21:0] pp [0:10];
  genvar gi;
  generate
    for (gi = 0; gi < 11; gi = gi + 1) begin : g_pp
      assign pp[gi] = ({11{ma[gi]}} & {11'b0, mb}) << gi;
    end
  endgenerate

  wire [21:0] s11 = pp[0] ^ pp[1] ^ pp[2];
  wire [21:0] c11 = ((pp[0] & pp[1]) | (pp[0] & pp[2]) | (pp[1] & pp[2])) << 1;
  wire [21:0] s12 = pp[3] ^ pp[4] ^ pp[5];
  wire [21:0] c12 = ((pp[3] & pp[4]) | (pp[3] & pp[5]) | (pp[4] & pp[5])) << 1;
  wire [21:0] s13 = pp[6] ^ pp[7] ^ pp[8];
  wire [21:0] c13 = ((pp[6] & pp[7]) | (pp[6] & pp[8]) | (pp[7] & pp[8])) << 1;
  wire [21:0] s14 = pp[9] ^ pp[10];
  wire [21:0] c14 = (pp[9] & pp[10]) << 1;

  wire [21:0] s21 = s11 ^ s12 ^ s13;
  wire [21:0] c21 = ((s11 & s12) | (s11 & s13) | (s12 & s13)) << 1;
  wire [21:0] s22 = c11 ^ c12 ^ c13;
  wire [21:0] c22 = ((c11 & c12) | (c11 & c13) | (c12 & c13)) << 1;
  wire [21:0] s23 = s14 ^ c14;
  wire [21:0] c23 = (s14 & c14) << 1;

  wire [21:0] s31 = s21 ^ s22 ^ s23;
  wire [21:0] c31 = ((s21 & s22) | (s21 & s23) | (s22 & s23)) << 1;
  wire [21:0] s32 = c21 ^ c22 ^ c23;
  wire [21:0] c32 = ((c21 & c22) | (c21 & c23) | (c22 & c23)) << 1;

  wire [21:0] s41 = s31 ^ c31 ^ s32;
  wire [21:0] c41 = ((s31 & c31) | (s31 & s32) | (c31 & s32)) << 1;
  wire [21:0] U   = s41 ^ c41 ^ c32;
  wire [21:0] V   = ((s41 & c41) | (s41 & c32) | (c41 & c32)) << 1;

  // Kogge-Stone prefix adder: prod = U + V (< 2^22, no carry out)
  wire [21:0] p0 = U ^ V;
  wire [21:0] g0 = U & V;
  wire [21:0] p1 = p0 & {p0[20:0], 1'b0};
  wire [21:0] g1 = g0 | (p0 & {g0[20:0], 1'b0});
  wire [21:0] p2 = p1 & {p1[19:0], 2'b0};
  wire [21:0] g2 = g1 | (p1 & {g1[19:0], 2'b0});
  wire [21:0] p3 = p2 & {p2[15:0], 4'b0};
  wire [21:0] g3 = g2 | (p2 & {g2[15:0], 4'b0});
  wire [21:0] p4 = p3 & {p3[7:0], 8'b0};
  wire [21:0] g4 = g3 | (p3 & {g3[7:0], 8'b0});
  wire [21:0] g5 = g4 | (p4 & {g4[5:0], 16'b0});
  wire [21:0] prod = p0 ^ {g5[20:0], 1'b0};

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

  // ------------------------------------------------ leading-one position / exponent
  wire [4:0] P;
  assign P[0] = h[0]|h[2]|h[4]|h[6]|h[8]|h[10]|h[12]|h[14]|h[16]|h[18]|h[20];
  assign P[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign P[2] = h[0]|h[1]|h[6]|h[7]|h[8]|h[9]|h[14]|h[15]|h[16]|h[17];
  assign P[3] = h[6]|h[7]|h[8]|h[9]|h[10]|h[11]|h[12]|h[13];
  assign P[4] = h[0]|h[1]|h[2]|h[3]|h[4]|h[5];

  wire [6:0] Va = Am + {2'b00, P};               // T - 35 (7-bit two's complement)
  wire under_sel = Va[6] | ~|Va[6:0];            // T <= 35 -> subnormal path
  wire under34   = Va[6];                        // T <= 34 -> exp 0
  wire over_T    = ~Va[6] & (Va[5] | &Va[4:0]);  // T >= 66 -> inf
  wire eq65      = ~Va[6] & Va[5] & &Va[4:1] & ~Va[0];
  wire over      = over_T | (eq65 & c);          // overflow to inf
  wire [4:0] field = c ? (Va[4:0] + 5'd1) : Va[4:0];  // E = T + c - 35

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
                     is_sub            ? {4'd0, carry_s} :
                     under34           ? 5'd0  : field;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      is_inf ? 10'd0 :
                      is_sub ? frac_sb : frac_r;
  wire       sg_y  = is_nan ? 1'b0 : sy;

  assign y = {sg_y, exp_y, frac_y};
endmodule

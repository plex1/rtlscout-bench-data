// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign / 5 exp / 10 frac), subnormals supported.
//
// Notation: prod = ma*mb (exact 11x11), P = position of its leading one,
// S = ea_eff + eb_eff.  Product value = 1.frac * 2^(P+S-50), so the biased
// result exponent is E = P + S - 35 + c (c = rounding carry).
//
// Timing structure:
//  * Am = S - 35 is formed in parallel with the multiplier; Va = Am + P is a
//    single 7-bit adder feeding exponent decode (under/over/field).
//  * The subnormal result path is anchored at 35 - S (early) and rounded once;
//    it is selected for Va <= 0 (T = P+S <= 35), which covers T = 35 as well --
//    there the 2^-24 grid of the subnormal path coincides with the normal
//    grid, so the select does not depend on the late rounding carry c.
//  * Zero products (only possible with S <= 31) fall through the subnormal
//    path and produce signed zeros.

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
  wire [6:0]  Am = S - 7'd35;                    // Va = Am + P = T - 35

  // ------------------------------------------------ product + prefix nets
  wire [21:0] prod = ma * mb;

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

  wire [10:0] fr     = {1'b0, N[9:0]} + {10'd0, r};
  wire        c      = fr[10];                 // rounding carry
  wire [9:0]  frac_r = fr[9:0];

  // ------------------------------------------------ leading-one position P
  wire [4:0] P;
  assign P[0] = h[0]|h[2]|h[4]|h[6]|h[8]|h[10]|h[12]|h[14]|h[16]|h[18]|h[20];
  assign P[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign P[2] = h[0]|h[1]|h[6]|h[7]|h[8]|h[9]|h[14]|h[15]|h[16]|h[17];
  assign P[3] = h[6]|h[7]|h[8]|h[9]|h[10]|h[11]|h[12]|h[13];
  assign P[4] = h[0]|h[1]|h[2]|h[3]|h[4]|h[5];

  // ------------------------------------------------ exponent decode
  wire [6:0] Va = Am + {2'b00, P};               // T - 35 (7-bit two's compl.)
  wire under_sel = Va[6] | ~|Va[6:0];            // T <= 35 -> subnormal path
  wire over_T    = ~Va[6] & (Va[5] | &Va[4:0]);  // T >= 66 -> inf
  wire eq65      = ~Va[6] & Va[5] & &Va[4:1] & ~Va[0];
  wire over      = over_T | (eq65 & c);          // T = 65 with rounding carry
  wire [4:0] Vb    = Va[4:0] + 5'd1;
  wire [4:0] field = c ? Vb : Va[4:0];           // E = T + c - 35

  // ------------------------------------------------ subnormal result path
  // value*2^24 = prod * 2^(S-26): right shift by 26-S (>=0 whenever used)
  wire [5:0]  sfull = 6'd26 - S;
  wire [4:0]  s5    = sfull[4:0];
  wire [21:0] shs   = prod >> s5;
  wire [9:0]  fsp   = shs[9:0];

  wire [31:0] pg32 = {9'd0, prod, 1'b0};         // pg32[m] = prod[m-1]
  wire [31:0] bt2  = {8'd0, bot, 2'b00};         // bt2[m]  = bot[m-2]
  wire guard_s  = pg32[s5];
  wire sticky_s = bt2[s5];
  wire r_s      = guard_s & (sticky_s | fsp[0]);

  wire [10:0] frs     = {1'b0, fsp} + {10'd0, r_s};
  wire        carry_s = frs[10];                 // round-up to smallest normal
  wire [9:0]  frac_sb = frs[9:0];

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

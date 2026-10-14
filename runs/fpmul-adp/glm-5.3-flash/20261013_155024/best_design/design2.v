// fp_mul_e5f10 -- binary16 multiplier (1 sign / 5 exp / 10 frac), subnormals supported.
//
// Structure:
//   * significand product (11x11) with a leading-zero counter built from
//     top-prefix ORs; the count selects an 11-bit window of the product
//     (N[j] = prod[j+11-lz]) realised as ONE right shifter on a zero padded
//     vector, so left/right normalisation cases share the same hardware.
//   * rounding is done in the packed domain: round_up is added to the 10-bit
//     fraction and its carry-out flows straight into the exponent, which
//     removes the mantissa post-shift of the classic datapath.
//   * underflowing products are re-aligned straight from the raw product
//     (shift = 26 - exp_sum) and rounded once, giving correct subnormal
//     results including the round-up to the smallest normal.

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

  wire nan_a  = eao & ~faz;
  wire nan_b  = ebo & ~fbz;
  wire inf_a  = eao &  faz;
  wire inf_b  = ebo &  fbz;
  wire zero_a = eaz &  faz;
  wire zero_b = ebz &  fbz;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ------------------------------------------------ significands / exponent sum
  wire [10:0] ma = {~eaz, fa};
  wire [10:0] mb = {~ebz, fb};
  wire [5:0]  S  = {1'b0, (eaz ? 5'd1 : ea)} + {1'b0, (ebz ? 5'd1 : eb)};

  // ------------------------------------------------ product + leading-zero count
  wire [21:0] prod = ma * mb;
  wire        pn   = |prod;

  wire [21:1] nz;                       // nz[j] = |prod[21:j]
  wire [21:0] h;                        // h[k] : leading one is prod[21-k]
  genvar j, k;
  generate
    for (j = 1; j <= 21; j = j + 1) begin : g_nz
      assign nz[j] = |prod[21:j];
    end
    for (k = 1; k <= 21; k = k + 1) begin : g_h
      assign h[k] = prod[21-k] & ~nz[22-k];
    end
  endgenerate
  assign h[0] = prod[21];

  wire [4:0] lz;
  assign lz[0] = h[1]|h[3]|h[5]|h[7]|h[9]|h[11]|h[13]|h[15]|h[17]|h[19]|h[21];
  assign lz[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign lz[2] = h[4]|h[5]|h[6]|h[7]|h[12]|h[13]|h[14]|h[15]|h[20]|h[21];
  assign lz[3] = h[8]|h[9]|h[10]|h[11]|h[12]|h[13]|h[14]|h[15];
  assign lz[4] = h[16]|h[17]|h[18]|h[19]|h[20]|h[21];

  wire [31:0] bot;                      // bot[m] = |prod[m:0]  (zero padded)
  assign bot[0] = prod[0];
  generate
    for (k = 1; k <= 21; k = k + 1) begin : g_bot
      assign bot[k] = |prod[k:0];
    end
  endgenerate
  assign bot[31:22] = 10'b0;

  // ------------------------------------------------ normal path
  // window of the product: N[j] = prod[j+11-lz], out-of-range reads are 0
  wire [31:0] pe   = {prod, 10'd0};
  wire [4:0]  kamt = 5'd21 - lz;
  wire [31:0] shw  = pe >> kamt;
  wire [10:0] N    = shw[10:0];         // N[10] is the leading one

  wire pick_right = (lz <= 5'd11);
  wire guard_n  = pick_right & pe[5'd20 - lz];          // prod[10-lz], 0 if none
  wire sticky_n = pick_right & bot[5'd9 - lz];          // |prod[9-lz:0]
  wire round_n  = pick_right & guard_n & (sticky_n | N[0]); // round to nearest even

  wire [10:0] fr     = {1'b0, N[9:0]} + {10'd0, round_n};
  wire        c      = fr[10];          // rounding carried into the exponent
  wire [9:0]  frac_r = fr[9:0];

  // ------------------------------------------------ exponent
  // biased exponent of the rounded normal result (8-bit two's complement)
  wire [7:0] E8    = {2'b00, S} + {7'b0, c} - {3'b0, lz} - 8'd14;
  wire       under = pn & (E8[7] | ~|(E8[6:0]));        // subnormal (or zero) result
  wire       over  = pn & ~E8[7] & (E8[5] | (&E8[4:0])); // overflow to inf
  wire [4:0] field = E8[4:0];

  // ------------------------------------------------ subnormal path
  // value*2^24 = prod * 2^(S-26)  ->  right shift by 26-S (0..24 when under)
  wire [5:0]  sfull = 6'd26 - S;
  wire [4:0]  s5    = sfull[4:0];
  wire [21:0] shs   = prod >> s5;
  wire [9:0]  fsp   = shs[9:0];

  wire [31:0] pe3 = {9'd0, prod, 1'b0};                 // pe3[t] = prod[t-1]
  wire guard_s  = under & pe3[s5];                      // prod[s-1]
  wire sticky_s = under & bot[s5 - 5'd2];   // |prod[s-2:0]
  wire round_s  = pn & guard_s & (sticky_s | fsp[0]);

  wire [10:0] frs     = {1'b0, fsp} + {10'd0, round_s};
  wire        carry_s = frs[10];                        // rounded up to smallest normal
  wire [9:0]  frac_sb = frs[9:0];

  // ------------------------------------------------ pack
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | over);
  wire is_sub  = ~is_nan & ~is_inf & under;
  wire is_zero = ~is_nan & ~is_inf & (zero_in | ~pn | (is_sub & ~carry_s & ~|frac_sb));

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_zero           ? 5'd0  :
                     is_sub            ? {4'd0, carry_s} : field;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 :
                      is_sub ? frac_sb : frac_r;
  wire       sg_y  = is_nan ? 1'b0 : sy;

  assign y = {sg_y, exp_y, frac_y};
endmodule

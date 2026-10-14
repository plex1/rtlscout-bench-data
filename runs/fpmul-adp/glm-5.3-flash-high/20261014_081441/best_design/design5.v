// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
//
// Strategy: pre-normalise subnormal inputs with a small 10-bit CLZ + shifter so that both
// significands always have their leading one at bit 10.  The 11x11 product then has its
// leading one at bit 20 or 21, so:
//   - normal results need only a FIXED shift (10+c): pure slice selects, no shifter
//   - subnormal results use one 22-bit barrel shifter q = p >> t (t = 15 - esum);
//     fraction = q[20:11], guard = q[11], sticky = |q[10:0]| + bits shifted out
//   - sticky comes from a 4-level prefix-OR tree over the product
//   - rounding uses a parallel-prefix conditional increment (no ripple adder)
//   - zero inputs fold into the effective exponent (no extra output muxes)
//   - overflow is split into an early part (muxes) and a late part (exponent only)

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // 1. unpack
  wire       sa = a[15];
  wire       sb = b[15];
  wire [4:0] ea = a[14:10];
  wire [4:0] eb = b[14:10];
  wire [9:0] fa = a[9:0];
  wire [9:0] fb = b[9:0];

  // 2. classify
  wire ea_z = ~|ea, eb_z = ~|eb;
  wire ea_o =  &ea, eb_o =  &eb;
  wire fa_z = ~|fa, fb_z = ~|fb;

  wire za   = ea_z & fa_z;
  wire zb   = eb_z & fb_z;
  wire na   = ea_o & ~fa_z;
  wire nb   = eb_o & ~fb_z;
  wire ia   = ea_o &  fa_z;
  wire ib   = eb_o &  fb_z;
  wire suba = ea_z & ~fa_z;
  wire subb = eb_z & ~fb_z;

  wire nan_in  = na | nb | (ia & zb) | (za & ib);
  wire inf_in  = ia | ib;
  wire sy      = sa ^ sb;

  // 3. normalise subnormal inputs: leading-zero count of the 10-bit fraction
  function [3:0] clz10;
    input [9:0] f;
    begin
      clz10 = f[9] ? 4'd0 :
              f[8] ? 4'd1 :
              f[7] ? 4'd2 :
              f[6] ? 4'd3 :
              f[5] ? 4'd4 :
              f[4] ? 4'd5 :
              f[3] ? 4'd6 :
              f[2] ? 4'd7 :
              f[1] ? 4'd8 :
                     4'd9;
    end
  endfunction

  wire [3:0] lz_a = clz10(fa);
  wire [3:0] lz_b = clz10(fb);

  // 11-bit significands, leading one always at bit 10 (zero inputs -> 0)
  wire [10:0] ma = za   ? 11'd0 :
                   suba ? ({fa, 1'b0} << lz_a) : {1'b1, fa};
  wire [10:0] mb = zb   ? 11'd0 :
                   subb ? ({fb, 1'b0} << lz_b) : {1'b1, fb};

  // 4. significand product, leading one always at bit 20 or 21 (nonzero case)
  wire [21:0] p = ma * mb;
  wire        c = p[21];

  // effective exponents; subnormal: e' = -lz; zero input: very negative so the
  // subnormal path (p = 0) naturally yields a zero result with the right sign
  wire signed [6:0] ea_eff = za   ? -7'sd24 :
                             suba ? -$signed({3'b0, lz_a}) : $signed({2'b0, ea});
  wire signed [6:0] eb_eff = zb   ? -7'sd24 :
                             subb ? -$signed({3'b0, lz_b}) : $signed({2'b0, eb});

  wire signed [6:0] esum = ea_eff + eb_eff;            // [-48, 60]
  wire signed [6:0] E    = esum - 7'sd15 + {6'b0, c};  // [-63, 46]

  wire is_sub = ~nan_in & ~inf_in & (E <= 7'sd0);

  // subnormal shift amount: value = p * 2^(esum-50); fraction int = p >> (26 - esum)
  wire signed [6:0] tw = 7'sd15 - esum;
  wire [5:0] t    = tw[5:0];        // in [0, 63] whenever is_sub
  wire      tbig  = t[5] | t[4] | (t[3] & t[2]);   // t >= 12: everything out

  // 5. prefix-OR over the product: b4[k] = |p[k:0] (4 stages cover bits 0..15)
  wire [21:0] b1 = p  | (p  <<  1);
  wire [21:0] b2 = b1 | (b1 <<  2);
  wire [21:0] b3 = b2 | (b2 <<  4);
  wire [21:0] b4 = b3 | (b3 <<  8);

  // 6a. normal result: fixed right shift by 10+c, rounding bits are a fixed slice.
  //     Rounding field F = c ? p[20:11] : p[19:10] (10 bits, leading one excluded):
  //     F[9:1] = c ? p[20:12] : p[19:11], F[0] = c ? p[11] : p[10].
  wire [8:0]  F_hi   = c ? p[20:12] : p[19:11];
  wire        F_lo   = c ? p[11]   : p[10];
  wire        rnd_n  = c ? p[10]   : p[9];
  wire        stky_n = c ? b4[9]   : b4[8];
  wire        ru_n   = rnd_n & (stky_n | F_lo);

  // 6b. subnormal result: one 22-bit barrel shifter q = p >> t yields
  //     fraction = q[20:11], guard = q[10], sticky = |q[9:0]| + bits shifted out.
  wire [21:0] q0 = t[3] ? { 8'b0, p[21:8]} : p;
  wire [21:0] q1 = t[2] ? { 4'b0, q0[21:4]} : q0;
  wire [21:0] q2 = t[1] ? { 2'b0, q1[21:2]} : q1;
  wire [21:0] qs = t[0] ? { 1'b0, q2[21:1]} : q2;
  wire [21:0] q  = tbig ? 22'd0 : qs;   // t >= 12: everything shifted out

  wire [9:0]  F_sub   = q[20:11];
  wire        guard_s = q[10];
  wire        out3    = t[3] & (|p[7:0]);
  wire        out2    = t[2] & (|q0[3:0]);
  wire        out1    = t[1] & (|q1[1:0]);
  wire        out0    = t[0] & q2[0];
  wire        stky_s  = (|q[9:0]) | out3 | out2 | out1 | out0;

  wire        ru_s    = guard_s & (stky_s | F_sub[0]);

  // round-to-nearest-even conditional increment: sum = F + r (10-bit), plus carry.
  // prefix-AND: A[i] = &F[i:0];  sum[i] = F[i] ^ (r & A[i-1]), sum[0] = F[0] ^ r,
  // carry = r & A[9].  AND-gating (no mux) keeps the increment shallow.
  function [10:0] rnd_inc;   // {carry_out, sum[9:0]}
    input        r;
    input [9:0]  F;
    reg [9:0] A1, A2, A3, A4, G;
    begin
      A1 = F & {F[8:0], 1'b1};
      A2 = A1 & {A1[7:0], 2'b11};
      A3 = A2 & {A2[5:0], 4'b1111};
      A4 = A3 & {A3[1:0], 8'b11111111};
      G  = {10{r}} & A4;
      rnd_inc = {(r & A4[9]), F ^ {G[8:0], r}};
    end
  endfunction

  wire [10:0] res_n = rnd_inc(ru_n, {F_hi, F_lo});
  wire [9:0]  frac_n  = res_n[9:0];
  wire        carry_n = res_n[10];   // rounded past 0x7FF -> exponent +1

  wire [10:0] res_s = rnd_inc(ru_s, F_sub);
  wire [9:0]  frac_s  = res_s[9:0];
  wire        carry_s = res_s[10];   // rounded up to the smallest normal

  // 7. pack.  Overflow: E >= 31 (early), or E == 30 with a rounding carry -- in the
  //    latter case E[4:0] + carry_n = 31 already yields the right exponent and the
  //    fraction wraps to zero by itself, so no late muxing is needed.
  wire is_nan  = nan_in;
  wire inf_e   = inf_in | (E > 7'sd30);
  wire is_inf  = ~is_nan & (inf_e | ((&E[4:0]) & carry_n & (E[6:5] == 2'b00)));
  wire [4:0] E_out = E[4:0] + {4'd0, carry_n};

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_sub            ? {4'd0, carry_s} :
                                         E_out;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      inf_e  ? 10'd0   :
                      is_sub ? frac_s  : frac_n;

  assign y = {(is_nan ? 1'b0 : sy), exp_y, frac_y};
endmodule

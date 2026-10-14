// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
//
// Strategy: pre-normalise subnormal inputs with a small 10-bit CLZ + shifter so that both
// significands always have their leading one at bit 10.  The 11x11 product then has its
// leading one at bit 20 or 21, so:
//   - normal results need only a FIXED shift (10 or 11): pure slice selects, no shifter
//   - subnormal results shift the top 11 product bits by t = 15 - esum (small 11-bit barrel)
//   - sticky comes from a 5-level prefix-OR tree over the product

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
  wire zero_in = za | zb;
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

  // 4. significand product, always with leading one at bit 20 or 21 (nonzero case)
  wire [21:0] p = ma * mb;
  wire        c = p[21];

  // effective (normalised) exponents; subnormal: e' = -lz
  wire signed [6:0] ea_eff = suba ? $signed(-$signed({3'b0, lz_a}))
                                  : $signed({2'b0, ea});
  wire signed [6:0] eb_eff = subb ? $signed(-$signed({3'b0, lz_b}))
                                  : $signed({2'b0, eb});

  wire signed [6:0] esum = ea_eff + eb_eff;          // [-18, 60]
  wire signed [6:0] E    = esum - 7'sd15 + {6'b0, c}; // [-34, 46]

  wire is_sub = ~nan_in & ~inf_in & ~zero_in & (E <= 7'sd0);

  // subnormal shift amount: value = p * 2^(esum-50); fraction int = p >> (26 - esum)
  wire signed [6:0] tw = 7'sd15 - esum;
  wire [5:0] t     = tw[5:0];        // in [0, 33] whenever is_sub
  wire      tge12  = (t >= 6'd12);   // everything shifted out: result is zero

  // 5. prefix-OR over the product: b5[k] = |p[k:0]
  wire [21:0] b1 = p  | (p  <<  1);
  wire [21:0] b2 = b1 | (b1 <<  2);
  wire [21:0] b3 = b2 | (b2 <<  4);
  wire [21:0] b4 = b3 | (b3 <<  8);
  wire [21:0] b5 = b4 | (b4 << 16);

  // 6a. normal result: fixed right shift by 10+c, rounding bits are a fixed slice
  wire [10:0] F_pre  = c ? p[21:11] : p[20:10];
  wire        rnd_n  = c ? p[10]    : p[9];
  wire        stky_n = c ? b5[9]    : b5[8];
  wire        ru_n   = rnd_n & (stky_n | F_pre[0]);
  wire [11:0] F12n   = {1'b0, F_pre} + ru_n;
  wire        carry_n = F12n[11];
  wire [9:0]  frac_n  = carry_n ? 10'd0 : F12n[9:0];

  // 6b. subnormal result: one 22-bit barrel shifter q = p >> t yields
  //     fraction, guard and sticky at once:  F_sub = q[21:11], guard = q[10],
  //     sticky = |q[9:0].  For t >= 12 the guard is naturally 0 -> zero result.
  wire [21:0] q0 = t[3] ? { 8'b0, p[21:8]} : p;
  wire [21:0] q1 = t[2] ? { 4'b0, q0[21:4]} : q0;
  wire [21:0] q2 = t[1] ? { 2'b0, q1[21:2]} : q1;
  wire [21:0] qs = t[0] ? { 1'b0, q2[21:1]} : q2;
  wire [21:0] q  = tge12 ? 22'd0 : qs;   // t >= 12: everything shifted out

  wire [10:0] F_sub   = q[21:11];
  wire        guard_s = q[10];
  // sticky: bits remaining below the guard plus the bits shifted out at each stage
  wire        out3    = t[3] & (|p[7:0]);
  wire        out2    = t[2] & (|q0[3:0]);
  wire        out1    = t[1] & (|q1[1:0]);
  wire        out0    = t[0] & q2[0];
  wire        stky_s  = (|q[9:0]) | out3 | out2 | out1 | out0;

  wire        ru_s    = guard_s & (stky_s | F_sub[0]);
  wire [11:0] F12s    = {1'b0, F_sub} + ru_s;
  wire        carry_s = F12s[10];   // F_sub <= 1023, rounding up to 1024 = min normal
  wire [9:0]  frac_s  = carry_s ? 10'd0 : F12s[9:0];

  // 7. pack
  wire ovf      = (E > 7'sd30) | (E == 7'sd30 & carry_n);
  wire is_nan   = nan_in;
  wire is_inf   = ~is_nan & (inf_in | ovf);
  wire is_zero  = ~is_nan & ~is_inf & zero_in;

  wire [6:0] E_out  = E + {6'b0, carry_n};
  wire [4:0] exp_n  = E_out[4:0];

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_zero           ? 5'd0  :
                     is_sub            ? {4'd0, carry_s} :
                                         exp_n;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 :
                      is_sub ? frac_s : frac_n;

  assign y = {(is_nan ? 1'b0 : sy), exp_y, frac_y};
endmodule

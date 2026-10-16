// fp_mul_e5f10 -- binary16 FP multiplier.
// v3: subnormal inputs pre-normalized (10-bit LZC + small left shift), so the
// 11x11 significand product is always normalized. One shared right shifter
// (13-bit window V = prod[21:9]) handles both normal mantissa extraction and
// subnormal-result alignment; one shared rounding incrementer; 7-bit exponent
// arithmetic only.
//
// Derivation: value(prod) = prod * 2^(Ea+Eb-50); er = Ea+Eb-15+p21 is the
// biased exponent as-if-normal. With V = prod[21:9] and
//   s = subnormal-path ? (2+p21-er) : (1+p21)
// we have mant = (V>>s)[9:0], guard = V[s-1], sticky = |V[s-2:0] | |prod[8:0].
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack & classify ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire ea_zero = (ea == 5'd0),  eb_zero = (eb == 5'd0);
  wire ea_ones = &ea,           eb_ones = &eb;
  wire fa_zero = (fa == 10'd0), fb_zero = (fb == 10'd0);
  wire zero_a  = ea_zero & fa_zero,  zero_b = eb_zero & fb_zero;
  wire nan_a   = ea_ones & ~fa_zero, nan_b  = eb_ones & ~fb_zero;
  wire inf_a   = ea_ones & fa_zero,  inf_b  = eb_ones & fb_zero;

  wire any_nan  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire any_inf  = inf_a | inf_b;
  wire any_zero = zero_a | zero_b;
  wire sign_p   = sa ^ sb;

  // ---- leading-zero count of the 10-bit fractions ----
  reg [3:0] lza, lzb;
  integer i;
  always @* begin
    lza = 4'd9;
    for (i = 0; i < 10; i = i + 1) if (fa[i]) lza = 4'd9 - i[3:0];
  end
  always @* begin
    lzb = 4'd9;
    for (i = 0; i < 10; i = i + 1) if (fb[i]) lzb = 4'd9 - i[3:0];
  end

  // ---- normalized significands; effective exponents (7-bit signed) ----
  wire [10:0] fa_sh = {1'b0, fa} << (lza + 4'd1);
  wire [10:0] fb_sh = {1'b0, fb} << (lzb + 4'd1);
  wire [10:0] ma = ea_zero ? {1'b1, fa_sh[9:0]} : {1'b1, fa};
  wire [10:0] mb = eb_zero ? {1'b1, fb_sh[9:0]} : {1'b1, fb};

  wire signed [6:0] Ea = ea_zero ? -$signed({3'd0, lza}) : $signed({2'd0, ea});
  wire signed [6:0] Eb = eb_zero ? -$signed({3'd0, lzb}) : $signed({2'd0, eb});

  // ---- significand product: always normalized (leading one at bit 21/20) ----
  wire [21:0] prod = ma * mb;
  wire        p21  = prod[21];

  wire signed [6:0] er = Ea + Eb - 7'sd15 + {6'd0, p21};
  wire              sub_path = (er <= 7'sd0);

  // ---- unified right shift on the 13-bit window V = prod[21:9] ----
  wire signed [6:0] s_sub = 7'sd2 + {6'd0, p21} - er;      // >= 2 on subnormal path
  wire [6:0]        s     = sub_path ? s_sub : (7'd1 + {6'd0, p21});
  wire              big   = (s > 7'd14);
  wire [3:0]        sc    = big ? 4'd14 : s[3:0];

  wire [15:0] V16    = {3'd0, prod[21:9]};
  wire [15:0] V_sh   = V16 >> sc;
  wire [9:0]  mant   = V_sh[9:0];
  wire        guard = (sc == 4'd0) ? 1'b0 : V16[sc - 4'd1];

  // sticky: |V[sc-2:0] | |prod[8:0], prefix OR over V
  wire [12:0] belowV;
  assign belowV[0] = prod[9];
  genvar k;
  generate
    for (k = 1; k < 13; k = k + 1) assign belowV[k] = belowV[k-1] | prod[9+k];
  endgenerate
  wire below8 = |prod[8:0];
  wire sticky = below8 | ((sc >= 4'd2) ? belowV[(sc >= 4'd14) ? 4'd12 : sc - 4'd2] : 1'b0);

  // ---- rounding (RTNE), single shared incrementer ----
  wire        rup      = guard & (sticky | mant[0]);
  wire [10:0] frac_rnd = {1'b0, mant} + {10'd0, rup};
  wire        carry    = frac_rnd[10];
  wire [9:0]  frac_out = frac_rnd[9:0];

  // ---- result exponent ----
  wire [4:0] er5   = er[4:0];
  wire       ovf   = ~sub_path & ((er + {6'd0, carry}) >= 7'sd31);
  wire [4:0] e_y   = ovf ? 5'd31 : (sub_path ? {4'd0, carry} : (er5 + {4'd0, carry}));

  // signed zero: zero input, or subnormal-path result that rounded to zero
  wire is_inf   = any_inf | ovf;
  wire res_zero = ~any_nan & ~is_inf &
                  (any_zero | (sub_path & (frac_out == 10'd0) & ~carry));

  // ---- pack ----
  wire [4:0] exp_y  = any_nan ? 5'd31 : is_inf ? 5'd31 : res_zero ? 5'd0 : e_y;
  wire [9:0] frac_y = any_nan ? 10'h200 : (is_inf | res_zero) ? 10'd0 : frac_out;
  wire sign_y = any_nan ? 1'b0 : sign_p;

  assign y = {sign_y, exp_y, frac_y};
endmodule

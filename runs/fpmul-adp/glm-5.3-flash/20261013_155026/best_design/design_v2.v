// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// v2: per-operand normalization.  Every operand significand is normalized to have its
// leading one at bit 10 (subnormal fractions are left-justified, with the exponent
// decremented accordingly), so the 11x11 product always lies in [2^20, 2^22).
// One unified signed exponent U2 (= ea'+eb'+18+adj, with E' = U2-33) drives:
//   - normal results (U2 >= 34): fixed slices of the product, standard RNE rounding.
//   - subnormal results (U2 <= 33): one variable right shift, guard/lsb from the shifted
//     vector, sticky from a fixed low OR plus a variable-mask OR.
// U2 <= 22 makes the value strictly below half an ulp -> rounds to zero.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // 1. unpack
  wire        s_a = a[15];
  wire        s_b = b[15];
  wire [4:0]  e_a = a[14:10];
  wire [4:0]  e_b = b[14:10];
  wire [9:0]  f_a = a[9:0];
  wire [9:0]  f_b = b[9:0];

  // 2. classify
  wire e_a_zero = ~|e_a;
  wire e_b_zero = ~|e_b;
  wire e_a_ones =  &e_a;
  wire e_b_ones =  &e_b;
  wire f_a_zero = ~|f_a;
  wire f_b_zero = ~|f_b;
  wire zero_a   = e_a_zero &  f_a_zero;
  wire zero_b   = e_b_zero &  f_b_zero;
  wire nan_a    = e_a_ones & ~f_a_zero;
  wire nan_b    = e_b_ones & ~f_b_zero;
  wire inf_a    = e_a_ones &  f_a_zero;
  wire inf_b    = e_b_ones &  f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // 3. per-operand normalize: index of the leading one of the fraction
  function [3:0] l1idx;
    input [9:0] f;
    integer i;
    begin
      l1idx = 4'd0;
      for (i = 0; i <= 9; i = i + 1)
        if (f[i]) l1idx = i[3:0];
    end
  endfunction

  wire [3:0] i_a = l1idx(f_a);
  wire [3:0] i_b = l1idx(f_b);

  // normalized 11-bit significands (leading one at bit 10); zero inputs give 0
  wire [15:0] ma_ext = {6'd0, f_a} << (5'd10 - {1'b0, i_a});
  wire [15:0] mb_ext = {6'd0, f_b} << (5'd10 - {1'b0, i_b});
  wire [10:0] m_a    = ma_ext[10:0];
  wire [10:0] m_b    = mb_ext[10:0];

  // shifted biased exponents: EA = e' + 9   (e' = i-9 for subnormals, e for normals)
  wire [5:0] EA = e_a_zero ? {2'b00, i_a} : (6'd9 + e_a);
  wire [5:0] EB = e_b_zero ? {2'b00, i_b} : (6'd9 + e_b);

  // 4. significand product and one-step normalization
  wire [21:0] prod = m_a * m_b;
  wire        adj  = prod[21];
  wire [21:0] P    = adj ? prod : {prod[20:0], 1'b0};
  wire [10:0] P11  = P[21:11];

  // combined exponent: U2 = ea' + eb' + 18 + adj ; result biased exponent E' = U2 - 33
  wire [6:0] U2 = {1'b0, EA} + {1'b0, EB} + adj;

  // 5. normal path (U2 >= 34): frac = P[20:11], guard = P[10], sticky = |P[9:0]
  wire [6:0]  exp_pre = U2 - 7'd33;
  wire        n_round = P[10] & (|P[9:0] | P[11]);
  wire [10:0] n_sum   = {1'b0, P[20:11]} + n_round;
  wire [6:0]  n_exp   = exp_pre + n_sum[10];
  wire        n_inf   = (n_exp >= 7'd31);

  // 6. subnormal path (U2 <= 33); U2 <= 22 => value < half ulp => rounded zero
  wire        sub_is   = (U2 <= 7'd33);
  wire        sub_zero = (U2 <= 7'd22);
  wire [5:0]  qv       = 6'd33 - {1'b0, U2[5:0]};   // = t-11, guard/frac shift
  wire [10:0] W        = P11 >> qv;                 // W[0]=guard, W[10:1]=frac
  wire [10:0] smask    = 11'h7FF >> (6'd11 - qv);   // bits [qv-1:0]
  wire        s_sticky = |P[10:0] | (|(P11 & smask));
  wire        s_round  = W[0] & (s_sticky | W[1]);
  wire [10:0] s_sum    = {1'b0, W[10:1]} + s_round; // s_sum[10]: up to smallest normal

  // 7. select / pack
  wire        res_inf  = ~sub_is & n_inf;
  wire        res_zero = sub_is & sub_zero;
  wire [4:0]  exp_g    = sub_is ? (sub_zero ? 5'd0 : {4'd0, s_sum[10]})
                                : (n_inf ? 5'd31 : n_exp[4:0]);
  wire [9:0]  frac_g   = sub_is ? (sub_zero ? 10'd0 : s_sum[9:0]) : n_sum[9:0];

  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | res_inf);
  wire is_zero = ~is_nan & ~is_inf & (zero_in | res_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_g;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_g;

  assign y = is_nan ? 16'h7E00 : {s_y, exp_y, frac_y};
endmodule

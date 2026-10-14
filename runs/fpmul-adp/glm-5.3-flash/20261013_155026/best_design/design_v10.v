// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// v10: fused per-operand normalizer (one-hot leading-one + AND-OR left-justify),
// sticky accumulation through the subnormal shift stages, early exp/exp+1 muxed by the
// rounding carry.
// Each operand significand is normalized to have its leading one at bit 10 (normals:
// {1,f}; subnormals: fraction left-justified with the exponent decremented), so the
// 11x11 product always lies in [2^20, 2^22).  Unified exponent U2 = EA + EB + adj gives
// biased result exponent U2-33 and drives:
//   - normal results (U2 >= 34): fixed slices of P, standard RNE rounding.
//   - subnormal results (U2 <= 33): variable right shift of P11, guard/lsb from the
//     shifted vector, sticky from the fixed low OR plus the bits discarded by the
//     shift stages.  U2 <= 22 (clamped shift 11) => rounded zero.

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

  // 3a. one-hot of the leading one of an 11-bit vector (all zero if v == 0)
  function [10:0] lead1oh;
    input [10:0] v;
    reg [10:0] g, oh;
    integer l, k;
    begin
      g = v;
      for (l = 1; l <= 4; l = l + 1)
        for (k = 10; k >= 1; k = k - 1)
          if (k + (1 << (l-1)) <= 10)
            g[k] = g[k] | g[k + (1 << (l-1))];   // g[k] = |v[10:k]
      oh[10] = v[10];
      for (k = 9; k >= 0; k = k - 1)
        oh[k] = v[k] & ~g[k+1];
      lead1oh = oh;
    end
  endfunction

  // 3b. left-justify v's leading one onto bit 10 (identity if it is already there)
  function [10:0] ljust;
    input [10:0] v;
    input [10:0] oh;
    reg [10:0] m, sh;
    integer k;
    begin
      m  = 11'd0;
      sh = v;
      for (k = 10; k >= 0; k = k - 1) begin
        m  = m | ({11{oh[k]}} & sh);
        sh = {sh[9:0], 1'b0};
      end
      ljust = m;
    end
  endfunction

  // 3c. binary encode of the leading-one index from the one-hot
  function [3:0] enc11;
    input [10:0] oh;
    begin
      enc11 = {|oh[10:8], |oh[7:4], (|oh[7:6]) | (|oh[3:2]),
               oh[9] | oh[7] | oh[5] | oh[3] | oh[1]};
    end
  endfunction

  wire [10:0] v_a = {~e_a_zero, f_a};
  wire [10:0] v_b = {~e_b_zero, f_b};
  wire [10:0] oh_a = lead1oh(v_a);
  wire [10:0] oh_b = lead1oh(v_b);
  wire [10:0] m_a  = ljust(v_a, oh_a);
  wire [10:0] m_b  = ljust(v_b, oh_b);

  // shifted biased exponents: EA = e' + 9  (e' = i-9 for subnormals, e for normals)
  wire [5:0] EA = e_a_zero ? {2'b00, enc11(oh_a)} : (6'd9 + e_a);
  wire [5:0] EB = e_b_zero ? {2'b00, enc11(oh_b)} : (6'd9 + e_b);

  // 4. significand product and one-step normalization (leading one onto bit 21)
  wire [21:0] prod = m_a * m_b;
  wire        adj  = prod[21];
  wire [10:0] gP11  = adj ? prod[21:11] : prod[20:10];   // P[21:11]
  wire        gGd   = adj ? prod[10]    : prod[9];       // P[10]  (guard)
  wire        stLo  = adj ? |prod[9:0]  : |prod[8:0];    // |P[9:0]

  // 5. normal path (U2 >= 34): frac = gP11[9:0], guard = gGd, lsb = gP11[0]
  wire        n_round = gGd & (stLo | gP11[0]);
  wire [9:0]  frac_n  = gP11[9:0];
  wire [9:0]  fA      = andpre10(frac_n);
  wire [9:0]  frac_r;
  genvar g;
  generate
    assign frac_r[0] = frac_n[0] ^ n_round;
    for (g = 1; g <= 9; g = g + 1)
      assign frac_r[g] = frac_n[g] ^ (n_round & fA[g-1]);
  endgenerate
  wire c = n_round & fA[9];   // rounding carry into the exponent

  // 6. exponent: U2 = EA + EB + adj; biased result exponent = U2 - 33 (+ carry)
  wire [6:0] U2      = {1'b0, EA} + {1'b0, EB} + adj;
  wire       sub_is  = (U2 <= 7'd33);
  wire [6:0] exp_pre = U2 - 7'd33;                  // 1..48 on the normal path
  wire [6:0] exp_hi  = exp_pre + 7'd1;
  wire [5:0] n_exp   = c ? exp_hi[5:0] : exp_pre[5:0];
  wire       n_inf   = exp_pre[5] | (&exp_pre[4:0])
                      | (c & (exp_pre[4:0] == 5'd30) & ~exp_pre[5]);
  wire [4:0] exp_n   = n_inf ? 5'd31 : n_exp[4:0];

  // 7. subnormal path (U2 <= 33); clamped qv = 11 yields the rounded zero
  wire        qclamp = (U2 <= 7'd22);
  wire [5:0]  qv0    = 6'd33 - U2[5:0];
  wire [3:0]  qv     = qclamp ? 4'd11 : qv0[3:0];
  wire [10:0] w1 = qv[3] ? {8'b00000000, gP11[10:8]} : gP11;
  wire        q1 = qv[3] & |gP11[7:0];
  wire [10:0] w2 = qv[2] ? {4'b0000, w1[10:4]} : w1;
  wire        q2 = q1 | (qv[2] & |w1[3:0]);
  wire [10:0] w3 = qv[1] ? {2'b00, w2[10:2]} : w2;
  wire        q3 = q2 | (qv[1] & |w2[1:0]);
  wire [10:0] W  = qv[0] ? {1'b0, w3[10:1]} : w3;
  wire        q4 = q3 | (qv[0] & w3[0]);
  wire        stLo2 = |prod[9:0] | (adj & prod[10]);   // |P[10:0]
  wire        s_sticky = q4 | stLo2;
  wire        s_round  = W[0] & (s_sticky | W[1]);
  wire [9:0]  frac_s0  = W[10:1];
  wire [9:0]  sA       = andpre10(frac_s0);
  wire [9:0]  frac_sr;
  generate
    assign frac_sr[0] = frac_s0[0] ^ s_round;
    for (g = 1; g <= 9; g = g + 1)
      assign frac_sr[g] = frac_s0[g] ^ (s_round & sA[g-1]);
  endgenerate
  wire csub = s_round & sA[9];   // rounded up to the smallest normal

  // 8. select / pack
  wire [4:0] exp_g  = sub_is ? {4'd0, csub} : exp_n;
  wire [9:0] frac_g = sub_is ? frac_sr : frac_r;

  wire res_inf  = ~sub_is & n_inf;
  wire is_nan   = nan_in;
  wire is_inf   = ~is_nan & (inf_in | res_inf);
  wire is_zero  = ~is_nan & ~is_inf & zero_in;

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_g;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_g;

  assign y = is_nan ? 16'h7E00 : {s_y, exp_y, frac_y};

  // parallel-prefix AND (c[j] = &x[j:0], depth 4)
  function [9:0] andpre10;
    input [9:0] x;
    integer j, l;
    reg [9:0] c;
    begin
      c = x;
      for (l = 0; l < 4; l = l + 1)
        for (j = 9; j >= (1 << l); j = j - 1)
          c[j] = c[j] & c[j - (1 << l)];
      andpre10 = c;
    end
  endfunction
endmodule

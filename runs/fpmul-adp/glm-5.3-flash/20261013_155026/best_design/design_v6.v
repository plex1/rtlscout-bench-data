// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// v6: fast front end (one-hot leading-one decode + 4-stage barrel left shift), per-operand
// normalization, prefix-AND rounding increments, clamped subnormal shift.
// Every nonzero operand significand is normalized to have its leading one at bit 10
// (normals: {1,f}; subnormals: fraction left-justified so its leading one lands on bit 10,
// with the exponent decremented), so the 11x11 product always lies in [2^20, 2^22).
// Unified exponent U2 = EA + EB + adj gives biased result exponent U2-33 and drives:
//   - normal results (U2 >= 34): fixed slices of P, standard RNE rounding.
//   - subnormal results (U2 <= 33): variable right shift of P11 (clamped to 11, which
//     yields the rounded-to-zero result for U2 <= 22).
// U2 <= 22 makes the value strictly below half an ulp -> rounded zero.

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

  // 3a. leading-one one-hot of the fraction (all zero if f == 0)
  function [9:0] lead1;
    input [9:0] f;
    integer k;
    reg hi;
    reg [9:0] o;
    begin
      hi = 1'b1;
      for (k = 9; k >= 0; k = k - 1) begin
        o[k] = f[k] & hi;
        hi   = hi & ~f[k];
      end
      lead1 = o;
    end
  endfunction

  wire [9:0] oh_a = lead1(f_a);
  wire [9:0] oh_b = lead1(f_b);

  // 3b. binary encode of the leading-one index (parallel, from the one-hot)
  wire [3:0] i_a = {|oh_a[9:8], |oh_a[7:4], (|oh_a[7:6]) | (|oh_a[3:2]),
                    oh_a[9] | oh_a[7] | oh_a[5] | oh_a[3] | oh_a[1]};
  wire [3:0] i_b = {|oh_b[9:8], |oh_b[7:4], (|oh_b[7:6]) | (|oh_b[3:2]),
                    oh_b[9] | oh_b[7] | oh_b[5] | oh_b[3] | oh_b[1]};

  // 3c. left-justify the fraction so its leading one lands on bit 10 (subnormals only)
  //     m = f << (10 - i): leading one at bit 10, value preserved (no hidden bit!)
  wire [3:0]  d_a = 4'd10 - {1'b0, i_a};
  wire [3:0]  d_b = 4'd10 - {1'b0, i_b};
  wire [10:0] sa0 = {1'b0, f_a};
  wire [10:0] sa1 = d_a[0] ? {sa0[9:0], 1'b0} : sa0;
  wire [10:0] sa2 = d_a[1] ? {sa1[9:0], 2'b0} : sa1;
  wire [10:0] sa3 = d_a[2] ? {sa2[7:0], 4'b0} : sa2;
  wire [10:0] sa4 = d_a[3] ? {sa3[3:0], 8'b0} : sa3;
  wire [10:0] sb0 = {1'b0, f_b};
  wire [10:0] sb1 = d_b[0] ? {sb0[9:0], 1'b0} : sb0;
  wire [10:0] sb2 = d_b[1] ? {sb1[9:0], 2'b0} : sb1;
  wire [10:0] sb3 = d_b[2] ? {sb2[7:0], 4'b0} : sb2;
  wire [10:0] sb4 = d_b[3] ? {sb3[3:0], 8'b0} : sb3;

  // normalized 11-bit significands (leading one at bit 10); zero inputs give 0
  wire [10:0] m_a = e_a_zero ? sa4 : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? sb4 : {1'b1, f_b};

  // shifted biased exponents: EA = e' + 9  (e' = i-9 for subnormals, e for normals)
  wire [5:0] EA = e_a_zero ? {2'b00, i_a} : (6'd9 + e_a);
  wire [5:0] EB = e_b_zero ? {2'b00, i_b} : (6'd9 + e_b);

  // 4. significand product and one-step normalization (leading one onto bit 21)
  wire [21:0] prod = m_a * m_b;
  wire        adj  = prod[21];
  wire [21:0] P    = adj ? prod : {prod[20:0], 1'b0};
  wire [10:0] P11  = P[21:11];

  // combined exponent: biased result exponent = U2 - 33
  wire [6:0] U2 = {1'b0, EA} + {1'b0, EB} + adj;

  // 5. normal path (U2 >= 34): frac = P[20:11], guard = P[10], sticky = |P[9:0]
  wire        stickn  = |P[9:0];
  wire        n_round = P[10] & (stickn | P[11]);
  wire [10:0] nx      = {1'b0, P[20:11]};
  wire [10:0] nA      = andpre(nx);
  wire [11:0] n_sumx;
  genvar g;
  generate
    assign n_sumx[0] = nx[0] ^ n_round;
    for (g = 1; g <= 10; g = g + 1)
      assign n_sumx[g] = nx[g] ^ (n_round & nA[g-1]);
  endgenerate
  assign n_sumx[11] = n_round & nA[10];
  wire [6:0] exp_pre = U2 - 7'd33;
  wire [6:0] n_exp   = exp_pre + n_sumx[11];
  wire       n_inf   = n_exp[5] | (&n_exp[4:0]);

  // 6. subnormal path (23 <= U2 <= 33); clamped qv = 11 yields the rounded zero
  wire        sub_is = (U2 <= 7'd33);
  wire        qclamp = (U2 <= 7'd22);
  wire [5:0]  qv0    = 6'd33 - U2[5:0];
  wire [3:0]  qv     = qclamp ? 4'd11 : qv0[3:0];
  wire [10:0] w1 = qv[3] ? {8'b00000000, P11[10:8]} : P11;
  wire [10:0] w2 = qv[2] ? {4'b0000, w1[10:4]} : w1;
  wire [10:0] w3 = qv[1] ? {2'b00, w2[10:2]} : w2;
  wire [10:0] W  = qv[0] ? {1'b0, w3[10:1]} : w3;
  wire [10:0] smask;
  generate
    for (g = 0; g <= 10; g = g + 1)
      assign smask[g] = (qv > g[3:0]);
  endgenerate
  wire        mo       = |(P11 & smask);
  wire        s_sticky = |P[10:0] | mo;
  wire        s_round  = W[0] & (s_sticky | W[1]);
  wire [10:0] sx       = {1'b0, W[10:1]};
  wire [10:0] sA       = andpre(sx);
  wire [11:0] s_sumx;
  generate
    assign s_sumx[0] = sx[0] ^ s_round;
    for (g = 1; g <= 10; g = g + 1)
      assign s_sumx[g] = sx[g] ^ (s_round & sA[g-1]);
  endgenerate
  assign s_sumx[11] = s_round & sA[10];

  // 7. select / pack
  wire [4:0] exp_n  = n_inf ? 5'd31 : n_exp[4:0];
  wire [4:0] exp_g  = sub_is ? {4'd0, s_sumx[11]} : exp_n;
  wire [9:0] frac_g = sub_is ? s_sumx[9:0] : n_sumx[9:0];

  wire res_inf  = ~sub_is & n_inf;
  wire is_nan   = nan_in;
  wire is_inf   = ~is_nan & (inf_in | res_inf);
  wire is_zero  = ~is_nan & ~is_inf & zero_in;

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_g;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_g;

  assign y = is_nan ? 16'h7E00 : {s_y, exp_y, frac_y};

  // parallel-prefix AND (c[j] = &x[j:0], depth 4)
  function [10:0] andpre;
    input [10:0] x;
    integer j, l;
    reg [10:0] c;
    begin
      c = x;
      for (l = 0; l < 4; l = l + 1)
        for (j = 10; j >= (1 << l); j = j - 1)
          c[j] = c[j] & c[j - (1 << l)];
      andpre = c;
    end
  endfunction
endmodule

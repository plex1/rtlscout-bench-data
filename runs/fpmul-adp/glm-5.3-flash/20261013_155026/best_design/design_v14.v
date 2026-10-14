// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// v14: per-operand normalization + truncated multiplier + adj-free subnormal shift.
// Each nonzero operand significand m is normalized to have its leading one at bit 10
// (normals: {1,f}; subnormals: fraction left-justified with the exponent decremented),
// so m = 2^10 + t and the product prod = 2^20 + 2^10*(ta+tb) + ta*tb.
// Only prod[21:10] (prod_hi) and |prod[9:0] (sticky_low) are computed in full:
//   prod_hi = 2^10 + (ta+tb) + (ta*tb)[19:10],  sticky_low = |(ta*tb)[9:0]|.
// With P = prod<<(1-adj) (adj = prod[21], so P's leading one is at bit 21) and
// E' = U2-33 (U2 = EA+EB+adj):  value = P * 2^(E'-36).
//   - normal results (E' >= 1): frac = P[20:11], guard = P[10], sticky = |P[9:0].
//   - subnormal results (E' <= 0): k = value*2^24 = prod*2^(EA+EB-44)  [adj-free!],
//     so with qn = 33-(EA+EB): guard = prod[10+qn], frac[j] = prod[11+j+qn],
//     sticky = |prod[9+qn:0] (accumulated through the shift stages + sticky_low).
//     qn >= 12 => guard = 0 and frac = 0 => rounded zero automatically (clamped shift).

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

  // 4. truncated product: prod = 2^20 + 2^10*(ta+tb) + ta*tb
  wire [9:0]  ta = m_a[9:0];
  wire [9:0]  tb = m_b[9:0];
  wire [19:0] TT = ta * tb;
  wire [10:0] S1 = {1'b0, ta} + {1'b0, tb};
  wire [11:0] prod_hi = 12'd1024 + {1'b0, S1} + {2'b0, TT[19:10]};
  wire        adj        = prod_hi[11];
  wire        sticky_low = |TT[9:0];
  wire        sticky_lw0 = |TT[8:0];

  // 5. unified exponent
  wire [6:0] U2n     = {1'b0, EA} + {1'b0, EB};
  wire [6:0] U2      = U2n + adj;
  wire       sub_is  = (U2 <= 7'd33);

  // subnormal shift amount (adj-free, early): qn = 33 - (EA+EB), clamp >= 12 to 15
  wire [6:0] qn7 = 7'd33 - {1'b0, EA} - {1'b0, EB};
  wire       qbig = qn7[6] | qn7[5] | qn7[4] | (qn7[3] & qn7[2]);
  wire [3:0] qn_c = qbig ? 4'd15 : qn7[3:0];

  // 6. normal path (U2 >= 34): frac = P[20:11], guard = P[10], sticky = |P[9:0]
  wire [9:0] frac_n  = adj ? prod_hi[10:1] : prod_hi[9:0];
  wire       guard_n = adj ? prod_hi[0]   : TT[9];
  wire       stLo    = adj ? sticky_low   : sticky_lw0;
  wire       n_round = guard_n & (stLo | frac_n[0]);
  wire [9:0] fA      = andpre10(frac_n);
  wire [9:0] frac_r;
  genvar g;
  generate
    assign frac_r[0] = frac_n[0] ^ n_round;
    for (g = 1; g <= 9; g = g + 1)
      assign frac_r[g] = frac_n[g] ^ (n_round & fA[g-1]);
  endgenerate
  wire c = n_round & fA[9];   // rounding carry into the exponent

  // exponent of the normal result
  wire [6:0] exp_pre = U2 - 7'd33;
  wire [6:0] exp_hi  = exp_pre + 7'd1;
  wire [5:0] n_exp   = c ? exp_hi[5:0] : exp_pre[5:0];
  wire       n_inf   = exp_pre[5] | (&exp_pre[4:0])
                      | (c & (exp_pre[4:0] == 5'd30) & ~exp_pre[5]);
  wire [4:0] exp_n   = n_inf ? 5'd31 : n_exp[4:0];

  // 7. subnormal path: guard = prod[10+qn], frac[j] = prod[11+j+qn] via a 12-bit
  //    window V = prod[21:10] >> qn; sticky = discarded bits + sticky_low
  wire [15:0] V0 = {4'b0000, prod_hi};
  wire [15:0] V1 = qn_c[0] ? {1'b0, V0[15:1]} : V0;
  wire        q1 = qn_c[0] & V0[0];
  wire [15:0] V2 = qn_c[1] ? {2'b00, V1[15:2]} : V1;
  wire        q2 = q1 | (qn_c[1] & (|V1[1:0]));
  wire [15:0] V3 = qn_c[2] ? {4'b0000, V2[15:4]} : V2;
  wire        q3 = q2 | (qn_c[2] & (|V2[3:0]));
  wire [15:0] V4 = qn_c[3] ? {8'b00000000, V3[15:8]} : V3;
  wire        q4 = q3 | (qn_c[3] & (|V3[7:0]));
  wire [10:0] W  = V4[10:0];
  wire        s_sticky = q4 | sticky_low;
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

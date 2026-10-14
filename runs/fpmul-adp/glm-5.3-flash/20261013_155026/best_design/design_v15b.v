// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// v15b: per-operand normalization + truncated multiplier + adj-free subnormal shift +
// early (adj-free) exponent bookkeeping with a late 2-select increment mux.
// Each nonzero operand significand m is normalized to have its leading one at bit 10
// (normals: {1,f}; subnormals: fraction left-justified with the exponent decremented),
// so m = 2^10 + t and prod = m_a*m_b = 2^20 + 2^10*(ta+tb) + ta*tb, with
//   prod_hi = prod[21:10] = 2^10 + S1 + TT[19:10]   (S1 = ta+tb, TT = ta*tb),
//   adj = prod_hi[11] (carry: leading one of prod at bit 21 vs 20).
// P = prod << (1-adj) has its leading one at bit 21; result exponent E' = EA+EB+adj-33
// (+1 if rounding carries).  With S = EA+EB (early): E' = E0 + adj + c, E0 = S-33.
//   - normal results (E' >= 1): frac[j] = prod[10+j+adj], guard = prod[9+adj],
//     sticky = |prod[8+adj:0]; X = E'+c <= 29 so exp = X[4:0] and inf <=> X == 31.
//   - subnormal results (E' <= 0): qn = 33-S (early, adj-free): guard = prod[10+qn],
//     frac[j] = prod[11+j+qn], sticky = discarded bits + |prod[9:0].
//     qn >= 12 => guard = frac = 0 => rounded zero automatically (clamped shift).

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

  // 4. product: prod = 2^20 + 2^10*S1 + TT ; prod_hi = prod[21:10]
  wire [9:0]  ta = m_a[9:0];
  wire [9:0]  tb = m_b[9:0];
  wire [19:0] TT = ta * tb;
  wire [10:0] S1 = {1'b0, ta} + {1'b0, tb};
  wire [11:0] prod_hi = 12'd1024 + {1'b0, S1} + {2'b0, TT[19:10]};
  wire        adj        = prod_hi[11];
  wire        sticky_low = |TT[9:0];
  wire        sticky_lw0 = |TT[8:0];

  // 5. early (adj-free) exponent bookkeeping: S = EA+EB, E0 = S-33
  wire [6:0] S    = {1'b0, EA} + {1'b0, EB};
  wire [6:0] E0   = S - 7'd33;
  wire [6:0] E0p1 = E0 + 7'd1;
  wire [6:0] E0p2 = E0 + 7'd2;
  wire       ge29 = ~E0[6] & (E0[5] | (E0[4] & E0[3] & E0[2] & (E0[1] | E0[0])));
  wire       ge30 = ~E0[6] & (E0[5] | (E0[4] & E0[3] & E0[2] & E0[1]));
  wire       ge31 = ~E0[6] & (E0[5] | (&E0[4:0]));
  wire       le32 = ~S[5] | (~|S[4:0]);
  wire       le33 = ~S[5] | (~|S[4:1]);
  wire [6:0] qn7  = 7'd33 - S;
  wire       qbig = qn7[6] | qn7[5] | qn7[4] | (qn7[3] & qn7[2]);
  wire [3:0] qn_c = qbig ? 4'd15 : qn7[3:0];

  // 6. normal path: frac_n[j] = prod_hi[j+adj], guard = prod[9+adj], sticky below
  wire [9:0] frac_n  = adj ? prod_hi[10:1] : prod_hi[9:0];
  wire       guard_n = adj ? prod_hi[0]    : TT[9];
  wire       stLo    = adj ? sticky_low    : sticky_lw0;
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

  // 7. subnormal path: guard = prod[10+qn], frac[j] = prod[11+j+qn] via a 12-bit
  //    window V = prod_hi >> qn; sticky = discarded bits + sticky_low
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

  // 8. late exponent assembly: X = E0 + adj + c  (<= 29, so inf <=> X == 31)
  wire       sel1   = adj | c;
  wire       sel2   = adj & c;
  wire [4:0] X5     = sel2 ? E0p2[4:0] : sel1 ? E0p1[4:0] : E0[4:0];
  wire       n_inf  = sel2 ? ge29 : sel1 ? ge30 : ge31;
  wire       sub_is = adj ? le32 : le33;

  // 9. select / pack
  wire [4:0] exp_g  = sub_is ? {4'd0, csub} : X5;
  wire [9:0] frac_g = sub_is ? frac_sr : frac_r;

  wire res_inf  = ~sub_is & n_inf;
  wire is_nan   = nan_in;
  wire is_inf   = ~is_nan & (inf_in | res_inf);
  wire is_zero  = ~is_nan & ~is_inf & zero_in;

  wire [4:0] exp_y  = is_zero ? 5'd0 : (is_nan | inf_in) ? 5'd31 : exp_g;
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

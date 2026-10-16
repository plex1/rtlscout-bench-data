// fp_mul_e5f10 -- IEEE binary16 (1+5+10) multiplier, subnormals supported, RNE.
// v9: v6 structure (pre-normalized operands, de-serialized subnormal shifter)
// with flattened predicates, balanced incrementer, flat packing.

module lzc10 (input [9:0] f, output [3:0] c);
  wire [15:0] F  = {f, 6'b0};
  wire        c3 = ~|F[15:8];
  wire [7:0]  F2 = c3 ? F[7:0] : F[15:8];
  wire        c2 = ~|F2[7:4];
  wire [3:0]  F3 = c2 ? F2[3:0] : F2[7:4];
  wire        c1 = ~|F3[3:2];
  wire [1:0]  F4 = c1 ? F3[1:0] : F3[3:2];
  wire        c0 = ~F4[1];
  assign c = {c3, c2, c1, c0};
endmodule

// balanced carry-lookahead incrementer: S = A + Cin
module inc10 (input [9:0] A, input Cin, output [9:0] S, output Cout);
  wire [9:0] p = A;
  wire [9:0] G [0:4];
  wire [9:0] P [0:4];
  assign G[0] = 10'd0;
  assign P[0] = p;
  genvar i, l;
  generate
    for (l = 1; l <= 4; l = l + 1) begin : lvl
      for (i = 0; i < 10; i = i + 1) begin : b
        if (i >= (1 << (l-1))) begin : hi
          assign G[l][i] = G[l-1][i] | (P[l-1][i] & G[l-1][i - (1 << (l-1))]);
          assign P[l][i] = P[l-1][i] & P[l-1][i - (1 << (l-1))];
        end else begin : lo
          assign G[l][i] = G[l-1][i];
          assign P[l][i] = P[l-1][i];
        end
      end
    end
  endgenerate
  wire [9:0] carry;
  assign carry[0] = Cin;
  generate
    for (i = 1; i < 10; i = i + 1) begin : c
      assign carry[i] = G[4][i-1] | (P[4][i-1] & Cin);
    end
  endgenerate
  assign S = p ^ carry;
  assign Cout = G[4][9] | (P[4][9] & Cin);
endmodule

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire        s_a = a[15];
  wire        s_b = b[15];
  wire [4:0]  e_a = a[14:10];
  wire [4:0]  e_b = b[14:10];
  wire [9:0]  f_a = a[9:0];
  wire [9:0]  f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'd31);
  wire e_b_ones = (e_b == 5'd31);
  wire zero_a = e_a_zero & (f_a == 10'd0);
  wire zero_b = e_b_zero & (f_b == 10'd0);
  wire nan_a  = e_a_ones & (f_a != 10'd0);
  wire nan_b  = e_b_ones & (f_b != 10'd0);
  wire inf_a  = e_a_ones & (f_a == 10'd0);
  wire inf_b  = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire s_y     = s_a ^ s_b;

  wire [3:0] lzc_a, lzc_b;
  lzc10 u_lzc_a (.f(f_a), .c(lzc_a));
  lzc10 u_lzc_b (.f(f_b), .c(lzc_b));

  // normalized significands (hidden one at bit 10) and signed exponents.
  // Zero operands get a very negative exponent so the arithmetic path
  // naturally produces a signed zero (product is 0 anyway).
  wire [10:0] m_a = e_a_zero ? ({f_a, 1'b0} << lzc_a) : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? ({f_b, 1'b0} << lzc_b) : {1'b1, f_b};
  wire signed [6:0] E_a = e_a_zero ? (zero_a ? -7'sd20 : -$signed({3'd0, lzc_a}))
                                   : $signed({2'd0, e_a});
  wire signed [6:0] E_b = e_b_zero ? (zero_b ? -7'sd20 : -$signed({3'd0, lzc_b}))
                                   : $signed({2'd0, e_b});

  wire [21:0] prod = m_a * m_b;
  wire norm = prod[21];

  // E_sum = E_a + E_b ; W0 = E_a + E_b + norm
  wire signed [6:0] E_sum = E_a + E_b;
  wire signed [6:0] W0    = E_sum + {6'd0, norm};

  // below[j] = |prod[j:0], saturated above bit 21
  wire [24:0] below;
  genvar k;
  generate
    for (k = 0; k < 25; k = k + 1) begin : g_below
      assign below[k] = (k < 22) ? |prod[k:0] : |prod[21:0];
    end
  endgenerate

  // --- normal result path -------------------------------------------------
  wire [10:0] mant_pre = norm ? prod[21:11] : prod[20:10];
  wire        G_n      = norm ? prod[10] : prod[9];
  wire        S_n      = norm ? below[9] : below[8];
  wire        round_n  = G_n & (S_n | mant_pre[0]);
  wire [11:0] mant_sum = {1'b0, mant_pre} + {11'd0, round_n};
  wire        cy       = mant_sum[11];
  wire [9:0]  frac_n   = cy ? mant_sum[10:1] : mant_sum[9:0];

  // predicates from W0 = E_a + E_b + norm (pre-rounding biased exponent + 15)
  wire sub_res = W0[6] | ~|W0[6:4];                    // W0 <= 15
  wire w0_hi   = ~W0[6] & W0[5];
  wire gt45    = w0_hi & (W0[4] | (W0[3] & W0[2] & W0[1]));
  wire eq45    = w0_hi & ~W0[4] & W0[3] & W0[2] & ~W0[1] & W0[0];
  wire ovf     = gt45 | (eq45 & cy);
  wire [4:0] ey5 = W0[4:0] + (cy ? 5'd18 : 5'd17);     // W0 + cy - 15

  // --- subnormal result path ---------------------------------------------
  // value = prod * 2^(E_sum-50); subnormal frac = prod >> (26 - E_sum)
  // 26 - E_sum = 11 + x  with  x = 15 - E_sum  (for E_sum in [1,15]).
  // The shifter select uses only the low sum bits; the sign of E_sum forces
  // a zero result in parallel with the shift (E_sum <= 0 => value tiny).
  wire        zforce  = E_sum[6] | ~|E_sum[5:0];       // E_sum <= 0
  wire [3:0]  x       = ~E_sum[3:0];
  wire [15:0] w_ext   = {4'd0, prod[21:10]};           // prod >> 10
  wire [15:0] w_sh    = w_ext >> ({1'b0, x} + 5'd1);
  wire [9:0]  frac_t  = zforce ? 10'd0 : w_sh[9:0];
  wire        G_s     = zforce ? 1'b0 : w_ext[x];
  wire        S_s     = below[{1'b0, x} + 5'd9];
  wire        round_s = G_s & (S_s | frac_t[0]);
  wire [9:0]  frac_s;
  wire        fs_co;
  inc10 u_fs (.A(frac_t), .Cin(round_s), .S(frac_s), .Cout(fs_co));

  // --- pack (flat) ---------------------------------------------------------
  wire spec_e = nan_in | inf_in | ovf;    // exponent field 31
  wire spec_f = inf_in | ovf;             // fraction field 0 (non-NaN specials)

  wire [4:0] exp_y;
  wire [9:0] frac_y;
  assign exp_y[4:1] = {4{spec_e}} | ({4{~sub_res}} & ey5[4:1]);
  assign exp_y[0]   = spec_e | (sub_res & fs_co) | (~sub_res & ey5[0]);
  assign frac_y[8:0] = {9{~nan_in & ~spec_f}} &
                       (sub_res ? frac_s[8:0] : frac_n[8:0]);
  assign frac_y[9]  = nan_in | (~spec_f &
                       (sub_res ? frac_s[9] : frac_n[9]));
  wire sign_y = ~nan_in & s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

// fp_mul_e5f10 -- IEEE binary16 (1+5+10) multiplier, subnormals supported, RNE.
// v14: v13 with joint window/sticky extraction: shifting {w_ext, prefixOR, 0}
// places frac/G/S at fixed positions, removing the wide sticky mux.

module lzc10 (input [9:0] f, output [3:0] c);
  wire t1  = ~f[9];
  wire t2  = ~|f[9:8];
  wire t3  = ~|f[9:7];
  wire t4  = ~|f[9:6];
  wire t5  = ~|f[9:5];
  wire t6  = ~|f[9:4];
  wire t7  = ~|f[9:3];
  wire t8  = ~|f[9:2];
  wire t9  = ~|f[9:1];
  wire t10 = ~|f[9:0];
  assign c[3] = t8;
  assign c[2] = t4 & ~t8;
  assign c[1] = (t2 & ~t4) | (t6 & ~t8) | t10;
  assign c[0] = (t1 & ~t2) | (t3 & ~t4) | (t5 & ~t6) | (t7 & ~t8) | (t9 & ~t10);
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

  // --- normal result path -------------------------------------------------
  wire [10:0] mant_pre = norm ? prod[21:11] : prod[20:10];
  wire        G_n      = norm ? prod[10] : prod[9];
  wire        S_n      = norm ? (|prod[9:0]) : (|prod[8:0]);
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
  // frac = prod >> (26 - E_sum) = (prod >> 11) >> x,  x = 15 - E_sum.
  // Joint extraction: shift {w_ext, prefixOR(w_ext), 1'b0} by x so that frac,
  // G and S land at fixed positions (no wide sticky mux).
  wire [15:0] w_ext = {4'd0, prod[21:10]};             // prod >> 10
  wire        S0    = |prod[9:0];                      // bits below w_ext's base
  wire [15:0] bext;                                    // bext[j] = |prod[0..9+j]
  genvar j;
  generate
    for (j = 0; j < 16; j = j + 1) begin : gB
      if (j == 0) assign bext[j] = S0;
      else        assign bext[j] = S0 | |w_ext[j-1:0];
    end
  endgenerate
  wire [32:0] Pv = {w_ext, bext, 1'b0};
  wire [3:0]  x  = ~E_sum[3:0];
  wire [32:0] Q  = Pv >> x;                            // shift by x
  wire        zforce = E_sum[6] | ~|E_sum[5:0];        // E_sum <= 0
  wire [9:0]  frac_t = zforce ? 10'd0 : Q[27:18];
  wire        G_s    = zforce ? 1'b0 : Q[17];
  wire        S_s    = Q[1];
  wire        round_s = G_s & (S_s | frac_t[0]);
  wire [10:0] frac_sum = {1'b0, frac_t} + {10'd0, round_s};
  wire        fs_co    = frac_sum[10];

  // --- pack ---------------------------------------------------------------
  wire spec_e = nan_in | inf_in | ovf;    // exponent field 31
  wire spec_f = inf_in | ovf;             // fraction field 0 (non-NaN specials)

  wire [4:0] exp_y  = spec_e ? 5'd31 :
                      sub_res ? {4'd0, fs_co} : ey5;
  wire [9:0] frac_y = nan_in  ? 10'h200 :
                      spec_f  ? 10'd0 :
                      sub_res ? frac_sum[9:0] : frac_n;
  wire sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

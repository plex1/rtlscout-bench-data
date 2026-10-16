// fp_mul_e5f10 -- IEEE binary16 (1+5+10) multiplier, subnormals supported, RNE.
// Strategy: pre-normalize subnormal operands (LZC + 11-bit left shift),
// multiply normalized 11-bit significands (product in [2^20, 2^22)),
// single-bit normalisation, one rounding, compact subnormal right shift.

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
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [3:0] lzc_a, lzc_b;
  lzc10 u_lzc_a (.f(f_a), .c(lzc_a));
  lzc10 u_lzc_b (.f(f_b), .c(lzc_b));

  // normalized significands (hidden one at bit 10) and signed exponents
  wire [10:0] m_a = e_a_zero ? ({f_a, 1'b0} << lzc_a) : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? ({f_b, 1'b0} << lzc_b) : {1'b1, f_b};
  wire signed [7:0] E_a   = e_a_zero ? -$signed({4'd0, lzc_a}) : $signed({3'd0, e_a});
  wire signed [7:0] E_b   = e_b_zero ? -$signed({4'd0, lzc_b}) : $signed({3'd0, e_b});
  wire signed [7:0] E_sum = E_a + E_b;

  wire [21:0] prod = m_a * m_b;
  wire norm = prod[21];

  // below[j] = |prod[j:0], saturated above bit 21
  wire [24:0] below;
  genvar k;
  generate
    for (k = 0; k < 25; k = k + 1) begin : g_below
      assign below[k] = (k < 22) ? |prod[k:0] : |prod[21:0];
    end
  endgenerate

  // --- normal result path -------------------------------------------------
  wire [10:0] mant_pre  = norm ? prod[21:11] : prod[20:10];
  wire        G_n       = norm ? prod[10] : prod[9];
  wire        S_n       = norm ? below[9] : below[8];
  wire        round_n   = G_n & (S_n | mant_pre[0]);
  wire [11:0] mant_sum  = {1'b0, mant_pre} + {11'd0, round_n};
  wire        cy        = mant_sum[11];
  wire [9:0]  frac_n    = cy ? mant_sum[10:1] : mant_sum[9:0];

  wire signed [7:0] exp12    = E_sum - 8'sd15 + {7'd0, norm};
  wire signed [7:0] E_y      = exp12 + {7'd0, cy};
  wire              sub_res  = (exp12 <= 8'sd0);
  wire              overflow = (E_y > 8'sd30);

  // --- subnormal result path ---------------------------------------------
  // value = prod * 2^(E_sum-50); subnormal frac = prod >> (26 - E_sum)
  // 26 - E_sum = 11 + x  with  x = 15 - E_sum  (clamped at 15)
  wire [3:0]  x      = (E_sum <= 8'sd0) ? 4'd15 : (4'd15 - E_sum[3:0]);
  wire [15:0] w_ext  = {4'd0, prod[21:10]};        // prod >> 10
  wire [15:0] w_sh   = w_ext >> ({1'b0, x} + 5'd1);
  wire [9:0]  frac_t = w_sh[9:0];
  wire        G_s    = w_ext[x];
  wire        S_s    = below[{1'b0, x} + 5'd9];
  wire        round_s  = G_s & (S_s | frac_t[0]);
  wire [10:0] frac_sum = {1'b0, frac_t} + {10'd0, round_s};

  // --- pack ---------------------------------------------------------------
  wire [15:0] y_norm = {s_y, E_y[4:0], frac_n};
  wire [15:0] y_sub  = {s_y, 4'd0, frac_sum};

  assign y = nan_in    ? 16'h7E00 :
             inf_in    ? {s_y, 5'd31, 10'd0} :
             zero_in   ? {s_y, 15'd0} :
             overflow  ? {s_y, 5'd31, 10'd0} :
             sub_res   ? y_sub :
                         y_norm;
endmodule

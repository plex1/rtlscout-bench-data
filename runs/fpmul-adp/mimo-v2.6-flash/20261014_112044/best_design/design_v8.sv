// fp_mul_e5f10 v8: v6 with only the multiply operand order swapped
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // 1. unpack + classify
  wire        s_a = a[15];
  wire        s_b = b[15];
  wire [4:0]  e_a = a[14:10];
  wire [4:0]  e_b = b[14:10];
  wire [9:0]  f_a = a[9:0];
  wire [9:0]  f_b = b[9:0];

  wire e_a_zero = ~|e_a;
  wire e_b_zero = ~|e_b;
  wire e_a_ones =  &e_a;
  wire e_b_ones =  &e_b;
  wire zero_a = e_a_zero & ~|f_a;
  wire zero_b = e_b_zero & ~|f_b;
  wire nan_a  = e_a_ones &  |f_a;
  wire nan_b  = e_b_ones &  |f_b;
  wire inf_a  = e_a_ones & ~|f_a;
  wire inf_b  = e_b_ones & ~|f_b;

  wire nan_in = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire s_y    = s_a ^ s_b;

  // 2. effective significands / exponents (subnormals: no hidden one, eff exp 1)
  wire [10:0] m_a = {~e_a_zero, f_a};
  wire [10:0] m_b = {~e_b_zero, f_b};
  wire [5:0]  exp_sum = {1'b0, (e_a_zero ? 5'd1 : e_a)} + {1'b0, (e_b_zero ? 5'd1 : e_b)};

  // 3. significand product (22 bits) and balanced leading-zero count
  wire [21:0] prod = m_b * m_a;

  wire g_hi = |prod[21:16];
  wire g_md = |prod[15:8];

  wire h01 = |prod[21:20];
  wire h23 = |prod[19:18];
  wire [2:0] czh0 = prod[21] ? 3'd0 : 3'd1;
  wire [2:0] czh1 = prod[19] ? 3'd2 : 3'd3;
  wire [2:0] czh2 = prod[17] ? 3'd4 : 3'd5;
  wire [2:0] clz_hi = h01 ? czh0 : h23 ? czh1 : czh2;

  wire m01 = |prod[15:12];
  wire [3:0] czm0 = prod[15] ? 4'd0 : prod[14] ? 4'd1 : prod[13] ? 4'd2 : 4'd3;
  wire [3:0] czm1 = prod[11] ? 4'd4 : prod[10] ? 4'd5 : prod[9] ? 4'd6 : 4'd7;
  wire [3:0] clz_md = m01 ? czm0 : czm1;

  wire l01 = |prod[7:4];
  wire [3:0] czl0 = prod[7] ? 4'd0 : prod[6] ? 4'd1 : prod[5] ? 4'd2 : 4'd3;
  wire [3:0] czl1 = prod[3] ? 4'd4 : prod[2] ? 4'd5 : prod[1] ? 4'd6 : prod[0] ? 4'd7 : 4'd8;
  wire [3:0] clz_lo = l01 ? czl0 : czl1;

  wire [4:0] lz = g_hi ? {2'b00, clz_hi} :
                  g_md ? (5'd6 + {1'b0, clz_md}) :
                         (5'd14 - {1'b0, clz_lo});

  // 4. prefix reduction: below[k] = |prod[k:0] (shared sticky tree)
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_red
      assign below[k] = |prod[k:0];
    end
  endgenerate

  wire [63:0] prod_e  = {42'd0, prod};
  wire [31:0] below_x = {10'd0, below};

  // 5. regime selection (carry-free):  E_raw = d - 14 ; sub <=> E_raw <= 0
  wire signed [6:0] d    = $signed({1'b0, exp_sum}) - $signed({2'b0, lz});
  wire              sub  = (d <= 7'sd14);

  // 6. normal path: sr = 11 - lz (lz <= 11 whenever ~sub for valid inputs)
  wire [4:0]  sr         = 5'd11 - lz;
  wire [21:0] shifted_n  = prod >> sr;
  wire [10:0] mant_n     = shifted_n[10:0];
  wire [4:0]  idx_gn     = sr - 5'd1;
  wire [4:0]  idx_sn     = sr - 5'd2;
  wire guard_n  = prod_e[idx_gn];
  wire sticky_n = below_x[idx_sn];
  wire round_n  = guard_n & (sticky_n | mant_n[0]);
  wire [9:0] frac_n  = mant_n[9:0] + {9'd0, round_n};
  wire       carry_n = guard_n & (&mant_n[10:0]);

  // 7. subnormal path: ts = 26 - exp_sum
  wire [4:0]  ts         = 6'd26 - exp_sum;
  wire [21:0] shifted_s  = prod >> ts;
  wire [10:0] mant_s     = shifted_s[10:0];
  wire [4:0]  idx_gs     = ts - 5'd1;
  wire [4:0]  idx_ss     = ts - 5'd2;
  wire guard_s  = prod_e[idx_gs];
  wire sticky_s = below_x[idx_ss];
  wire round_s  = guard_s & (sticky_s | mant_s[0]);
  wire [9:0] frac_s  = mant_s[9:0] + {9'd0, round_s};
  wire       carry_s = guard_s & (&mant_s[9:0]);

  // 8. final exponent (normal path): E = d - 14 + carry_n
  wire signed [6:0] e_norm = d - 7'sd14 + {6'd0, carry_n};
  wire              ovf    = (e_norm > 7'sd30);

  // 9. pack
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_a | inf_b | ovf);
  wire is_sub   = ~nan_in & ~is_inf & sub;
  wire zero_in  = zero_a | zero_b;
  wire sub_zero = is_sub & ~carry_s & ~|frac_s;
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 :
                      is_sub ? {4'd0, carry_s} : e_norm[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 :
                      sub ? frac_s : frac_n;
  wire       sign_y = ~is_nan & s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

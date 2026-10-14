// fp_mul_e5f10 v4: v3 + split multiplier (parallel partial products)
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

  // 3. significand product (22 bits), split into parallel partial products
  wire [15:0] p_hi = m_a[10:6] * m_b;
  wire [16:0] p_lo = m_a[5:0]  * m_b;
  wire [21:0] prod = {6'd0, p_lo} + {p_hi, 6'd0};

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

  // 4. prefix reduction: below[k] = |prod[k:0] (sticky)
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_red
      assign below[k] = |prod[k:0];
    end
  endgenerate

  // 5. exponent bookkeeping (signed):
  //    E_raw = exp_sum + p - 35 = d - 14,  d = exp_sum - lz
  //    subnormal result  <=>  E_raw <= 0  <=>  d <= 14
  wire signed [6:0] d    = $signed({1'b0, exp_sum}) - $signed({2'b0, lz});
  wire              sub  = (d <= 7'sd14);

  // 6. shift amount: S = max(11 - lz, 26 - exp_sum, 0)
  wire signed [6:0] s_n  = 7'sd11 - $signed({2'b00, lz});
  wire signed [6:0] s_s  = 7'sd26 - $signed({1'b0, exp_sum});
  wire signed [6:0] s_m  = (s_n > s_s) ? s_n : s_s;
  wire [4:0]        S    = s_m[6] ? 5'd0 : s_m[4:0];

  // 7. shift + round (single rounding stage, ties-to-even)
  wire [21:0] shifted   = prod >> S;
  wire [10:0] mant_pre  = shifted[10:0];

  wire [63:0] prod_e = {42'd0, prod};
  wire [31:0] below_x = {10'd0, below};
  wire [4:0]  idx_g   = S - 5'd1;

  wire guard    = prod_e[idx_g];               // first bit shifted out
  wire sticky   = below_x[S - 5'd2];           // any bit below the guard
  wire round_up = guard & (sticky | mant_pre[0]);

  // 8. rounding carry-out without the incrementer:
  //    carry <=> guard & mantissa window all ones (guard is ANDed in anyway)
  wire carry  = guard & (&mant_pre[9:0]) & (sub | mant_pre[10]);

  wire [9:0]  frac_r     = mant_pre[9:0] + {9'd0, round_up};

  // 9. final exponent: normal E = d - 14 + carry ; subnormal field = carry
  wire signed [6:0] e_fin = d - 7'sd14 + {6'd0, carry};
  wire              ovf   = (e_fin > 7'sd30);   // only reachable on normal path
  wire [4:0]        exp_n = e_fin[4:0];         // valid when 1..30

  // 10. pack
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_a | inf_b | ovf);
  wire is_sub   = ~nan_in & ~is_inf & sub;
  wire zero_in  = zero_a | zero_b;
  wire sub_zero = is_sub & ~carry & ~|frac_r;
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 :
                      is_sub ? {4'd0, carry} : exp_n;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_r;
  wire       sign_y = ~is_nan & s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

// fp_mul_e5f10 -- optimized v2: barrel shifter for guard/sticky
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS         = 15;
  localparam MAX_E        = 31;
  localparam MAX_FINITE_E = 30;
  localparam PROD_W       = 22;

  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == MAX_E);
  wire e_b_ones = (e_b == MAX_E);
  wire zero_a   = e_a_zero & (f_a == 10'd0);
  wire zero_b   = e_b_zero & (f_b == 10'd0);
  wire nan_a    = e_a_ones & (f_a != 10'd0);
  wire nan_b    = e_b_ones & (f_b != 10'd0);
  wire inf_a    = e_a_ones & (f_a == 10'd0);
  wire inf_b    = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  wire [PROD_W-1:0] prod = m_a * m_b;

  // priority encoder for leading-zero count
  wire [4:0] lz =     prod[21] ? 5'd0  :
                      prod[20] ? 5'd1  :
                      prod[19] ? 5'd2  :
                      prod[18] ? 5'd3  :
                      prod[17] ? 5'd4  :
                      prod[16] ? 5'd5  :
                      prod[15] ? 5'd6  :
                      prod[14] ? 5'd7  :
                      prod[13] ? 5'd8  :
                      prod[12] ? 5'd9  :
                      prod[11] ? 5'd10 :
                      prod[10] ? 5'd11 :
                      prod[9]  ? 5'd12 :
                      prod[8]  ? 5'd13 :
                      prod[7]  ? 5'd14 :
                      prod[6]  ? 5'd15 :
                      prod[5]  ? 5'd16 :
                      prod[4]  ? 5'd17 :
                      prod[3]  ? 5'd18 :
                      prod[2]  ? 5'd19 :
                      prod[1]  ? 5'd20 :
                      prod[0]  ? 5'd21 : 5'd0;

  // Normal path: shift prod so leading 1 is at bit 10
  wire              need_right = (lz <= 5'd11);
  wire [5:0]        sr         = 6'd11 - {1'b0, lz};
  wire [5:0]        sl         = {1'b0, lz} - 6'd11;
  wire [PROD_W-1:0] shifted    = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0]       mant_pre   = shifted[10:0];

  // Guard/sticky for normal path using left-shift approach
  // prod << (PROD_W - sr) puts guard at MSB, sticky in lower bits
  // Only valid when need_right; when need_left, guard=sticky=0
  wire [5:0]        lshift_n   = 6'd22 - sr;
  wire [PROD_W-1:0] prod_l_n   = prod << lshift_n;
  wire              guard_n    = need_right ? prod_l_n[PROD_W-1] : 1'b0;
  wire              sticky_n   = need_right ? |prod_l_n[PROD_W-2:0] : 1'b0;
  wire              round_up  = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // Exponent of normal result
  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // Subnormal path
  wire signed [7:0] total_shift_s = $signed({1'b0, 8'd26}) - $signed({1'b0, exp_sum});
  wire [7:0]  total_shift = total_shift_s[7:0];
  wire        ts_valid   = (total_shift_s >= 0) & (total_shift_s < PROD_W);
  wire [PROD_W-1:0] sub_shifted = ts_valid ? (prod >> total_shift) : {PROD_W{1'b0}};
  wire [9:0]  frac_trunc = sub_shifted[9:0];

  // Guard/sticky for subnormal path using left-shift
  wire [7:0]        lshift_s  = 8'd22 - total_shift; // can wrap if total_shift > 22
  wire [PROD_W-1:0] prod_l_s  = (total_shift_s >= 1 && total_shift_s <= PROD_W) ? (prod << (PROD_W - total_shift)) : {PROD_W{1'b0}};
  wire              guard_s   = prod_l_s[PROD_W-1];
  wire              sticky_s  = |prod_l_s[PROD_W-2:0];
  wire              round_s   = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];
  wire [4:0]  exp_sub    = {4'd0, carry_s};

  // Pack
  wire is_nan      = nan_in;
  wire is_inf      = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_sub      = ~is_nan & ~is_inf & underflow;
  wire sub_is_zero = is_sub & (exp_sub == 5'd0) & (frac_sub == 10'd0);
  wire is_zero     = ~is_nan & ~is_inf & (zero_in | sub_is_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : is_sub ? exp_sub : exp_norm;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : is_sub ? frac_sub : frac_norm;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

// fp_mul_e5f10 -- optimized variant 4: avoid below[] array, use shift-based sticky
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

  wire e_a_zero = ~|e_a;
  wire e_b_zero = ~|e_b;
  wire e_a_ones = &e_a;
  wire e_b_ones = &e_b;
  wire zero_a   = e_a_zero & ~|f_a;
  wire zero_b   = e_b_zero & ~|f_b;
  wire nan_a    = e_a_ones & |f_a;
  wire nan_b    = e_b_ones & |f_b;
  wire inf_a    = e_a_ones & ~|f_a;
  wire inf_b    = e_b_ones & ~|f_b;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  wire [PROD_W-1:0] prod = m_a * m_b;

  // leading-zero count
  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < PROD_W; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // sticky: OR all bits below position idx
  // For normal path: sticky = |(prod >> (sr-1))  ... actually below[sr-2]
  // Instead, compute sticky by shifting prod right and checking if any bits remain

  // normalise and round
  wire              need_right = (lz <= 5'd11);
  wire [4:0]        sr         = 5'd11 - lz;
  wire [4:0]        sl         = lz - 5'd11;
  wire [PROD_W-1:0] shifted    = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0]       mant_pre   = shifted[10:0];

  // guard bit: prod[sr-1] when need_right
  // sticky bit: |prod[sr-2:0] when need_right
  // Compute sticky by OR-ing bits below guard position
  wire [PROD_W-1:0] guard_mask = (PROD_W'(1) << (sr - 1));
  wire [PROD_W-1:0] sticky_mask = guard_mask - 1'b1;  // bits below guard
  wire guard_n  = need_right & |(prod & guard_mask);
  wire sticky_n = need_right & |(prod & sticky_mask);
  wire round_up = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // subnormal result
  integer           total_shift;
  reg [PROD_W-1:0]  sub_shifted;
  always @* begin
    total_shift = (BIAS + 10 + 1) - $signed({1'b0, exp_sum});
    sub_shifted = (total_shift >= 0 && total_shift < PROD_W) ? (prod >> total_shift) : {PROD_W{1'b0}};
  end
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  // guard/sticky for subnormal
  wire [PROD_W-1:0] guard_mask_s = (PROD_W'(1) << (total_shift - 1));
  wire [PROD_W-1:0] sticky_mask_s = guard_mask_s - 1'b1;
  wire guard_s    = |(prod & guard_mask_s);
  wire sticky_s   = |(prod & sticky_mask_s);
  wire round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];
  wire [4:0]  exp_sub    = {4'd0, carry_s};

  // pack
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

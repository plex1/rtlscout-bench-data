// fp_mul_e5f10 -- optimized IEEE binary16 multiplier (barrel-shift guard/sticky)
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a, b;
  output [15:0] y;

  localparam BIAS = 15;

  // unpack
  wire       s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0],   f_b = b[9:0];

  // classify
  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'd31);
  wire e_b_ones = (e_b == 5'd31);
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

  // effective significands and exponents
  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [5:0]  ea_eff = e_a_zero ? 6'd1 : {1'b0, e_a};
  wire [5:0]  eb_eff = e_b_zero ? 6'd1 : {1'b0, e_b};

  // significand product
  wire [21:0] prod = m_a * m_b;

  // leading-zero count (priority encoder via for loop)
  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // suffix OR: below[k] = |prod[k:0]
  wire [21:0] below;
  assign below[0] = prod[0];
  genvar k;
  generate
    for (k = 1; k < 22; k = k + 1) begin : g_below
      assign below[k] = below[k-1] | prod[k];
    end
  endgenerate

  // ---- normal path ----
  wire        need_right = (lz <= 5'd11);
  wire [4:0]  sr = 5'd11 - lz;
  wire [4:0]  sl = lz - 5'd11;
  wire [21:0] shifted_n = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0] mant_pre  = shifted_n[10:0];

  // guard/sticky via barrel shift (sr-1 wraps to large value when sr=0, giving 0)
  wire [5:0] gs_shift_n = {1'b0, sr} - 6'd1;
  wire [5:0] st_shift_n = {1'b0, sr} - 6'd2;
  wire [21:0] prod_gs_n = prod >> gs_shift_n[4:0];
  wire [21:0] below_st_n = below >> st_shift_n[4:0];
  wire guard_n  = need_right ? prod_gs_n[0] : 1'b0;
  wire sticky_n = need_right ? below_st_n[0] : 1'b0;
  wire round_n  = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_n};
  wire        carry_n    = mant_round[11];
  wire [9:0]  frac_norm  = carry_n ? mant_round[10:1] : mant_round[9:0];

  // exponent (7-bit)
  wire [6:0] exp_sum   = {1'b0, ea_eff} + {1'b0, eb_eff};
  wire [6:0] lhs       = exp_sum + 7'd1 + {6'd0, carry_n};
  wire [6:0] limit_und = BIAS + lz;
  wire [6:0] limit_ovr = BIAS + 7'd30 + lz;
  wire       underflow = (lhs <= limit_und);
  wire       overflow  = (lhs >  limit_ovr);
  wire [6:0] e_norm    = lhs - limit_und;
  wire [4:0] exp_norm  = e_norm[4:0];

  // ---- subnormal path ----
  wire [6:0] total_shift = 7'd26 - exp_sum;
  wire       ts_valid    = (total_shift < 7'd22);
  wire [21:0] sub_shifted = ts_valid ? (prod >> total_shift[4:0]) : 22'd0;
  wire [9:0]  frac_trunc = sub_shifted[9:0];

  wire [6:0] gs_shift_s = total_shift - 7'd1;
  wire [6:0] st_shift_s = total_shift - 7'd2;
  wire [21:0] prod_gs_s = prod >> gs_shift_s[4:0];
  wire [21:0] below_st_s = below >> st_shift_s[4:0];
  wire        guard_s   = prod_gs_s[0];
  wire        sticky_s  = below_st_s[0];
  wire        round_s   = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum  = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s   = frac_sum[10];
  wire [9:0]  frac_sub  = frac_sum[9:0];
  wire [4:0]  exp_sub   = {4'd0, carry_s};

  // ---- pack ----
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

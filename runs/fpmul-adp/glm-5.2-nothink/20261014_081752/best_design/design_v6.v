// fp_mul_e5f10 -- optimized IEEE binary16 multiplier v6
// Normal path: left-shift prod by lz, fixed-position extraction.
// Subnormal path: explicit guard/sticky without bit_at function.

module fp_mul_e5f10 (input  [15:0] a, input [15:0] b, output [15:0] y);

  localparam BIAS = 15;
  localparam MAX_E = 31;
  localparam MAX_FINITE_E = 30;
  localparam PROD_W = 22;

  // 1. unpack
  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  // 2. classify
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

  // 3. effective significands and exponents
  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  // 4. significand product
  wire [PROD_W-1:0] prod = m_a * m_b;

  // leading-zero count
  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < PROD_W; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // below[k] = |prod[k:0] (prefix OR) -- for subnormal sticky
  wire [PROD_W-1:0] below;
  assign below[0] = prod[0];
  generate
    for (genvar k = 1; k < PROD_W; k = k + 1)
      assign below[k] = below[k-1] | prod[k];
  endgenerate

  // 5. NORMAL PATH: left-shift prod by lz to align MSB at bit 21
  wire [PROD_W-1:0] prod_aligned = prod << lz;
  wire [10:0]       mant_pre   = prod_aligned[21:11];
  wire              guard      = prod_aligned[10];
  wire              sticky     = |prod_aligned[9:0];
  wire              round_up  = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // 6. exponent of a normal result
  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // 7. SUBNORMAL PATH
  // total_shift = 26 - exp_sum (right shift of prod to subnormal grid)
  // Use 6-bit to handle range -4..26
  wire [6:0] ts7 = 7'd26 - {1'b0, exp_sum};  // 26 - exp_sum (7-bit, can be negative via wrap)
  wire       ts_neg = exp_sum > 7'd26;  // total_shift would be negative
  wire [4:0] total_shift = ts7[4:0];
  wire       ts_ge0      = ~ts_neg;
  wire       ts_lt22     = ts_ge0 & (total_shift < 5'd22);

  // Subnormal fraction: prod >> total_shift
  wire [PROD_W-1:0] sub_shifted = ts_lt22 ? (prod >> total_shift) : {PROD_W{1'b0}};
  wire [9:0]  frac_trunc = sub_shifted[9:0];

  // Guard and sticky using explicit indexing with range checks
  // guard = prod[total_shift - 1] if total_shift >= 1 and total_shift <= 21
  // sticky = below[total_shift - 2] if total_shift >= 2 and total_shift <= 21
  // For total_shift = 0: no guard (shift by 0, nothing lost)
  // For total_shift = 1: guard = prod[0], sticky = 0
  // For total_shift >= 22: guard = prod[21] if total_shift = 22, else 0

  // Use a 32-bit extended prod/below to avoid x from out-of-bounds indexing
  wire [31:0] prod_ext  = prod;
  wire [31:0] below_ext = below;
  wire [4:0]  ts_m1     = total_shift - 5'd1;
  wire [4:0]  ts_m2     = total_shift - 5'd2;

  wire guard_s  = ts_ge0 & (total_shift >= 5'd1) & prod_ext[ts_m1];
  wire sticky_s = ts_ge0 & (total_shift >= 5'd2) & below_ext[ts_m2];

  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];
  wire [4:0]  exp_sub    = {4'd0, carry_s};

  // 8. pack
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

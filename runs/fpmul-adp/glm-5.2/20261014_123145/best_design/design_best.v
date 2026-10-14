// fp_mul_e5f10 -- optimized: simplified output mux + parallel exp
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

  // leading-zero count
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

  // helper: extract bit at variable index (0 outside range)
  function bit_at;
    input [21:0] v;
    input integer idx;
    begin
      bit_at = (idx >= 0 && idx < 22) ? v[idx] : 1'b0;
    end
  endfunction

  // ---- normal path ----
  wire        need_right = (lz <= 5'd11);
  wire [4:0]  sr = 5'd11 - lz;
  wire [4:0]  sl = lz - 5'd11;
  wire [21:0] shifted_n = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0] mant_pre  = shifted_n[10:0];

  wire guard_n  = need_right ? bit_at(prod,  sr - 1) : 1'b0;
  wire sticky_n = need_right ? bit_at(below, sr - 2) : 1'b0;
  wire round_n  = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_n};
  wire        carry_n    = mant_round[11];
  wire [9:0]  frac_norm  = carry_n ? mant_round[10:1] : mant_round[9:0];

  // exponent: pre-compute for carry=0 and carry=1
  wire [6:0] exp_sum    = {1'b0, ea_eff} + {1'b0, eb_eff};
  wire [6:0] lhs_0      = exp_sum + 7'd1;
  wire [6:0] lhs_1      = exp_sum + 7'd2;
  wire [6:0] limit_und  = BIAS + lz;
  wire [6:0] limit_ovr  = BIAS + 7'd30 + lz;

  wire underflow_0 = (lhs_0 <= limit_und);
  wire underflow_1 = (lhs_1 <= limit_und);
  wire overflow_0  = (lhs_0 >  limit_ovr);
  wire overflow_1  = (lhs_1 >  limit_ovr);

  wire underflow = carry_n ? underflow_1 : underflow_0;
  wire overflow  = carry_n ? overflow_1  : overflow_0;

  wire [6:0] e_norm_0 = lhs_0 - limit_und;
  wire [6:0] e_norm_1 = lhs_1 - limit_und;
  wire [6:0] e_norm   = carry_n ? e_norm_1 : e_norm_0;
  wire [4:0] exp_norm = e_norm[4:0];

  // ---- subnormal path ----
  wire [6:0] total_shift = 7'd26 - exp_sum;
  wire       ts_valid    = (total_shift < 7'd22);
  wire [21:0] sub_shifted = ts_valid ? (prod >> total_shift) : 22'd0;
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s   = bit_at(prod,  total_shift - 1);
  wire        sticky_s  = bit_at(below, total_shift - 2);
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

  wire is_special = is_nan | is_inf | is_zero | is_sub;
  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_zero ? 5'd0 :
                     is_sub ? exp_sub : exp_norm;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 :
                      is_sub ? frac_sub : frac_norm;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

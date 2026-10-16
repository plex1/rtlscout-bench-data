// fp_mul_e5f10 -- starting-point datapath, but guard/sticky are produced by
// a logarithmic right-shifter with sticky accumulation rather than per-bit
// prefix-OR/mux trees.
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

  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < PROD_W; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // Normal-path logarithmic right shifter with sticky accumulation for the
  // normal result (guard is bit 10-lz of the product).
  wire [7:0]  exp_sum      = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire need_right = (lz <= 5'd11);
  wire [4:0] sl   = lz - 5'd11;
  wire [4:0] sh_amt = 5'd10 - lz;

  wire [PROD_W-1:0] d0 = sh_amt[0] ? (prod >> 5'd1)  : prod;
  wire              s0 = sh_amt[0] ? prod[0] : 1'b0;
  wire [PROD_W-1:0] d1 = sh_amt[1] ? (d0 >> 5'd2)    : d0;
  wire              s1 = s0 | (sh_amt[1] ? |d0[1:0] : 1'b0);
  wire [PROD_W-1:0] d2 = sh_amt[2] ? (d1 >> 5'd4)    : d1;
  wire              s2 = s1 | (sh_amt[2] ? |d1[3:0] : 1'b0);
  wire [PROD_W-1:0] d3 = sh_amt[3] ? (d2 >> 5'd8)    : d2;
  wire              s3 = s2 | (sh_amt[3] ? |d2[7:0] : 1'b0);
  wire [PROD_W-1:0] d4 = sh_amt[4] ? (d3 >> 5'd16)   : d3;
  wire              s4 = s3 | (sh_amt[4] ? |d3[15:0] : 1'b0);

  wire [PROD_W-1:0] rshifted  = d4;
  wire              rstick    = s4;

  // Normal result significand before round (split left/right paths)
  wire [PROD_W-1:0] lshifted  = prod << sl;
  wire              lz_eq_11  = (lz == 5'd11);
  wire [10:0]       mant_pre  = need_right ? (lz_eq_11 ? prod[10:0] : rshifted[11:1]) : lshifted[10:0];
  wire              guard     = need_right ? (lz_eq_11 ? 1'b0 : rshifted[0]) : 1'b0;
  wire              sticky    = need_right ? (lz_eq_11 ? 1'b0 : rstick) : 1'b0;
  wire              round_up  = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // Subnormal path: shift by q-1, where q = 26-exp_sum (1..22).  The first
  // discarded (guard) bit is then qfinal[0]; the carried fraction is qfinal[10:1].
  wire [7:0]  total_shift8 = 8'd26 - {3'd0, exp_sum};
  wire [4:0]  q_m1         = total_shift8[4:0] - 5'd1;
  wire        q_gt0        = (total_shift8 != 8'd0);

  wire [PROD_W-1:0] q0 = q_m1[0] ? (prod >> 5'd1)  : prod;
  wire              u0 = q_m1[0] ? prod[0] : 1'b0;
  wire [PROD_W-1:0] q1 = q_m1[1] ? (q0 >> 5'd2)    : q0;
  wire              u1 = u0 | (q_m1[1] ? |q0[1:0] : 1'b0);
  wire [PROD_W-1:0] q2 = q_m1[2] ? (q1 >> 5'd4)    : q1;
  wire              u2 = u1 | (q_m1[2] ? |q1[3:0] : 1'b0);
  wire [PROD_W-1:0] q3 = q_m1[3] ? (q2 >> 5'd8)    : q2;
  wire              u3 = u2 | (q_m1[3] ? |q2[7:0] : 1'b0);
  wire [PROD_W-1:0] q4 = q_m1[4] ? (q3 >> 5'd16)   : q3;
  wire              u4 = u3 | (q_m1[4] ? |q3[15:0] : 1'b0);

  wire [9:0] frac_trunc_s = q4[10:1];
  wire       guard_s    = q4[0];
  wire       sticky_s   = u4;
  wire       round_s    = guard_s & (sticky_s | frac_trunc_s[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc_s} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];
  wire [4:0]  exp_sub    = {4'd0, carry_s};

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
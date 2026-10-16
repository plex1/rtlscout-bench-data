// fp_mul_e5f10 -- v11: v8 with only the normal shifter truncated to its
// 12 useful output bits.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS         = 15;
  localparam MAX_E        = 31;
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

  wire [4:0]  lz;
  wire [11:0] top = prod[21:10];
  assign lz = top[11] ? 5'd0  :
              top[10] ? 5'd1  :
              top[9]  ? 5'd2  :
              top[8]  ? 5'd3  :
              top[7]  ? 5'd4  :
              top[6]  ? 5'd5  :
              top[5]  ? 5'd6  :
              top[4]  ? 5'd7  :
              top[3]  ? 5'd8  :
              top[2]  ? 5'd9  :
              top[1]  ? 5'd10 :
              top[0]  ? 5'd11 : 5'd12;

  wire [PROD_W-1:0] below;
  genvar k;
  generate
    for (k = 0; k < PROD_W; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  wire [5:0] exp_sum = {1'b0, ea_eff} + {1'b0, eb_eff};

  // Normal path: shift by sr-1, so shifted_n[0] is the guard bit and
  // shifted_n[11:1] is the 11-bit pre-round significand.
  wire [4:0] sr = 5'd11 - lz;
  wire [11:0] shifted_n;
  assign shifted_n = (sr == 5'd0) ? {prod[10:0], 1'b0} :
                                    (prod >> (sr - 5'd1));

  wire [10:0] mant_pre = shifted_n[11:1];
  wire        guard    = shifted_n[0];
  wire        sticky   = (sr <= 5'd1) ? 1'b0 : below[sr - 5'd2];
  wire        round_up = guard & (sticky | mant_pre[0]);

  wire [9:0] frac_inc  = mant_pre[9:0] + {9'd0, round_up};
  wire       carry     = round_up & (&mant_pre[9:0]);
  wire [9:0] frac_norm = carry ? 10'd0 : frac_inc;

  // exp_c = exp_sum + carry; normal unbiased exponent is exp_c-(14+lz).
  wire [5:0] exp_c        = exp_sum + {5'd0, carry};
  wire [5:0] limit_under  = 6'd14 + {1'b0, lz};
  wire [5:0] limit_over   = 6'd44 + {1'b0, lz};
  wire       underflow    = (exp_c <= limit_under);
  wire       overflow     = (exp_c >  limit_over);
  wire [5:0] e_norm       = exp_c - limit_under;
  wire [4:0] exp_norm     = e_norm[4:0];

  // Subnormal path: total_shift = 26 - exp_sum.
  wire [6:0] total_shift = 7'd26 - {1'b0, exp_sum};
  wire       ts_valid    = ({1'b0, exp_sum} <= 7'd25);
  wire [4:0] ts_m1       = total_shift[4:0] - 5'd1;
  wire [PROD_W-1:0] shifted_s = ts_valid ? (prod >> ts_m1) : {PROD_W{1'b0}};

  wire [9:0]  frac_trunc = shifted_s[10:1];
  wire        guard_s    = shifted_s[0];
  wire        sticky_s   = (total_shift >= 7'd2 && total_shift < 7'd23) ? below[total_shift[4:0] - 5'd2] : 1'b0;
  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [9:0]  frac_sub_inc = frac_trunc + {9'd0, round_s};
  wire        carry_s    = round_s & (&frac_trunc);
  wire [9:0]  frac_sub   = carry_s ? 10'd0 : frac_sub_inc;

  wire is_nan      = nan_in;
  wire is_inf      = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_sub      = ~is_nan & ~is_inf & underflow;
  wire sub_is_zero = is_sub & ~carry_s & (frac_sub_inc == 10'd0);
  wire is_zero     = ~is_nan & ~is_inf & (zero_in | sub_is_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : is_sub ? {4'd0, carry_s} : exp_norm;
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : is_sub ? frac_sub : frac_norm;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule
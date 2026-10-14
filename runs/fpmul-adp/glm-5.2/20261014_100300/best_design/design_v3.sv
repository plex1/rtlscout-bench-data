// fp_mul_e5f10 -- optimized variant 3: narrower exp arithmetic, tree LZ
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

  // classify
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

  // Tree-based leading-zero count
  wire [3:0] lz_lo4;  // LZ of lower 8 bits (0-8)
  wire [3:0] lz_hi4;  // LZ of bits [21:14] (0-8)
  wire [2:0] lz_mid;  // LZ of bits [13:8] (0-6)
  wire       mid_nz = |prod[13:8];
  wire       hi_nz  = |prod[21:14];
  // lower 8 bits [7:0]
  assign lz_lo4 =
    prod[7] ? 4'd0 : prod[6] ? 4'd1 : prod[5] ? 4'd2 : prod[4] ? 4'd3 :
    prod[3] ? 4'd4 : prod[2] ? 4'd5 : prod[1] ? 4'd6 : prod[0] ? 4'd7 : 4'd8;
  // mid 6 bits [13:8]
  assign lz_mid =
    prod[13] ? 3'd0 : prod[12] ? 3'd1 : prod[11] ? 3'd2 :
    prod[10] ? 3'd3 : prod[9]  ? 3'd4 : prod[8]  ? 3'd5 : 3'd6;
  // upper 8 bits [21:14]
  assign lz_hi4 =
    prod[21] ? 4'd0 : prod[20] ? 4'd1 : prod[19] ? 4'd2 : prod[18] ? 4'd3 :
    prod[17] ? 4'd4 : prod[16] ? 4'd5 : prod[15] ? 4'd6 : prod[14] ? 4'd7 : 4'd8;

  wire [4:0] lz =
    hi_nz  ? {1'b0, lz_hi4}            :  // 0-8
    mid_nz ? {2'b01, lz_mid}          :  // 8+6 = 14 max -> 8 to 13
    (|prod[7:0]) ? (5'd14 + {1'b0, lz_lo4}) :  // 14 to 21
    5'd22;

  wire [PROD_W-1:0] below;
  genvar k;
  generate
    for (k = 0; k < PROD_W; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  function bit_at;
    input [PROD_W-1:0] v;
    input integer      idx;
    begin
      bit_at = (idx >= 0 && idx < PROD_W) ? v[idx] : 1'b0;
    end
  endfunction

  // normalise and round
  wire              need_right = (lz <= 5'd11);
  wire [4:0]        sr         = 5'd11 - lz;
  wire [4:0]        sl         = lz - 5'd11;
  wire [PROD_W-1:0] shifted    = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0]       mant_pre   = shifted[10:0];

  wire guard    = need_right ? bit_at(prod,  sr - 1) : 1'b0;
  wire sticky   = need_right ? bit_at(below, sr - 2) : 1'b0;
  wire round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // exponent - use 7-bit arithmetic (max value is 67)
  wire [6:0] exp_sum     = {2'd0, ea_eff} + {2'd0, eb_eff};
  wire [6:0] lhs         = exp_sum + 7'd1 + {6'd0, carry};
  wire [6:0] limit_under = 7'd15 + {2'd0, lz};
  wire [6:0] limit_over  = 7'd45 + {2'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [6:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // subnormal result
  integer           total_shift;
  reg [PROD_W-1:0]  sub_shifted;
  always @* begin
    total_shift = (BIAS + 10 + 1) - $signed({1'b0, exp_sum});
    sub_shifted = (total_shift >= 0 && total_shift < PROD_W) ? (prod >> total_shift) : {PROD_W{1'b0}};
  end
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s    = bit_at(prod,  total_shift - 1);
  wire        sticky_s   = bit_at(below, total_shift - 2);
  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
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

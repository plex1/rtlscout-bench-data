// fp_mul_e5f10 -- optimized IEEE binary16 multiplier

module fp_mul_e5f10 (input  [15:0] a, input [15:0] b, output [15:0] y);

  localparam BIAS = 15;

  // unpack
  wire       s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0], f_b = b[9:0];

  // classify
  wire e_a_zero = (e_a == 5'd0), e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'd31), e_b_ones = (e_b == 5'd31);
  wire zero_a = e_a_zero & ~|f_a, zero_b = e_b_zero & ~|f_b;
  wire nan_a  = e_a_ones &  |f_a, nan_b  = e_b_ones &  |f_b;
  wire inf_a  = e_a_ones & ~|f_a, inf_b  = e_b_ones & ~|f_b;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // effective significands and exponents
  wire [10:0] m_a = {~e_a_zero, f_a};
  wire [10:0] m_b = {~e_b_zero, f_b};
  wire [5:0]  ea_eff = e_a_zero ? 6'd1 : {1'b0, e_a};
  wire [5:0]  eb_eff = e_b_zero ? 6'd1 : {1'b0, e_b};

  // significand product (22 bits)
  wire [21:0] prod = m_a * m_b;

  // leading-zero count of prod (priority encoder using OR-reductions)
  wire [21:0] lz_in = prod;
  wire [4:0] lz;
  wire v21 = lz_in[21];
  wire v20 = |lz_in[21:20];
  wire v19 = |lz_in[21:19];
  wire v18 = |lz_in[21:18];
  wire v17 = |lz_in[21:17];
  wire v16 = |lz_in[21:16];
  wire v15 = |lz_in[21:15];
  wire v14 = |lz_in[21:14];
  wire v13 = |lz_in[21:13];
  wire v12 = |lz_in[21:12];
  wire v11 = |lz_in[21:11];
  wire v10 = |lz_in[21:10];
  wire v9  = |lz_in[21:9];
  wire v8  = |lz_in[21:8];
  wire v7  = |lz_in[21:7];
  wire v6  = |lz_in[21:6];
  wire v5  = |lz_in[21:5];
  wire v4  = |lz_in[21:4];
  wire v3  = |lz_in[21:3];
  wire v2  = |lz_in[21:2];
  wire v1  = |lz_in[21:1];
  wire v0  = |lz_in[21:0];

  assign lz = v21 ? 5'd0 :
              v20 ? 5'd1 :
              v19 ? 5'd2 :
              v18 ? 5'd3 :
              v17 ? 5'd4 :
              v16 ? 5'd5 :
              v15 ? 5'd6 :
              v14 ? 5'd7 :
              v13 ? 5'd8 :
              v12 ? 5'd9 :
              v11 ? 5'd10 :
              v10 ? 5'd11 :
              v9  ? 5'd12 :
              v8  ? 5'd13 :
              v7  ? 5'd14 :
              v6  ? 5'd15 :
              v5  ? 5'd16 :
              v4  ? 5'd17 :
              v3  ? 5'd18 :
              v2  ? 5'd19 :
              v1  ? 5'd20 : 5'd21;

  // sticky: below[k] = |prod[k:0] (prefix OR)
  wire [21:0] below;
  assign below[0] = prod[0];
  generate
    for (genvar k = 1; k < 22; k = k + 1)
      assign below[k] = below[k-1] | prod[k];
  endgenerate

  // bit_at: safely read prod[idx] or below[idx], returning 0 if idx out of [0,21]
  function bit_at;
    input [21:0] v;
    input [4:0] idx;
    begin
      bit_at = (idx <= 5'd21) ? v[idx] : 1'b0;
    end
  endfunction

  // === NORMAL PATH ===
  wire need_right = (lz <= 5'd11);
  wire [4:0] sr = 5'd11 - lz;
  wire [4:0] sl = lz - 5'd11;

  wire [21:0] shifted_r = prod >> sr;
  wire [21:0] shifted_l = prod << sl;
  wire [21:0] shifted = need_right ? shifted_r : shifted_l;
  wire [10:0] mant_pre = shifted[10:0];

  wire guard_n    = need_right ? bit_at(prod, sr - 5'd1) : 1'b0;
  wire sticky_n   = need_right ? bit_at(below, sr - 5'd2) : 1'b0;
  wire round_up   = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // exponent for normal result
  wire [7:0] exp_sum     = {2'd0, ea_eff} + {2'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = 8'd15 + {3'd0, lz};
  wire [7:0] limit_over  = 8'd45 + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // === SUBNORMAL PATH ===
  wire [7:0] total_shift_8 = 8'd26 - exp_sum;
  wire [4:0] total_shift   = total_shift_8[4:0];
  wire       ts_lt_22      = (total_shift_8[7] == 1'b0) && (total_shift < 5'd22);
  wire [21:0] sub_shifted  = ts_lt_22 ? (prod >> total_shift) : 22'd0;
  wire [9:0]  frac_trunc   = sub_shifted[9:0];
  wire        guard_s      = bit_at(prod,  total_shift - 5'd1);
  wire        sticky_s     = bit_at(below, total_shift - 5'd2);
  wire        round_s      = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum     = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s      = frac_sum[10];
  wire [9:0]  frac_sub     = frac_sum[9:0];
  wire [4:0]  exp_sub      = {4'd0, carry_s};

  // === PACK ===
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

// fp_mul_e5f10 -- optimized IEEE binary16 multiplier (1+5+10)
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS = 15;
  localparam PROD_W = 22;

  // unpack
  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

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
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  // significand product
  wire [PROD_W-1:0] prod = m_a * m_b;

  // OR-below tree (needed for sticky bits)
  wire [21:0] below;
  assign below[0] = prod[0];
  assign below[1] = prod[1] | below[0];
  assign below[2] = prod[2] | below[1];
  assign below[3] = prod[3] | below[2];
  assign below[4] = prod[4] | below[3];
  assign below[5] = prod[5] | below[4];
  assign below[6] = prod[6] | below[5];
  assign below[7] = prod[7] | below[6];
  assign below[8] = prod[8] | below[7];
  assign below[9] = prod[9] | below[8];
  assign below[10] = prod[10] | below[9];
  assign below[11] = prod[11] | below[10];
  assign below[12] = prod[12] | below[11];
  assign below[13] = prod[13] | below[12];
  assign below[14] = prod[14] | below[13];
  assign below[15] = prod[15] | below[14];
  assign below[16] = prod[16] | below[15];
  assign below[17] = prod[17] | below[16];
  assign below[18] = prod[18] | below[17];
  assign below[19] = prod[19] | below[18];
  assign below[20] = prod[20] | below[19];
  assign below[21] = prod[21] | below[20];

  function [2:0] lead_grp6;
    input [5:0] np; // inverted prod bits (1 = zero)
    begin
      if (~np[5])      lead_grp6 = 3'd0;
      else if (~np[4]) lead_grp6 = 3'd1;
      else if (~np[3]) lead_grp6 = 3'd2;
      else if (~np[2]) lead_grp6 = 3'd3;
      else if (~np[1]) lead_grp6 = 3'd4;
      else             lead_grp6 = 3'd5;
    end
  endfunction

  function [1:0] lead_grp4;
    input [3:0] np;
    begin
      if (~np[3])      lead_grp4 = 2'd0;
      else if (~np[2]) lead_grp4 = 2'd1;
      else if (~np[1]) lead_grp4 = 2'd2;
      else             lead_grp4 = 2'd3;
    end
  endfunction

  // leading-zero count via priority encoder
  wire [21:0] not_prod = ~prod;
  wire g3 = |not_prod[21:18];
  wire g2 = |not_prod[17:12];
  wire g1 = |not_prod[11:6];
  wire g0 = |not_prod[5:0];

  reg [4:0] lz;
  always @(*) begin
    if (prod == 0)    lz = 5'd0;
    else if (g3)      lz = {2'b0, lead_grp4(not_prod[21:18])};
    else if (g2)      lz = 5'd4 + {1'b0, lead_grp6(not_prod[17:12])};
    else if (g1)      lz = 5'd10 + {1'b0, lead_grp6(not_prod[11:6])};
    else              lz = 5'd16 + {1'b0, lead_grp6(not_prod[5:0])};
  end

  function bit_at;
    input [PROD_W-1:0] v;
    input integer      idx;
    begin
      bit_at = (idx >= 0 && idx < PROD_W) ? v[idx] : 1'b0;
    end
  endfunction

  // normalise and round a normal result
  wire        need_right = (lz <= 5'd11);
  wire [4:0]  sr         = 5'd11 - lz;
  wire [4:0]  sl         = lz - 5'd11;
  wire [PROD_W-1:0] shifted = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0] mant_pre  = shifted[10:0];

  wire guard    = need_right ? bit_at(prod,  sr - 1) : 1'b0;
  wire sticky   = need_right ? bit_at(below, sr - 2) : 1'b0;
  wire round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // exponent of a normal result
  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = 8'd15 + {3'd0, lz};
  wire [7:0] limit_over  = 8'd45 + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // subnormal result
  integer           total_shift;
  reg [PROD_W-1:0]  sub_shifted;
  always @* begin
    total_shift = (15 + 10 + 1) - $signed({1'b0, exp_sum});
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

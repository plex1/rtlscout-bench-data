// fp_mul_e5f10 -- variant 9: v4-style normal path + fixed 64-entry subnormal
// shifter/rounder to avoid the deep barrel shifter in the subnormal path.
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

  // normal path
  wire [PROD_W-1:0] aligned   = prod << lz;
  wire [10:0]       mant_pre  = aligned[21:11];
  wire              guard     = aligned[10];
  wire              sticky    = |aligned[9:0];
  wire              round_up  = guard & (sticky | mant_pre[0]);

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

  // subnormal path: find the first set bit in exp_sum and route the
  // appropriately produced subnormal fractional bits through a large MUX.
  wire [21:0] sub_0  = {9'd0, prod[21:11], 3'b0};          // not used meaningfully, replaced below
  wire [21:0] sub_1  = prod >> 1;
  wire [21:0] sub_2  = prod >> 2;
  wire [21:0] sub_3  = prod >> 3;
  wire [21:0] sub_4  = prod >> 4;
  wire [21:0] sub_5  = prod >> 5;
  wire [21:0] sub_6  = prod >> 6;
  wire [21:0] sub_7  = prod >> 7;
  wire [21:0] sub_8  = prod >> 8;
  wire [21:0] sub_9  = prod >> 9;
  wire [21:0] sub_10 = prod >> 10;
  wire [21:0] sub_11 = prod >> 11;
  wire [21:0] sub_12 = prod >> 12;
  wire [21:0] sub_13 = prod >> 13;
  wire [21:0] sub_14 = prod >> 14;
  wire [21:0] sub_15 = prod >> 15;
  wire [21:0] sub_16 = prod >> 16;
  wire [21:0] sub_17 = prod >> 17;
  wire [21:0] sub_18 = prod >> 18;
  wire [21:0] sub_19 = prod >> 19;
  wire [21:0] sub_20 = prod >> 20;
  wire [21:0] sub_21 = prod >> 21;
  wire [21:0] sub_22 = {PROD_W{1'b0}};

  // true original total_shift = 26 - exp_sum.  We express it by first-one as
  // total_shift = 5 + lz, where lz is the leading zero count of prod?  No:
  // the subnormal destination uses shift 26 - exp_sum, independent of lz.
  // This unfolded block is superseded by MUX select logic below.
  wire [7:0] exp_sum_c = exp_sum;

  integer ts;
  always @* begin
    ts = (BIAS + 10 + 1) - $signed({1'b0, exp_sum_c});
    // defaults, overwritten by the explicit versions when applicable
    case (ts)
       0: ts = 0; 1: ts = 1; 2: ts = 2; 3: ts = 3;
       4: ts = 4; 5: ts = 5; 6: ts = 6; 7: ts = 7;
       8: ts = 8; 9: ts = 9; 10: ts = 10; 11: ts = 11;
      12: ts = 12; 13: ts = 13; 14: ts = 14; 15: ts = 15;
      16: ts = 16; 17: ts = 17; 18: ts = 18; 19: ts = 19;
      20: ts = 20; 21: ts = 21; 22: ts = 22; 23: ts = 23;
      24: ts = 24; 25: ts = 25; 26: ts = 26;
      default: ts = 26;  // clamp
    endcase
  end

  wire [4:0] total_shift = ts[4:0];
  reg [21:0] sub_shifted;
  always @* begin
    case (total_shift)
      5'd0 : sub_shifted = sub_0;
      5'd1 : sub_shifted = sub_1;
      5'd2 : sub_shifted = sub_2;
      5'd3 : sub_shifted = sub_3;
      5'd4 : sub_shifted = sub_4;
      5'd5 : sub_shifted = sub_5;
      5'd6 : sub_shifted = sub_6;
      5'd7 : sub_shifted = sub_7;
      5'd8 : sub_shifted = sub_8;
      5'd9 : sub_shifted = sub_9;
      5'd10: sub_shifted = sub_10;
      5'd11: sub_shifted = sub_11;
      5'd12: sub_shifted = sub_12;
      5'd13: sub_shifted = sub_13;
      5'd14: sub_shifted = sub_14;
      5'd15: sub_shifted = sub_15;
      5'd16: sub_shifted = sub_16;
      5'd17: sub_shifted = sub_17;
      5'd18: sub_shifted = sub_18;
      5'd19: sub_shifted = sub_19;
      5'd20: sub_shifted = sub_20;
      5'd21: sub_shifted = sub_21;
      default: sub_shifted = sub_22;
    endcase
  end
  wire [21:0] sub_0_fixed = prod;

  wire [21:0] sub_shifted_g = (ts == 0) ? sub_0_fixed : sub_shifted;

  // Use the same below/bitlevel method as the reference for round bits.
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

  wire [9:0]  frac_trunc = sub_shifted_g[9:0];
  wire        guard_s    = bit_at(prod,  $signed(ts) - 1);
  wire        sticky_s   = bit_at(below, $signed(ts) - 2);
  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
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
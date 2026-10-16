// fp_mul_e5f10 -- binary16 multiplier with input pre-normalization.
// Each input significand is normalized before the multiply, so the product
// leading bit is either bit 20 or 21 (no wide product LZ encoder is needed).
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

  wire [4:0] ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0] eb_eff = e_b_zero ? 5'd1 : e_b;

  // Leading-zero count of the raw fraction (only meaningful for subnormals).
  function [3:0] lzc10(input [9:0] v);
    begin
      lzc10 = (v[9]) ? 4'd0 :
              (v[8]) ? 4'd1 :
              (v[7]) ? 4'd2 :
              (v[6]) ? 4'd3 :
              (v[5]) ? 4'd4 :
              (v[4]) ? 4'd5 :
              (v[3]) ? 4'd6 :
              (v[2]) ? 4'd7 :
              (v[1]) ? 4'd8 :
              (v[0]) ? 4'd9 : 4'd10;
    end
  endfunction

  wire [3:0] lzc_a = lzc10(f_a);
  wire [3:0] lzc_b = lzc10(f_b);

  // Shift needed to make a subnormal significand's leading one sit on bit 10.
  wire [4:0] extra_a = e_a_zero ? ({1'b0, lzc_a} + 5'd1) : 5'd0;
  wire [4:0] extra_b = e_b_zero ? ({1'b0, lzc_b} + 5'd1) : 5'd0;

  // Pre-normalized 11-bit significands.
  wire [10:0] raw_a = {1'b0, f_a};
  wire [10:0] raw_b = {1'b0, f_b};
  wire [10:0] m_a_n = e_a_zero ? (raw_a << extra_a) : {1'b1, f_a};
  wire [10:0] m_b_n = e_b_zero ? (raw_b << extra_b) : {1'b1, f_b};

  wire [PROD_W-1:0] prod = m_a_n * m_b_n;

  // Product of two normalized significands is in [2^20, 2^22).
  wire [4:0] extrasum = {1'b0, extra_a} + {1'b0, extra_b};
  wire       lz       = prod[21] ? 5'd0 : 5'd1;
  wire [4:0] lz_orig  = zero_in ? 5'd0 : (extrasum + lz);

  function bit_at;
    input [PROD_W-1:0] v;
    input integer      idx;
    begin
      bit_at = (idx >= 0 && idx < PROD_W) ? v[idx] : 1'b0;
    end
  endfunction

  wire [PROD_W-1:0] below;
  genvar k;
  generate
    for (k = 0; k < PROD_W; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  // Normal-product rounding: shift right by 11 if bit21 set, else by 10.
  wire              lz_is0    = (lz == 5'd0);
  wire [PROD_W-1:0] shifted   = lz_is0 ? (prod >> 5'd11) : (prod >> 5'd10);
  wire [10:0]       mant_pre  = shifted[10:0];
  wire              guard     = lz_is0 ? bit_at(prod, 5'd10) : bit_at(prod, 5'd9);
  wire              sticky    = lz_is0 ? bit_at(below, 5'd9) : bit_at(below, 5'd8);
  wire              round_up  = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz_orig};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz_orig};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // Subnormal path: shift normalized product by q = 26 - exp_sum + extrasum.
  integer          q;
  reg [PROD_W-1:0] sub_shifted;
  always @* begin
    q = (BIAS + 10 + 1) - $signed({1'b0, exp_sum}) + $signed({1'b0, extrasum});
    sub_shifted = (q >= 0 && q < PROD_W) ? (prod >> q) : {PROD_W{1'b0}};
  end
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s    = bit_at(prod, q - 1);
  wire        sticky_s   = bit_at(below, q - 2);
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
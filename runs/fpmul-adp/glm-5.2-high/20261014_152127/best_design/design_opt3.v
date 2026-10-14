// fp_mul_e5f10 -- optimized v3: left-shift normalization for normal path
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

  // Leading-zero count: priority encoder from MSB
  wire [4:0] lz =     prod[21] ? 5'd0  :
                      prod[20] ? 5'd1  :
                      prod[19] ? 5'd2  :
                      prod[18] ? 5'd3  :
                      prod[17] ? 5'd4  :
                      prod[16] ? 5'd5  :
                      prod[15] ? 5'd6  :
                      prod[14] ? 5'd7  :
                      prod[13] ? 5'd8  :
                      prod[12] ? 5'd9  :
                      prod[11] ? 5'd10 :
                      prod[10] ? 5'd11 :
                      prod[9]  ? 5'd12 :
                      prod[8]  ? 5'd13 :
                      prod[7]  ? 5'd14 :
                      prod[6]  ? 5'd15 :
                      prod[5]  ? 5'd16 :
                      prod[4]  ? 5'd17 :
                      prod[3]  ? 5'd18 :
                      prod[2]  ? 5'd19 :
                      prod[1]  ? 5'd20 :
                      prod[0]  ? 5'd21 : 5'd0;

  // Normal path: left-shift by lz to put leading 1 at bit 21 (MSB)
  // mant_pre = bits [21:11], guard = bit 10, sticky = |bits[9:0]
  wire [PROD_W-1:0] prod_norm = prod << lz;
  wire [10:0]       mant_pre  = prod_norm[21:11];
  wire              guard_n   = prod_norm[10];
  wire              sticky_n  = |prod_norm[9:0];
  wire              round_up = guard_n & (sticky_n | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // Exponent of normal result
  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = BIAS + {3'd0, lz};
  wire [7:0] limit_over  = BIAS + MAX_FINITE_E + {3'd0, lz};
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // Subnormal path: right-shift prod by total_shift
  // below array for sticky computation
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

  // Pack
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

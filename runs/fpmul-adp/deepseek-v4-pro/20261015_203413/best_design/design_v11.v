// fp_mul_e5f10 -- variant 11: normalise subnormal operands first.  The
// 11x11 product then only needs a 1-bit normal shift; the complex product
// leading-zero count and the left barrel shifter disappear.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam MAX_E = 31;
  localparam PROD_W = 22;

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

  // Operand leading-zero counts (0 .. 9).  Zero fractions are handled by the
  // zero path; their lz value is arbitrary (we leave it at 9 through the loop
  // default only for nonzero; zero gives 0), but the mantissa is zero anyway.
  reg [3:0] lz_a;
  reg [3:0] lz_b;
  integer   ni;
  always @* begin
    lz_a = 4'd0;
    lz_b = 4'd0;
    for (ni = 0; ni < 10; ni = ni + 1) begin
      if (f_a[ni]) lz_a = 4'd9 - ni[3:0];
      if (f_b[ni]) lz_b = 4'd9 - ni[3:0];
    end
  end

  // Normalised significands: always in [1024,2047] for non-zero operands.
  wire [10:0] m_a = e_a_zero ? (f_a << (lz_a + 4'd1)) : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? (f_b << (lz_b + 4'd1)) : {1'b1, f_b};

  // Value = m * 2^q.
  wire signed [8:0] q_a = e_a_zero ? (-25 - $signed({1'b0, lz_a}))
                                   : ($signed({1'b0, e_a}) - 25);
  wire signed [8:0] q_b = e_b_zero ? (-25 - $signed({1'b0, lz_b}))
                                   : ($signed({1'b0, e_b}) - 25);
  wire signed [8:0] q_sum = q_a + q_b;

  wire [PROD_W-1:0] prod = m_a * m_b;

  // p: leading one at bit 21 or bit 20 for a non-zero product.
  wire        p_hi     = prod[21];
  wire [10:0] mant_pre = p_hi ? prod[21:11] : prod[20:10];
  wire        guard    = p_hi ? prod[10]    : prod[9];
  wire        sticky   = p_hi ? |prod[9:0]  : |prod[8:0];
  wire        round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // Biased exponent for a normal rounded result.
  // E = q_sum + h + 15 + carry.
  reg signed [8:0] Eres;
  always @* begin
    Eres = q_sum + (p_hi ? 9'sd21 : 9'sd20) + 9'sd15 + {8'sd0, carry};
  end

  wire overflow = ~zero_in & (Eres > 30);
  wire is_sub   = ~zero_in & (Eres < 1) & ~overflow;

  // Subnormal path.  Because prod = original_prod << LZtotal_eff and
  // q_sum is the normalised-product value exponent, the right shift that
  // lands the product on the subnormal grid is -q_sum - 24.
  wire signed [8:0] total_shift = -q_sum - 9'sd24;

  wire       shift_ok = (total_shift >= 0) && (total_shift < PROD_W);
  wire [4:0] sh_amt   = total_shift[4:0];

  wire [PROD_W-1:0] sub_shifted = shift_ok ? (prod >> sh_amt) : {PROD_W{1'b0}};

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

  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s    = bit_at(prod,  total_shift - 1);
  wire        sticky_s   = bit_at(below, total_shift - 2);
  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];
  wire [4:0]  exp_sub    = {4'd0, carry_s};

  wire is_nan      = nan_in;
  wire is_inf      = ~is_nan & (inf_in | overflow);
  wire sub_is_zero = is_sub & (exp_sub == 5'd0) & (frac_sub == 10'd0);
  wire is_zero     = ~is_nan & ~is_inf & (zero_in | sub_is_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 :
                      is_sub ? exp_sub : Eres[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 :
                      is_sub ? frac_sub : frac_norm;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule
// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
// Optimized: one significand product, one barrel shift, one rounding adder shared by
// the normal and the subnormal paths.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam MAX_E = 31;

  // 1. unpack
  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  // 2. classify
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

  // 3. effective significands and exponents
  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  // 4. significand product, leading-zero count, prefix ORs below[k] = |prod[k:0]
  wire [21:0] prod = m_a * m_b;

  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = (k == 0) ? prod[0] : (below[k-1] | prod[k]);
    end
  endgenerate

  // 5. exponent bookkeeping
  //    true result = prod * 2^(exp_sum-50) = 1.p * 2^(exp_sum-29-lz)
  wire [7:0] exp_sum     = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] limit_under = 8'd15 + {3'd0, lz};   // result < 2^-14  <=>  exp_sum < limit_under
  wire [7:0] limit_over  = 8'd45 + {3'd0, lz};

  wire underflow = (exp_sum < limit_under);
  wire [7:0] e_norm0 = exp_sum + 8'd1 - limit_under;   // >= 1 when !underflow

  // 6. single normalising shift
  //    !underflow : lz <= 11, shift = 11 - lz   (leading one lands on bit 10)
  //    underflow  : exp_sum <= 26, shift = 26 - exp_sum (result lands on bits 9:0)
  wire [4:0] sh = underflow ? (5'd26 - exp_sum[4:0]) : (5'd11 - lz);

  wire [10:0] mant_pre = prod >> sh;
  wire        guard    = (sh >= 5'd1) && (sh <= 5'd22) ? prod[sh - 5'd1] : 1'b0;
  wire        sticky   = (sh >= 5'd2) && (sh <= 5'd23) ? below[sh - 5'd2] : 1'b0;
  wire        round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] rnd      = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry    = underflow ? rnd[10] : rnd[11];

  wire [9:0] frac_n = (carry & ~underflow) ? rnd[10:1] : rnd[9:0];

  // 7. exponent / overflow
  wire [7:0] exp_raw = underflow ? {7'd0, carry} : (e_norm0 + {7'd0, carry});
  wire       overflow = (exp_raw > 8'd30);

  // 8. pack
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_zero = ~is_nan & ~is_inf & (zero_in | (underflow & ~carry & (frac_n == 10'd0)));

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_raw[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_n;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

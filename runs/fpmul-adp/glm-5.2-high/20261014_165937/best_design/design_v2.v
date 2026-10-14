// fp_mul_e5f10 -- optimized IEEE binary16 multiplier
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS = 15;

  // 1. unpack
  wire        sa = a[15];
  wire        sb = b[15];
  wire [4:0]  ea = a[14:10];
  wire [4:0]  eb = b[14:10];
  wire [9:0]  fa = a[9:0];
  wire [9:0]  fb = b[9:0];

  // 2. classify
  wire ea_zero = ~|ea;
  wire eb_zero = ~|eb;
  wire ea_max  = &ea;
  wire eb_max  = &eb;
  wire a_zero  = ea_zero & ~|fa;
  wire b_zero  = eb_zero & ~|fb;
  wire a_nan   = ea_max  & |fa;
  wire b_nan   = eb_max  & |fb;
  wire a_inf   = ea_max  & ~|fa;
  wire b_inf   = eb_max  & ~|fb;

  wire nan_in  = a_nan | b_nan | (a_inf & b_zero) | (a_zero & b_inf);
  wire inf_in  = a_inf | b_inf;
  wire zero_in = a_zero | b_zero;
  wire s_y     = sa ^ sb;

  // 3. effective significands and exponents
  wire [10:0] m_a     = {~ea_zero, fa};
  wire [10:0] m_b     = {~eb_zero, fb};
  wire [5:0]  ea_eff  = ea_zero ? 6'd1 : {1'b0, ea};
  wire [5:0]  eb_eff  = eb_zero ? 6'd1 : {1'b0, eb};

  // 4. product
  wire [21:0] prod = m_a * m_b;

  // 5. leading-zero count via priority encoder
  reg [4:0] lz;
  always @* begin
    casez (prod)
      22'b1????????????????????: lz = 5'd0;
      22'b01???????????????????: lz = 5'd1;
      22'b001??????????????????: lz = 5'd2;
      22'b0001?????????????????: lz = 5'd3;
      22'b00001????????????????: lz = 5'd4;
      22'b000001???????????????: lz = 5'd5;
      22'b0000001??????????????: lz = 5'd6;
      22'b00000001?????????????: lz = 5'd7;
      22'b000000001????????????: lz = 5'd8;
      22'b0000000001???????????: lz = 5'd9;
      22'b00000000001??????????: lz = 5'd10;
      22'b000000000001?????????: lz = 5'd11;
      22'b0000000000001????????: lz = 5'd12;
      22'b00000000000001???????: lz = 5'd13;
      22'b000000000000001??????: lz = 5'd14;
      22'b0000000000000001?????: lz = 5'd15;
      22'b00000000000000001????: lz = 5'd16;
      22'b000000000000000001???: lz = 5'd17;
      22'b0000000000000000001??: lz = 5'd18;
      22'b00000000000000000001?: lz = 5'd19;
      22'b000000000000000000001: lz = 5'd20;
      22'b000000000000000000000: lz = 5'd21;
      default:                    lz = 5'd21;
    endcase
  end

  // 6. sticky tree: below[k] = |prod[k:0]
  wire [21:0] below;
  assign below[0] = prod[0];
  assign below[1] = prod[1] | prod[0];
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

  // 7. normalize: place leading 1 at bit 10
  wire        need_right = (lz <= 5'd11);
  wire [5:0]  sr = 6'd11 - {1'b0, lz};
  wire [5:0]  sl = {1'b0, lz} - 6'd11;
  wire [21:0] shifted_r = prod >> sr;
  wire [21:0] shifted_l = prod << sl;
  wire [21:0] shifted   = need_right ? shifted_r : shifted_l;
  wire [10:0] mant_pre  = shifted[10:0];

  // guard & sticky from raw product
  wire guard    = need_right ? (sr >= 6'd1  ? prod[sr-1 +: 1] : 1'b0) : 1'b0;
  wire sticky   = need_right ? (sr >= 6'd2  ? below[sr-2] : 1'b0) : 1'b0;
  wire round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [10:0] mant_post  = carry ? mant_round[11:1] : mant_round[10:0];
  wire [9:0]  frac_norm  = mant_post[9:0];

  // 8. exponent of normal result
  wire [7:0] exp_sum     = {2'd0, ea_eff} + {2'd0, eb_eff};
  wire [7:0] lhs         = exp_sum + 8'd1 + {7'd0, carry};
  wire [7:0] limit_under = 8'd15 + {3'd0, lz};
  wire [7:0] limit_over  = 8'd45 + {3'd0, lz};  // 15 + 30 + lz
  wire       underflow   = (lhs <= limit_under);
  wire       overflow    = (lhs >  limit_over);
  wire [7:0] e_norm      = lhs - limit_under;
  wire [4:0] exp_norm    = e_norm[4:0];

  // 9. subnormal result
  wire signed [8:0] total_shift = $signed({1'b0, 8'd26}) - $signed({1'b0, exp_sum});
  wire [21:0] sub_shifted = (total_shift >= 0 && total_shift < 22) ? (prod >> total_shift) : 22'd0;
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s   = (total_shift >= 1 && total_shift < 22) ? prod[total_shift - 1] : 1'b0;
  wire        sticky_s  = (total_shift >= 2 && total_shift < 22) ? below[total_shift - 2] : 1'b0;
  wire        round_s   = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum  = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s   = frac_sum[10];
  wire [9:0]  frac_sub  = frac_sum[9:0];
  wire [4:0]  exp_sub   = {4'd0, carry_s};

  // 10. pack
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

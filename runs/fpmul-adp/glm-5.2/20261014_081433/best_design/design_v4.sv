// fp_mul_e5f10 -- behavioral, optimized for synthesis
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS = 15;

  wire        sa = a[15];
  wire        sb = b[15];
  wire [4:0]  ea = a[14:10];
  wire [4:0]  eb = b[14:10];
  wire [9:0]  fa = a[9:0];
  wire [9:0]  fb = b[9:0];

  // classify
  wire a_is_zero = (ea == 0) && (fa == 0);
  wire b_is_zero = (eb == 0) && (fb == 0);
  wire a_is_inf  = (ea == 31) && (fa == 0);
  wire b_is_inf  = (eb == 31) && (fb == 0);
  wire a_is_nan  = (ea == 31) && (fa != 0);
  wire b_is_nan  = (eb == 31) && (fb != 0);
  wire a_sub    = (ea == 0) && (fa != 0);
  wire b_sub    = (eb == 0) && (fb != 0);

  wire is_nan = a_is_nan | b_is_nan | (a_is_inf & b_is_zero) | (a_is_zero & b_is_inf);
  wire is_inf_in = a_is_inf | b_is_inf;
  wire is_zero_in = a_is_zero | b_is_zero;
  wire sign_out = sa ^ sb;

  // effective significands (with hidden bit for normal, without for subnormal)
  wire [10:0] ma = {(ea != 0), fa};
  wire [10:0] mb = {(eb != 0), fb};
  // effective exponents
  wire [5:0] ea_eff = ea[4:0] == 0 ? 6'd1 : {1'b0, ea};
  wire [5:0] eb_eff = eb[4:0] == 0 ? 6'd1 : {1'b0, eb};

  wire [21:0] prod = ma * mb;

  // leading zero count
  wire [4:0] lz;
  assign lz = (prod == 0) ? 5'd22 :
    prod[21] ? 5'd0 : prod[20] ? 5'd1 : prod[19] ? 5'd2 : prod[18] ? 5'd3 :
    prod[17] ? 5'd4 : prod[16] ? 5'd5 : prod[15] ? 5'd6 : prod[14] ? 5'd7 :
    prod[13] ? 5'd8 : prod[12] ? 5'd9 : prod[11] ? 5'd10 : prod[10] ? 5'd11 :
    prod[9] ? 5'd12 : prod[8] ? 5'd13 : prod[7] ? 5'd14 : prod[6] ? 5'd15 :
    prod[5] ? 5'd16 : prod[4] ? 5'd17 : prod[3] ? 5'd18 : prod[2] ? 5'd19 :
    prod[1] ? 5'd20 : 5'd21;

  // normal path: normalize to leading 1 at bit 10
  wire        need_right = (lz <= 5'd11);
  wire [5:0]  sr = 6'd11 - {1'b0, lz};
  wire [5:0]  sl = {1'b0, lz} - 6'd11;
  wire [21:0] shifted_n = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0] mant_pre = shifted_n[10:0];
  wire        guard_n  = need_right ? ((sr >= 1 && sr <= 21) ? prod[sr-1] : 1'b0) : 1'b0;
  wire        sticky_n = need_right ? ((sr >= 2) ? |prod[sr-2:0] : 1'b0) : 1'b0;
  wire        round_n  = guard_n & (sticky_n | mant_pre[0]);
  wire [11:0] mant_rnd = {1'b0, mant_pre} + {11'd0, round_n};
  wire        carry_n  = mant_rnd[11];
  wire [9:0]  frac_n   = carry_n ? mant_rnd[10:1] : mant_rnd[9:0];

  // normal exponent
  wire [6:0] exp_sum = {1'b0, ea_eff} + {1'b0, eb_eff};
  wire [7:0] exp_raw = {1'b0, exp_sum} + 8'd1 + {7'd0, carry_n} - {3'd0, lz} - 8'd15;
  wire       overflow = (exp_raw > 8'd30) || (exp_raw[7]);
  wire       underflow = (exp_raw == 0) || (exp_raw[7]);
  wire [4:0] exp_n = exp_raw[4:0];

  // subnormal path
  wire [6:0] tshift = 7'd26 - {1'b0, exp_sum};
  wire       tvalid = (tshift < 7'd22) && ~tshift[6];
  wire [21:0] shifted_s = tvalid ? (prod >> tshift) : 22'd0;
  wire [9:0]  frac_s_trunc = shifted_s[9:0];
  wire        guard_s  = (tshift >= 1 && tshift < 22) ? prod[tshift-1] : 1'b0;
  wire        sticky_s = (tshift >= 2 && tshift < 22) ? |prod[tshift-2:0] : 1'b0;
  wire        round_s  = guard_s & (sticky_s | frac_s_trunc[0]);
  wire [10:0] frac_s_rnd = {1'b0, frac_s_trunc} + {10'd0, round_s};
  wire        carry_s  = frac_s_rnd[10];
  wire [9:0]  frac_s   = frac_s_rnd[9:0];

  // select
  wire is_inf_out = ~is_nan & (is_inf_in | (overflow & ~is_zero_in));
  wire is_sub_out = ~is_nan & ~is_inf_out & underflow;
  wire sub_zero   = is_sub_out & ~carry_s & (frac_s == 0);
  wire is_zero_out = ~is_nan & ~is_inf_out & (is_zero_in | sub_zero);

  wire [4:0] exp_out = (is_nan | is_inf_out) ? 5'd31 : is_zero_out ? 5'd0 : is_sub_out ? {4'd0, carry_s} : exp_n;
  wire [9:0] frac_out = is_nan ? 10'h200 : (is_inf_out | is_zero_out) ? 10'd0 : is_sub_out ? frac_s : frac_n;

  assign y = {(is_nan ? 1'b0 : sign_out), exp_out, frac_out};
endmodule

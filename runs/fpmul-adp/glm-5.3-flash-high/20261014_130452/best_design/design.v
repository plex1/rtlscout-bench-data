// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits)
// with full subnormal support, round-to-nearest ties-to-even.
//
// Optimisations vs. a plain implementation:
//  * each significand is pre-normalised (leading one onto bit 10) with a small
//    11-bit CLZ + shift and its exponent adjusted; the product's leading one
//    is then always at bit 21 or 20, so no wide priority encoder is needed.
//  * with ma_n = 2^10 + fa, mb_n = 2^10 + fb:
//      prod = 2^20 + (fa+fb)<<10 + fa*fb = ((fa+fb+2^10)<<10) + fa*fb
//    i.e. only a 10x10 multiplier plus one 20-bit add.
//  * a single barrel shifter (prod >> sm1) yields mantissa and guard; sticky
//    comes from a constant mask shifted by the same amount.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack / classify ----
  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = ~|e_a;
  wire e_b_zero = ~|e_b;
  wire e_a_max  =  &e_a;
  wire e_b_max  =  &e_b;
  wire f_a_zero = ~|f_a;
  wire f_b_zero = ~|f_b;

  wire zero_a = e_a_zero & f_a_zero;
  wire zero_b = e_b_zero & f_b_zero;
  wire nan_a  = e_a_max  & ~f_a_zero;
  wire nan_b  = e_b_max  & ~f_b_zero;
  wire inf_a  = e_a_max  &  f_a_zero;
  wire inf_b  = e_b_max  &  f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // ---- pre-normalise significands ----
  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  function [3:0] clz11;
    input [10:0] m;
    begin
      clz11 = m[10] ? 4'd0  :
              m[9]  ? 4'd1  :
              m[8]  ? 4'd2  :
              m[7]  ? 4'd3  :
              m[6]  ? 4'd4  :
              m[5]  ? 4'd5  :
              m[4]  ? 4'd6  :
              m[3]  ? 4'd7  :
              m[2]  ? 4'd8  :
              m[1]  ? 4'd9  :
              m[0]  ? 4'd10 : 4'd11;
    end
  endfunction

  wire [3:0]  la   = clz11(m_a);
  wire [3:0]  lb   = clz11(m_b);
  wire [10:0] ma_n = m_a << la;
  wire [10:0] mb_n = m_b << lb;
  wire [9:0]  fa   = ma_n[9:0];
  wire [9:0]  fb   = mb_n[9:0];

  // ---- significand product: ((fa+fb+2^10)<<10) + fa*fb ----
  wire [10:0] fsum   = fa + fb;
  wire [11:0] fsum12 = {1'b0, fsum} + 12'd1024;
  wire [19:0] fx     = fa * fb;
  wire [21:0] prod   = {fsum12, 10'd0} + {2'd0, fx};
  wire        hi     = prod[21];      // leading one at bit 21 rather than 20

  // ---- exponents (signed) ----
  wire [7:0] ea_n = {3'd0, ea_eff} - {4'd0, la};
  wire [7:0] eb_n = {3'd0, eb_eff} - {4'd0, lb};
  wire [7:0] S    = ea_n + eb_n;
  wire [7:0] Sx   = S + {7'd0, hi};

  // result is subnormal (or zero) iff value < 2^-14, i.e. signed Sx <= 15
  wire under = Sx[7] | ~|Sx[7:4];

  // alignment: subnormal grid -> prod >> (26 - S); normal -> prod >> (10 + hi)
  wire [7:0] s8  = 8'd26 - S;
  wire [5:0] s   = under ? s8[5:0] : (6'd10 + {5'd0, hi});
  wire [5:0] sm1 = s - 6'd1;
  // clamp: S <= 3 (s >= 23) -> result is 0; shift 22 already yields
  // guard = prod[21] = 0-equivalent and sticky = |prod, so 5-bit control suffices
  wire       sbig = S[7] | ~|(S[6:2]);   // signed S <= 3  ->  result is zero
  wire [4:0] sm1c = sbig ? 5'd22 : sm1[4:0];

  wire [21:0] shg      = prod >> sm1c;
  wire [10:0] mant_pre = shg[11:1];
  wire        guard    = shg[0];
  wire        sticky   = |(prod & ~(22'h3FFFFF << sm1c));

  // ---- round to nearest, ties to even ----
  wire        round_up = guard & (sticky | mant_pre[0]);
  wire [11:0] mant_r   = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry    = under ? mant_r[10] : mant_r[11];

  wire [7:0] expn = Sx + {7'd0, carry} - 8'd15;
  wire       ovf  = ~under & (expn >= 8'd31);

  // ---- pack ----
  assign y = nan_in         ? 16'h7E00 :
             (inf_in | ovf) ? {s_y, 5'd31, 10'd0} :
             zero_in        ? {s_y, 15'd0} :
             under          ? (carry ? {s_y, 5'd1, 10'd0} : {s_y, 5'd0, mant_r[9:0]}) :
                              {s_y, expn[4:0], carry ? mant_r[10:1] : mant_r[9:0]};
endmodule

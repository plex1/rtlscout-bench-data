// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits)
// with full subnormal support, round-to-nearest ties-to-even.
//
// Structure:
//  * each significand is pre-normalised (leading one onto bit 10) with a small
//    11-bit CLZ + shift and its exponent adjusted; the product's leading one
//    is then always at bit 21 or 20.
//  * product decomposition: with ma_n = 2^10+fa, mb_n = 2^10+fb,
//      prod = 2^20 + (fa+fb)<<10 + fa*fb = (A << 10) + lo10
//    where A = fa + fb + 2^10 + fx[19:10] (one 12-bit add, max 4093)
//    and lo10 = fx[9:0].  prod[21:10] = A, prod[9:0] = lo10.
//  * Because lo10 < 2^10, prod >> k is carry-free for every k used:
//      normal result (k = 9+hi):      pure muxes on A/lo10
//      subnormal result (k = 25-S):   a 12-bit shift of A only
//    so no wide adder and no wide barrel shifter are needed.

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

  // ---- product decomposition: prod = (A << 10) + lo ----
  wire [10:0] fs = fa + fb;
  wire [19:0] fx = fa * fb;
  wire [11:0] A  = {1'b0, fs} + 12'd1024 + {2'b0, fx[19:10]};
  wire [9:0]  lo = fx[9:0];
  wire        hi = A[11];              // leading one at bit 21 rather than 20

  // ---- exponents (signed) ----
  wire [7:0] ea_n = {3'd0, ea_eff} - {4'd0, la};
  wire [7:0] eb_n = {3'd0, eb_eff} - {4'd0, lb};
  wire [7:0] S    = ea_n + eb_n;
  wire [7:0] Sx   = S + {7'd0, hi};

  // result is subnormal (or zero) iff value < 2^-14, i.e. signed Sx <= 15
  wire under = Sx[7] | ~|Sx[7:4];

  // ---- normal result path: shg = prod >> (9+hi), carry-free decomposition ----
  wire [10:0] mant_n  = hi ? A[11:1] : A[10:0];
  wire        guard_n = hi ? A[0]    : lo[9];
  wire        sticky_n = hi ? (|lo) : (|lo[8:0]);

  // ---- subnormal result path ----
  // subnormal grid value: y*2^24 = prod >> (26-S); with prod[j>=10] = A[j-10]:
  //   mant[j] = A[j+16-S] (0 when out of range), guard = A[15-S],
  //   sticky = |lo | |A[15-S:0]
  wire [5:0]  r      = 6'd15 - S[5:0];        // r = 15-S: sub_sh[0] = guard
  wire [11:0] sub_sh = A >> r;                // >> >= 12 gives 0 automatically
  wire [10:0] mant_s  = sub_sh[11:1];
  wire        guard_s = sub_sh[0];

  wire [5:0]  bits = 6'd15 - S[5:0];          // number of low A bits below the grid
  wire [11:0] mask = (12'h1 << bits) - 12'h1; // saturates to all-ones for bits > 12
  wire        sticky_s = (|lo) | (|(A & mask));

  // ---- round to nearest, ties to even ----
  wire [10:0] mant_pre = under ? mant_s : mant_n;
  wire        guard    = under ? guard_s : guard_n;
  wire        sticky   = under ? sticky_s : sticky_n;

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

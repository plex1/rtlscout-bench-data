// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits)
// with full subnormal support, round-to-nearest ties-to-even.
//
// Structure:
//  * each significand is pre-normalised (leading one onto bit 10) with a
//    balanced-tree leading-one detector + small shift; its exponent is
//    adjusted, so the product's leading one is always at bit 21 or 20.
//  * with ma_n = 2^10+fa, mb_n = 2^10+fb:
//      prod = 2^20 + (fa+fb)<<10 + fa*fb = (A << 10) + lo
//    where A = fa + fb + 2^10 + fx[19:10] (one 12-bit add) and lo = fx[9:0].
//    Since lo < 2^10, prod >> k is carry-free for every shift used:
//      normal result:      pure muxes on A / lo
//      subnormal result:   a 12-bit shift of A only (4-bit control)
//  * sticky comes from a constant mask shifted by the same 4-bit amount.
//  * exponent arithmetic is done early (T = ea+eb-la-lb-15); only a tiny
//    increment by hi/carry remains at the end.

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

  // ---- significands / effective exponents ----
  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  // ---- leading-one index, balanced tree (11 when m == 0) ----
  // operand a
  // ---- leading-one index: 4-bit mux tree (11 when m == 0) ----
  // {10:8} -> idx 8..10, {7:4} -> idx 4..7, {3:0} -> idx 0..3
  wire g8a = |m_a[10:8], g4a = |m_a[7:4], g0a = |m_a[3:0];
  wire [3:0] idxa = g8a ? {2'b10, m_a[10] | m_a[9], ~m_a[10] & m_a[9]} :
                     g4a ? {2'b01, m_a[7] | m_a[6], m_a[7] | m_a[5]} :
                     g0a ? {2'b00, m_a[3] | m_a[2], m_a[3] | m_a[1]} : 4'd11;
  // ---- leading-one index: 4-bit mux tree (11 when m == 0) ----
  // {10:8} -> idx 8..10, {7:4} -> idx 4..7, {3:0} -> idx 0..3
  wire g8b = |m_b[10:8], g4b = |m_b[7:4], g0b = |m_b[3:0];
  wire [3:0] idxb = g8b ? {2'b10, m_b[10] | m_b[9], ~m_b[10] & m_b[9]} :
                     g4b ? {2'b01, m_b[7] | m_b[6], m_b[7] | m_b[5]} :
                     g0b ? {2'b00, m_b[3] | m_b[2], m_b[3] | m_b[1]} : 4'd11;

  // normalise: shift left by (10 - idx) == clz (garbage for m == 0, overridden)
  wire [3:0]  la   = 4'd10 - idxa;
  wire [3:0]  lb   = 4'd10 - idxb;
  wire [10:0] ma_n = m_a << la;
  wire [10:0] mb_n = m_b << lb;
  wire [9:0]  fa   = ma_n[9:0];
  wire [9:0]  fb   = mb_n[9:0];

  // ---- product decomposition: prod = (A << 10) + lo ----
  wire [10:0] fs = fa + fb;
  wire [19:0] fx = fa * fb;
  wire [11:0] A  = {1'b0, fs} + 12'd1024 + {2'b0, fx[19:10]};
  wire [9:0]  lo = fx[9:0];
  wire        hi = A[11];          // leading one at bit 21 rather than 20

  // ---- exponents ----
  wire [5:0] p  = ea_eff + eb_eff;              // 2..60
  wire [4:0] q  = la + lb;                      // 0..24 (garbage if an input is 0)
  wire [6:0] T  = {1'b0, p} - {2'b0, q};        // signed
  wire [6:0] Tu = T + {6'd0, hi};
  wire under = Tu[6] | (~Tu[5] & ~Tu[4]);       // Sx <= 15  ->  subnormal result
  wire [6:0] Tm = T - 7'd15;                    // exponent bias correction, early

  // ---- normal result path (always computed) ----
  // leading one at bit 20+hi: mant = prod[20+hi:10+hi], guard = prod[9+hi]
  wire [10:0] mant_n   = hi ? A[11:1] : A[10:0];
  wire        guard_n  = hi ? A[0]    : lo[9];
  wire        sticky_n = hi ? |lo      : |lo[8:0];

  // ---- subnormal result path ----
  // frac = prod >> (26-T); mant[j] = A[j+16-T], guard = A[15-T],
  // sticky = |lo | |A[14-T:0].  c = 15-T fits 4 bits for every T in range.
  wire [3:0]  c      = T[6] ? 4'hF : (4'd15 - T[3:0]);
  wire [11:0] sub_sh = A >> c;                  // c >= 12 gives 0 automatically
  wire [10:0] mant_s  = sub_sh[11:1];
  wire        guard_s = sub_sh[0];
  wire [11:0] mask    = (12'h1 << c) - 12'd1;   // saturates to all-ones
  wire        sticky_s = |lo | (|(A & mask));

  // ---- select path ----
  wire [10:0] mant_pre = under ? mant_s  : mant_n;
  wire        guard    = under ? guard_s : guard_n;
  wire        sticky   = under ? sticky_s : sticky_n;

  // ---- round to nearest, ties to even ----
  wire        round_up = guard & (sticky | mant_pre[0]);
  wire [9:0]  frac_r   = mant_pre[9:0] + {9'd0, round_up};
  wire        ov9      = round_up & &mant_pre[9:0];
  wire        carry_n  = round_up & &mant_pre[10:0];   // mantissa overflow (normal)
  wire        carry_s  = mant_pre[10] ^ ov9;           // mantissa overflow (subnormal)

  wire [6:0] expn7 = Tm + {6'd0, hi} + {6'd0, carry_n};  wire       ovf   = ~under & (expn7[5] | (&expn7[4:0]));

  wire       kill  = under & carry_s;   // rounded up to smallest normal
  wire [9:0] frac  = kill ? 10'd0 : frac_r;

  // ---- pack ----
  assign y = nan_in         ? 16'h7E00 :
             (inf_in | ovf) ? {s_y, 5'd31, 10'd0} :
             zero_in        ? {s_y, 15'd0} :
             under          ? {s_y, 4'd0, carry_s, frac} :
                              {s_y, expn7[4:0], frac};
endmodule

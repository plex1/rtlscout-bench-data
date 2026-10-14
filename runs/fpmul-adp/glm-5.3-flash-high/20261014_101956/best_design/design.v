// fp_mul_e5f10 -- IEEE binary16 (E5M10) multiplier with full subnormal support.
//
// Architecture:
//   * Subnormal fractions are pre-normalised: the 10-bit fraction is shifted
//     left so its leading one sits at bit 10, with the exponent adjusted by
//     (-nlz).  Every non-zero significand therefore has its leading one in
//     bit 10, so the 22-bit product always has its leading one in bit 21 or
//     bit 20.  The normal path needs no leading-zero counter and only a 2:1
//     shift, with guard/sticky taken from two fixed bit positions.
//   * A subnormal result is produced by right-shifting the raw product by
//     t = 26 - (ea_adj + eb_adj) and rounding once (round-to-nearest-even),
//     so no double rounding occurs.
//   value = (m_a * m_b) * 2^(ea_adj + eb_adj - 50)

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ------------------------------------------------------------- unpack
  wire        s_a = a[15];
  wire        s_b = b[15];
  wire [4:0]  e_a = a[14:10];
  wire [4:0]  e_b = b[14:10];
  wire [9:0]  f_a = a[9:0];
  wire [9:0]  f_b = b[9:0];

  // ------------------------------------------------------------- classify
  wire e_a_zero = ~|e_a;
  wire e_b_zero = ~|e_b;
  wire e_a_max  =  &e_a;
  wire e_b_max  =  &e_b;
  wire f_a_zero = ~|f_a;
  wire f_b_zero = ~|f_b;

  wire zero_a = e_a_zero &  f_a_zero;
  wire zero_b = e_b_zero &  f_b_zero;
  wire nan_a  = e_a_max  & ~f_a_zero;
  wire nan_b  = e_b_max  & ~f_b_zero;
  wire inf_a  = e_a_max  &  f_a_zero;
  wire inf_b  = e_b_max  &  f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // ------------------------------------- normalise subnormal fractions
  // leading-zero count within the 10-bit fraction field
  function [3:0] nlz10;
    input [9:0] f;
    integer i;
    reg [3:0] c;
    begin
      c = 4'd10;
      for (i = 0; i <= 9; i = i + 1)
        if (f[i]) c = 4'd9 - i[3:0];
      nlz10 = c;
    end
  endfunction

  wire [3:0] nlz_a = e_a_zero ? nlz10(f_a) : 4'd0;
  wire [3:0] nlz_b = e_b_zero ? nlz10(f_b) : 4'd0;

  // 11-bit significands, leading one in bit 10; value = m * 2^(e_adj - 25)
  wire [10:0] m_a = e_a_zero ? ({1'b0, f_a} << (nlz_a + 4'd1)) : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? ({1'b0, f_b} << (nlz_b + 4'd1)) : {1'b1, f_b};

  wire [5:0] ea_adj = e_a_zero ? (6'd0 - {2'b0, nlz_a}) : {1'b0, e_a};
  wire [5:0] eb_adj = e_b_zero ? (6'd0 - {2'b0, nlz_b}) : {1'b0, e_b};

  wire signed [7:0] s_adj = $signed({{2{ea_adj[5]}}, ea_adj})
                          + $signed({{2{eb_adj[5]}}, eb_adj});

  // ------------------------------------------------- significand product
  wire [21:0] prod = m_a * m_b;
  wire        p21  = prod[21];

  // prefix-OR: below[k] = |prod[k:0]  (5 OR2 levels)
  wire [21:0] or1  = prod | {prod[20:0], 1'b0};
  wire [21:0] or2  = or1  | {or1[19:0],  2'b0};
  wire [21:0] or4  = or2  | {or2[17:0],  4'b0};
  wire [21:0] or8  = or4  | {or4[13:0],  8'b0};
  wire [21:0] below = or8 | {or8[5:0],  16'b0};

  // ---------------------------------------------- early shift machinery
  // subnormal result iff e_pre = s_adj - 15 + p21 <= 0
  wire        sub_res = p21 ? ($signed(s_adj) <= 8'sd14)
                            : ($signed(s_adj) <= 8'sd15);
  wire signed [7:0] t_full = 8'sd26 - $signed(s_adj);   // >= 11 if sub_res
  wire        sub23  = (t_full > 8'sd22);   // t >= 23: product rounds to zero
  wire [5:0]  t_sel  = sub23 ? 6'd22 : t_full[5:0];

  // ------------------------- normal path (result >= 2^-14), 2:1 shift
  wire [10:0] mant_pre = p21 ? prod[21:11] : prod[20:10];
  wire        guard_n  = p21 ? prod[10]    : prod[9];
  wire        sticky_n = p21 ? below[9]    : below[8];
  wire        round_n  = guard_n & (sticky_n | mant_pre[0]);
  wire [11:0] mant_rnd = {1'b0, mant_pre} + {11'd0, round_n};
  wire        carry_n  = mant_rnd[11];
  wire [9:0]  frac_nrm = carry_n ? mant_rnd[10:1] : mant_rnd[9:0];

  wire signed [7:0] e_pre = $signed(s_adj) - 8'sd15 + $signed({7'b0, p21});
  wire signed [7:0] e_nrm = e_pre + $signed({7'b0, carry_n});
  wire        overflow = e_nrm > 8'sd30;

  // -------------------- subnormal path: fraction = RNE(prod * 2^(s_adj-26))
  // (prod << (22 - t))[31:22] = prod[t+9 : t]
  wire [5:0]  sh_amt  = 6'd22 - t_sel;                  // 0..11 when sub_res
  wire [32:0] pshift  = {11'b0, prod} << sh_amt;
  wire [9:0]  frac_pre = sub23 ? 10'd0 : pshift[31:22];
  wire        guard_s  = prod[t_sel - 6'd1];            // = prod[t-1]
  wire        sticky_s = below[t_sel - 6'd2];           // = |prod[t-2:0]
  wire        round_s  = guard_s & ~sub23 & (sticky_s | frac_pre[0]);
  wire [10:0] frac_sum = {1'b0, frac_pre} + {10'd0, round_s};
  wire        carry_s  = frac_sum[10];
  wire [9:0]  frac_sub = carry_s ? 10'd0 : frac_sum[9:0];
  wire [4:0]  exp_sub  = {4'd0, carry_s};

  // ------------------------------------------------------------- pack
  wire sub_zero = sub_res & ~carry_s & (frac_sum[9:0] == 10'd0);  // incl. t >= 23
  wire is_nan   = nan_in;
  wire is_inf   = inf_in | (overflow & ~zero_in);
  wire is_zero  = zero_in | sub_zero;

  wire [4:0] exp_y = is_nan  ? 5'd31 :
                     is_inf  ? 5'd31 :
                     is_zero ? 5'd0  :
                     sub_res ? exp_sub :
                               e_nrm[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 :
                      sub_res ? frac_sub :
                                frac_nrm;

  assign y = {(is_nan ? 1'b0 : s_y), exp_y, frac_y};
endmodule

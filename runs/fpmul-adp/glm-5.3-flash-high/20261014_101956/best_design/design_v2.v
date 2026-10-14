// fp_mul_e5f10 -- IEEE binary16 (E5M10) multiplier with full subnormal support.
//
// Architecture:
//   * Subnormal input fractions are pre-normalised (balanced-tree leading-zero
//     encoder + small barrel shift) so every non-zero significand has its
//     leading one in bit 10 of an 11-bit significand m; the exponent is
//     adjusted by (-nlz).  value(m) = m * 2^(e_adj - 25).
//   * With that, the 22-bit product always has its leading one in bit 21 or 20.
//     A single alignment position
//         S = clamp(max(26 - s_adj, 10 + p21), 10, 23)
//     (s_adj = ea_adj + eb_adj, p21 = prod[21]) serves both result classes:
//       - normal results  (S = 10 + p21):  mantissa = {1'b1, F}
//       - subnormal results (S = 26 - s_adj >= 12): mantissa = {1'b0, F}
//     where F = (prod >> S)[9:0].  One round-to-nearest-even step (guard =
//     prod[S-1], sticky = |prod[S-2:0]) is applied -- no double rounding.
//   * F, guard and sticky come from a one-hot AND-OR plane over S in [10,23];
//     the one-hot codes are constant comparisons of s_adj, so they are ready
//     long before the multiplier finishes.
//   * exp = clamp(e_pre, 0) + carry with e_pre = s_adj - 15 + p21; this yields
//     0 for subnormal/zero results and 1 when a subnormal rounds up to the
//     smallest normal.  Overflow (exp > 30) gives inf.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ------------------------------------------------------------- unpack
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];
  wire       s_y = a[15] ^ b[15];

  // ------------------------------------------------------------- classify
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

  wire is_nan  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire is_inf  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  // ------------------------------------- normalise subnormal fractions
  // nlz5: number of leading zeros of a 5-bit vector (0..5), balanced logic
  function [2:0] nlz5;
    input [4:0] x;
    begin
      nlz5[2] = ~x[4] & ~x[3] & ~x[2] & ~x[1];
      nlz5[1] = ~x[4] & ~x[3] & (x[2] | x[1]);
      nlz5[0] = (~x[4] & (x[3] | (~x[2] & x[1]))) | ~|x;
    end
  endfunction

  // nlz10: leading zeros of a 10-bit vector (0..10), two-level structure
  function [3:0] nlz10;
    input [9:0] f;
    begin
      nlz10 = (|f[9:5]) ? {1'b0, nlz5(f[9:5])}
                        : (4'd5 + {1'b0, nlz5(f[4:0])});
    end
  endfunction

  wire [3:0] nlz_a = e_a_zero ? nlz10(f_a) : 4'd0;
  wire [3:0] nlz_b = e_b_zero ? nlz10(f_b) : 4'd0;

  // 11-bit significands with the leading one in bit 10
  wire [10:0] m_a = e_a_zero ? ({1'b0, f_a} << (nlz_a + 4'd1)) : {1'b1, f_a};
  wire [10:0] m_b = e_b_zero ? ({1'b0, f_b} << (nlz_b + 4'd1)) : {1'b1, f_b};

  wire [5:0] ea_adj = e_a_zero ? (6'd0 - {2'b0, nlz_a}) : {1'b0, e_a};
  wire [5:0] eb_adj = e_b_zero ? (6'd0 - {2'b0, nlz_b}) : {1'b0, e_b};

  wire signed [6:0] s_adj = {{1{ea_adj[5]}}, ea_adj}
                          + {{1{eb_adj[5]}}, eb_adj};

  // ------------------------------------------------- significand product
  wire [21:0] prod = m_a * m_b;
  wire        p21  = prod[21];

  // prefix-OR: below[k] = |prod[k:0]  (5 OR2 levels)
  wire [21:0] or1   = prod | {prod[20:0], 1'b0};
  wire [21:0] or2   = or1  | {or1[19:0],  2'b0};
  wire [21:0] or4   = or2  | {or2[17:0],  4'b0};
  wire [21:0] or8   = or4  | {or4[13:0],  8'b0};
  wire [21:0] below = or8  | {or8[5:0],  16'b0};

  // --------------------------------------- alignment position, one-hot
  // S = clamp(max(26 - s_adj, 10 + p21), 10, 23)
  wire oh10 = (s_adj >= 7'sd16) & ~p21;                       // S = 10
  wire oh11 = (s_adj == 7'sd15) | (p21 & (s_adj >= 7'sd16));  // S = 11
  wire oh23 = (s_adj <= 7'sd3);                               // S = 23

  wire [23:10] oh;
  assign oh[10] = oh10;
  assign oh[11] = oh11;
  assign oh[23] = oh23;

  genvar gk;
  generate
    for (gk = 12; gk <= 22; gk = gk + 1) begin : g_oh
      assign oh[gk] = (s_adj == (26 - gk));
    end
  endgenerate

  // implicit leading one (normal result) and exponent clamp flag
  wire imp = oh[10] | (oh[11] & p21);

  // ------------------------------------ one-hot AND-OR extraction plane
  wire [9:0] F;
  wire       guard, sticky;

  genvar gj;
  generate
    for (gj = 0; gj < 10; gj = gj + 1) begin : g_f
      wire [23:10] terms;
      for (gk = 10; gk <= 23; gk = gk + 1) begin : g_t
        if (gk + gj <= 21)
          assign terms[gk] = oh[gk] & prod[gk + gj];   // F[j] = prod[S+j]
        else
          assign terms[gk] = 1'b0;
      end
      assign F[gj] = |terms;
    end

    wire [23:10] tg, ts;
    for (gk = 10; gk <= 23; gk = gk + 1) begin : g_gs
      if (gk <= 22)
        assign tg[gk] = oh[gk] & prod[gk - 1];         // guard = prod[S-1]
      else
        assign tg[gk] = 1'b0;                          // S = 23: guard = 0
      assign ts[gk] = oh[gk] & below[gk - 2];          // sticky = |prod[S-2:0]
    end
    assign guard  = |tg;
    assign sticky = |ts;
  endgenerate

  // ---------------------------------------- round to nearest, ties even
  wire round = guard & (sticky | F[0]);

  // increment F by `round` via a 4-level prefix-AND (mask[j] = &F[j-1:0]):
  // F+1 toggles exactly those bits below (and including) the first zero
  wire [9:0] q1 = F  & {F[8:0],   1'b1};
  wire [9:0] q2 = q1 & {q1[7:0],  2'b11};
  wire [9:0] q4 = q2 & {q2[5:0],  4'b1111};
  wire [9:0] q8 = q4 & {q4[1:0],  8'b11111111};   // q8[j] = &F[j:0]

  // carry happens exactly when the 11-bit significand field is all ones
  wire        carry    = round & q8[9];
  wire [9:0]  frac_pre = F ^ (round ? {q8[8:0], 1'b1} : 10'b0);

  // ------------------------------------------------------- exponent path
  wire signed [6:0] ec7 = s_adj - 7'sd15 + $signed({6'b0, p21});
  wire [5:0] e_c   = imp ? ec7[5:0] : 6'd0;      // clamp(e_pre, 0)
  wire [5:0] e_c_p1 = e_c + 6'd1;
  wire [5:0] exp6   = carry ? e_c_p1 : e_c;
  wire       ovf    = (ec7 >= 7'sd31) | (carry & (e_c == 6'd30));

  // ------------------------------------------------------------- pack
  wire is_inf_y = is_inf | ovf;

  wire [4:0] exp_y = (is_nan | is_inf_y) ? 5'd31 :
                     zero_in             ? 5'd0  :
                                           exp6[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      (is_inf_y | zero_in) ? 10'd0 :
                                             frac_pre;

  assign y = {(is_nan ? 1'b0 : s_y), exp_y, frac_y};
endmodule

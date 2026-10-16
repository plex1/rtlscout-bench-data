// fp_mul_e5f10 v2 -- unified round-then-pack architecture.
// value = P * 2^(T-50), P = m_a*m_b (22 bits), T = ea_eff + eb_eff.
// Alignment position L (bit of P that becomes the guard bit):
//   Normal path   : L = 10 - lz   (lz = leading zeros of P; L in [-1..10])
//   Subnormal path: L = 25 - T    (L in [-10..23])
// U[m] = P[m+L] (zero outside): guard = U[0], R[j] = U[j+1] = P[L+1+j].
// sticky = |P[L-1:0]|.  One right/left shift, one rounding increment, then pack.
// Rounding carry promotes subnormal->smallest-normal and normal->E+1 identically.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack / classify ----
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = &e_a;
  wire e_b_ones = &e_b;
  wire f_a_zero = (f_a == 10'd0);
  wire f_b_zero = (f_b == 10'd0);

  wire zero_a = e_a_zero & f_a_zero;
  wire zero_b = e_b_zero & f_b_zero;
  wire nan_a  = e_a_ones & ~f_a_zero;
  wire nan_b  = e_b_ones & ~f_b_zero;
  wire inf_a  = e_a_ones & f_a_zero;
  wire inf_b  = e_b_ones & f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = a[15] ^ b[15];

  // ---- significands and exponent sum ----
  wire [10:0] m_a = {~e_a_zero, f_a};
  wire [10:0] m_b = {~e_b_zero, f_b};
  wire [5:0] T = {1'b0, e_a_zero ? 5'd1 : e_a} + {1'b0, e_b_zero ? 5'd1 : e_b};

  wire [21:0] P = m_a * m_b;

  // below[k] = |P[k:0]  (sticky helper)
  wire [21:0] below;
  genvar j;
  generate
    for (j = 0; j < 22; j = j + 1) begin : g_below
      assign below[j] = |P[j:0];
    end
  endgenerate

  // leading-zero count of P (lz = 0 when P == 0)
  reg [4:0] lz;
  integer i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (P[i]) lz = 5'd21 - i[4:0];
  end

  // ---- alignment position L (signed) and single alignment shift ----
  wire is_sub = ({1'b0, lz} + 6'd14 >= T);   // unrounded exponent field <= 0
  wire signed [5:0] Ls = is_sub ? ($signed(6'sd25) - $signed({1'b0, T}))
                                : ($signed(6'sd10) - $signed({2'b00, lz}));
  wire       lneg   = Ls[5];
  wire [4:0] shamt  = lneg ? (5'd0 - Ls[4:0]) : Ls[4:0];
  wire [21:0] U     = lneg ? (P << shamt) : (P >> shamt);
  wire        guard = U[0];
  wire [9:0]  R     = U[10:1];

  // sticky = |P[L-1:0]| = Bs[L+1] for L >= 0, else 0
  wire [24:0] Bs = {below[21], below, 2'b00};
  wire        sticky = ~lneg & Bs[Ls[4:0] + 5'd1];

  // ---- single rounding ----
  wire        round_up = guard & (sticky | R[0]);
  wire [10:0] Rr       = {1'b0, R} + {10'd0, round_up};
  wire        C        = Rr[10];
  wire [9:0]  Rf       = Rr[9:0];

  // ---- normal exponent field ----
  wire [5:0] E        = T - {1'b0, lz} - 6'd14 + {5'd0, C};
  wire       overflow = ~is_sub & (E >= 6'd31);

  // ---- pack ----
  assign y = nan_in              ? 16'h7E00 :
             (inf_in | overflow) ? {s_y, 5'h1F, 10'd0} :
             zero_in             ? {s_y, 15'd0} :
             is_sub              ? {s_y, 4'd0, C, Rf} :
                                   {s_y, E[4:0], Rf};
endmodule

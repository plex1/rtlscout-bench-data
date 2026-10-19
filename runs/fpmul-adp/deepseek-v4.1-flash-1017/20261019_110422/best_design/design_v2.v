// fp_mul_e5f10 -- binary16 (e5f10) multiplier, single-shift unified rounding.
//
// value(a) = m_a * 2^(ea_eff - 25)   (subnormal: hidden bit 0, ea_eff = 1)
// prod P = m_a*m_b (22b) ; value = P * 2^(S-50),  S = ea_eff + eb_eff
// Normal result:  k = msb(P)-10 ,  E = msb(P)+S-35+carry , frac from P>>k
// Subnormal  :    k = 26-S      ,  E = 0 (or 1 on round-up to smallest normal)
// Both cases share one barrel shift of P and one rounding add.

module fp_mul_e5f10 (
  input  wire [15:0] a,
  input  wire [15:0] b,
  output wire [15:0] y
);

  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];

  // ---- classification -------------------------------------------------
  wire ea0 = ~(|ea);
  wire eb0 = ~(|eb);
  wire nan_a = (ea == 5'h1f) & (|fa);
  wire nan_b = (eb == 5'h1f) & (|fb);
  wire inf_a = (ea == 5'h1f) & ~(|fa);
  wire inf_b = (eb == 5'h1f) & ~(|fb);
  wire zero_a = ea0 & ~(|fa);
  wire zero_b = eb0 & ~(|fb);
  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  // ---- significand product -------------------------------------------
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [21:0] P  = ma * mb;

  // leading-zero count of P (0 when P == 0)
  reg [4:0] lz;
  integer i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (P[i]) lz = 5'd21 - i[4:0];
  end
  wire [4:0] msb = 5'd21 - lz;

  wire [4:0] eae = ea0 ? 5'd1 : ea;
  wire [4:0] ebe = eb0 ? 5'd1 : eb;
  wire [5:0] S   = {1'b0, eae} + {1'b0, ebe};

  // normal iff msb + S >= 36
  wire [6:0] Nsum   = {1'b0, S} + {2'b0, msb};
  wire       is_norm = (Nsum >= 7'd36);

  wire [5:0] An = {1'b0, msb} - 6'd10;    // msb-10 (signed)
  wire [5:0] Bn = 6'd26 - {1'b0, S};      // 26-S   (signed)
  wire [4:0] k  = is_norm ? An[4:0] : Bn[4:0];

  // ---- one barrel shift: Q = P >> (k-1); guard = Q[0], M = Q[11:1] ---
  wire [24:0] PX = {3'b0, P};
  wire [4:0]  d  = (k == 5'd0) ? 5'd0 : (k - 5'd1);
  wire [24:0] Q  = PX >> d;

  // sticky: any bit of P below the guard
  wire [24:0] msk = (25'd1 << d) - 1;
  wire sticky = |(PX & msk);

  wire [10:0] M      = (k == 5'd0) ? P[10:0] : Q[11:1];
  wire        guard  = (k != 5'd0) & Q[0];
  wire        rnd    = guard & (sticky | M[0]);

  wire [11:0] Msum  = {1'b0, M} + {11'd0, rnd};
  wire        carry = Msum[11];

  // ---- exponent ------------------------------------------------------
  wire signed [6:0] En = {1'b0, S} - {2'b0, lz} - 7'sd14 + {6'd0, carry};
  wire ovf = (En > 7'sd30);

  wire [4:0] e_fin = is_norm ? ((zero_in) ? 5'd0 : En[4:0]) : {4'd0, Msum[10]};
  wire [9:0] f_fin = Msum[9:0];

  wire nan_out = nan_in;
  wire inf_out = ~nan_out & (inf_in | (ovf & ~zero_in));

  wire [4:0] e_y = (nan_out | inf_out) ? 5'h1f : e_fin;
  wire [9:0] f_y = nan_out ? 10'h200 : (inf_out ? 10'd0 : f_fin);
  wire       s_y = nan_out ? 1'b0 : (sa ^ sb);

  assign y = {s_y, e_y, f_y};

endmodule

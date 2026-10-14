// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
//
// Delay-oriented structure:
//   - one-hot leading-one positions of the significands (prefix-OR), combined
//     into a one-hot of ij = i_a + i_b; the normalisation window is a one-hot
//     AND-OR funnel V = X >> ij, giving p1 = V[12] (product leading one).
//   - single shared rounding increment: the unrounded fraction fields
//     (normal window vs subnormal window) are muxed first, then rounded once
//     with the muxed round condition (nearest, ties to even).
//   - exponent: H0 = ea_eff + eb_eff + ij, compares in parallel, three
//     candidates (+29/+30/+31) selected by s = p1 + carry.
//   - subnormal results: one-hot funnel prod >> ts, ts decoded early.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- 1. unpack / classify ----
  wire       s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0],  f_b = b[9:0];

  wire e_a_zero = ~|e_a, e_b_zero = ~|e_b;
  wire e_a_ones =  &e_a, e_b_ones =  &e_b;
  wire f_a_zero = ~|f_a, f_b_zero = ~|f_b;
  wire zero_a = e_a_zero & f_a_zero, zero_b = e_b_zero & f_b_zero;
  wire nan_a  = e_a_ones & ~f_a_zero, nan_b = e_b_ones & ~f_b_zero;
  wire inf_a  = e_a_ones &  f_a_zero, inf_b = e_b_ones &  f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // ---- 2. significand product ----
  wire h_a = ~e_a_zero, h_b = ~e_b_zero;
  wire [10:0] m_a = {h_a, f_a}, m_b = {h_b, f_b};
  wire [21:0] prod = m_a * m_b;

  // ---- 3. one-hot leading-one position of each significand ----
  genvar gi, gj, gp, gk, gq;
  wire [10:0] oh_a, oh_b;
  generate
    for (gi = 0; gi <= 10; gi = gi + 1) begin : g_oha
      if (gi == 10) assign oh_a[gi] = m_a[10];
      else          assign oh_a[gi] = m_a[gi] & ~|m_a[10:gi+1];
    end
    for (gj = 0; gj <= 10; gj = gj + 1) begin : g_ohb
      if (gj == 10) assign oh_b[gj] = m_b[10];
      else          assign oh_b[gj] = m_b[gj] & ~|m_b[10:gj+1];
    end
  endgenerate

  wire [120:0] pair;
  generate
    for (gi = 0; gi <= 10; gi = gi + 1)
      for (gj = 0; gj <= 10; gj = gj + 1)
        assign pair[gi*11+gj] = oh_a[gi] & oh_b[gj];
  endgenerate

  wire [20:0] ohij;
  generate
    for (gp = 0; gp <= 20; gp = gp + 1) begin : g_ohij
      wire [10:0] o;
      for (gi = 0; gi <= 10; gi = gi + 1) begin : g_i
        if (gp >= gi && gp - gi <= 10) assign o[gi] = pair[gi*11 + (gp-gi)];
        else                           assign o[gi] = 1'b0;
      end
      assign ohij[gp] = |o;
    end
  endgenerate

  wire [4:0] ij;
  generate
    for (gk = 0; gk < 5; gk = gk + 1) begin : g_ij
      wire [20:0] t;
      for (gp = 0; gp <= 20; gp = gp + 1)
        assign t[gp] = ohij[gp] & ((gp >> gk) & 1);
      assign ij[gk] = |t;
    end
  endgenerate

  // ---- 4. normalisation window: V[k] = prod[ij-11+k], k = 0..12 ----
  wire [32:0] X = {prod, 11'b0};
  wire [12:0] V;
  generate
    for (gk = 0; gk <= 12; gk = gk + 1) begin : g_V
      wire [20:0] t;
      for (gp = 0; gp <= 20; gp = gp + 1)
        assign t[gp] = ohij[gp] & X[gp + gk];
      assign V[gk] = |t;
    end
  endgenerate

  wire p1 = V[12];

  wire [21:0] below;
  generate
    for (gk = 0; gk < 22; gk = gk + 1) begin : g_below
      assign below[gk] = |prod[gk:0];
    end
  endgenerate

  wire sticky0, sticky1;
  generate
    begin : g_stick
      wire [20:0] t0, t1;
      for (gp = 0; gp <= 20; gp = gp + 1) begin : g_p
        assign t0[gp] = (gp >= 12) ? ohij[gp] & below[gp-12] : 1'b0;
        assign t1[gp] = (gp >= 11) ? ohij[gp] & below[gp-11] : 1'b0;
      end
      assign sticky0 = |t0;
      assign sticky1 = |t1;
    end
  endgenerate

  wire [10:0] mant_pre;
  generate
    for (gk = 0; gk < 11; gk = gk + 1) begin : g_mant
      assign mant_pre[gk] = p1 ? V[gk+2] : V[gk+1];
    end
  endgenerate
  wire guard  = p1 ? V[1] : V[0];
  wire sticky = p1 ? sticky1 : sticky0;

  // ---- 5. rounding conditions ----
  wire r = guard & (sticky | mant_pre[0]);

  // rounding carry into the exponent: all-ones mantissa => mant[0]=1 => r=guard
  wire all1_hi = &V[12:2];
  wire all1_lo = &V[11:1];
  wire carry   = guard & (p1 ? all1_hi : all1_lo);

  // ---- 6. exponent ----
  wire [4:0] ea_eff = e_a | {4'd0, e_a_zero};
  wire [4:0] eb_eff = e_b | {4'd0, e_b_zero};
  wire [6:0] H0 = ea_eff + eb_eff + {2'b0, ij};   // e_norm = H0 + s - 35, s = p1+carry

  wire [4:0] EN0 = H0[4:0] + 5'd29;               // (H0 - 35) mod 32
  wire [4:0] EN1 = H0[4:0] + 5'd30;
  wire [4:0] EN2 = H0[4:0] + 5'd31;

  wire u0 = (H0 <= 7'd35), u1 = (H0 <= 7'd34), u2 = (H0 <= 7'd33);
  wire o0 = (H0 >= 7'd66), o1 = (H0 >= 7'd65), o2 = (H0 >= 7'd64);

  wire s1 = p1 & carry, s0 = p1 ^ carry;
  wire unf = s1 ? u2 : s0 ? u1 : u0;
  wire ovf = s1 ? o2 : s0 ? o1 : o0;
  wire [4:0] EN_sel = s1 ? EN2 : s0 ? EN1 : EN0;

  // ---- 7. subnormal-result path: frac = prod >> (26 - exp_sum), rounded once ----
  // ts = 26 - exp_sum in [-34..24]; negative ts -> all-zero (out-of-range shift)
  wire [6:0] ts7 = 7'd26 - {1'b0, ea_eff} - {1'b0, eb_eff};
  wire       ts_valid = ~ts7[6];
  wire [4:0] ts1   = ts7[4:0];                   // meaningful only when ts_valid

  // one-hot decode of ts, then guard = prod[ts-1], frac[k] = prod[ts+k], sticky = |prod[ts-2:0]|
  wire [24:0] oh_ts;
  wire [9:0]  frac_trunc;
  wire        guard_s, sticky_s;
  generate
    for (gq = 0; gq <= 24; gq = gq + 1) begin : g_ts
      assign oh_ts[gq] = ts_valid & (ts1 == gq[4:0]);
    end
    for (gk = 0; gk < 10; gk = gk + 1) begin : g_F
      wire [24:0] t;
      for (gq = 0; gq <= 24; gq = gq + 1) begin : g_q
        wire [7:0] n = gq + gk;
        assign t[gq] = oh_ts[gq] & ((n < 22) ? prod[n[4:0]] : 1'b0);
      end
      assign frac_trunc[gk] = |t;
    end
    begin : g_gs
      wire [24:0] tg, ts;
      for (gq = 0; gq <= 24; gq = gq + 1) begin : g_q
        wire [7:0] n = gq - 1;
        wire [7:0] m = gq - 2;
        assign tg[gq] = oh_ts[gq] & ((gq >= 1) ? ((n < 22) ? prod[n[4:0]] : 1'b0) : 1'b0);
        assign ts[gq] = oh_ts[gq] & ((gq >= 2) ? below[m[4:0]] : 1'b0);
      end
      assign guard_s  = |tg;
      assign sticky_s = |ts;
    end
  endgenerate

  wire round_s = guard_s & (sticky_s | frac_trunc[0]);
  wire carry_s = round_s & &frac_trunc;

  // ---- 8. shared round + flat pack ----
  wire FINF = nan_in | inf_in | ovf;
  wire FR   = ~FINF;

  wire [9:0] fmux = unf ? frac_trunc : mant_pre[9:0];
  wire rsel  = unf ? round_s : r;

  wire [9:0] cinc;
  generate
    for (gk = 0; gk < 10; gk = gk + 1) begin : g_cinc
      if (gk == 0) assign cinc[gk] = rsel;
      else         assign cinc[gk] = rsel & &fmux[gk-1:0];
    end
  endgenerate
  wire [9:0] frac_y;
  generate
    for (gk = 0; gk < 9; gk = gk + 1) begin : g_fy
      assign frac_y[gk] = FR & (fmux[gk] ^ cinc[gk]);
    end
  endgenerate
  assign frac_y[9] = nan_in | (FR & (fmux[9] ^ cinc[9]));

  wire [4:0] exp_y;
  generate
    for (gk = 1; gk < 5; gk = gk + 1) begin : g_ey
      assign exp_y[gk] = FINF | (~unf & EN_sel[gk]);
    end
  endgenerate
  assign exp_y[0] = FINF | (unf & carry_s) | (~unf & EN_sel[0]);

  wire sign_y = ~nan_in & s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

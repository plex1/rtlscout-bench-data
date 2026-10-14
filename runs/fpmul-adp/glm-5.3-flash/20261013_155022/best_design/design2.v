// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
//
// Delay-oriented structure:
//   - leading-one positions of the significands as one-hot vectors (prefix-OR),
//     combined into a one-hot of ij = i_a + i_b; the normalisation window is a
//     one-hot AND-OR funnel V = X >> ij, giving p1 = V[12] (product leading one)
//     without a separate 22:1 mux.
//   - mantissa = 2:1 mux of two adjacent 11-bit windows keyed by p1; single
//     rounding increment with parallel-prefix AND trees.
//   - exponent: H = exp_sum + ij computed early; underflow/overflow and the
//     exponent field are three early candidates selected by s = p1 + carry.
//   - subnormal results: separate funnel prod >> ts, rounded once.
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
  genvar gi, gj, gp, gk;
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

  // one-hot of ij = i + j  (0..20)
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

  // binary ij (for the exponent sum)
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

  wire p1 = V[12];   // product's leading one is one above ij

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

  // ---- 5. round (nearest, ties to even) ----
  wire r = guard & (sticky | mant_pre[0]);

  wire [10:0] cinc;
  generate
    for (gk = 0; gk < 11; gk = gk + 1) begin : g_cinc
      if (gk == 0) assign cinc[gk] = r;
      else         assign cinc[gk] = r & &mant_pre[gk-1:0];
    end
  endgenerate
  wire [10:0] mant_post = mant_pre ^ cinc;
  wire carry  = r & &mant_pre[10:0];

  // ---- 6. exponent ----
  wire [4:0] ea_eff = e_a | {4'd0, e_a_zero};
  wire [4:0] eb_eff = e_b | {4'd0, e_b_zero};
  wire [5:0] exp_sum = ea_eff + eb_eff;
  wire [6:0] H = {1'b0, exp_sum} + {2'b0, ij};   // e_norm = H + s - 35, s = p1+carry

  wire [4:0] EN0 = H[4:0] + 5'd29;               // (H - 35) mod 32
  wire [4:0] EN1 = EN0 + 5'd1;
  wire [4:0] EN2 = EN0 + 5'd2;

  wire u0 = (H <= 7'd35), u1 = (H <= 7'd34), u2 = (H <= 7'd33);
  wire o0 = (H >= 7'd66), o1 = (H >= 7'd65), o2 = (H >= 7'd64);

  wire [1:0] s = {1'b0, p1} + carry;
  wire unf = s[1] ? u2 : s[0] ? u1 : u0;
  wire ovf = s[1] ? o2 : s[0] ? o1 : o0;
  wire [4:0] EN_sel = s[1] ? EN2 : s[0] ? EN1 : EN0;

  // ---- 7. subnormal-result path ----
  wire [7:0] ts8   = 8'd26 - {2'b0, exp_sum};
  wire [7:0] ts1_8 = ts8 - 8'd1;
  wire [4:0] ts1   = ts1_8[7] ? 5'd0 : ts1_8[4:0];

  wire [21:0] v = prod >> ts1;
  wire [9:0]  frac_trunc = v[10:1];
  wire        guard_s    = v[0];
  wire [4:0]  ss_idx     = ts1 - 5'd1;
  wire        sticky_s   = (ts1 > 5'd0 && ts1 < 5'd23) ? below[ss_idx] : 1'b0;

  wire round_s = guard_s & (sticky_s | frac_trunc[0]);
  wire carry_s = round_s & &frac_trunc;

  wire [9:0] csub;
  generate
    for (gk = 0; gk < 10; gk = gk + 1) begin : g_csub
      if (gk == 0) assign csub[gk] = round_s;
      else         assign csub[gk] = round_s & &frac_trunc[gk-1:0];
    end
  endgenerate
  wire [9:0] frac_sub = frac_trunc ^ csub;
  wire       sub_zero = ~carry_s & ~|frac_sub;

  // ---- 8. pack ----
  wire c_nan  = nan_in;
  wire c_inf  = ~c_nan & (inf_in | ovf);
  wire c_zero = ~c_nan & ~c_inf & (zero_in | (unf & sub_zero));
  wire c_sub  = ~c_nan & ~c_inf & ~c_zero & unf;

  wire [4:0] exp_y = (c_nan | c_inf) ? 5'd31 :
                      c_zero         ? 5'd0  :
                      c_sub          ? {4'd0, carry_s} : EN_sel;
  wire [9:0] frac_y = c_nan ? 10'h200 :
                      (c_inf | c_zero) ? 10'd0 :
                      c_sub ? frac_sub : mant_post[9:0];
  wire sign_y = ~c_nan & s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

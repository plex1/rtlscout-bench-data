// fp_mul_e5f10 -- IEEE binary16 multiplier.
// Key facts: m_a,m_b in [2^10,2^11) => product p in [2^20,2^22) => lz(p) in {0,1}.
// And the result shift is always >= 9 to the right, so only p[21:9] (13 bits)
// plus a constant OR of p[8:0] are needed.
module fp_mul_e5f10 (a, b, y);
  input  wire [15:0] a, b;
  output wire [15:0] y;

  wire s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0],   f_b = b[9:0];

  wire eaz = ~|e_a, ebz = ~|e_b;
  wire eao =  &e_a, ebo =  &e_b;
  wire faz = ~|f_a, fbz = ~|f_b;
  wire z_a = eaz & faz, z_b = ebz & fbz;
  wire n_a = eao & ~faz, n_b = ebo & ~fbz;
  wire i_a = eao & faz,  i_b = ebo & fbz;

  wire nan_in = n_a | n_b | (i_a & z_b) | (z_a & i_b);
  wire inf_in = i_a | i_b;
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a = {~eaz, f_a};
  wire [10:0] m_b = {~ebz, f_b};
  wire [5:0]  ea  = eaz ? 6'd1 : {1'b0, e_a};
  wire [5:0]  eb  = ebz ? 6'd1 : {1'b0, e_b};
  wire [5:0]  E   = ea + eb;

  wire [21:0] p = m_a * m_b;

  wire        lz   = ~p[21];                 // 0 or 1
  wire signed [6:0] ebase = $signed({1'b0, E}) - 7'sd14 - $signed({2'b0, lz});
  wire        normal = (ebase >= 7'sd1);

  // x = p >> 9 (13 bits).  Result shift sx = normal ? lz : (16-E).
  wire [12:0] x   = p[21:9];
  wire [4:0]  sx  = normal ? {4'd0, lz} : (5'd16 - E[4:0]);
  wire [12:0] shr = x >> sx;

  wire [10:0] M_pre = shr[11:1];
  wire        guard = shr[0];
  wire        lsb   = shr[1];

  // sticky = |p[8:0] (always below the guard)  |  |x[sx-1:0]
  wire low9 = |p[8:0];
  wire [12:0] xp;
  assign xp[0] = x[0];
  genvar k;
  generate for (k = 1; k < 13; k = k + 1) begin : g_xp
    assign xp[k] = xp[k-1] | x[k];
  end endgenerate
  wire [13:0] xpe = {xp[12], xp};            // clamp index 13 -> xp[12]
  wire sticky = low9 | (|sx ? xpe[sx-1] : 1'b0);

  wire rup = guard & (sticky | lsb);

  wire [11:0] Mfull = {1'b0, M_pre} + {11'd0, rup};
  wire        carry  = Mfull[11];
  wire signed [6:0] ef = ebase + carry;

  wire [4:0] field = normal ? ef[4:0] : {4'd0, Mfull[10]};
  wire [9:0] frac  = normal ? (carry ? Mfull[10:1] : Mfull[9:0]) : Mfull[9:0];

  wire is_inf   = ~nan_in & (inf_in | (ef > 7'sd30));
  wire zero_res = ~nan_in & ~is_inf & (field == 5'd0) & (frac == 10'd0);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : zero_res ? 5'd0 : field;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | zero_res) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

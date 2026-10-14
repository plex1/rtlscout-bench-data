// fp_mul_e5f10 -- IEEE binary16 multiplier (1/5/10), subnormal support.
// p = m_a*m_b in [2^20,2^22)  =>  lz(p) in {0,1}  =>  lz = ~p[21].
// Result shift S = normal ? (11-lz) : (26-E) is always >= 10, so with
// y = p >> 9  (13 bits: y[12]=0, y[11:0]=p[21:10]) the mantissa window is
// simply  sh = y >> u ,  u = normal ? lz : (16-E)  in [0,14].
// Guard = sh[0], LSB = sh[1], sticky = |p[9:1] | |y[u-1:0].
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
  wire i_a = eao &  faz, i_b = ebo &  fbz;

  wire nan_in = n_a | n_b | (i_a & z_b) | (z_a & i_b);
  wire inf_in = i_a | i_b;
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a = {~eaz, f_a};
  wire [10:0] m_b = {~ebz, f_b};
  wire [5:0]  E   = (eaz ? 6'd1 : {1'b0, e_a}) + (ebz ? 6'd1 : {1'b0, e_b});

  wire [21:0] p  = m_a * m_b;
  wire        lz = ~p[21];                     // 0 or 1
  wire [5:0]  thr = lz ? 6'd16 : 6'd15;
  wire        normal = (E >= thr);
  wire [4:0]  u  = normal ? {4'd0, lz} : (5'd16 - E[4:0]);   // 0..14

  // y = p>>9 : 13 significant bits
  wire [12:0] y0 = {1'b0, p[21:10]};

  // prefix OR of y0 so a shift yields the sticky bit of the shifted range
  wire [12:0] yp;
  assign yp[0] = y0[0];
  genvar g;
  generate for (g = 1; g < 13; g = g + 1) begin : g_yp
    assign yp[g] = yp[g-1] | y0[g];
  end endgenerate

  wire st_hi = (u == 5'd0)  ? 1'b0 :
               (u == 5'd1)  ? yp[0]  :
               (u == 5'd2)  ? yp[1]  :
               (u == 5'd3)  ? yp[2]  :
               (u == 5'd4)  ? yp[3]  :
               (u == 5'd5)  ? yp[4]  :
               (u == 5'd6)  ? yp[5]  :
               (u == 5'd7)  ? yp[6]  :
               (u == 5'd8)  ? yp[7]  :
               (u == 5'd9)  ? yp[8]  :
               (u == 5'd10) ? yp[9]  :
               (u == 5'd11) ? yp[10] :
               (u == 5'd12) ? yp[11] : yp[12];
  wire sticky = st_hi | (|p[9:1]);

  wire [12:0] sh = y0 >> u[3:0];
  wire [10:0] M_pre = sh[11:1];
  wire        guard = sh[0];
  wire        lsb   = sh[1];
  wire        rup   = guard & (sticky | lsb);

  wire [11:0] Msum  = {1'b0, M_pre} + {11'd0, rup};
  wire        carry = Msum[11];
  wire signed [6:0] ef = $signed({1'b0, E}) - 7'sd14 - $signed({2'b0, lz}) + carry;

  wire [4:0] field = normal ? ef[4:0] : {4'd0, Msum[10]};
  wire [9:0] frac  = normal ? (carry ? Msum[10:1] : Msum[9:0]) : Msum[9:0];

  wire zero_in = z_a | z_b;
  wire is_inf  = ~nan_in & (inf_in | (ef > 7'sd30));
  wire is_zres = ~nan_in & ~is_inf & (zero_in | ((field == 5'd0) & (frac == 10'd0)));

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_zres ? 5'd0 : field;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_zres) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

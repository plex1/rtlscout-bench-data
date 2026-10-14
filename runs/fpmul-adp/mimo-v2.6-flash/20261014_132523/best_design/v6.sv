// fp_mul_e5f10 -- binary16 multiplier, v6.
// v6 vs v5: s extracted early as prod[21-ZS] (small 22:1 mux right after the
// multiply) instead of u[21] after the normalisation shifter, so jv/jj/thermo
// are ready before the shifter output w arrives.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- classify ----
  wire        sa = a[15],    sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire ea0 = (ea == 5'd0),  eb0 = (eb == 5'd0);
  wire eap = (ea == 5'd31), ebp = (eb == 5'd31);

  wire zero_a = ea0 & (fa == 10'd0), zero_b = eb0 & (fb == 10'd0);
  wire inf_a  = eap & (fa == 10'd0), inf_b  = ebp & (fb == 10'd0);
  wire nan_in = (eap & (fa != 10'd0)) | (ebp & (fb != 10'd0))
              | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ---- raw significands ----
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};

  // ---- tree leading-zero counts (parallel with multiply) ----
  wire        g0nz = ma[10] | ma[9],  g1nz = ma[8] | ma[7];
  wire        g2nz = ma[6] | ma[5],   g3nz = ma[4] | ma[3];
  wire        g4nz = ma[2] | ma[1];
  wire        g0z = ~ma[10], g1z = ~ma[8], g2z = ~ma[6], g3z = ~ma[4], g4z = ~ma[2];
  wire        h0nz = g0nz | g1nz;  wire [1:0] h0z = g0nz ? {1'b0, g0z} : {1'b1, g1z};
  wire        h1nz = g2nz | g3nz;  wire [1:0] h1z = g2nz ? {1'b0, g2z} : {1'b1, g3z};
  wire        h2nz = g4nz;         wire [1:0] h2z = g4nz ? {1'b0, g4z} : 2'd2;
  wire        f0nz = h0nz | h1nz;  wire [2:0] f0z = h0nz ? {1'b0, h0z} : {1'b1, h1z};
  wire [3:0]  za = f0nz ? {1'b0, f0z} : (4'd8 + {1'b0, h2z});

  wire        k0nz = mb[10] | mb[9],  k1nz = mb[8] | mb[7];
  wire        k2nz = mb[6] | mb[5],   k3nz = mb[4] | mb[3];
  wire        k4nz = mb[2] | mb[1];
  wire        k0z = ~mb[10], k1z = ~mb[8], k2z = ~mb[6], k3z = ~mb[4], k4z = ~mb[2];
  wire        l0nz = k0nz | k1nz;  wire [1:0] l0z = k0nz ? {1'b0, k0z} : {1'b1, k1z};
  wire        l1nz = k2nz | k3nz;  wire [1:0] l1z = k2nz ? {1'b0, k2z} : {1'b1, k3z};
  wire        l2nz = k4nz;         wire [1:0] l2z = k4nz ? {1'b0, k4z} : 2'd2;
  wire        m0nz = l0nz | l1nz;  wire [2:0] m0z = l0nz ? {1'b0, l0z} : {1'b1, l1z};
  wire [3:0]  zb = m0nz ? {1'b0, m0z} : (4'd8 + {1'b0, l2z});

  wire [4:0]  ZS  = {1'b0, za} + {1'b0, zb};      // 0..22

  // ---- exponent base: jv = -E = 15 - E0 + ZS - s ----
  wire [5:0]  E0  = {1'b0, (ea0 ? 5'd1 : ea)} + {1'b0, (eb0 ? 5'd1 : eb)};
  wire [6:0]  jbase = 7'd15 - {1'b0, E0} + {2'd0, ZS};

  // ---- raw product ----
  wire [21:0] prod = ma * mb;

  // ---- early s:  s = u[21] = prod[21-ZS]  (ZS==22 -> prod==0 -> s==0) ----
  reg s_pre;
  integer ki;
  always @* begin
    s_pre = 1'b0;
    for (ki = 0; ki <= 21; ki = ki + 1)
      if (ZS == ki[4:0]) s_pre = prod[21 - ki];
  end
  wire s = s_pre;

  // ---- normalise: u = prod << ZS (lossless: prod < 2^(22-ZS)) ----
  wire [21:0] u1 = ZS[0] ? {prod[20:0], 1'b0} : prod;
  wire [21:0] u2 = ZS[1] ? {u1[19:0], 2'b0}   : u1;
  wire [21:0] u4 = ZS[2] ? {u2[17:0], 4'b0}   : u2;
  wire [21:0] u8 = ZS[3] ? {u4[13:0], 8'b0}   : u4;
  wire [21:0] u  = ZS[4] ? {u8[5:0], 16'b0}   : u8;

  wire [20:0] w    = s ? u[21:1] : u[20:0];
  wire        tail = s & u[0];

  wire [6:0]  jv   = jbase - {6'd0, s};   // jv = -E ; >=0 -> subnormal result

  // ---- normal path: round w at 11-bit precision ----
  wire [10:0] rpre   = w[20:10];
  wire        gn     = w[9];
  wire        stn    = (|w[8:0]) | tail;
  wire        rup    = gn & (stn | rpre[0]);
  wire [11:0] rres   = {1'b0, rpre} + {11'd0, rup};
  wire        cr     = rres[11];
  wire [9:0]  frac_n = rres[9:0];

  // ---- shared exponent / selection from jv ----
  wire signed [6:0] expfs = $signed({6'b0, cr}) - $signed(jv);  // E + cr
  wire sub  = ~jv[6];                       // E <= 0
  wire ovf  = (expfs >= 7'sd31);            // E + cr >= 31 (only reachable when normal)
  wire juse = (jv <= 7'd10);                // 0 <= j <= 10
  wire [3:0] jj = jv[3:0];

  // ---- subnormal path: frac = round(w / 2^(11+j)) ----
  wire [31:0] wx = {11'd0, w};
  reg  [11:0] win;
  always @* begin
    win = 12'd0;
    if (juse)
      case (jj)
        4'd0:  win = wx[10 +: 12];
        4'd1:  win = wx[11 +: 12];
        4'd2:  win = wx[12 +: 12];
        4'd3:  win = wx[13 +: 12];
        4'd4:  win = wx[14 +: 12];
        4'd5:  win = wx[15 +: 12];
        4'd6:  win = wx[16 +: 12];
        4'd7:  win = wx[17 +: 12];
        4'd8:  win = wx[18 +: 12];
        4'd9:  win = wx[19 +: 12];
        4'd10: win = wx[20 +: 12];
        default: win = 12'd0;
      endcase
  end
  wire        gs   = win[0];
  wire [9:0]  t10  = win[11:1];
  // sticky = |w[9+j:0] = |w[9:0] | (|(w[19:10] & thermo)), thermo[t] = (j >= t+1)
  wire [9:0]  th   = {jj >= 4'd10, jj >= 4'd9, jj >= 4'd8, jj >= 4'd7, jj >= 4'd6,
                      jj >= 4'd5,  jj >= 4'd4, jj >= 4'd3, jj >= 4'd2, jj >= 4'd1};
  wire        sts  = (|w[9:0]) | (|(w[19:10] & th)) | tail;
  wire        rs   = gs & (sts | t10[0]);
  wire [10:0] fres = {1'b0, t10} + {10'd0, rs};
  wire        cs   = fres[10];

  // ---- pack: priority N > F > (sub|zero) > normal ----
  wire        v_nan  = nan_in;
  wire        v_inf  = inf_in | ovf;
  wire        selA   = ~(v_nan | v_inf);       // not nan, not inf
  wire        zs     = sub | zero_in;          // subnormal/zero result path
  wire        selZS  = selA & zs;
  wire        selNRM = selA & ~zs;

  assign y = v_nan  ? 16'h7E00 :
             v_inf  ? {sy, 5'd31, 10'd0} :
             selZS  ? {sy, 4'd0, cs, fres[9:0]} :
                      {sy, expfs[4:0], frac_n};
endmodule

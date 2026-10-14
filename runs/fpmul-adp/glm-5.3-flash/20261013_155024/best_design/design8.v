// fp_mul_e5f10 -- bisection variant: explicit concatenation partial products,
// CSA tree, plain final adder.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire       sa = a[15],      sb = b[15];
  wire [4:0] ea = a[14:10],   eb = b[14:10];
  wire [9:0] fa = a[9:0],     fb = b[9:0];

  wire eaz = ~|ea,  ebz = ~|eb;
  wire eao =  &ea,  ebo =  &eb;
  wire faz = ~|fa,  fbz = ~|fb;

  wire nan_in = (eao & ~faz) | (ebo & ~fbz)
              | ((eao & faz) & (ebz & fbz))
              | ((eaz & faz) & (ebo & fbz));
  wire inf_in = (eao & faz) | (ebo & fbz);
  wire sy     = sa ^ sb;

  wire [10:0] ma = {~eaz, fa};
  wire [10:0] mb = {~ebz, fb};
  wire [5:0]  S  = (eaz ? 6'd1 : {1'b0, ea}) + (ebz ? 6'd1 : {1'b0, eb});
  wire [6:0]  Am = S - 7'd35;

  // ------------------------------------------------ partial products
  wire [10:0] m0  = {11{ma[0]}}  & mb;
  wire [10:0] m1  = {11{ma[1]}}  & mb;
  wire [10:0] m2  = {11{ma[2]}}  & mb;
  wire [10:0] m3  = {11{ma[3]}}  & mb;
  wire [10:0] m4  = {11{ma[4]}}  & mb;
  wire [10:0] m5  = {11{ma[5]}}  & mb;
  wire [10:0] m6  = {11{ma[6]}}  & mb;
  wire [10:0] m7  = {11{ma[7]}}  & mb;
  wire [10:0] m8  = {11{ma[8]}}  & mb;
  wire [10:0] m9  = {11{ma[9]}}  & mb;
  wire [10:0] m10 = {11{ma[10]}} & mb;

  wire [21:0] pp0  = {11'b0,               m0};
  wire [21:0] pp1  = {10'b0,               m1,  1'b0};
  wire [21:0] pp2  = {9'b0,                m2,  2'b0};
  wire [21:0] pp3  = {8'b0,                m3,  3'b0};
  wire [21:0] pp4  = {7'b0,                m4,  4'b0};
  wire [21:0] pp5  = {6'b0,                m5,  5'b0};
  wire [21:0] pp6  = {5'b0,                m6,  6'b0};
  wire [21:0] pp7  = {4'b0,                m7,  7'b0};
  wire [21:0] pp8  = {3'b0,                m8,  8'b0};
  wire [21:0] pp9  = {2'b0,                m9,  9'b0};
  wire [21:0] pp10 = {1'b0,                m10, 10'b0};

  // ------------------------------------------------ carry-save tree
  wire [21:0] s11 = pp0 ^ pp1 ^ pp2;
  wire [21:0] c11 = ((pp0 & pp1) | (pp0 & pp2) | (pp1 & pp2)) << 1;
  wire [21:0] s12 = pp3 ^ pp4 ^ pp5;
  wire [21:0] c12 = ((pp3 & pp4) | (pp3 & pp5) | (pp4 & pp5)) << 1;
  wire [21:0] s13 = pp6 ^ pp7 ^ pp8;
  wire [21:0] c13 = ((pp6 & pp7) | (pp6 & pp8) | (pp7 & pp8)) << 1;
  wire [21:0] s14 = pp9 ^ pp10;
  wire [21:0] c14 = (pp9 & pp10) << 1;

  wire [21:0] s21 = s11 ^ s12 ^ s13;
  wire [21:0] c21 = ((s11 & s12) | (s11 & s13) | (s12 & s13)) << 1;
  wire [21:0] s22 = c11 ^ c12 ^ c13;
  wire [21:0] c22 = ((c11 & c12) | (c11 & c13) | (c12 & c13)) << 1;
  wire [21:0] s23 = s14 ^ c14;
  wire [21:0] c23 = (s14 & c14) << 1;

  wire [21:0] s31 = s21 ^ s22 ^ s23;
  wire [21:0] c31 = ((s21 & s22) | (s21 & s23) | (s22 & s23)) << 1;
  wire [21:0] s32 = c21 ^ c22 ^ c23;
  wire [21:0] c32 = ((c21 & c22) | (c21 & c23) | (c22 & c23)) << 1;

  wire [21:0] s41 = s31 ^ c31 ^ s32;
  wire [21:0] c41 = ((s31 & c31) | (s31 & s32) | (c31 & s32)) << 1;
  wire [21:0] U   = s41 ^ c41 ^ c32;
  wire [21:0] V   = ((s41 & c41) | (s41 & c32) | (c41 & c32)) << 1;

  wire [22:0] prodf = U + V;                     // plain adder for bisection
  wire [21:0] prod  = prodf[21:0];  // ------------------------------------------------ prefix nets
  wire [21:0] nz, bot, h;
  genvar i;
  generate
    for (i = 1; i < 22; i = i + 1) begin : g_pf
      assign nz[i]  = |prod[21:i];
      assign bot[i] = |prod[i:0];
      assign h[i]   = prod[21-i] & ~nz[22-i];
    end
  endgenerate
  assign nz[0]  = |prod;
  assign bot[0] = prod[0];
  assign h[0]   = prod[21];

  // ------------------------------------------------ normal window
  wire [10:0] N;
  genvar j, k;
  generate
    for (j = 0; j < 11; j = j + 1) begin : g_win
      wire [11:0] t;
      for (k = 0; k < 12; k = k + 1) begin : g_t
        assign t[k] = h[k] & prod[j+11-k];
      end
      assign N[j] = |t;
    end
  endgenerate

  wire [10:0] gt;
  wire [9:0]  st;
  generate
    for (k = 0; k <= 10; k = k + 1) begin : g_g
      assign gt[k] = h[k] & prod[10-k];
    end
    for (k = 0; k <= 9; k = k + 1) begin : g_s
      assign st[k] = h[k] & bot[9-k];
    end
  endgenerate
  wire guard_n  = |gt;
  wire sticky_n = |st;
  wire r        = guard_n & (sticky_n | N[0]);

  wire [9:0] frac_r = N[9:0] + r;
  wire       c      = r & &N[9:0];

  // ------------------------------------------------ leading-one position / exponent
  wire [4:0] P;
  assign P[0] = h[0]|h[2]|h[4]|h[6]|h[8]|h[10]|h[12]|h[14]|h[16]|h[18]|h[20];
  assign P[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign P[2] = h[0]|h[1]|h[6]|h[7]|h[8]|h[9]|h[14]|h[15]|h[16]|h[17];
  assign P[3] = h[6]|h[7]|h[8]|h[9]|h[10]|h[11]|h[12]|h[13];
  assign P[4] = h[0]|h[1]|h[2]|h[3]|h[4]|h[5];

  wire [6:0] Va = Am + {2'b00, P};
  wire under_sel = Va[6] | ~|Va[6:0];
  wire under34   = Va[6];
  wire over_T    = ~Va[6] & (Va[5] | &Va[4:0]);
  wire eq65      = ~Va[6] & Va[5] & &Va[4:1] & ~Va[0];
  wire over      = over_T | (eq65 & c);
  wire [4:0] field = c ? (Va[4:0] + 5'd1) : Va[4:0];

  // ------------------------------------------------ subnormal result path
  wire [5:0]  sfull = 6'd26 - S;
  wire [4:0]  s5    = sfull[4:0];
  wire [21:0] shs   = prod >> s5;
  wire [9:0]  fsp   = shs[9:0];

  wire [31:0] pg32 = {9'd0, prod, 1'b0};
  wire [31:0] bt2  = {8'd0, bot, 2'b00};
  wire guard_s  = pg32[s5];
  wire sticky_s = bt2[s5];
  wire r_s      = guard_s & (sticky_s | fsp[0]);

  wire [9:0] frac_sb = fsp + r_s;
  wire       carry_s = r_s & &fsp[9:0];

  // ------------------------------------------------ pack
  wire is_nan = nan_in;
  wire is_inf = inf_in | over;
  wire is_sub = under_sel & ~is_nan & ~inf_in;

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_sub            ? {4'd0, carry_s} :
                     under34           ? 5'd0  : field;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      is_inf ? 10'd0 :
                      is_sub ? frac_sb : frac_r;
  wire       sg_y  = is_nan ? 1'b0 : sy;

  assign y = {sg_y, exp_y, frac_y};
endmodule

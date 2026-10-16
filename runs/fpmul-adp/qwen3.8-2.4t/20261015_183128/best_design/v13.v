// fp_mul_e5f10 v13 -- v8 core, but the subnormal alignment + normal/subnormal
// frac/guard/sticky selection is a single AND-OR plane off the multiplier
// product bits (one-hot shift decode, mode selects embedded), replacing the
// funnel shifter and its muxes.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  // ---- classify ----
  wire ez_a = ~|ea, ez_b = ~|eb;
  wire eo_a = &ea,  eo_b = &eb;
  wire fz_a = ~|fa, fz_b = ~|fb;
  wire zero_a = ez_a & fz_a,  zero_b = ez_b & fz_b;
  wire inf_a  = eo_a & fz_a,  inf_b  = eo_b & fz_b;
  wire nan_a  = eo_a & ~fz_a, nan_b  = eo_b & ~fz_b;
  wire sub_a  = ez_a & ~fz_a, sub_b  = ez_b & ~fz_b;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ---- leading-zero count of the 10-bit fractions ----
  wire [3:0] lz_a = fa[9] ? 4'd0 : fa[8] ? 4'd1 : fa[7] ? 4'd2 : fa[6] ? 4'd3 :
                    fa[5] ? 4'd4 : fa[4] ? 4'd5 : fa[3] ? 4'd6 : fa[2] ? 4'd7 :
                    fa[1] ? 4'd8 : 4'd9;
  wire [3:0] lz_b = fb[9] ? 4'd0 : fb[8] ? 4'd1 : fb[7] ? 4'd2 : fb[6] ? 4'd3 :
                    fb[5] ? 4'd4 : fb[4] ? 4'd5 : fb[3] ? 4'd6 : fb[2] ? 4'd7 :
                    fb[1] ? 4'd8 : 4'd9;

  // ---- normalized significands (leading 1 at bit 10) ----
  wire [10:0] ma = sub_a ? ({fa, 1'b0} << lz_a) : {1'b1, fa};
  wire [10:0] mb = sub_b ? ({fb, 1'b0} << lz_b) : {1'b1, fb};

  // ---- exponent sum (early) ----
  wire signed [6:0] ea_e = sub_a ? -$signed({3'b000, lz_a}) : $signed({2'b00, ea});
  wire signed [6:0] eb_e = sub_b ? -$signed({3'b000, lz_b}) : $signed({2'b00, eb});
  wire signed [6:0] t_base = ea_e + eb_e;

  wire tb_le14 = (t_base <= 7'sd14);
  wire tb_le15 = (t_base <= 7'sd15);
  wire tb_ge45 = (t_base >= 7'sd45);
  wire tb_ge46 = (t_base >= 7'sd46);
  wire tb_eq44 = (t_base == 7'sd44);
  wire tb_eq45 = (t_base == 7'sd45);

  wire [4:0] em0 = t_base[4:0] - 5'd15;
  wire [4:0] em1 = t_base[4:0] - 5'd14;
  wire [4:0] em2 = t_base[4:0] - 5'd13;

  wire [3:0] R = (t_base >= 7'sd0) ? (4'd15 - t_base[3:0]) : 4'd15;

  // ---- one-hot decode of R ----
  wire [15:0] Req;
  genvar k;
  generate
    for (k = 0; k < 16; k = k + 1)
      assign Req[k] = (R == k[3:0]);
  endgenerate

  // ---- multiply ----
  wire [21:0] p    = ma * mb;
  wire        p_hi = p[21];

  wire is_sub  = p_hi ? tb_le14 : tb_le15;
  wire ovf_cmp = p_hi ? tb_ge45 : tb_ge46;
  wire t_eq45  = p_hi ? tb_eq44 : tb_eq45;

  // ---- mode selects (late: depend on p_hi) ----
  wire is_nor = ~is_sub;
  wire sel_n0 = is_nor & ~p_hi;    // normal, product leading 1 at bit 20
  wire sel_n1 = is_nor &  p_hi;    // normal, product leading 1 at bit 21
  wire [15:0] sel_s;               // subnormal, shift by R
  generate
    for (k = 0; k < 16; k = k + 1)
      assign sel_s[k] = is_sub & Req[k];
  endgenerate

  // ---- AND-OR plane: frac_o[i] = p[11+k+i] selected by sel_s[k],
  //      or p[10+i] / p[11+i] for normal ----
  wire [9:0] frac_o;
  genvar i;
  generate
    for (i = 0; i < 10; i = i + 1) begin : gfrac
      wire [12:0] term;
      genvar kk;
      for (kk = 0; kk < 11; kk = kk + 1) begin : gs
        if (11 + kk + i <= 21)
          assign term[kk] = sel_s[kk] & p[11+kk+i];
        else
          assign term[kk] = 1'b0;
      end
      assign term[11] = sel_n0 & p[10+i];
      assign term[12] = sel_n1 & p[11+i];
      assign frac_o[i] = |term;
    end
  endgenerate

  // ---- guard ----
  // subnormal: p[10+R] for R <= 11, else 0 ; normal: p[9] or p[10]
  wire [11:0] gterm;
  generate
    for (k = 0; k < 12; k = k + 1)
      assign gterm[k] = sel_s[k] & p[10+k];
  endgenerate
  wire guard = (|gterm) | (sel_n0 & p[9]) | (sel_n1 & p[10]);

  // ---- sticky: prefix OR of p gives st[k] = |p[9+k:0] ----
  wire [21:0] pfx;
  assign pfx[0] = p[0];
  genvar q;
  generate
    for (q = 1; q < 22; q = q + 1) assign pfx[q] = pfx[q-1] | p[q];
  endgenerate
  wire [15:0] sterm;
  generate
    for (k = 0; k < 16; k = k + 1) begin : gst
      if (9 + k <= 21)
        assign sterm[k] = sel_s[k] & pfx[9+k];
      else
        assign sterm[k] = sel_s[k] & pfx[21];
    end
  endgenerate
  wire sticky = (|sterm) | (sel_n0 & pfx[8]) | (sel_n1 & pfx[9]);

  // ---- rounding ----
  wire round_up = guard & (sticky | frac_o[0]);

  // ripple prefix-AND of frac_o for the increment
  wire [9:0] pa;
  assign pa[0] = frac_o[0];
  generate
    for (i = 1; i < 10; i = i + 1) assign pa[i] = pa[i-1] & frac_o[i];
  endgenerate
  wire [9:0] frac_r;
  assign frac_r[0] = frac_o[0] ^ round_up;
  generate
    for (i = 1; i < 10; i = i + 1)
      assign frac_r[i] = frac_o[i] ^ (round_up & pa[i-1]);
  endgenerate
  wire carry = round_up & pa[9];

  // ---- final exponent ----
  wire [1:0] ksel = {1'b0, p_hi} + {1'b0, carry};
  wire [4:0] exp_norm = (ksel == 2'd0) ? em0 : (ksel == 2'd1) ? em1 : em2;
  wire [4:0] exp_o    = is_sub ? {4'd0, carry} : exp_norm;
  wire       ovf      = ovf_cmp | (t_eq45 & carry);

  // ---- flattened packing ----
  wire not_special = ~nan_in & ~inf_in & ~zero_in;
  wire gate_f      = not_special & ~ovf;

  wire [9:0] frac_y;
  assign frac_y[8:0] = frac_r[8:0] & {9{gate_f}};
  assign frac_y[9]   = nan_in | (frac_r[9] & gate_f);

  wire [4:0] exp_y = (nan_in | inf_in) ? 5'd31 :
                     zero_in ? 5'd0 :
                     ovf ? 5'd31 : exp_o;

  wire sign_y = sy & ~nan_in;

  assign y = {sign_y, exp_y, frac_y};
endmodule

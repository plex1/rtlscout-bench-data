// fp_mul_e5f10 v18 -- v14 with balanced prefix OR (sticky) and balanced carry
// computation: pfx uses Kogge-Stone prefix OR (5 levels instead of 21-level
// ripple); carries use &frac computed straight from p (bypassing the ripple).

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // balanced prefix OR: out[i] = |in[i:0] (Kogge-Stone, 22-bit, 5 stages)
  function automatic [21:0] pfx_or22;
    input [21:0] x;
    reg [21:0] g;
    integer i;
    begin
      g = x;
      for (i=21;i>=1;i=i-1) g[i] = g[i] | g[i-1];    // d=1
      for (i=21;i>=2;i=i-1) g[i] = g[i] | g[i-2];    // d=2
      for (i=21;i>=4;i=i-1) g[i] = g[i] | g[i-4];    // d=4
      for (i=21;i>=8;i=i-1) g[i] = g[i] | g[i-8];    // d=8
      for (i=21;i>=16;i=i-1) g[i] = g[i] | g[i-16];  // d=16
      pfx_or22 = g;
    end
  endfunction

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

  // ---- one-hot decode of R (early) ----
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

  // ---- balanced prefix OR of p for sticky ----
  wire [21:0] pfx = pfx_or22(p);

  // ================= subnormal path: AND-OR plane =================
  wire [9:0] frac_s;
  genvar i;
  generate
    for (i = 0; i < 10; i = i + 1) begin : gfs
      wire [11:0] t;
      genvar kk;
      for (kk = 0; kk < 12; kk = kk + 1) begin : gt
        if (11 + kk + i <= 21)
          assign t[kk] = Req[kk] & p[11+kk+i];
        else
          assign t[kk] = 1'b0;
      end
      assign frac_s[i] = |t;
    end
  endgenerate

  wire [11:0] gt;
  generate
    for (k = 0; k < 12; k = k + 1)
      assign gt[k] = Req[k] & p[10+k];
  endgenerate
  wire guard_s = |gt;

  // sticky_s = OR_k Req[k] & pfx[9+k];  k=0..11 valid; k>=12 -> pfx[21]
  wire [11:0] st;
  generate
    for (k = 0; k < 12; k = k + 1)
      assign st[k] = Req[k] & pfx[9+k];
  endgenerate
  wire sticky_s = (|st) | ((Req[12]|Req[13]|Req[14]|Req[15]) & pfx[21]);

  wire round_s = guard_s & (sticky_s | frac_s[0]);

  // ripple prefixes for the rounded fraction bits (needed for XOR outputs)
  wire [9:0] pa_s;
  assign pa_s[0] = frac_s[0];
  generate
    for (i = 1; i < 10; i = i + 1) assign pa_s[i] = pa_s[i-1] & frac_s[i];
  endgenerate
  wire [9:0] frac_s_o;
  assign frac_s_o[0] = frac_s[0] ^ round_s;
  generate
    for (i = 1; i < 10; i = i + 1)
      assign frac_s_o[i] = frac_s[i] ^ (round_s & pa_s[i-1]);
  endgenerate
  // balanced carry for the exponent path (&frac_s = pa_s[9], 10-bit AND tree)
  wire carry_s = round_s & (&frac_s);

  // ================= normal rounding path =================
  wire [20:0] mant = p_hi ? p[21:1] : p[20:0];
  wire [9:0] frac_n  = mant[19:10];
  wire       guard_n = mant[9];
  wire       sticky_n = (|mant[8:0]) | (p_hi & p[0]);
  wire       round_n = guard_n & (sticky_n | frac_n[0]);

  wire [9:0] pa_n;
  assign pa_n[0] = frac_n[0];
  generate
    for (i = 1; i < 10; i = i + 1) assign pa_n[i] = pa_n[i-1] & frac_n[i];
  endgenerate
  wire [9:0] frac_n_o;
  assign frac_n_o[0] = frac_n[0] ^ round_n;
  generate
    for (i = 1; i < 10; i = i + 1)
      assign frac_n_o[i] = frac_n[i] ^ (round_n & pa_n[i-1]);
  endgenerate
  // balanced carry straight from p (no mux, no ripple):
  // (&frac_n) = p_hi ? (&p[20:11]) : (&p[19:10])
  wire carry_n = round_n & (p_hi ? (&p[20:11]) : (&p[19:10]));

  // ================= merge =================
  wire [9:0] frac_o = is_sub ? frac_s_o : frac_n_o;

  // ---- final exponent ----
  wire [1:0] ksel = {1'b0, p_hi} + {1'b0, carry_n};
  wire [4:0] exp_norm = (ksel == 2'd0) ? em0 : (ksel == 2'd1) ? em1 : em2;
  wire [4:0] exp_o    = is_sub ? {4'd0, carry_s} : exp_norm;
  wire       ovf      = ovf_cmp | (t_eq45 & carry_n);

  // ---- flattened packing ----
  wire not_special = ~nan_in & ~inf_in & ~zero_in;
  wire gate_f      = not_special & ~ovf;

  wire [9:0] frac_y;
  assign frac_y[8:0] = frac_o[8:0] & {9{gate_f}};
  assign frac_y[9]   = nan_in | (frac_o[9] & gate_f);

  wire [4:0] exp_y = (nan_in | inf_in) ? 5'd31 :
                     zero_in ? 5'd0 :
                     ovf ? 5'd31 : exp_o;

  wire sign_y = sy & ~nan_in;

  assign y = {sign_y, exp_y, frac_y};
endmodule

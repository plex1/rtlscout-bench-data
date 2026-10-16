// fp_mul_e5f10 v10 -- v8 with balanced (Kogge-Stone) parallel-prefix AND for
// the rounding increment, replacing the linear ripple prefix chain. Every
// prefix pa[i]=&frac[i:0] (and carry) is ~log2(10) levels, computed in
// parallel with the late round signal.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- parallel-prefix AND: pfx[i] = &x[i:0] ----
  // (Kogge-Stone; 10-bit, 4 stages)
  function automatic [9:0] pfx_and10;
    input [9:0] x;
    reg [9:0] g;
    integer i;
    begin
      g = x;                                   // stage d=0
      // d=1
      for (i=9;i>=1;i=i-1) g[i] = g[i] & g[i-1];
      // d=2
      for (i=9;i>=2;i=i-1) g[i] = g[i] & g[i-2];
      // d=4
      for (i=9;i>=4;i=i-1) g[i] = g[i] & g[i-4];
      // d=8
      for (i=9;i>=8;i=i-1) g[i] = g[i] & g[i-8];
      pfx_and10 = g;
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

  // ---- multiply ----
  wire [21:0] p    = ma * mb;
  wire        p_hi = p[21];

  wire is_sub  = p_hi ? tb_le14 : tb_le15;
  wire ovf_cmp = p_hi ? tb_ge45 : tb_ge46;
  wire t_eq45  = p_hi ? tb_eq44 : tb_eq45;

  // ================= subnormal rounding path =================
  wire [11:0] e0 = p[21:10];
  wire [11:0] e1 = R[0] ? (e0 >> 1) : e0;
  wire        o1 = R[0] & e0[0];
  wire [11:0] e2 = R[1] ? (e1 >> 2) : e1;
  wire        o2 = R[1] & (|e1[1:0]);
  wire [11:0] e4 = R[2] ? (e2 >> 4) : e2;
  wire        o4 = R[2] & (|e2[3:0]);
  wire [11:0] e8 = R[3] ? (e4 >> 8) : e4;
  wire        o8 = R[3] & (|e4[7:0]);

  wire [9:0] frac_s  = e8[10:1];
  wire       guard_s = e8[0];
  wire       sticky_s = (|p[9:0]) | o1 | o2 | o4 | o8;
  wire       round_s = guard_s & (sticky_s | frac_s[0]);

  wire [9:0] pa_s = pfx_and10(frac_s);
  wire [9:0] frac_s_o;
  assign frac_s_o[0] = frac_s[0] ^ round_s;
  genvar i;
  generate
    for (i = 1; i < 10; i = i + 1)
      assign frac_s_o[i] = frac_s[i] ^ (round_s & pa_s[i-1]);
  endgenerate
  wire carry_s = round_s & pa_s[9];

  // ================= normal rounding path =================
  wire [9:0] frac_na = p[19:10];
  wire [9:0] frac_nb = p[20:11];
  wire [9:0] frac_n  = p_hi ? frac_nb : frac_na;
  wire       guard_n = p_hi ? p[10] : p[9];
  wire       sticky_n = p_hi ? (|p[9:0]) : (|p[8:0]);
  wire       round_n = guard_n & (sticky_n | frac_n[0]);

  wire [9:0] pa_na = pfx_and10(frac_na);
  wire [9:0] pa_nb = pfx_and10(frac_nb);
  wire [9:0] pa_n  = p_hi ? pa_nb : pa_na;
  wire [9:0] frac_n_o;
  assign frac_n_o[0] = frac_n[0] ^ round_n;
  generate
    for (i = 1; i < 10; i = i + 1)
      assign frac_n_o[i] = frac_n[i] ^ (round_n & pa_n[i-1]);
  endgenerate
  wire carry_n = round_n & pa_n[9];

  // ================= merge =================
  wire [9:0] frac_o = is_sub ? frac_s_o : frac_n_o;
  wire       carry  = is_sub ? carry_s  : carry_n;

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

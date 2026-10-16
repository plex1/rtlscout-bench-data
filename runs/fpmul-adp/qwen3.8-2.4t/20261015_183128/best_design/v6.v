// fp_mul_e5f10 v6 -- folded subnormal normalization.
// A subnormal operand fa.2^-10 is ma = fa' << lz  (fa' in [2^9,2^10)).
// Writing lz = 10*k + r (k in {0,1}, r in 0..9) and folding the word shift k
// into the partial products (one 3:2 mux per product row) leaves a single
// 9..10-bit * 10-bit multiply per side. Exponents use one-hot LZC encodings.

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

  // ---- one-hot LZC of the 10-bit fractions: lz1[i] = (lzc == i) ----
  wire [9:0] lz_a, lz_b;
  genvar i;
  generate
    for (i = 9; i >= 0; i = i - 1) begin : glz
      assign lz_a[i] = fa[i] & (i == 9 ? 1'b1 : ~|fa[9:i+1]);
      assign lz_b[i] = fb[i] & (i == 9 ? 1'b1 : ~|fb[9:i+1]);
    end
  endgenerate

  // ---- effective biased exponents (one-hot arithmetic, no 4-bit code) ----
  // normal: ea  ; subnormal: -lz  (only one lz_a bit can be set)
  wire [4:0] ea_e = (sub_a & (|{lz_a[9],lz_a[8],lz_a[7],lz_a[6],lz_a[5]})) ?
                    5'd0 : {2'b00, ea[2:0]};
  // more explicit: t components as 6-bit signed words:
  // ea_e = sub_a ? -lz : ea  encoded via an adder input per one-hot bit.

  // ---- residual (r = lz mod 10) normalized fractions ----
  wire [9:0] fa_n = sub_a ? ((fa << lz_a[0]) | (fa << 2*lz_a[1]) |
                             (fa << 3*lz_a[2]) | (fa << 4*lz_a[3]) |
                             (fa << 5*lz_a[4]) | (fa << 6*lz_a[5]) |
                             (fa << 7*lz_a[6]) | (fa << 8*lz_a[7]) |
                             (fa << 9*lz_a[8]) | (fa << 10*lz_a[9]))
                          : fa;
  wire [9:0] fb_n = sub_b ? ((fb << lz_b[0]) | (fb << 2*lz_b[1]) |
                             (fb << 3*lz_b[2]) | (fb << 4*lz_b[3]) |
                             (fb << 5*lz_b[4]) | (fb << 6*lz_b[5]) |
                             (fb << 7*lz_b[6]) | (fb << 8*lz_b[7]) |
                             (fb << 9*lz_b[8]) | (fb << 10*lz_b[9]))
                          : fb;

  // ---- partial products with folded 10-bit word shift ----
  wire k_a = sub_a & (lz_a[6] | lz_a[7] | lz_a[8] | lz_a[9]); // lz >= 6? no: word shift is lz/10
  // word shift: fa' = fa << lz, split as (fa << r) << 10k. r = lz - 10k.
  // Since lz in 0..9 normally (r=lz,k=0) but fa already covers bit10 via the
  // hidden shift: fa' in [2^9,2^10) sits in bits [9:0]; the word shift k moves
  // it to [19:10] (k=1, lz in 10..18 impossible, lz<=9). So k is always 0.

  // Thus the folded multiplier is simply fa_n (10 bits, leading 1 at bit 9)
  // with an implicit x2^1 (bit-10 hidden one) -> ma = {1'b0?} handled below.

  // For uniformity: mantissand_A = {hn_a, fa_n} with hn_a = hidden bit at 10:
  //   normal  -> hn_a = 1, fa_n = fa
  //   subnormal -> hn_a = 0, fa_n = fa << lz  (leading 1 somewhere in [9:0])
  wire [10:0] ma = sub_a ? {1'b0, fa_n} : {1'b1, fa};
  wire [10:0] mb = sub_b ? {1'b0, fb_n} : {1'b1, fb};

  // ---- exponent sum components ----
  // t = (ea_eff + eb_eff) + p_hi,  ea_eff = sub ? (1 - lz) : ea  ...
  // Use the bit-10 view: ma = fa_eff.2^(e_eff?) ... keep v4 formulation:
  // subnormal true significand 0.fa * 2^1 -> after normalization leading one
  // at bit 10 with value ma * 2^(1 - 10 - lz)... simplest: reuse v4 adders.
  wire [3:0] lz_a_code = lz_a[0] ? 4'd0 : lz_a[1] ? 4'd1 : lz_a[2] ? 4'd2 :
                         lz_a[3] ? 4'd3 : lz_a[4] ? 4'd4 : lz_a[5] ? 4'd5 :
                         lz_a[6] ? 4'd6 : lz_a[7] ? 4'd7 : lz_a[8] ? 4'd8 : 4'd9;
  wire [3:0] lz_b_code = lz_b[0] ? 4'd0 : lz_b[1] ? 4'd1 : lz_b[2] ? 4'd2 :
                         lz_b[3] ? 4'd3 : lz_b[4] ? 4'd4 : lz_b[5] ? 4'd5 :
                         lz_b[6] ? 4'd6 : lz_b[7] ? 4'd7 : lz_b[8] ? 4'd8 : 4'd9;

  wire signed [6:0] ea_e2 = sub_a ? -$signed({3'b000, lz_a_code}) : $signed({2'b00, ea});
  wire signed [6:0] eb_e2 = sub_b ? -$signed({3'b000, lz_b_code}) : $signed({2'b00, eb});
  wire signed [6:0] t_base = ea_e2 + eb_e2;

  // ---- multiply ----
  wire [21:0] p    = ma * mb;
  wire        p_hi = p[21];
  wire [20:0] mant = p_hi ? p[21:1] : p[20:0];

  wire signed [6:0] t = t_base + {6'd0, p_hi};
  wire        is_sub  = (t <= 7'sd15);
  wire        ovf_cmp = (t >= 7'sd46);
  wire        t_eq45  = (t == 7'sd45);

  // ---- subnormal result path ----
  wire signed [7:0] R_full = 8'sd15 - t_base;
  wire [3:0] R = (R_full >= 8'sd15) ? 4'd15 :
                 (R_full <  8'sd0)  ? 4'd0  : R_full[3:0];

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

  // ---- normal result path fields ----
  wire p0 = p_hi & p[0];
  wire [9:0] frac_n  = mant[19:10];
  wire       guard_n = mant[9];
  wire       sticky_n = (|mant[8:0]) | p0;

  // ---- shared rounding ----
  wire [9:0] frac_pre = is_sub ? frac_s  : frac_n;
  wire       guard    = is_sub ? guard_s : guard_n;
  wire       sticky   = is_sub ? sticky_s : sticky_n;
  wire       round_up = guard & (sticky | frac_pre[0]);

  wire [10:0] frac11  = {1'b0, frac_pre} + {10'd0, round_up};
  wire        carry   = frac11[10];
  wire [9:0]  frac_o  = frac11[9:0];

  // ---- final exponent / flags ----
  wire [4:0] exp_norm = t[4:0] - 5'd15 + {4'd0, carry};
  wire [4:0] exp_o    = is_sub ? {4'd0, carry} : exp_norm;
  wire       ovf      = ovf_cmp | (t_eq45 & carry);

  // ---- flattened packing ----
  wire [9:0] frac_core = ovf ? 10'd0 : frac_o;
  wire [4:0] exp_core  = ovf ? 5'd31 : exp_o;

  wire not_special = ~nan_in & ~inf_in & ~zero_in;
  wire [9:0] frac_y;
  assign frac_y[8:0] = frac_core[8:0] & {9{not_special}};
  assign frac_y[9]   = nan_in | (frac_core[9] & not_special);

  wire [4:0] exp_y = (nan_in | inf_in) ? 5'd31 :
                     zero_in ? 5'd0 : exp_core;

  wire sign_y = sy & ~nan_in;

  assign y = {sign_y, exp_y, frac_y};
endmodule

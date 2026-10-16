// fp_mul_e5f10 v4 -- input-side subnormal normalization (like v2) but the
// subnormal-result alignment shift is computed from t_base = ea_e + eb_e only,
// so the product MSB (p_hi, from the multiplier) is removed from that path.
//
// Derivation: product value = p * 2^(t_base-50), p = ma*mb in [2^20,2^22).
// subnormal value = frac * 2^-24  =>  frac_real = p * 2^(t_base-26) = p >> S,
// S = 26 - t_base.  Equivalently shift p' = p[21:10] right by R = S-11 = 15-t_base,
// with |p[9:0] folded into sticky.  R and the shift are independent of p_hi.

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

  // ---- exponent sums ----
  wire signed [6:0] ea_e = sub_a ? -$signed({3'b000, lz_a}) : $signed({2'b00, ea});
  wire signed [6:0] eb_e = sub_b ? -$signed({3'b000, lz_b}) : $signed({2'b00, eb});
  wire signed [6:0] t_base = ea_e + eb_e;              // pre-multiplier

  // ---- multiply ----
  wire [21:0] p    = ma * mb;
  wire        p_hi = p[21];
  wire [20:0] mant = p_hi ? p[21:1] : p[20:0];         // leading 1 at bit 20

  wire signed [6:0] t = t_base + {6'd0, p_hi};
  wire        is_sub  = (t <= 7'sd15);
  wire        ovf_cmp = (t >= 7'sd46);
  wire        t_eq45  = (t == 7'sd45);

  // ---- subnormal result path (p_hi-independent shift) ----
  wire signed [7:0] R_full = 8'sd15 - t_base;          // = S - 11
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

  // ---- pack ----
  wire [15:0] core = ovf ? {sy, 5'h1F, 10'd0} : {sy, exp_o, frac_o};

  assign y = nan_in  ? 16'h7E00 :
             inf_in  ? {sy, 5'h1F, 10'd0} :
             zero_in ? {sy, 15'd0} :
             core;
endmodule

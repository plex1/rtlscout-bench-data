// fp_mul_e5f10 -- binary16 FP multiplier.
// v7 = v4 + decomposed multiplier (ma=2^10+A, mb=2^10+B:
//      prod = 2^20 + (A+B)*2^10 + A*B, 10x10 core instead of 11x11)
//      + OR-form output packing (flat 2-level logic instead of mux chains).
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack & classify ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire ea_zero = (ea == 5'd0),  eb_zero = (eb == 5'd0);
  wire ea_ones = &ea,           eb_ones = &eb;
  wire fa_zero = (fa == 10'd0), fb_zero = (fb == 10'd0);
  wire zero_a  = ea_zero & fa_zero,  zero_b = eb_zero & fb_zero;
  wire nan_a   = ea_ones & ~fa_zero, nan_b  = eb_ones & ~fb_zero;
  wire inf_a   = ea_ones & fa_zero,  inf_b  = eb_ones & fb_zero;

  wire any_nan  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire any_inf  = inf_a | inf_b;
  wire any_zero = zero_a | zero_b;
  wire sign_p   = sa ^ sb;

  // ---- leading-zero count of the 10-bit fractions ----
  reg [3:0] lza, lzb;
  integer i;
  always @* begin
    lza = 4'd9;
    for (i = 0; i < 10; i = i + 1) if (fa[i]) lza = 4'd9 - i[3:0];
  end
  always @* begin
    lzb = 4'd9;
    for (i = 0; i < 10; i = i + 1) if (fb[i]) lzb = 4'd9 - i[3:0];
  end

  // ---- normalized significands ----
  wire [10:0] fa_sh = {1'b0, fa} << (lza + 4'd1);
  wire [10:0] fb_sh = {1'b0, fb} << (lzb + 4'd1);
  wire [9:0]  A = ea_zero ? fa_sh[9:0] : fa;   // ma = {1'b1, A}
  wire [9:0]  B = eb_zero ? fb_sh[9:0] : fb;   // mb = {1'b1, B}

  // ---- effective exponents & everything derived before the multiply ----
  wire signed [6:0] Ea = ea_zero ? -$signed({3'd0, lza}) : $signed({2'd0, ea});
  wire signed [6:0] Eb = eb_zero ? -$signed({3'd0, lzb}) : $signed({2'd0, eb});
  wire signed [6:0] pre = Ea + Eb - 7'sd15;

  wire signed [6:0] t7 = 7'sd11 - pre;       // subnormal-path shift (p21 cancels)
  wire       pre_neg = pre[6];
  wire       pre_eq0 = (pre[5:0] == 6'd0);
  wire [4:0] pre5    = pre[4:0];
  wire       ge31_0  = ~pre[6] & (pre[5:0] >= 6'd31);
  wire       eq30_0  = (pre[6:0] == 7'd30);
  wire       eq29_0  = (pre[6:0] == 7'd29);
  wire [4:0] e00 = pre5;
  wire [4:0] e01 = pre5 + 5'd1;
  wire [4:0] e11 = pre5 + 5'd2;

  // ---- significand product: (2^10+A)(2^10+B) = 2^20 + (A+B)2^10 + AB ----
  wire [10:0] S    = A + B;
  wire [19:0] P    = A * B;
  wire [21:0] prod = {1'b0, S, 10'd0} + {2'b00, P} + 22'h100000;
  wire        p21  = prod[21];

  wire under = pre_neg | (pre_eq0 & ~p21);     // er <= 0
  wire ge31  = ge31_0 | (eq30_0 & p21);        // er >= 31
  wire eq30  = (eq30_0 & ~p21) | (eq29_0 & p21);

  // ---- normal path: fraction fields straight from the product ----
  wire [9:0] frac_raw = p21 ? prod[20:11] : prod[19:10];
  wire       g_n      = p21 ? prod[10]    : prod[9];
  wire       s_n      = p21 ? |prod[9:0]  : |prod[8:0];
  wire       rup      = g_n & (s_n | frac_raw[0]);
  wire [10:0] frac_n  = {1'b0, frac_raw} + {10'd0, rup};
  wire        carry   = frac_n[10];

  // ---- subnormal path: right shift prod by t7 ----
  reg [9:0] f_tr;
  reg       g_s, st_s;
  always @* begin
    case (t7)
      7'sd11: begin f_tr = prod[20:11];         g_s = prod[10]; st_s = |prod[9:0];  end
      7'sd12: begin f_tr = prod[21:12];         g_s = prod[11]; st_s = |prod[10:0]; end
      7'sd13: begin f_tr = {1'b0, prod[21:13]}; g_s = prod[12]; st_s = |prod[11:0]; end
      7'sd14: begin f_tr = {2'd0, prod[21:14]}; g_s = prod[13]; st_s = |prod[12:0]; end
      7'sd15: begin f_tr = {3'd0, prod[21:15]}; g_s = prod[14]; st_s = |prod[13:0]; end
      7'sd16: begin f_tr = {4'd0, prod[21:16]}; g_s = prod[15]; st_s = |prod[14:0]; end
      7'sd17: begin f_tr = {5'd0, prod[21:17]}; g_s = prod[16]; st_s = |prod[15:0]; end
      7'sd18: begin f_tr = {6'd0, prod[21:18]}; g_s = prod[17]; st_s = |prod[16:0]; end
      7'sd19: begin f_tr = {7'd0, prod[21:19]}; g_s = prod[18]; st_s = |prod[17:0]; end
      7'sd20: begin f_tr = {8'd0, prod[21:20]}; g_s = prod[19]; st_s = |prod[18:0]; end
      7'sd21: begin f_tr = {9'd0, prod[21]};    g_s = prod[20]; st_s = |prod[19:0]; end
      7'sd22: begin f_tr = 10'd0;               g_s = prod[21]; st_s = |prod[20:0]; end
      default: begin f_tr = 10'd0;              g_s = 1'b0;     st_s = 1'b0;        end
    endcase
  end
  wire        rup_s   = g_s & (st_s | f_tr[0]);
  wire [10:0] frac_s  = {1'b0, f_tr} + {10'd0, rup_s};
  wire        carry_s = frac_s[10];

  // ---- result exponent & overflow ----
  wire       ovf  = ~under & (ge31 | (eq30 & carry));
  wire [4:0] t1   = p21 ? e01 : e00;      // er mod 32
  wire [4:0] t2   = p21 ? e11 : e01;      // er + 1 mod 32
  wire [4:0] e_n  = carry ? t2 : t1;

  // valid (finite nonzero) result term; all class terms are mutually exclusive
  wire is_inf   = any_inf | ovf;
  wire res_zero = ~any_nan & ~is_inf &
                  (any_zero | (under & ~carry_s & (frac_s[9:0] == 10'd0)));
  wire valid    = ~any_nan & ~is_inf & ~res_zero;

  // selected fraction before packing
  wire [9:0] f_sel = under ? (carry_s ? 10'd0 : frac_s[9:0]) : frac_n[9:0];

  // ---- OR-form packing (conditions mutually exclusive) ----
  wire spc = any_nan | is_inf;                  // special exponent 31
  wire ex0 = spc | (valid & (under ? carry_s : e_n[0]));  // exp[0]
  assign y[15] = ~any_nan & sign_p;
  assign y[14] = spc | (valid & ~under & e_n[4]);
  assign y[13] = spc | (valid & ~under & e_n[3]);
  assign y[12] = spc | (valid & ~under & e_n[2]);
  assign y[11] = spc | (valid & ~under & e_n[1]);
  assign y[10] = ex0;
  assign y[9]  = any_nan | (valid & f_sel[9]);
  assign y[8]  = valid & f_sel[8];
  assign y[7]  = valid & f_sel[7];
  assign y[6]  = valid & f_sel[6];
  assign y[5]  = valid & f_sel[5];
  assign y[4]  = valid & f_sel[4];
  assign y[3]  = valid & f_sel[3];
  assign y[2]  = valid & f_sel[2];
  assign y[1]  = valid & f_sel[1];
  assign y[0]  = valid & f_sel[0];
endmodule

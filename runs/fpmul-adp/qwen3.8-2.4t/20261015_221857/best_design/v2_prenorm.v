// fp_mul_e5f10 -- binary16 FP multiplier, subnormals pre-normalized before multiply.
// Architecture: 10-bit LZC + small left shift normalize subnormal operands so the
// 11x11 significand product is always normalized (leading one at bit 21 or 20).
// Normal-result normalization is then a 2:1 mux; the subnormal-result path is a
// small right shifter on the product.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  // ---- classify ----
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

  // ---- leading-zero count of the 10-bit fractions (lz = 9 - MSB index) ----
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

  // ---- normalized significands (hidden one restored for subnormals) ----
  wire [10:0] fa_sh = {1'b0, fa} << (lza + 4'd1);
  wire [10:0] fb_sh = {1'b0, fb} << (lzb + 4'd1);
  wire [10:0] ma = ea_zero ? {1'b1, fa_sh[9:0]} : {1'b1, fa};
  wire [10:0] mb = eb_zero ? {1'b1, fb_sh[9:0]} : {1'b1, fb};

  // ---- effective biased exponents: normal e, subnormal -lz ----
  wire signed [6:0] Ea = ea_zero ? -$signed({3'd0, lza}) : $signed({2'd0, ea});
  wire signed [6:0] Eb = eb_zero ? -$signed({3'd0, lzb}) : $signed({2'd0, eb});

  // ---- significand product (always in [2^20, 2^22)) ----
  wire [21:0] prod = ma * mb;
  wire        p21  = prod[21];

  // ---- biased exponent of the exact product, as if normal ----
  wire signed [6:0] er = Ea + Eb - 7'sd15 + {6'd0, p21};

  // ---- mantissa fields of the exact product ----
  wire [9:0] frac_raw = p21 ? prod[20:11] : prod[19:10];
  wire       g_n      = p21 ? prod[10]    : prod[9];
  wire       s_n      = p21 ? |prod[9:0]  : |prod[8:0];

  // ---- normal-path rounding (RTNE) ----
  wire        rup    = g_n & (s_n | frac_raw[0]);
  wire [10:0] frac_r = {1'b0, frac_raw} + {10'd0, rup};
  wire        carry  = frac_r[10];

  wire signed [6:0] exp_n = er + {6'd0, carry};
  wire              under = (er <= 7'sd0);                // subnormal-result path
  wire              over  = ~under & (exp_n >= 7'sd31);   // overflow -> inf

  // ---- subnormal-result path: shift exact product onto the subnormal grid ----
  // frac16_exact = prod >> t, t = 11 + p21 - er ; t in [11, 44], >=23 gives 0
  wire [5:0] t = 6'd11 + {5'd0, p21} - er[5:0];
  reg [9:0] f_tr;
  reg       g_s, st_s;
  always @* begin
    case (t)
      6'd11: begin f_tr = prod[20:11];         g_s = prod[10]; st_s = |prod[9:0];  end
      6'd12: begin f_tr = prod[21:12];         g_s = prod[11]; st_s = |prod[10:0]; end
      6'd13: begin f_tr = {1'b0, prod[21:13]}; g_s = prod[12]; st_s = |prod[11:0]; end
      6'd14: begin f_tr = {2'd0, prod[21:14]}; g_s = prod[13]; st_s = |prod[12:0]; end
      6'd15: begin f_tr = {3'd0, prod[21:15]}; g_s = prod[14]; st_s = |prod[13:0]; end
      6'd16: begin f_tr = {4'd0, prod[21:16]}; g_s = prod[15]; st_s = |prod[14:0]; end
      6'd17: begin f_tr = {5'd0, prod[21:17]}; g_s = prod[16]; st_s = |prod[15:0]; end
      6'd18: begin f_tr = {6'd0, prod[21:18]}; g_s = prod[17]; st_s = |prod[16:0]; end
      6'd19: begin f_tr = {7'd0, prod[21:19]}; g_s = prod[18]; st_s = |prod[17:0]; end
      6'd20: begin f_tr = {8'd0, prod[21:20]}; g_s = prod[19]; st_s = |prod[18:0]; end
      6'd21: begin f_tr = {9'd0, prod[21]};    g_s = prod[20]; st_s = |prod[19:0]; end
      6'd22: begin f_tr = 10'd0;               g_s = prod[21]; st_s = |prod[20:0]; end
      default: begin f_tr = 10'd0;             g_s = 1'b0;     st_s = 1'b0;        end
    endcase
  end

  wire        rup_s   = g_s & (st_s | f_tr[0]);
  wire [10:0] frac_s  = {1'b0, f_tr} + {10'd0, rup_s};
  wire        carry_s = frac_s[10];   // rounded up into the smallest normal
  wire        sub_zero = ~carry_s & (frac_s[9:0] == 10'd0);

  // ---- pack ----
  wire is_nan  = any_nan;
  wire is_inf  = ~is_nan & (any_inf | over);
  wire is_zero = ~is_nan & ~is_inf & (any_zero | (under & sub_zero));
  wire is_sub  = ~is_nan & ~is_inf & under & ~sub_zero;

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 :
                      is_zero             ? 5'd0 :
                      is_sub              ? {4'd0, carry_s} :
                                            exp_n[4:0];
  wire [9:0] frac_y = is_nan  ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 :
                      is_sub  ? (carry_s ? 10'd0 : frac_s[9:0]) :
                                frac_r[9:0];
  wire sign_y = is_nan ? 1'b0 : sign_p;

  assign y = {sign_y, exp_y, frac_y};
endmodule

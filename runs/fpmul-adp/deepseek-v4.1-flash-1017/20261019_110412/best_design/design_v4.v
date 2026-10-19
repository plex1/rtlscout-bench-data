// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
//
// Structure
//   * raw significand product p = ma*mb (22b); its leading one sits at bit 20-lz
//     or 21-lz, where lz = lz(ma)+lz(mb) is derived from the operands (so the
//     leading-zero logic runs in parallel with the multiplier).
//   * a single 23-bit right shifter produces Y = (2*p) >> s, so that Y[0] is the
//     round bit, Y[10:1] the fraction and Y[12:2] the mantissa window.
//     s = 10-lz for a normal result, 26-ea_eff-eb_eff for a subnormal one.
//   * one extra right shift handles the single ambiguous exponent case.
module fp_mul_e5f10 (input [15:0] a, input [15:0] b, output [15:0] y);

  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];

  wire ea_z = ~|ea, eb_z = ~|eb, ea_o = &ea, eb_o = &eb, fa_z = ~|fa, fb_z = ~|fb;
  wire inf_a  = ea_o &  fa_z;
  wire inf_b  = eb_o &  fb_z;
  wire zero_a = ea_z &  fa_z;
  wire zero_b = eb_z &  fb_z;
  wire nan_in = (ea_o & ~fa_z) | (eb_o & ~fb_z) | (inf_a & zero_b) | (zero_a & inf_b);

  // ---------------- leading-zero count of an 11-bit significand -------------
  function [3:0] lz4(input [3:0] v);
    lz4 = v[3] ? 4'd0 : v[2] ? 4'd1 : v[1] ? 4'd2 : v[0] ? 4'd3 : 4'd4;
  endfunction

  function [3:0] lz11(input [10:0] v);
    lz11 = |v[10:7] ? lz4(v[10:7]) :
           |v[6:3]  ? (4'd4 + lz4(v[6:3])) :
           v[2] ? 4'd8 : v[1] ? 4'd9 : v[0] ? 4'd10 : 4'd11;
  endfunction

  wire [10:0] ma = {~ea_z, fa};
  wire [10:0] mb = {~eb_z, fb};
  wire [4:0]  lz = {1'b0, lz11(ma)} + {1'b0, lz11(mb)};

  wire [21:0] p = ma * mb;

  // ---------------- exponent bookkeeping (parallel to the multiplier) --------
  wire [5:0] eae = ea_z ? 6'd1 : {1'b0, ea};
  wire [5:0] ebe = eb_z ? 6'd1 : {1'b0, eb};
  // Z = ea_eff + eb_eff - lz - 14 ; result exponent (pre-round) = Z + e - 1
  wire [6:0] Z = {1'b0, eae} + {1'b0, ebe} - 7'd14 - {2'd0, lz};

  wire [4:0] s_n   = 5'd10 - lz;                  // normal-result shift
  wire [4:0] s_sub = 5'd26 - eae[4:0] - ebe[4:0]; // subnormal-result shift
  // exponent of the normalised interpretation (0 means subnormal encoding)
  wire [6:0] Zp    = Z + {6'd0, p[21]};
  wire       sub   = Zp[6] | ~|Zp[6:1];           // Zp <= 1 : subnormal result
  wire [4:0] s     = sub ? s_sub : s_n;

  // ---------------- shifter with sticky accumulation ------------------------
  wire [22:0] X  = {p, 1'b0};
  wire [22:0] t1 = s[0] ? (X  >> 1)  : X;
  wire        q0 = s[0] & X[0];
  wire [22:0] t2 = s[1] ? (t1 >> 2)  : t1;
  wire        q1 = s[1] & (|t1[1:0]);
  wire [22:0] t3 = s[2] ? (t2 >> 4)  : t2;
  wire        q2 = s[2] & (|t2[3:0]);
  wire [22:0] t4 = s[3] ? (t3 >> 8)  : t3;
  wire        q3 = s[3] & (|t3[7:0]);
  wire [22:0] t5 = s[4] ? (t4 >> 16) : t4;
  wire        q4 = s[4] & (|t4[15:0]);

  wire [12:0] Y  = t5[12:0];
  wire        st = q0 | q1 | q2 | q3 | q4;

  // ---------------- alignment / rounding ------------------------------------
  wire        e     = Y[12];                       // mantissa needs one more shift
  wire        sh1   = sub ? 1'b0 : e;

  wire [10:0] M   = sub ? Y[10:0] : (e ? Y[12:2] : Y[11:1]);
  wire        g   = sh1 ? Y[1] : Y[0];
  wire        stk = st | (sh1 & Y[0]);
  wire        rup = g & (stk | M[0]);

  wire [11:0] S   = {1'b0, M} + {11'd0, rup};
  wire [6:0]  En  = Z + {6'd0, e} - 7'd1;
  wire [6:0]  exp_n = En + {6'd0, S[11]};

  wire [6:0]  exp_pre = sub ? {6'd0, S[10]} : exp_n;
  wire [9:0]  frac    = sub ? S[9:0] : S[10:0];

  // ---------------- packing --------------------------------------------------
  wire is_inf  = ~nan_in & (inf_a | inf_b | (exp_pre > 7'd30));
  wire is_z    = ~nan_in & ~is_inf & (zero_a | zero_b);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_z ? 5'd0 : exp_pre[4:0];
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_z) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {sign_y, exp_y, frac_y};
endmodule

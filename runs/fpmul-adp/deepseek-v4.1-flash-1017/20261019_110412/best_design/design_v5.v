// fp_mul_e5f10 -- IEEE binary16 multiplier (1s,5e,10f) with subnormal support.
//
// * significands are normalised with a leading-zero count before the multiply
//   (identity for normal operands), so the 22-bit product always has its leading
//   one at bit 20 or 21 -> post-multiply normalisation is a single bit select.
// * Y0 = (2*prod)>>10 is free wiring; Y0[0] is the round bit for the normal
//   case and a 4-stage 13-bit shifter realigns the subnormal case.
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

  // ---------------- leading-zero count -------------------------------------
  function [3:0] lz4(input [3:0] v);
    lz4 = v[3] ? 4'd0 : v[2] ? 4'd1 : v[1] ? 4'd2 : v[0] ? 4'd3 : 4'd4;
  endfunction

  function [3:0] lz11(input [10:0] v);
    lz11 = |v[10:7] ? lz4(v[10:7]) :
           |v[6:3]  ? (4'd4 + lz4(v[6:3])) :
           v[2] ? 4'd8 : v[1] ? 4'd9 : v[0] ? 4'd10 : 4'd11;
  endfunction

  function [10:0] shl11(input [10:0] v, input [3:0] s);
    reg [10:0] t;
    begin
      t = s[0] ? {v[9:0], 1'b0} : v;
      t = s[1] ? {t[8:0], 2'b0} : t;
      t = s[2] ? {t[6:0], 4'b0} : t;
      t = s[3] ? {t[2:0], 8'b0} : t;
      shl11 = t;
    end
  endfunction

  wire [10:0] ma = {~ea_z, fa};
  wire [10:0] mb = {~eb_z, fb};
  wire [3:0]  lza = lz11(ma), lzb = lz11(mb);
  wire [10:0] na = shl11(ma, lza);
  wire [10:0] nb = shl11(mb, lzb);

  wire [21:0] p = na * nb;
  wire        e = p[21];

  // ---------------- alignment window ---------------------------------------
  wire [12:0] Y0 = {p[21], p[20:9]};
  wire        st0 = |p[8:0];

  // ---------------- exponent bookkeeping -----------------------------------
  wire [5:0] esum = (ea_z ? 6'd1 : {1'b0, ea}) + (eb_z ? 6'd1 : {1'b0, eb});
  wire [5:0] lzs  = {2'd0, lza} + {2'd0, lzb};
  wire [6:0] ebi  = {1'b0, esum} - {1'b0, lzs} + {6'd0, e} - 7'd15;
  wire       sub  = ebi[6] | ~|ebi;

  wire [6:0] ss = 7'd16 + {1'b0, lzs} - {1'b0, esum};
  wire [3:0] sh = sub ? ((|ss[6:4]) ? 4'd15 : ss[3:0]) : {3'd0, e};

  // ---------------- 13-bit right shifter with sticky ------------------------
  wire [12:0] u1 = sh[0] ? {1'b0, Y0[12:1]} : Y0;
  wire        q0 = sh[0] & Y0[0];
  wire [12:0] u2 = sh[1] ? {2'b0, u1[12:2]} : u1;
  wire        q1 = sh[1] & (|u1[1:0]);
  wire [12:0] u3 = sh[2] ? {4'b0, u2[12:4]} : u2;
  wire        q2 = sh[2] & (|u2[3:0]);
  wire [12:0] u4 = sh[3] ? {8'b0, u3[12:8]} : u3;
  wire        q3 = sh[3] & (|u3[7:0]);
  wire        st = q0 | q1 | q2 | q3 | st0;

  // ---------------- round to nearest, ties to even -------------------------
  wire [10:0] M   = u4[11:1];
  wire        rup = u4[0] & (st | M[0]);
  wire [11:0] S   = {1'b0, M} + {11'd0, rup};

  wire [6:0]  exp_pre = sub ? {6'd0, S[10]} : (ebi + {6'd0, S[11]});
  wire [9:0]  frac    = S[9:0];

  // ---------------- packing ------------------------------------------------
  wire is_inf = ~nan_in & (inf_a | inf_b | (exp_pre > 7'd30));
  wire is_z   = ~nan_in & ~is_inf & (zero_a | zero_b);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_z ? 5'd0 : exp_pre[4:0];
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_z) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {sign_y, exp_y, frac_y};
endmodule

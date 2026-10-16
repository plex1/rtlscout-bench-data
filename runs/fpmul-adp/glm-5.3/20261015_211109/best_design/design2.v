// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
// Optimized: single unified right shifter, balanced priority encoder, or-into-sticky.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack / classify ----
  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],  fb = b[9:0];

  wire ea0 = (ea == 5'd0),  eb0 = (eb == 5'd0);
  wire eaM = (ea == 5'd31), ebM = (eb == 5'd31);
  wire fza = (fa == 10'd0), fzb = (fb == 10'd0);

  wire zero_a = ea0 & fza,  zero_b = eb0 & fzb;
  wire inf_a  = eaM & fza,  inf_b  = ebM & fzb;
  wire nan_a  = eaM & ~fza, nan_b  = ebM & ~fzb;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  // ---- effective significands and exponents (subnormal: e=1, no hidden bit) ----
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [5:0]  es = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};
  wire [21:0] prod = ma * mb;

  // ---- leading-one position of prod, balanced tree ----
  function [1:0] prio3;
    input [2:0] v;
    prio3 = v[2] ? 2'd2 : {1'b0, v[1]};
  endfunction
  function [2:0] prio6;
    input [5:0] v;
    prio6 = (|v[5:3]) ? (3'd3 + {1'b0, prio3(v[5:3])}) : {1'b0, prio3(v[2:0])};
  endfunction
  function [3:0] prio11;
    input [10:0] v;
    prio11 = (|v[10:6]) ? (4'd6 + {1'b0, prio6(v[10:6])}) : {1'b0, prio6(v[5:0])};
  endfunction
  function [4:0] prio22;
    input [21:0] v;
    prio22 = (|v[21:11]) ? (5'd11 + {1'b0, prio11(v[21:11])}) : {1'b0, prio11(v[10:0])};
  endfunction

  wire [4:0] p = prio22(prod);          // leading-one index (0 if prod==0)

  // result is normal iff p + es >= 36  (unbiased exp = p+es-50 >= -14)
  wire [6:0] t = {2'b0, p} + {1'b0, es};
  wire       is_norm = (t >= 7'd36);

  // unified shift: SH = max(p-10, 26-es); >= 0 always
  wire [4:0] sh_n = p - 5'd10;
  wire [4:0] sh_s = 6'd26 - {1'b0, es};
  wire [4:0] SH   = is_norm ? sh_n : sh_s;

  // ---- sticky-or right shifter: win[i] = prod[i+SH-1] ----
  wire [22:0] v0 = {prod, 1'b0};
  wire [22:0] v1 = SH[0] ? {1'b0,  v0[22:1]}  : v0;
  wire [22:0] v2 = SH[1] ? {2'b0,  v1[22:2]}  : v1;
  wire [22:0] v3 = SH[2] ? {4'b0,  v2[22:4]}  : v2;
  wire [22:0] v4 = SH[3] ? {8'b0,  v3[22:8]}  : v3;
  wire [22:0] v5 = SH[4] ? {16'b0, v4[22:16]} : v4;
  wire sticky = (SH[0] & v0[0])    | (SH[1] & |v1[1:0]) | (SH[2] & |v2[3:0])
              | (SH[3] & |v3[7:0]) | (SH[4] & |v4[15:0]);
  wire [11:0] win   = v5[11:0];
  wire [10:0] mant  = win[11:1];       // leading one at bit 10 (normal case)
  wire        guard = win[0];

  // ---- round to nearest, ties to even ----
  wire        rnd  = guard & (sticky | mant[0]);
  wire [11:0] msum = {1'b0, mant} + {11'd0, rnd};

  wire [6:0]  field_n = t - 7'd35 + {6'd0, msum[11]};
  wire        ovf     = is_norm & (field_n > 7'd30);

  wire is_inf    = ~nan_in & (inf_in | (ovf & ~zero_in));
  wire sub_zero  = ~is_norm & (msum[10:0] == 11'd0);
  wire is_zero   = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 :
                      is_zero ? 5'd0 :
                      is_norm ? field_n[4:0] : {4'd0, msum[10]};
  wire [9:0] frac_y = nan_in ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 : msum[9:0];

  assign y = {(nan_in ? 1'b0 : (sa ^ sb)), exp_y, frac_y};
endmodule

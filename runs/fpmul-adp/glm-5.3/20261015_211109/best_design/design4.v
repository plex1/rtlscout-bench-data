// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
// 12-bit CLZ + unified sticky shifter + flat output tail.

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

  wire zero_a  = ea0 & fza,  zero_b  = eb0 & fzb;
  wire inf_a   = eaM & fza,  inf_b   = ebM & fzb;
  wire nan_a   = eaM & ~fza, nan_b   = ebM & ~fzb;

  wire n_in    = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  // ---- effective significands / exponent sum (subnormal: e=1, no hidden bit) ----
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [5:0]  es = ea + eb + ea0 + eb0;
  wire [21:0] prod = ma * mb;

  // ---- leading zeros of the top 12 product bits ----
  function [1:0] clz3;  input [2:0] v; clz3  = v[2] ? 2'd0 : v[1] ? 2'd1 : 2'd2; endfunction
  function [2:0] clz6;  input [5:0] v; clz6  = (|v[5:3]) ? {1'b0, clz3(v[5:3])}
                                                         : (3'd3 + {1'b0, clz3(v[2:0])}); endfunction
  function [3:0] clz12; input [11:0] v; clz12 = (|v[11:6]) ? {1'b0, clz6(v[11:6])}
                                                           : (4'd6 + {1'b0, clz6(v[5:0])}); endfunction

  wire [3:0] s = clz12(prod[21:10]);

  // R = es - 14 - s : >= 1  <=>  result is normal
  wire [6:0] esm14 = {3'b0, es} - 7'd14;
  wire [6:0] R     = esm14 - {3'd0, s};
  wire       is_norm = ~zero_in & ~R[6] & (|R[5:0]);

  // ---- unified shift ----
  wire [4:0] sh_n = 5'd11 - {1'b0, s};
  wire [4:0] sh_s = 6'd26 - {1'b0, es};
  wire [4:0] SH   = is_norm ? sh_n : sh_s;

  // ---- sticky-or right shifter ----
  wire [22:0] v0 = {prod, 1'b0};
  wire [22:0] v1 = SH[0] ? {1'b0,  v0[22:1]}  : v0;
  wire [22:0] v2 = SH[1] ? {2'b0,  v1[22:2]}  : v1;
  wire [22:0] v3 = SH[2] ? {4'b0,  v2[22:4]}  : v2;
  wire [22:0] v4 = SH[3] ? {8'b0,  v3[22:8]}  : v3;
  wire [22:0] v5 = SH[4] ? {16'b0, v4[22:16]} : v4;
  wire sticky = (SH[0] & v0[0])    | (SH[1] & |v1[1:0]) | (SH[2] & |v2[3:0])
              | (SH[3] & |v3[7:0]) | (SH[4] & |v4[15:0]);
  wire [11:0] win   = v5[11:0];
  wire [10:0] mant  = win[11:1];
  wire        guard = win[0];

  // ---- round to nearest, ties to even ----
  wire        rnd  = guard & (sticky | mant[0]);
  wire [11:0] msum = {1'b0, mant} + {11'd0, rnd};
  wire        carry = msum[11];

  // ---- exponent field and overflow ----
  wire [6:0] F   = R + {6'd0, carry};
  wire       ovf = ~F[6] & (F[5] | &F[4:0]);
  wire       i_o = ~n_in & (inf_in | ovf);

  wire [4:0] exp_y  = (n_in | i_o) ? 5'd31 :
                      is_norm ? F[4:0] : msum[10];
  wire [9:0] frac_y = n_in ? 10'h200 : (i_o ? 10'd0 : msum[9:0]);

  assign y = {(n_in ? 1'b0 : (sa ^ sb)), exp_y, frac_y};
endmodule

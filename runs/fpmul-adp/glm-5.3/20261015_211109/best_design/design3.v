// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
// Wallace-tree multiplier + 12-bit CLZ normalisation + unified sticky shifter.

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

  // ---- effective significands / exponents ----
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [5:0]  es = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};

  // ---- Wallace tree (11x11 -> carry-save) ----
  function [43:0] csa;
    input [21:0] x, y, z;
    csa = {x ^ y ^ z, ((x & y) | (x & z) | (y & z)) << 1};
  endfunction

  wire [241:0] R;
  genvar i;
  generate
    for (i = 0; i < 11; i = i + 1) begin : gpp
      assign R[i*22 +: 22] = {11'b0, ma & {11{mb[i]}}} << i;
    end
  endgenerate

  wire [21:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5, s6, c6, s7, c7, s8, c8, PS, PC;
  assign {s1, c1} = csa(R[21:0],    R[43:22],   R[65:44]);
  assign {s2, c2} = csa(R[87:66],   R[109:88],  R[131:110]);
  assign {s3, c3} = csa(R[153:132], R[175:154], R[197:176]);
  assign {s4, c4} = csa(s1, c1, s2);
  assign {s5, c5} = csa(c2, s3, c3);
  assign {s6, c6} = csa(s4, c4, s5);
  assign {s7, c7} = csa(c5, R[219:198], R[241:220]);
  assign {s8, c8} = csa(s6, c6, s7);
  assign {PS, PC} = csa(s8, c8, c7);

  wire [21:0] prod = PS + PC;

  // ---- leading zeros of the top 12 product bits (normal-path shift) ----
  // prod < 2^10 implies both operands subnormal (es=2) => subnormal result.
  function [1:0] clz3;  input [2:0] v; clz3  = v[2] ? 2'd0 : v[1] ? 2'd1 : 2'd2; endfunction
  function [2:0] clz6;  input [5:0] v; clz6  = (|v[5:3]) ? {1'b0, clz3(v[5:3])}
                                                         : (3'd3 + {1'b0, clz3(v[2:0])}); endfunction
  function [3:0] clz12; input [11:0] v; clz12 = (|v[11:6]) ? {1'b0, clz6(v[11:6])}
                                                           : (4'd6 + {1'b0, clz6(v[5:0])}); endfunction

  wire [11:0] V   = prod[21:10];
  wire        Vnz = |V;
  wire [3:0]  s   = clz12(V);
  wire [6:0]  t   = {1'b0, es} + (7'd21 - {3'd0, s});   // = p + es
  wire        is_norm = Vnz & (t >= 7'd36);

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
  wire [11:0] win  = v5[11:0];
  wire [10:0] mant = win[11:1];
  wire        guard = win[0];

  // ---- round to nearest, ties to even ----
  wire        rnd  = guard & (sticky | mant[0]);
  wire [11:0] msum = {1'b0, mant} + {11'd0, rnd};

  wire [6:0]  field_n = t - 7'd35 + {6'd0, msum[11]};
  wire        ovf     = is_norm & (field_n > 7'd30);

  wire is_inf   = ~nan_in & (inf_in | (ovf & ~zero_in));
  wire sub_zero = ~is_norm & (msum[10:0] == 11'd0);
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 :
                      is_zero ? 5'd0 :
                      is_norm ? field_n[4:0] : {4'd0, msum[10]};
  wire [9:0] frac_y = nan_in ? 10'h200 :
                      (is_inf | is_zero) ? 10'd0 : msum[9:0];

  assign y = {(nan_in ? 1'b0 : (sa ^ sb)), exp_y, frac_y};
endmodule

// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
// CSA-tree multiplier + Sklansky prefix CPA + prefix-AND rounding incrementer.

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

  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [5:0]  es = ea + eb + ea0 + eb0;

  // ================= 11x11 carry-save tree =================
  function [43:0] csa;
    input [21:0] x, y, z;
    csa = {x ^ y ^ z, ((x & y) | (x & z) | (y & z)) << 1};
  endfunction

  wire [241:0] R;
  genvar gi;
  generate
    for (gi = 0; gi < 11; gi = gi + 1) begin : gpp
      assign R[gi*22 +: 22] = {11'b0, ma & {11{mb[gi]}}} << gi;
    end
  endgenerate

  wire [21:0] q01,r01,q34,r34,q67,r67, w1,w2,w3,w4, x1,x2,x3,x4, z1,z2, PS,PC;
  assign {q01,r01} = csa(R[21:0],    R[43:22],  R[65:44]);
  assign {q34,r34} = csa(R[87:66],   R[109:88], R[131:110]);
  assign {q67,r67} = csa(R[153:132], R[175:154],R[197:176]);
  assign {w1,w2}   = csa(q01, r01, q34);
  assign {w3,w4}   = csa(r34, q67, r67);
  assign {x1,x2}   = csa(w1, w2, w3);
  assign {x3,x4}   = csa(w4, R[219:198], R[241:220]);
  assign {z1,z2}   = csa(x1, x2, x3);
  assign {PS,PC}   = csa(z1, z2, x4);

  // ================= Sklansky prefix CPA: prod = PS + PC =================
  wire [21:0] Ag = PS & PC;
  wire [21:0] Ap = PS ^ PC;
  wire [21:0] G1,P1,G2,P2,G3,P3,G4,P4,G5;
  wire [21:0] prod;
  genvar gb;
  generate
    for (gb = 0; gb < 22; gb = gb + 1) begin : lv1
      assign G1[gb] = ((gb % 2) == 1) ? (Ag[gb] | (Ap[gb] & Ag[gb-1])) : Ag[gb];
      assign P1[gb] = ((gb % 2) == 1) ? (Ap[gb] & Ap[gb-1]) : Ap[gb];
    end
    for (gb = 0; gb < 22; gb = gb + 1) begin : lv2
      assign G2[gb] = ((gb % 4) >= 2) ? (G1[gb] | (P1[gb] & G1[gb-(gb%4)+1])) : G1[gb];
      assign P2[gb] = ((gb % 4) >= 2) ? (P1[gb] & P1[gb-(gb%4)+1]) : P1[gb];
    end
    for (gb = 0; gb < 22; gb = gb + 1) begin : lv3
      assign G3[gb] = ((gb % 8) >= 4) ? (G2[gb] | (P2[gb] & G2[gb-(gb%8)+3])) : G2[gb];
      assign P3[gb] = ((gb % 8) >= 4) ? (P2[gb] & P2[gb-(gb%8)+3]) : P2[gb];
    end
    for (gb = 0; gb < 22; gb = gb + 1) begin : lv4
      assign G4[gb] = ((gb % 16) >= 8) ? (G3[gb] | (P3[gb] & G3[gb-(gb%16)+7])) : G3[gb];
      assign P4[gb] = ((gb % 16) >= 8) ? (P3[gb] & P3[gb-(gb%16)+7]) : P3[gb];
    end
    for (gb = 0; gb < 22; gb = gb + 1) begin : lv5
      assign G5[gb] = (gb >= 16) ? (G4[gb] | (P4[gb] & G4[15])) : G4[gb];
    end
    for (gb = 1; gb < 22; gb = gb + 1) begin : sbit
      assign prod[gb] = Ap[gb] ^ G5[gb-1];
    end
  endgenerate
  assign prod[0] = Ap[0];

  // ---- leading zeros of the top 12 product bits ----
  function [1:0] clz3;  input [2:0] v; clz3  = v[2] ? 2'd0 : v[1] ? 2'd1 : 2'd2; endfunction
  function [2:0] clz6;  input [5:0] v; clz6  = (|v[5:3]) ? {1'b0, clz3(v[5:3])}
                                                         : (3'd3 + {1'b0, clz3(v[2:0])}); endfunction
  function [3:0] clz12; input [11:0] v; clz12 = (|v[11:6]) ? {1'b0, clz6(v[11:6])}
                                                           : (4'd6 + {1'b0, clz6(v[5:0])}); endfunction

  wire [3:0] s = clz12(prod[21:10]);
  wire [6:0] t = {1'b0, es} + (7'd21 - {3'd0, s});   // = p + es
  wire       is_norm = ~zero_in & (t >= 7'd36);

  wire [4:0] sh_n = 5'd11 - {1'b0, s};
  wire [4:0] sh_s = 6'd26 - {1'b0, es};
  wire [4:0] SH   = is_norm ? sh_n : sh_s;

  // ---- sticky-or right shifter ----
  wire [22:0] d0 = {prod, 1'b0};
  wire [22:0] d1 = SH[0] ? {1'b0,  d0[22:1]}  : d0;
  wire [22:0] d2 = SH[1] ? {2'b0,  d1[22:2]}  : d1;
  wire [22:0] d3 = SH[2] ? {4'b0,  d2[22:4]}  : d2;
  wire [22:0] d4 = SH[3] ? {8'b0,  d3[22:8]}  : d3;
  wire [22:0] d5 = SH[4] ? {16'b0, d4[22:16]} : d4;
  wire sticky = (SH[0] & d0[0])    | (SH[1] & |d1[1:0]) | (SH[2] & |d2[3:0])
              | (SH[3] & |d3[7:0]) | (SH[4] & |d4[15:0]);
  wire [11:0] win   = d5[11:0];
  wire [10:0] mant  = win[11:1];
  wire        guard = win[0];

  // ---- round to nearest, ties to even: prefix-AND incrementer ----
  wire rnd = guard & (sticky | mant[0]);

  wire [10:0] A1, A2, A3, A4;
  genvar q;
  generate
    for (q = 0; q < 11; q = q + 1) begin : pa1
      assign A1[q] = ((q % 2) == 1) ? (mant[q] & mant[q-1]) : mant[q];
    end
    for (q = 0; q < 11; q = q + 1) begin : pa2
      assign A2[q] = ((q % 4) >= 2) ? (A1[q] & A1[q-(q%4)+1]) : A1[q];
    end
    for (q = 0; q < 11; q = q + 1) begin : pa3
      assign A3[q] = ((q % 8) >= 4) ? (A2[q] & A2[q-(q%8)+3]) : A2[q];
    end
    for (q = 0; q < 11; q = q + 1) begin : pa4
      assign A4[q] = (q >= 8) ? (A3[q] & A3[7]) : A3[q];
    end
  endgenerate

  wire [11:0] msum = {rnd & A4[10], mant ^ ({11{rnd}} & {A4[9:0], 1'b1})};
  wire        carry = msum[11];

  // ---- exponent field, overflow, pack ----
  wire [6:0]  F   = t - 7'd35 + {6'd0, carry};
  wire        ovf = is_norm & (F > 7'd30);
  wire        i_o = ~n_in & (inf_in | ovf);

  wire sub_zero = ~is_norm & (msum[10:0] == 11'd0);
  wire is_zero  = ~n_in & ~i_o & (zero_in | sub_zero);

  wire [4:0] exp_y  = (n_in | i_o) ? 5'd31 :
                      is_zero ? 5'd0 :
                      is_norm ? F[4:0] : {4'd0, msum[10]};
  wire [9:0] frac_y = n_in ? 10'h200 : (i_o ? 10'd0 : msum[9:0]);

  assign y = {(n_in ? 1'b0 : (sa ^ sb)), exp_y, frac_y};
endmodule

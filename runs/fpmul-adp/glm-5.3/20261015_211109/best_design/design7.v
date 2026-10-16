// fp_mul_e5f10 -- IEEE binary16 multiplier (1s, 5e, 10f), subnormals supported.
// Key idea: lz(prod) = la + lb + (1-d) with d = prod[21-L], L = la+lb (operand CLZs,
// computed in parallel with the multiplier). Both shift controls are therefore early.

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
  wire [21:0] prod = ma * mb;

  // ---- operand leading-zero counts (saturating at 10) ----
  function [1:0] clz2;  input [1:0] w; clz2  = w[1] ? 2'd0 : 2'd1; endfunction
  function [1:0] clz3;  input [2:0] w; clz3  = w[2] ? 2'd0 : (w[1] ? 2'd1 : 2'd2); endfunction
  function [2:0] clz4;  input [3:0] w; clz4  = (|w[3:2]) ? {1'b0, clz2(w[3:2])}
                                                         : (2'd2 + {1'b0, clz2(w[1:0])}); endfunction
  function [2:0] clz6;  input [5:0] w; clz6  = (|w[5:3]) ? {1'b0, clz3(w[5:3])}
                                                         : (3'd3 + {1'b0, clz3(w[2:0])}); endfunction
  function [3:0] clz10; input [9:0] w; clz10 = (|w[9:4]) ? {1'b0, clz6(w[9:4])}
                                                         : (4'd6 + {2'b00, clz4(w[3:0])}); endfunction

  wire [3:0] la = ma[10] ? 4'd0 : (4'd1 + clz10(ma[9:0]));
  wire [3:0] lb = mb[10] ? 4'd0 : (4'd1 + clz10(mb[9:0]));
  wire [4:0] L  = {1'b0, la} + {1'b0, lb};          // 0..20
  wire [4:0] shA = 5'd11 - L;                       // 1..11 (valid range)
  wire [6:0] K   = {1'b0, es} + 7'd32 - {2'b00, L}; // es - L + 32, in [14,94]

  // d = prod[21-L] : 1 iff the product's leading one is at bit 21-L
  wire d = prod[21 - L];

  wire is_norm = ~zero_in & ((K >= 7'd48) | ((K == 7'd47) & d));

  // ---- shifter A (normal path): {prod,00} >> (11-L) ----
  wire [23:0] VA = {prod, 2'b0};
  wire [23:0] A1 = shA[0] ? {1'b0,  VA[23:1]}  : VA;
  wire [23:0] A2 = shA[1] ? {2'b0,  A1[23:2]}  : A1;
  wire [23:0] A3 = shA[2] ? {4'b0,  A2[23:4]}  : A2;
  wire [23:0] A4 = shA[3] ? {8'b0,  A3[23:8]}  : A3;
  wire stA = (shA[0] & VA[0]) | (shA[1] & |A1[1:0])
           | (shA[2] & |A2[3:0]) | (shA[3] & |A3[7:0]);
  wire [12:0] T = A4[12:0];                      // T[k] = prod[k+9-L]

  wire [10:0] mant_n = d ? T[12:2] : T[11:1];
  wire        guard_n  = d ? T[1]  : T[0];
  wire        sticky_n = d ? (stA | T[0]) : stA;

  // ---- shifter B (subnormal path): {prod,0} >> (26-es) ----
  wire [4:0]  shB = 6'd26 - {1'b0, es};
  wire [22:0] VB = {prod, 1'b0};
  wire [22:0] B1 = shB[0] ? {1'b0,  VB[22:1]}  : VB;
  wire [22:0] B2 = shB[1] ? {2'b0,  B1[22:2]}  : B1;
  wire [22:0] B3 = shB[2] ? {4'b0,  B2[22:4]}  : B2;
  wire [22:0] B4 = shB[3] ? {8'b0,  B3[22:8]}  : B3;
  wire [22:0] B5 = shB[4] ? {16'b0, B4[22:16]} : B4;
  wire stB = (shB[0] & VB[0]) | (shB[1] & |B1[1:0])
           | (shB[2] & |B2[3:0]) | (shB[3] & |B3[7:0]) | (shB[4] & |B4[15:0]);
  wire [11:0] U = B5[11:0];

  // ---- select path, round to nearest / ties to even ----
  wire [10:0] mant   = is_norm ? mant_n   : U[11:1];
  wire        guard  = is_norm ? guard_n  : U[0];
  wire        sticky = is_norm ? sticky_n : stB;

  wire        rnd  = guard & (sticky | mant[0]);
  wire [11:0] msum = {1'b0, mant} + {11'd0, rnd};
  wire        carry = msum[11];

  // ---- exponent / overflow / pack ----
  wire [6:0] Kh = K + {6'd0, d};
  wire [6:0] En = Kh - 7'd47;                     // exp field minus rounding carry
  wire [4:0] exp_n = carry ? (En[4:0] + 5'd1) : En[4:0];
  wire       ovf = (Kh >= 7'd78) | ((Kh == 7'd77) & carry);
  wire       i_o = ~n_in & (inf_in | ovf);

  wire [4:0] exp_y  = (n_in | i_o) ? 5'd31 : (is_norm ? exp_n : {4'd0, msum[10]});
  wire [9:0] frac_y = n_in ? 10'h200 : (i_o ? 10'd0 : msum[9:0]);

  assign y = {(n_in ? 1'b0 : (sa ^ sb)), exp_y, frac_y};
endmodule

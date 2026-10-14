// fp_mul_e5f10 v5 -- same structure as v3, but sticky is obtained from the
// barrel shifters themselves: the OR of all bits shifted out at each stage
// equals |prod[S-2:0].  This removes the Kogge-Stone prefix array and the
// two 22:1 sticky index muxes, replacing them with small OR trees over
// existing stage outputs.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- classify ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire ea0 = ~|ea,  eb0 = ~|eb;
  wire ea1 =  &ea,  eb1 =  &eb;
  wire zero_a = ea0 & ~|fa, zero_b = eb0 & ~|fb;
  wire inf_a  = ea1 & ~|fa, inf_b  = eb1 & ~|fb;
  wire nan_a  = ea1 &  |fa, nan_b  = eb1 &  |fb;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ---- significands / effective exponents ----
  wire [10:0] ma   = {~ea0, fa};
  wire [10:0] mb   = {~eb0, fb};
  wire [4:0]  ea_e = ea0 ? 5'd1 : ea;
  wire [4:0]  eb_e = eb0 ? 5'd1 : eb;
  wire [5:0]  esum = {1'b0, ea_e} + {1'b0, eb_e};

  wire [21:0] prod = ma * mb;

  // ---- balanced msb index of prod (padded to 32) ----
  wire [31:0] v = {10'd0, prod};
  wire h4 = |v[31:16];
  wire h3 = h4 ? |v[31:24] : |v[15:8];
  wire h2 = h4 ? (h3 ? |v[31:28] : |v[23:20])
               : (h3 ? |v[15:12] : |v[7:4]);
  wire h1 = h4 ? (h3 ? (h2 ? |v[31:30] : |v[27:26])
                     : (h2 ? |v[23:22] : |v[19:18]))
               : (h3 ? (h2 ? |v[15:14] : |v[11:10])
                     : (h2 ? |v[7:6]   : |v[3:2]));
  wire h0 = h4 ? (h3 ? (h2 ? (h1 ? |v[31] : |v[29])
                           : (h1 ? |v[27] : |v[25]))
                     : (h2 ? (h1 ? |v[23] : |v[21])
                           : (h1 ? |v[19] : |v[17])))
               : (h3 ? (h2 ? (h1 ? |v[15] : |v[13])
                           : (h1 ? |v[11] : |v[9]))
                     : (h2 ? (h1 ? |v[7]  : |v[5])
                           : (h1 ? |v[3]  : |v[1])));
  wire [4:0] msb = {h4, h3, h2, h1, h0};

  // ---- shift amounts, exponent terms, branch decision ----
  wire [4:0] sr  = msb - 5'd10;        // >= 0 whenever branch = normal
  wire [5:0] t6  = 6'd26 - esum;
  wire [4:0] tsh = t6[4:0];            // in [1,24] whenever branch = sub

  wire [6:0] n0    = {1'b0, esum} + {2'b00, msb};
  wire       sub   = (n0 <= 7'd35);
  wire [6:0] npre  = n0 - 7'd35;
  wire       gt65  = (n0 > 7'd65);
  wire       eq65  = (n0 == 7'd65);

  // ---- barrel shifter with dropped-bit sticky (normal branch) ----
  wire [22:0] prodL = {prod, 1'b0};

  wire [22:0] n1 = sr[0] ? (prodL >> 1) : prodL;
  wire [22:0] n2 = sr[1] ? (n1   >> 2) : n1;
  wire [22:0] n3 = sr[2] ? (n2   >> 4) : n2;
  wire [22:0] n4 = sr[3] ? (n3   >> 8) : n3;
  wire [22:0] n5 = sr[4] ? (n4  >> 16) : n4;
  wire [11:0] shN = n5[11:0];
  // bits dropped at each stage: stage k drops input bits [2^k-1:0]
  wire stN = (sr[1] & (n1[0] | n1[1]))
           | (sr[2] & (|n2[3:0]))
           | (sr[3] & (|n3[7:0]))
           | (sr[4] & (|n4[15:0]));
  // sr[0] drops prodL[0] which is always 0 -- nothing to add.

  // ---- barrel shifter with dropped-bit sticky (subnormal branch) ----
  wire [22:0] t1 = tsh[0] ? (prodL >> 1) : prodL;
  wire [22:0] t2 = tsh[1] ? (t1   >> 2) : t1;
  wire [22:0] t3 = tsh[2] ? (t2   >> 4) : t2;
  wire [22:0] t4 = tsh[3] ? (t3   >> 8) : t3;
  wire [22:0] t5 = tsh[4] ? (t4  >> 16) : t4;
  wire [11:0] shS = t5[11:0];
  wire stS = (tsh[1] & (t1[0] | t1[1]))
           | (tsh[2] & (|t2[3:0]))
           | (tsh[3] & (|t3[7:0]))
           | (tsh[4] & (|t4[15:0]));

  // ---- select branch ----
  wire [11:0] sh     = sub ? shS  : shN;
  wire        sticky = sub ? stS  : stN;

  // ---- shared round-to-nearest-even ----
  wire guard    = sh[0];
  wire rlsb     = sh[1];
  wire round_up = guard & (sticky | rlsb);

  // carry out of the 10-bit window + round_up: only when window all ones
  wire        cs       = round_up & (&sh[10:1]);
  wire [10:0] sum10    = {1'b0, sh[10:1]} + {10'd0, round_up};
  wire [9:0]  frac_sub = sum10[9:0];
  wire        carry_n  = sh[11] & cs;
  wire [9:0]  frac_norm = carry_n ? 10'd0 : frac_sub;

  // ---- exponent / overflow ----
  wire [4:0] enm      = npre[4:0] + {4'd0, carry_n};
  wire       overflow = gt65 | (eq65 & carry_n);

  // ---- pack ----
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_subf = ~is_nan & ~is_inf & sub;

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31
                    : is_subf            ? {4'd0, cs}
                    :                      enm;
  wire [9:0] frac_y = is_nan  ? 10'h200
                    : is_inf  ? 10'd0
                    : is_subf ? frac_sub
                    :           frac_norm;
  wire       sign_y = is_nan ? 1'b0 : sy;

  assign y = {sign_y, exp_y, frac_y};
endmodule

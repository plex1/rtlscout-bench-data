// fp_mul_e5f10 v6 -- parallel two-branch datapath (like the starting point:
// each branch shifts/rounds independently, selection only at pack), plus:
//  * balanced mux-tree msb encoder (no 22-deep priority chain)
//  * Kogge-Stone prefix-OR for sticky (shared, log depth)
//  * rounding carry-out computed as round_up & window-all-ones (the sum
//    adder is then only needed for the fraction data, not the exponent)
//  * exponent pre-offset A = n0-35 computed early; final exponent is a
//    2-level selection between A and A+1 by the rounding carry
//  * overflow = gt65 | (eq65 & carry) with cheap bit tests on n0
//  * pack without the is_zero select tree (zero inputs are always in the
//    sub branch and round to cs=0, frac=0)
//  * dead left-shift path of the normal branch removed (normal results
//    always have msb >= 10; for msb < 10 the carry-out is 0 anyway, which
//    keeps the branch selection intact)

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS         = 15;
  localparam MAX_E        = 31;
  localparam MAX_FINITE_E = 30;

  // ---- unpack / classify ----
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

  // ---- Kogge-Stone prefix-OR for sticky ----
  localparam NW = 22;
  wire [6*NW-1:0] lv;
  assign lv[NW-1:0] = prod;
  genvar gi, gj;
  generate
    for (gi = 0; gi < 5; gi = gi + 1) begin : g_pref
      localparam integer SH = (1 << gi);
      for (gj = 0; gj < NW; gj = gj + 1) begin : g_bit
        if (gj >= SH)
          assign lv[(gi+1)*NW+gj] = lv[gi*NW+gj] | lv[gi*NW+gj-SH];
        else
          assign lv[(gi+1)*NW+gj] = lv[gi*NW+gj];
      end
    end
  endgenerate
  wire [NW-1:0] below   = lv[5*NW +: NW];
  wire [31:0]   below32 = {10'd0, below};

  // ---- branch decision and shift amounts (early) ----
  wire [4:0] sr  = msb - 5'd10;   // >= 0 whenever branch = normal
  wire [5:0] t6  = 6'd26 - esum;
  wire [4:0] tsh = t6[4:0];       // in [1,24] whenever branch = sub

  wire [6:0] n0    = {1'b0, esum} + {2'b00, msb};
  wire       sub   = ~n0[6] & (~n0[5] | ~(|n0[4:2]));   // n0 <= 35
  wire       gt65  = n0[6] & (|n0[5:1]);                // n0 >  65
  wire       eq65  = n0[6] & ~(|n0[5:1]) & n0[0];       // n0 == 65
  wire [4:0] npre  = n0[4:0] - 5'd35;                   // (n0-35) low bits
  wire [4:0] npre1 = npre + 5'd1;                       // early +1

  // ---- parallel shift + round: normal branch ----
  wire [22:0] prodL   = {prod, 1'b0};
  wire [11:0] shN     = prodL >> sr;
  wire        guardN  = shN[0];
  wire        lsbN    = shN[1];
  wire        stickyN = below32[sr - 5'd2];
  wire        ruN     = guardN & (stickyN | lsbN);
  wire [10:0] mantN   = shN[11:1];
  wire [10:0] sumN    = {1'b0, mantN} + {10'd0, ruN};
  wire        carry_n = ruN & (&mantN);          // == sumN[11]
  wire [9:0]  frac_norm = carry_n ? 10'd0 : sumN[9:0];

  // ---- parallel shift + round: subnormal branch ----
  wire [11:0] shS     = prodL >> tsh;
  wire        guardS  = shS[0];
  wire        lsbS    = shS[1];
  wire        stickyS = below32[tsh - 5'd2];
  wire        ruS     = guardS & (stickyS | lsbS);
  wire [9:0]  fracS   = shS[10:1];
  wire [10:0] sumS    = {1'b0, fracS} + {10'd0, ruS};
  wire        cs      = ruS & (&fracS);          // == sumS[10]
  wire [9:0]  frac_sub = sumS[9:0];

  // ---- exponent / overflow ----
  wire [4:0] enm      = carry_n ? npre1 : npre;
  wire       overflow = gt65 | (eq65 & carry_n);

  // ---- pack ----
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_subf = ~is_nan & ~is_inf & sub;
  wire M       = is_nan | is_inf;

  wire [4:0] exp_y  = M ? 5'd31 : (sub ? {4'd0, cs} : enm);
  wire [9:0] frac_y = is_nan ? 10'h200
                    : is_inf ? 10'd0
                    : is_subf ? frac_sub
                    :          frac_norm;
  wire       sign_y = is_nan ? 1'b0 : sy;

  assign y = {sign_y, exp_y, frac_y};
endmodule

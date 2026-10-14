// fp_mul_e5f10 v12 -- v8 with one change: flattened fraction pack.
// In the single-shared-shifter structure a selected sub branch always has
// carry_n = 0, so frac = nan ? payload : (is_inf | carry_n) ? 0 : frac_sub
// is equivalent to the nested 4:1 form but shallower.

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

  // ---- shift amounts, early exponent terms, branch decision ----
  wire [4:0] sr  = msb - 5'd10;        // >= 0 whenever branch = normal
  wire [5:0] t6  = 6'd26 - esum;
  wire [4:0] tsh = t6[4:0];            // in [1,24] whenever branch = sub

  wire [6:0] n0    = {1'b0, esum} + {2'b00, msb};
  wire       sub   = (n0 <= 7'd35);
  wire [6:0] npre  = n0 - 7'd35;
  wire [4:0] npre1 = npre[4:0] + 5'd1;     // early, parallel to shifters
  wire       gt65  = (n0 > 7'd65);
  wire       eq65  = (n0 == 7'd65);

  // ---- sticky: Kogge-Stone prefix-OR, one mux per branch ----
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

  wire sticky_n = below32[sr  - 5'd2];
  wire sticky_t = below32[tsh - 5'd2];
  wire sticky   = sub ? sticky_t : sticky_n;

  // ---- two right shifters, mux after shift ----
  wire [22:0] prodL = {prod, 1'b0};
  wire [11:0] shN   = prodL >> sr;
  wire [11:0] shS   = prodL >> tsh;
  wire [11:0] sh    = sub ? shS : shN;

  // ---- shared round-to-nearest-even ----
  wire guard    = sh[0];
  wire rlsb     = sh[1];
  wire round_up = guard & (sticky | rlsb);

  wire        cs       = round_up & (&sh[10:1]);
  wire [10:0] sum10    = {1'b0, sh[10:1]} + {10'd0, round_up};
  wire [9:0]  frac_sub = sum10[9:0];
  wire        carry_n  = sh[11] & cs;
  wire [9:0]  frac_norm = carry_n ? 10'd0 : frac_sub;

  // ---- exponent / overflow ----
  wire [4:0] enm      = carry_n ? npre1 : npre[4:0];
  wire       overflow = gt65 | (eq65 & carry_n);

  // ---- pack ----
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_subf = ~is_nan & ~is_inf & sub;

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31
                    : is_subf            ? {4'd0, cs}
                    :                      enm;
  wire [9:0] frac_y = is_nan ? 10'h200
                    : (is_inf | carry_n) ? 10'd0
                    : frac_sub;
  wire       sign_y = is_nan ? 1'b0 : sy;

  assign y = {sign_y, exp_y, frac_y};
endmodule

// fp_mul_e5f10 -- IEEE binary16 multiplier, subnormal support, unified single-rounding path.
//
// Derivation: value = prod * 2^(exp_sum-50), prod = ma*mb (11x11).
// With m = msb index of prod, T = m + exp_sum, the result exponent field is
// E = T - 35.  E <= 0  <=>  subnormal result, in which case the fixed subnormal
// grid gives M = round(prod >> (26-exp_sum)); otherwise M = round(prod >> (m-10)).
// Both shift amounts are provably in [0,24] (never negative), so ONE right-only
// barrel shifter and ONE rounder serve both paths.  Zero products fall out as M=0.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],  fb = b[9:0];

  wire eza = ~|ea, ezb = ~|eb;
  wire e1a = &ea,   e1b = &eb;
  wire za  = eza & ~|fa, zb = ezb & ~|fb;
  wire na  = e1a & |fa,  nb = e1b & |fb;
  wire ia  = e1a & ~|fa, ib = e1b & ~|fb;

  wire nan_in = na | nb | (ia & zb) | (za & ib);
  wire inf_in = ia | ib;
  wire sy     = sa ^ sb;

  // effective significands / exponents (subnormal: hidden 0, exponent 1)
  wire [10:0] ma = {~eza, fa}, mb = {~ezb, fb};
  wire [4:0]  ea_eff = eza ? 5'd1 : ea;
  wire [4:0]  eb_eff = ezb ? 5'd1 : eb;
  wire [5:0]  exp_sum = ea_eff + eb_eff;          // 2..62
  wire [21:0] prod = ma * mb;

  // msb index of prod (0 if prod == 0)
  reg [4:0] m;
  integer   i;
  always @* begin
    m = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) m = i[4:0];
  end

  // below[k] = |prod[k:0]
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  wire [31:0] prod_ext  = prod;    // zero-extend so dynamic index is safe
  wire [31:0] below_ext = below;

  wire [6:0] T   = {2'b0, m} + {1'b0, exp_sum};
  wire       sub = (T <= 7'd35);                   // subnormal (or zero) result
  wire [4:0] s   = sub ? (6'd26 - {1'b0, exp_sum}) : (m - 5'd10);
  wire [4:0] sm1 = s - 5'd1;                       // wraps to >=30 when s==0
  wire [4:0] sm2 = s - 5'd2;                       // wraps to >=30 when s<=1

  wire [10:0] Mpre   = prod >> s;                  // 11-bit mantissa window
  wire        guard  = prod_ext[sm1];
  wire        sticky = below_ext[sm2];
  wire        rnd    = guard & (sticky | Mpre[0]); // round to nearest, ties to even
  wire [11:0] M      = {1'b0, Mpre} + {11'd0, rnd};
  wire [6:0]  Epre   = T - 7'd35;
  wire [6:0]  Efin   = Epre + {6'b0, M[11]};
  wire        over   = ~sub & (Efin > 7'd30);       // only possible when ~sub

  assign y = nan_in              ? 16'h7E00 :
             (inf_in | over)     ? {sy, 5'h1F, 10'd0} :
             sub                 ? {sy, 4'b0000, M[10], M[9:0]} :
                                    {sy, Efin[4:0], M[9:0]};
endmodule

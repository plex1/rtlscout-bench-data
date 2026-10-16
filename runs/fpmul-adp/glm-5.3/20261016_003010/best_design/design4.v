// fp_mul_e5f10 -- IEEE binary16 multiplier, subnormal support.
// msb(ma*mb) = msb(ma)+msb(mb)+{0,1}: the normalising shift amount comes from the
// OPERANDS (parallel with the multiplier); a 1-bit fixup mux (c) finishes it.
// Subnormal results use an exact left shift onto the fixed subnormal grid.
// One shared rounder (increment-style) serves both paths.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],  fb = b[9:0];

  wire eza = ~|ea, ezb = ~|eb;
  wire e1a = &ea,   e1b = &eb;
  wire za  = eza & ~|fa, zb = ezb & ~|fb;
  wire na  = e1a & |fa,  nb = e1b & |fb;
  wire ia  = e1a & ~|fa, ib = e1b & ~|fb;

  wire nan_in  = na | nb | (ia & zb) | (za & ib);
  wire inf_in  = ia | ib;
  wire sy      = sa ^ sb;

  wire [10:0] ma = {~eza, fa}, mb = {~ezb, fb};
  wire [4:0]  ea_eff = eza ? 5'd1 : ea;
  wire [4:0]  eb_eff = ezb ? 5'd1 : eb;
  wire [6:0]  exp_sum = {2'b0, ea_eff} + {2'b0, eb_eff};
  wire [21:0] prod = ma * mb;

  // operand-side leading-bit positions (early, off the multiply path)
  reg [3:0] pf, qf;
  integer   i1, i2;
  always @* begin pf = 4'd0; for (i1 = 0; i1 < 10; i1 = i1 + 1) if (fa[i1]) pf = i1[3:0]; end
  always @* begin qf = 4'd0; for (i2 = 0; i2 < 10; i2 = i2 + 1) if (fb[i2]) qf = i2[3:0]; end
  wire [4:0] p = eza ? {1'b0, pf} : 5'd10;   // msb index of ma
  wire [4:0] q = ezb ? {1'b0, qf} : 5'd10;   // msb index of mb
  wire [4:0] pq = p + q;                     // 0..20, msb(prod) = pq + c

  // ---- normal path: align leading one of prod to bit 21 ----
  wire [22:0] Lfull = {prod, 1'b0} << (5'd20 - pq);  // prod[pq]->bit21, prod[pq+1]->bit22
  wire        c   = Lfull[22];                  // leading one is one higher
  wire [21:0] L   = c ? Lfull[22:1] : Lfull[21:0];   // 1-bit fixup shift
  wire [10:0] mant_n = L[21:11];               // 1.f (11 bits)
  wire        g_n    = L[10];
  wire        st_n   = |L[9:0];

  // ---- subnormal path: fixed subnormal grid via exact left shift ----
  // sub => exp_sum in [2,14] => X = exp_sum-2 in [0,12]; Z = prod<<X is EXACT.
  wire [4:0]  Xs  = exp_sum[4:0] - 5'd2;
  wire [33:0] Z   = prod << Xs;
  wire [9:0]  f_s = Z[33:24];
  wire        g_s = Z[23];
  wire        st_s = |Z[22:0];

  // W = msb(prod) + exp_sum: result exponent = W + round_carry - 35
  wire [6:0] W   = {2'b0, pq} + exp_sum + {6'b0, c};
  // subnormal (or zero) result: a zero operand forces the sub grid (prod = 0
  // there, so the sub branch naturally emits the signed zero).
  wire za_ = eza & ~|fa;
  wire zb_ = ezb & ~|fb;
  wire sub = (za_ | zb_) | (W <= 7'd35);

  wire [10:0] Msel  = sub ? {1'b0, f_s} : mant_n;
  wire        guard = sub ? g_s : g_n;
  wire        sticky= sub ? st_s : st_n;
  wire        round = guard & (sticky | Msel[0]); // round to nearest, ties to even

  // increment-style rounder: M = Msel + round  (round in {0,1})
  reg [10:0] pr;   // pr[j] = &Msel[j:0]
  reg [11:0] M;
  integer    jj;
  always @* begin
    pr[0] = Msel[0];
    for (jj = 1; jj < 11; jj = jj + 1) pr[jj] = pr[jj-1] & Msel[jj];
    M[0]  = Msel[0] ^ round;
    for (jj = 1; jj < 11; jj = jj + 1) M[jj] = Msel[jj] ^ (round & pr[jj-1]);
    M[11] = round & pr[10];
  end

  wire [6:0] Wm    = W + {6'b0, M[11]};
  wire       over  = ~sub & (Wm > 7'd65);      // E2 = Wm-35 > 30
  wire [6:0] E2    = Wm - 7'd35;
  wire [9:0] fnorm = M[11] ? M[10:1] : M[9:0];

  assign y = nan_in          ? 16'h7E00 :
             (inf_in | over) ? {sy, 5'h1F, 10'd0} :
             sub             ? {sy, 4'b0000, M[10], M[9:0]} :
                                {sy, E2[4:0], fnorm};
endmodule

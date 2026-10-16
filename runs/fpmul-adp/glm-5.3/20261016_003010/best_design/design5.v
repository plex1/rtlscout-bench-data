// fp_mul_e5f10 -- IEEE binary16 multiplier, subnormal support.
// msb(ma*mb) = msb(ma)+msb(mb)+{0,1}: the normalising shift comes from the
// OPERANDS (parallel with the multiplier); a 1-bit fixup mux (c) finishes it.
// Subnormal results use an exact left shift onto the fixed subnormal grid.
// One shared increment-style rounder serves both paths.
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
  wire zero_in = za | zb;
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
  wire [6:0] s  = {2'b0, pq} + exp_sum;      // = msb(prod) + exp_sum - c

  // ---- normal path: align leading one of prod to bit 21 ----
  // live normal cases have pq in [10,20] so the shift 20-pq fits 4 bits
  wire [4:0]  sn    = 5'd20 - pq;
  wire [22:0] Lfull = {prod, 1'b0} << sn[3:0]; // prod[pq]->bit21, prod[pq+1]->bit22
  wire        c     = Lfull[22];             // leading one is one higher
  // c-fixup: with c=1 the leading one is at bit 22, so the window moves up 1
  wire [9:0]  mant_n = c ? Lfull[21:12] : Lfull[20:11];  // fraction of 1.f
  wire        g_n    = c ? Lfull[11]    : Lfull[10];
  wire        st_n   = |Lfull[9:0] | (c & Lfull[10]);

  // ---- subnormal path: fixed subnormal grid via exact left shift ----
  // sub => one operand subnormal => exp_sum <= 31 => X = exp_sum-2 in [0,29]
  wire [4:0]  Xs  = exp_sum[4:0] - 5'd2;
  wire [33:0] Z   = prod << Xs;
  wire [9:0]  f_s = Z[33:24];
  wire        g_s = Z[23];
  wire        st_s = |Z[22:0];

  // T = s + c = msb(prod) + exp_sum;  result exponent = T + round_carry - 35
  wire sub = (s < 7'd35) | ((s == 7'd35) & ~c);   // T <= 35
  wire [6:0] E0 = s - 7'd35;                      // early, c-less exponent base
  wire [6:0] E1 = E0 + {6'b0, c};

  wire [10:0] Msel  = {~sub, sub ? f_s : mant_n}; // top bit: 1 iff normal
  wire        guard = sub ? g_s : g_n;
  wire        sticky= sub ? st_s : st_n;
  wire        round = guard & (sticky | Msel[0]); // round to nearest, ties to even

  // increment-style rounder: M = Msel + round
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

  wire [6:0] E2    = E1 + {6'b0, M[11]};
  wire       over  = ~sub & ((E1 > 7'd30) | ((E1 == 7'd30) & M[11]));
  wire [9:0] fnorm = M[11] ? M[10:1] : M[9:0];

  assign y = nan_in          ? 16'h7E00 :
             (inf_in | over) ? {sy, 5'h1F, 10'd0} :
             zero_in         ? {sy, 15'd0} :
             sub             ? {sy, 4'b0000, M[10], M[9:0]} :
                                {sy, E2[4:0], fnorm};
endmodule

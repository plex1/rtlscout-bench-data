// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign / 5 exp / 10 frac), subnormals supported.
//
// Datapath notes:
//  * prod = ma*mb (11x11).  Let L be the position of the leading one of prod
//    (L = 21 - lz).  A one-hot bus h[k] (leading one at bit 21-k) drives the
//    normalisation window muxes directly, skipping a binary lz conversion on
//    that path.
//  * Result biased exponent E = L + S + c - 35 (S = ea_eff+eb_eff, c = rounding
//    carry).  T = L + S is formed early; c only enters one mux level before the
//    output.  under: T+c <= 35, over: T+c >= 66 -- the T-side comparisons are
//    precomputed in parallel.
//  * Subnormal results are re-derived straight from prod and S (shift 26-S),
//    which are both early signals, so that path is off the critical path.
//  * Zero inputs/zero products fall out of the subnormal path naturally.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ------------------------------------------------ unpack / classify
  wire       sa = a[15],      sb = b[15];
  wire [4:0] ea = a[14:10],   eb = b[14:10];
  wire [9:0] fa = a[9:0],     fb = b[9:0];

  wire eaz = ~|ea,  ebz = ~|eb;
  wire eao =  &ea,  ebo =  &eb;
  wire faz = ~|fa,  fbz = ~|fb;

  wire nan_in = (eao & ~faz) | (ebo & ~fbz)
              | ((eao & faz) & (ebz & fbz))
              | ((eaz & faz) & (ebo & fbz));
  wire inf_in = (eao & faz) | (ebo & fbz);
  wire s_y    = sa ^ sb;

  // ------------------------------------------------ significands / exponent sum
  wire [10:0] ma = {~eaz, fa};
  wire [10:0] mb = {~ebz, fb};
  wire [5:0]  S  = (eaz ? 6'd1 : {1'b0, ea}) + (ebz ? 6'd1 : {1'b0, eb});

  // ------------------------------------------------ product, leading-one bus
  wire [21:0] prod = ma * mb;
  wire        pn   = |prod;

  wire [21:0] nz;    // nz[m] = |prod[21:m]
  wire [21:0] bot;   // bot[m] = |prod[m:0]
  wire [21:0] h;     // h[k] : leading one of prod at bit 21-k
  genvar i;
  generate
    for (i = 1; i < 22; i = i + 1) begin : g_pf
      assign nz[i]  = |prod[21:i];
      assign bot[i] = |prod[i:0];
      assign h[i]   = prod[21-i] & ~nz[22-i];
    end
  endgenerate
  assign nz[0]  = |prod;
  assign bot[0] = prod[0];
  assign h[0]   = prod[21];

  // ------------------------------------------------ normal window (one-hot mux)
  // N[j] = prod[j+11-lz]; with lz = k the term is h[k] & prod[j+11-k]
  wire [10:0] N;
  genvar j, k;
  generate
    for (j = 0; j < 11; j = j + 1) begin : g_win
      wire [11:0] t;
      for (k = 0; k < 12; k = k + 1) begin : g_t
        assign t[k] = h[k] & prod[j+11-k];
      end
      assign N[j] = |t;
    end
  endgenerate

  // guard (prod[10-lz]) and sticky (|prod[9-lz:0]) of the shift
  wire [10:0] gt;
  wire [9:0]  st;
  generate
    for (k = 0; k <= 10; k = k + 1) begin : g_g
      assign gt[k] = h[k] & prod[10-k];
    end
    for (k = 0; k <= 9; k = k + 1) begin : g_s
      assign st[k] = h[k] & bot[9-k];
    end
  endgenerate
  wire guard_n  = |gt;
  wire sticky_n = |st;
  wire r        = guard_n & (sticky_n | N[0]);   // round to nearest, ties to even

  wire [10:0] fr     = {1'b0, N[9:0]} + {10'd0, r};
  wire        c      = fr[10];                   // rounding carry into the exponent
  wire [9:0]  frac_r = fr[9:0];

  // ------------------------------------------------ result exponent
  wire [4:0] L;                                  // leading-one position
  assign L[0] = h[0]|h[2]|h[4]|h[6]|h[8]|h[10]|h[12]|h[14]|h[16]|h[18]|h[20];
  assign L[1] = h[2]|h[3]|h[6]|h[7]|h[10]|h[11]|h[14]|h[15]|h[18]|h[19];
  assign L[2] = h[0]|h[1]|h[6]|h[7]|h[8]|h[9]|h[14]|h[15]|h[16]|h[17];
  assign L[3] = h[6]|h[7]|h[8]|h[9]|h[10]|h[11]|h[12]|h[13];
  assign L[4] = h[0]|h[1]|h[2]|h[3]|h[4]|h[5];

  wire [6:0] T = {2'b00, L} + {1'b0, S};         // E = T + c - 35

  wire [6:0] Va = T - 7'd35;
  wire [6:0] Vb = T - 7'd34;
  wire [4:0] field = c ? Vb[4:0] : Va[4:0];

  wire le34 = (T <= 7'd34);
  wire eq35 = (T == 7'd35);
  wire ge66 = (T >= 7'd66);
  wire eq65 = (T == 7'd65);
  wire under = le34 | (eq35 & ~c);               // subnormal (or zero) result
  wire over  = ge66 | (eq65 &  c);               // overflow to inf

  // ------------------------------------------------ subnormal result path
  // value*2^24 = prod * 2^(S-26): right shift by 26-S (1..24 when under)
  wire [5:0]  sfull = 6'd26 - S;
  wire [4:0]  s5    = sfull[4:0];
  wire [21:0] shs   = prod >> s5;
  wire [9:0]  fsp   = shs[9:0];

  wire [31:0] pg32  = {9'd0, prod, 1'b0};        // pg32[m]  = prod[m-1]
  wire [31:0] bot2  = {{8{pn}}, bot, 2'b00};     // bot2[m]  = bot[m-2]
  wire guard_s  = pg32[s5];
  wire sticky_s = bot2[s5];
  wire r_s      = guard_s & (sticky_s | fsp[0]);

  wire [10:0] frs     = {1'b0, fsp} + {10'd0, r_s};
  wire        carry_s = frs[10];                 // rounded up to smallest normal
  wire [9:0]  frac_sb = frs[9:0];

  // ------------------------------------------------ pack
  wire is_nan = nan_in;
  wire is_inf = ~is_nan & (inf_in | over);
  wire is_sub = ~is_nan & ~is_inf & under;

  wire [4:0] exp_y = (is_nan | is_inf) ? 5'd31 :
                     is_sub            ? {4'd0, carry_s} : field;
  wire [9:0] frac_y = is_nan ? 10'h200 :
                      is_inf ? 10'd0 :
                      is_sub ? frac_sb : frac_r;
  wire       sg_y  = is_nan ? 1'b0 : s_y;

  assign y = {sg_y, exp_y, frac_y};
endmodule

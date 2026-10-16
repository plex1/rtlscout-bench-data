// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Structure: pre-normalise the (at most one) subnormal significand so the
// 11x11 product is already normalised (leading 1 at bit 21 or 20).
// The subnormal-result path shifts the product with an INPUT-ONLY shift amount
// s = 15 - (ea+eb) + lz, so no product bit ever feeds a shifter select.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire        sa = a[15];
  wire        sb = b[15];
  wire [4:0]  ea = a[14:10];
  wire [4:0]  eb = b[14:10];
  wire [9:0]  fa = a[9:0];
  wire [9:0]  fb = b[9:0];

  wire ea0 = (ea == 5'd0);
  wire eb0 = (eb == 5'd0);
  wire eaX = (ea == 5'd31);
  wire ebX = (eb == 5'd31);
  wire faz = (fa == 10'd0);
  wire fbz = (fb == 10'd0);
  wire za  = ea0 & faz;
  wire zb  = eb0 & fbz;
  wire na  = eaX & ~faz;
  wire nb  = ebX & ~fbz;
  wire ia  = eaX & faz;
  wire ib  = ebX & fbz;

  wire nan_in  = na | nb | (ia & zb) | (za & ib);
  wire inf_in  = ia | ib;
  wire zero_in = za | zb;
  wire sgn     = sa ^ sb;

  // ---- normalise the (at most one) subnormal significand ----
  wire [9:0] fsub = ea0 ? fa : fb;
  wire [3:0] clz =
      fsub[9]   ? 4'd0 :
      fsub[8]   ? 4'd1 :
      fsub[7]   ? 4'd2 :
      fsub[6]   ? 4'd3 :
      fsub[5]   ? 4'd4 :
      fsub[4]   ? 4'd5 :
      fsub[3]   ? 4'd6 :
      fsub[2]   ? 4'd7 :
      fsub[1]   ? 4'd8 : 4'd9;

  wire [9:0]  gn  = fsub << clz;                  // msb at bit 9
  wire [4:0]  lza = 5'd1 + {1'b0, clz};
  wire [4:0]  lz  = (ea0 ? lza : 5'd0) + (eb0 ? lza : 5'd0);

  wire [10:0] ma = ea0 ? {gn, 1'b0} : {1'b1, fa};
  wire [10:0] mb = eb0 ? {gn, 1'b0} : {1'b1, fb};

  wire [21:0] pp = ma * mb;                       // = raw product << lz

  // ---- exponent sums (input only) ----
  wire [5:0] exps = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};

  // ---- normal result path ----
  wire        n1    = pp[21];
  wire [22:0] nrm   = n1 ? {1'b0, pp} : {pp, 1'b0};
  wire [10:0] mpre  = nrm[21:11];
  wire        grdn  = nrm[10];
  wire        stkn  = |nrm[9:0];
  wire        rndn  = grdn & (stkn | mpre[0]);
  wire [10:0] minc  = {1'b0, mpre[9:0]} + {10'd0, rndn};
  wire        carry = minc[10];
  wire [9:0]  fracn = minc[9:0];

  wire [6:0]  enc   = {1'b0, exps} - 7'd15 - {2'b00, lz} + {6'd0, n1};
  wire [6:0]  eE    = enc + {6'd0, carry};
  wire        sub_res = eE[6] | (eE == 7'd0);
  wire        ovr     = ~eE[6] & (eE > 7'd30);

  // ---- subnormal result path: shift amount from inputs only ----
  // s = 15 - exps + lz   (subnormal result iff enc <= 0, then s >= 0)
  wire [6:0] sfull = 7'd15 + {2'b00, lz} - {1'b0, exps};
  wire       sbig  = (sfull > 7'd11);
  wire [3:0] sh    = sfull[3:0];

  wire [12:0] v13  = pp[21:9];                    // v13[i] = pp[i+9]
  wire [12:0] w13  = v13 >> sh;                   // 4-stage barrel
  wire [9:0]  f10  = sbig ? 10'd0 : w13[11:2];
  wire        gs   = sbig ? 1'b0  : w13[1];

  // sticky = |pp[s+9:0]  (mask generated from inputs only)
  wire [11:0] maskA = {12{1'b1}} >> (4'd12 - sh);
  wire        sts   = |(pp[21:10] & maskA) | |pp[9:0];

  wire        rnds  = gs & (sts | f10[0]);
  wire [10:0] finc  = {1'b0, f10} + {10'd0, rnds};
  wire        carry_s = finc[10];
  wire [9:0]  fracs = finc[9:0];

  // ---- pack ----
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_in | (ovr & ~zero_in));
  wire is_sub   = ~nan_in & ~is_inf & sub_res;
  wire sub_zero = is_sub & ~carry_s & (fracs == 10'd0);
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] ey = (is_nan | is_inf) ? 5'd31 :
                  is_zero ? 5'd0 :
                  is_sub  ? {4'd0, carry_s} : eE[4:0];
  wire [9:0] fy = is_nan ? 10'h200 :
                  (is_inf | is_zero) ? 10'd0 :
                  is_sub ? fracs : fracn;
  wire sy = is_nan ? 1'b0 : sgn;

  assign y = {sy, ey, fy};
endmodule

// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Hand-built Wallace tree + Kogge-Stone prefix adder for the 11x11 product.

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

  wire [9:0]  gn  = fsub << clz;
  wire [4:0]  lza = 5'd1 + {1'b0, clz};
  wire [4:0]  lz  = (ea0 ? lza : 5'd0) + (eb0 ? lza : 5'd0);

  wire [10:0] ma = ea0 ? {gn, 1'b0} : {1'b1, fa};
  wire [10:0] mb = eb0 ? {gn, 1'b0} : {1'b1, fb};

  // ---- Wallace tree: 11 partial products -> 2 rows ----
  wire [21:0] p0  = {11'd0, ma} & {22{mb[0]}};
  wire [21:0] p1  = ({11'd0, ma} << 1)  & {22{mb[1]}};
  wire [21:0] p2  = ({11'd0, ma} << 2)  & {22{mb[2]}};
  wire [21:0] p3  = ({11'd0, ma} << 3)  & {22{mb[3]}};
  wire [21:0] p4  = ({11'd0, ma} << 4)  & {22{mb[4]}};
  wire [21:0] p5  = ({11'd0, ma} << 5)  & {22{mb[5]}};
  wire [21:0] p6  = ({11'd0, ma} << 6)  & {22{mb[6]}};
  wire [21:0] p7  = ({11'd0, ma} << 7)  & {22{mb[7]}};
  wire [21:0] p8  = ({11'd0, ma} << 8)  & {22{mb[8]}};
  wire [21:0] p9  = ({11'd0, ma} << 9)  & {22{mb[9]}};
  wire [21:0] p10 = ({11'd0, ma} << 10) & {22{mb[10]}};

  wire [21:0] s1 = p0 ^ p1 ^ p2;
  wire [21:0] c1 = ((p0 & p1) | (p2 & (p0 ^ p1))) << 1;
  wire [21:0] s2 = p3 ^ p4 ^ p5;
  wire [21:0] c2 = ((p3 & p4) | (p5 & (p3 ^ p4))) << 1;
  wire [21:0] s3 = p6 ^ p7 ^ p8;
  wire [21:0] c3 = ((p6 & p7) | (p8 & (p6 ^ p7))) << 1;
  // rows: s1 c1 s2 c2 s3 c3 p9 p10
  wire [21:0] s4 = s1 ^ c1 ^ s2;
  wire [21:0] c4 = ((s1 & c1) | (s2 & (s1 ^ c1))) << 1;
  wire [21:0] s5 = s3 ^ c3 ^ p9;
  wire [21:0] c5 = ((s3 & c3) | (p9 & (s3 ^ c3))) << 1;
  // rows: s4 c4 s5 c5 c2 p10
  wire [21:0] s6 = s4 ^ c4 ^ s5;
  wire [21:0] c6 = ((s4 & c4) | (s5 & (s4 ^ c4))) << 1;
  wire [21:0] s7 = c5 ^ c2 ^ p10;
  wire [21:0] c7 = ((c5 & c2) | (p10 & (c5 ^ c2))) << 1;
  // rows: s6 c6 s7 c7
  wire [21:0] s8 = s6 ^ c6 ^ s7;
  wire [21:0] c8 = ((s6 & c6) | (s7 & (s6 ^ c6))) << 1;
  // rows: s8 c8 c7
  wire [21:0] S = s8 ^ c8 ^ c7;
  wire [21:0] C = ((s8 & c8) | (c7 & (s8 ^ c8))) << 1;

  // ---- Kogge-Stone prefix adder: Q = S + C ----
  wire [21:0] Pp [0:5];
  wire [21:0] Gg [0:5];
  assign Pp[0] = S ^ C;
  assign Gg[0] = S & C;
  genvar l, i;
  generate
    for (l = 1; l <= 5; l = l + 1) begin : gl
      for (i = 0; i < 22; i = i + 1) begin : gb
        if (i >= (1 << (l - 1))) begin : black
          assign Pp[l][i] = Pp[l-1][i] & Pp[l-1][i - (1 << (l - 1))];
          assign Gg[l][i] = Gg[l-1][i] | (Pp[l-1][i] & Gg[l-1][i - (1 << (l - 1))]);
        end else begin : pass
          assign Pp[l][i] = Pp[l-1][i];
          assign Gg[l][i] = Gg[l-1][i];
        end
      end
    end
  endgenerate

  wire [21:0] Q;
  assign Q[0]    = Pp[0][0];
  assign Q[21:1] = Pp[0][21:1] ^ Gg[5][20:0];

  // ---- exponent sums (input only) ----
  wire [5:0] exps = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};

  // ---- normal result path ----
  wire        n1    = Q[21];
  wire [10:0] mpre  = n1 ? Q[21:11] : Q[20:10];
  wire        grdn  = n1 ? Q[10]    : Q[9];
  wire        stkn  = n1 ? |Q[9:0]  : |Q[8:0];
  wire        rndn  = grdn & (stkn | mpre[0]);
  wire [10:0] minc  = {1'b0, mpre[9:0]} + {10'd0, rndn};
  wire        carry = minc[10];
  wire [9:0]  fracn = minc[9:0];

  wire [6:0]  enc   = {1'b0, exps} - 7'd15 - {2'b00, lz} + {6'd0, n1};
  wire [6:0]  eE    = enc + {6'd0, carry};
  wire        sub_res = eE[6] | (eE == 7'd0);
  wire        ovr     = ~eE[6] & (eE > 7'd30);

  // ---- subnormal result path: shift amount from inputs only ----
  wire [6:0] sfull = 7'd15 + {2'b00, lz} - {1'b0, exps};
  wire       sbig  = (sfull > 7'd11);
  wire [3:0] sh    = sfull[3:0];

  wire [12:0] v13  = Q[21:9];
  wire [12:0] w13  = v13 >> sh;
  wire [9:0]  f10  = sbig ? 10'd0 : w13[11:2];
  wire        gs   = sbig ? 1'b0  : w13[1];

  wire [11:0] maskA = {12{1'b1}} >> (4'd12 - sh);
  wire        sts   = |(Q[21:10] & maskA) | |Q[9:0];

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

// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Wallace/Sklansky multiply, then:
//  - normal path: left barrel by lz (input-only select)
//  - subnormal path: right barrel by h = 25-exps (input-only select)
//  - both sticky bits via input-precomputed masks (AND + OR tree, no shifter)
//  - exponent/flags computed in parallel, results merged by muxes

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

  // ---- leading zero counts (parallel with the multiplier) ----
  wire [3:0] clza = fa[9] ? 4'd0 : fa[8] ? 4'd1 : fa[7] ? 4'd2 : fa[6] ? 4'd3 :
                    fa[5] ? 4'd4 : fa[4] ? 4'd5 : fa[3] ? 4'd6 : fa[2] ? 4'd7 :
                    fa[1] ? 4'd8 : 4'd9;
  wire [3:0] clzb = fb[9] ? 4'd0 : fb[8] ? 4'd1 : fb[7] ? 4'd2 : fb[6] ? 4'd3 :
                    fb[5] ? 4'd4 : fb[4] ? 4'd5 : fb[3] ? 4'd6 : fb[2] ? 4'd7 :
                    fb[1] ? 4'd8 : 4'd9;

  wire [4:0] lz  = (ea0 ? (5'd1 + {1'b0, clza}) : 5'd0)
                 + (eb0 ? (5'd1 + {1'b0, clzb}) : 5'd0);      // 0..20
  wire [3:0] lz4 = lz[3:0];                                    // 0..15

  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};

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
  wire [21:0] s4 = s1 ^ c1 ^ s2;
  wire [21:0] c4 = ((s1 & c1) | (s2 & (s1 ^ c1))) << 1;
  wire [21:0] s5 = s3 ^ c3 ^ p9;
  wire [21:0] c5 = ((s3 & c3) | (p9 & (s3 ^ c3))) << 1;
  wire [21:0] s6 = s4 ^ c4 ^ s5;
  wire [21:0] c6 = ((s4 & c4) | (s5 & (s4 ^ c4))) << 1;
  wire [21:0] s7 = c5 ^ c2 ^ p10;
  wire [21:0] c7 = ((c5 & c2) | (p10 & (c5 ^ c2))) << 1;
  wire [21:0] s8 = s6 ^ c6 ^ s7;
  wire [21:0] c8 = ((s6 & c6) | (s7 & (s6 ^ c6))) << 1;
  wire [21:0] S = s8 ^ c8 ^ c7;
  wire [21:0] C = ((s8 & c8) | (c7 & (s8 ^ c8))) << 1;

  // ---- Sklansky prefix adder: T = S + C ----
  wire [21:0] Pp [0:5];
  wire [21:0] Gg [0:5];
  assign Pp[0] = S ^ C;
  assign Gg[0] = S & C;
  genvar l, i;
  generate
    for (l = 1; l <= 5; l = l + 1) begin : gl
      for (i = 0; i < 22; i = i + 1) begin : gb
        if ((i % (1 << l)) >= (1 << (l - 1))) begin : black
          localparam integer j = (i & ~((1 << l) - 1)) + (1 << (l - 1)) - 1;
          assign Pp[l][i] = Pp[l-1][i] & Pp[l-1][j];
          assign Gg[l][i] = Gg[l-1][i] | (Pp[l-1][i] & Gg[l-1][j]);
        end else begin : pass
          assign Pp[l][i] = Pp[l-1][i];
          assign Gg[l][i] = Gg[l-1][i];
        end
      end
    end
  endgenerate

  wire [21:0] T;
  assign T[0]    = Pp[0][0];
  assign T[21:1] = Pp[0][21:1] ^ Gg[5][20:0];

  // ---- exponent arithmetic (input only) ----
  wire [5:0] exps = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};
  wire [6:0] P    = {1'b0, exps} - 7'd15 - {2'b00, lz};

  // ---- normal result path ----
  wire [21:0] pp   = T << lz4;
  wire        n1   = pp[21];
  wire [10:0] mpre = n1 ? pp[21:11] : pp[20:10];
  wire        grdn = n1 ? pp[10]    : pp[9];
  // sticky = |pp[8:0] | (pp[9] & n1)  ->  |T[8-lz:0] | (T[9-lz] & n1)
  wire        lzge10 = (lz >= 5'd10);
  wire [3:0]  gidx   = 4'd9 - lz4;
  wire [8:0]  maskn  = lzge10 ? 9'd0 : ((9'd1 << gidx) - 9'd1);
  wire        stkn   = (|(T[8:0] & maskn)) | (lzge10 ? 1'b0 : (T[gidx] & n1));
  wire        rndn  = grdn & (stkn | mpre[0]);
  wire [10:0] incn  = {1'b0, mpre[9:0]} + 11'd1;
  wire [9:0]  fracn = rndn ? incn[9:0] : mpre[9:0];
  wire        carry = rndn & incn[10];

  wire [6:0] eE0 = P + {6'd0, n1};
  wire [6:0] eE1 = eE0 + 7'd1;
  wire        sub_res = eE0[6] | ((eE0 == 7'd0) & ~carry);
  wire        ovr     = (~eE0[6] & (eE0[5:0] > 6'd30)) | ((eE0 == 7'd30) & carry);
  wire [4:0]  eElo   = carry ? eE1[4:0] : eE0[4:0];

  // ---- subnormal result path: h = 25 - exps (input only) ----
  wire [4:0]  h  = 5'd25 - exps[4:0];
  wire [21:0] w  = T >> h;
  wire [9:0]  f10 = w[10:1];
  wire        grds = w[0];
  wire [21:0] masks = ~(22'h3FFFFF << h);
  wire        stks = |(T & masks);
  wire        rnds = grds & (stks | f10[0]);
  wire [10:0] incs = {1'b0, f10} + 11'd1;
  wire [9:0]  fracs = rnds ? incs[9:0] : f10;
  wire        carry_s = rnds & incs[10];
  wire        f10_zero = ~|f10;

  // ---- pack ----
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_in | (ovr & ~zero_in));
  wire is_sub   = ~nan_in & ~is_inf & sub_res;
  wire sub_zero = is_sub & ~rnds & f10_zero;
  wire is_zero  = ~nan_in & ~inf_in & (zero_in | sub_zero);

  wire [4:0] ey = (is_nan | is_inf) ? 5'd31 :
                  is_zero ? 5'd0 :
                  is_sub  ? {4'd0, carry_s} : eElo;
  wire [9:0] fy = is_nan ? 10'h200 :
                  (is_inf | is_zero) ? 10'd0 :
                  is_sub ? fracs : fracn;
  wire sy = is_nan ? 1'b0 : sgn;

  assign y = {sy, ey, fy};
endmodule

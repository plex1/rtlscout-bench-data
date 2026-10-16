// fp_mul_e5f10 -- variant of designN with a plain behavioural multiply
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

  wire [3:0] clza = fa[9] ? 4'd0 : fa[8] ? 4'd1 : fa[7] ? 4'd2 : fa[6] ? 4'd3 :
                    fa[5] ? 4'd4 : fa[4] ? 4'd5 : fa[3] ? 4'd6 : fa[2] ? 4'd7 :
                    fa[1] ? 4'd8 : 4'd9;
  wire [3:0] clzb = fb[9] ? 4'd0 : fb[8] ? 4'd1 : fb[7] ? 4'd2 : fb[6] ? 4'd3 :
                    fb[5] ? 4'd4 : fb[4] ? 4'd5 : fb[3] ? 4'd6 : fb[2] ? 4'd7 :
                    fb[1] ? 4'd8 : 4'd9;

  wire [4:0] lz  = (ea0 ? (5'd1 + {1'b0, clza}) : 5'd0)
                 + (eb0 ? (5'd1 + {1'b0, clzb}) : 5'd0);
  wire [3:0] lz4 = lz[3:0];

  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};

  wire [21:0] T = ma * mb;

  // ---- exponent arithmetic (input only) ----
  wire [5:0] exps = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};
  wire [6:0] P    = {1'b0, exps} - 7'd15 - {2'b00, lz};

  // trailing ones of P
  wire [6:0] q1, q2, q3, trP;
  genvar i;
  generate
    for (i = 0; i < 7; i = i + 1) begin : pq0
      assign q1[i] = (i == 0) ? 1'b1 : P[i-1];
    end
    for (i = 0; i < 7; i = i + 1) begin : pq1
      if ((i % 2) >= 1) begin : b1
        localparam integer j1 = (i & ~1);
        assign q2[i] = q1[i] & q1[j1];
      end else begin : p1
        assign q2[i] = q1[i];
      end
    end
    for (i = 0; i < 7; i = i + 1) begin : pq2
      if ((i % 4) >= 2) begin : b2
        localparam integer j2 = (i & ~3) + 1;
        assign q3[i] = q2[i] & q2[j2];
      end else begin : p2
        assign q3[i] = q2[i];
      end
    end
    for (i = 0; i < 7; i = i + 1) begin : pq3
      if ((i % 8) >= 4) begin : b3
        localparam integer j3 = (i & ~7) + 3;
        assign trP[i] = q3[i] & q3[j3];
      end else begin : p3
        assign trP[i] = q3[i];
      end
    end
  endgenerate

  // ---- normal result path ----
  wire [21:0] pp   = T << lz4;
  wire        n1   = pp[21];
  wire [10:0] mpre = n1 ? pp[21:11] : pp[20:10];
  wire        grdn = n1 ? pp[10]    : pp[9];
  wire        lzge10 = (lz >= 5'd10);
  wire [3:0]  gidx   = 4'd9 - lz4;
  wire [8:0]  maskn  = lzge10 ? 9'd0 : ((9'd1 << gidx) - 9'd1);
  wire        stkn   = (|(T[8:0] & maskn)) | (lzge10 ? 1'b0 : (T[gidx] & n1));
  wire        rndn   = grdn & (stkn | mpre[0]);

  wire [9:0] a1, a2, a3, trM;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fm0
      assign a1[i] = (i == 0) ? 1'b1 : mpre[i-1];
    end
    for (i = 0; i < 10; i = i + 1) begin : fm1
      if ((i % 2) >= 1) begin : c1b
        localparam integer k1 = (i & ~1);
        assign a2[i] = a1[i] & a1[k1];
      end else begin : c1p
        assign a2[i] = a1[i];
      end
    end
    for (i = 0; i < 10; i = i + 1) begin : fm2
      if ((i % 4) >= 2) begin : c2b
        localparam integer k2 = (i & ~3) + 1;
        assign a3[i] = a2[i] & a2[k2];
      end else begin : c2p
        assign a3[i] = a2[i];
      end
    end
    for (i = 0; i < 10; i = i + 1) begin : fm3
      if ((i % 8) >= 4) begin : c3b
        localparam integer k3 = (i & ~7) + 3;
        assign trM[i] = a3[i] & a3[k3];
      end else begin : c3p
        assign trM[i] = a3[i];
      end
    end
  endgenerate
  wire [9:0] trMx;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fm5
      if (i >= 8) begin : c5b
        assign trMx[i] = trM[i] & trM[7];
      end else begin : c5p
        assign trMx[i] = trM[i];
      end
    end
  endgenerate

  wire [9:0] fracn = mpre[9:0] ^ (trMx & {10{rndn}});
  wire       carry = rndn & &mpre[9:0];

  wire [6:0] A0 = P ^ (trP & {7{n1}});
  wire       sub_res = A0[6] | ~|A0[5:0];

  wire [6:0] b1w, b2w, b3w, trA;
  generate
    for (i = 0; i < 7; i = i + 1) begin : fa0
      assign b1w[i] = (i == 0) ? 1'b1 : A0[i-1];
    end
    for (i = 0; i < 7; i = i + 1) begin : fa1
      if ((i % 2) >= 1) begin : d1b
        localparam integer k1 = (i & ~1);
        assign b2w[i] = b1w[i] & b1w[k1];
      end else begin : d1p
        assign b2w[i] = b1w[i];
      end
    end
    for (i = 0; i < 7; i = i + 1) begin : fa2
      if ((i % 4) >= 2) begin : d2b
        localparam integer k2 = (i & ~3) + 1;
        assign b3w[i] = b2w[i] & b2w[k2];
      end else begin : d2p
        assign b3w[i] = b2w[i];
      end
    end
    for (i = 0; i < 7; i = i + 1) begin : fa3
      if ((i % 8) >= 4) begin : d3b
        localparam integer k3 = (i & ~7) + 3;
        assign trA[i] = b3w[i] & b3w[k3];
      end else begin : d3p
        assign trA[i] = b3w[i];
      end
    end
  endgenerate

  wire [6:0] eE = A0 ^ (trA & {7{carry}});
  wire       ovr = ~sub_res & (A0[5] | (&A0[4:0]) |
                               ((A0[5:0] == 6'd30) & carry));

  // ---- subnormal result path ----
  wire [4:0]  h  = 5'd25 - exps[4:0];
  wire [10:0] w  = T >> h;
  wire [9:0]  f10 = w[10:1];
  wire        grds = w[0];
  wire [21:0] masks = ~(22'h3FFFFF << h);
  wire        stks = |(T & masks);
  wire        rnds = grds & (stks | f10[0]);

  wire [9:0] e1w, e2w, e3w, trF;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fs0
      assign e1w[i] = (i == 0) ? 1'b1 : f10[i-1];
    end
    for (i = 0; i < 10; i = i + 1) begin : fs1
      if ((i % 2) >= 1) begin : g1b
        localparam integer k1 = (i & ~1);
        assign e2w[i] = e1w[i] & e1w[k1];
      end else begin : g1p
        assign e2w[i] = e1w[i];
      end
    end
    for (i = 0; i < 10; i = i + 1) begin : fs2
      if ((i % 4) >= 2) begin : g2b
        localparam integer k2 = (i & ~3) + 1;
        assign e3w[i] = e2w[i] & e2w[k2];
      end else begin : g2p
        assign e3w[i] = e2w[i];
      end
    end
    for (i = 0; i < 10; i = i + 1) begin : fs3
      if ((i % 8) >= 4) begin : g3b
        localparam integer k3 = (i & ~7) + 3;
        assign trF[i] = e3w[i] & e3w[k3];
      end else begin : g3p
        assign trF[i] = e3w[i];
      end
    end
  endgenerate
  wire [9:0] trFx;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fs4
      if (i >= 8) begin : g4b
        assign trFx[i] = trF[i] & trF[7];
      end else begin : g4p
        assign trFx[i] = trF[i];
      end
    end
  endgenerate

  wire [9:0] fracs = f10 ^ (trFx & {10{rnds}});
  wire       carry_s = rnds & &f10;

  // ---- flags and pack ----
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_in | (ovr & ~zero_in));
  wire is_sub   = ~nan_in & ~is_inf & sub_res;
  wire sub_zero = is_sub & ~rnds & ~|f10;
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

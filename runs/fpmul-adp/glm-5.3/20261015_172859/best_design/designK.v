// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Unified data path: the result fraction is always the 11-bit window
// T[s+10 : s] of the raw product, right-shifted by d = s-1 so that the
// guard lands in bit 0.  The barrel select is input-only plus n1 (=T[21],
// available early), so no product bit drives a shifter select.
//   normal case : s = 10 - lz + n1     (lz = subnormal normalisation shift)
//   subnormal   : s = 26 - exps        (subnormal grid 2^-24)
// Both cases round with one trail-XOR incrementer; the exponent is built
// with prefix-AND carry-in increments (no adders on the late path).

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

  // ---- leading zero counts of the fractions (parallel with the multiplier) ----
  wire [3:0] clza = fa[9] ? 4'd0 : fa[8] ? 4'd1 : fa[7] ? 4'd2 : fa[6] ? 4'd3 :
                    fa[5] ? 4'd4 : fa[4] ? 4'd5 : fa[3] ? 4'd6 : fa[2] ? 4'd7 :
                    fa[1] ? 4'd8 : 4'd9;
  wire [3:0] clzb = fb[9] ? 4'd0 : fb[8] ? 4'd1 : fb[7] ? 4'd2 : fb[6] ? 4'd3 :
                    fb[5] ? 4'd4 : fb[4] ? 4'd5 : fb[3] ? 4'd6 : fb[2] ? 4'd7 :
                    fb[1] ? 4'd8 : 4'd9;

  wire [4:0] lz = (ea0 ? (5'd1 + {1'b0, clza}) : 5'd0)
                + (eb0 ? (5'd1 + {1'b0, clzb}) : 5'd0);

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
  wire n1 = T[21];

  // ---- exponent arithmetic (input only) ----
  wire [5:0] exps = {1'b0, ea0 ? 5'd1 : ea} + {1'b0, eb0 ? 5'd1 : eb};
  wire [6:0] P    = {1'b0, exps} - 7'd15 - {2'b00, lz};

  // prefix-AND of P (early): trP[i] = &P[i-1:0]
  wire [6:0] q1, q2, q3;
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
  endgenerate
  wire [6:0] trP;
  generate
    for (i = 0; i < 7; i = i + 1) begin : pq3
      if ((i % 8) >= 4) begin : b3
        localparam integer j3 = (i & ~7) + 3;
        assign trP[i] = q3[i] & q3[j3];
      end else begin : p3
        assign trP[i] = q3[i];
      end
    end
  endgenerate

  // eE0 = P + n1  (biased exponent field before rounding)
  wire [6:0] eE0 = P ^ (trP & {7{n1}});
  wire sub_sel = eE0[6] | (eE0[5:0] == 6'd0);

  // ---- barrel shift amount d = s-1 (input-only + n1) ----
  wire [4:0] dsub  = 5'd25 - exps[4:0];
  wire [4:0] dnorm = 5'd9 + {1'b0, n1} - {1'b0, lz[3:0]};
  wire [4:0] d     = sub_sel ? dsub : dnorm;

  // ---- sticky-accumulating right barrel: W = T >> d ----
  wire [21:0] v1 = d[4] ? {16'd0, T[21:16]} : T;
  wire        g1 = d[4] & |T[15:0];
  wire [18:0] v2 = d[3] ? {5'd0, v1[21:8]} : v1[18:0];
  wire        g2 = g1 | (d[3] & |v1[7:0]);
  wire [14:0] v3 = d[2] ? {4'd0, v2[18:4]} : v2[14:0];
  wire        g3 = g2 | (d[2] & |v2[3:0]);
  wire [12:0] v4 = d[1] ? {2'd0, v3[14:2]} : v3[12:0];
  wire        g4 = g3 | (d[1] & |v3[1:0]);
  wire [11:0] v5 = d[0] ? {1'd0, v4[12:1]} : v4[11:0];
  wire        g5 = g4 | (d[0] & v4[0]);

  wire [10:0] F    = v5[11:1];
  wire        guard = v5[0];
  wire        sticky = g5;

  // ---- round to nearest, ties to even (trail-XOR increment) ----
  wire [9:0] x = F[9:0];
  wire [9:0] r1, r2, r3;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fr0
      assign r1[i] = (i == 0) ? 1'b1 : x[i-1];
    end
    for (i = 0; i < 10; i = i + 1) begin : fr1
      if ((i % 2) >= 1) begin : fb1
        localparam integer k1 = (i & ~1);
        assign r2[i] = r1[i] & r1[k1];
      end else begin : fp1
        assign r2[i] = r1[i];
      end
    end
    for (i = 0; i < 10; i = i + 1) begin : fr2
      if ((i % 4) >= 2) begin : fb2
        localparam integer k2 = (i & ~3) + 1;
        assign r3[i] = r2[i] & r2[k2];
      end else begin : fp2
        assign r3[i] = r2[i];
      end
    end
  endgenerate
  wire [9:0] trF;
  wire       trF10;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fr3
      if ((i % 8) >= 4) begin : fb3
        localparam integer k3 = (i & ~7) + 3;
        assign trF[i] = r3[i] & r3[k3];
      end else begin : fp3
        assign trF[i] = r3[i];
      end
    end
  endgenerate
  assign trF10 = &x;

  wire [9:0] trFx;
  generate
    for (i = 0; i < 10; i = i + 1) begin : fr4
      if (i >= 8) begin : fb4
        assign trFx[i] = trF[i] & trF[7];
      end else begin : fp4
        assign trFx[i] = trF[i];
      end
    end
  endgenerate

  wire rnd  = guard & (sticky | x[0]);
  wire [9:0] frac = x ^ (trFx & {10{rnd}});
  wire carry = rnd & trF10;

  // ---- exponent field after rounding ----
  wire [6:0] u1, u2, u3;
  generate
    for (i = 0; i < 7; i = i + 1) begin : eu0
      assign u1[i] = (i == 0) ? 1'b1 : eE0[i-1];
    end
    for (i = 0; i < 7; i = i + 1) begin : eu1
      if ((i % 2) >= 1) begin : eb1
        localparam integer m1 = (i & ~1);
        assign u2[i] = u1[i] & u1[m1];
      end else begin : ep1
        assign u2[i] = u1[i];
      end
    end
    for (i = 0; i < 7; i = i + 1) begin : eu2
      if ((i % 4) >= 2) begin : eb2
        localparam integer m2 = (i & ~3) + 1;
        assign u3[i] = u2[i] & u2[m2];
      end else begin : ep2
        assign u3[i] = u2[i];
      end
    end
  endgenerate
  wire [6:0] trE;
  generate
    for (i = 0; i < 7; i = i + 1) begin : eu3
      if ((i % 8) >= 4) begin : eb3
        localparam integer m3 = (i & ~7) + 3;
        assign trE[i] = u3[i] & u3[m3];
      end else begin : ep3
        assign trE[i] = u3[i];
      end
    end
  endgenerate

  wire [6:0] eE = eE0 ^ (trE & {7{carry}});
  wire [4:0] e_field = sub_sel ? {4'd0, carry} : eE[4:0];
  wire ovr = ~sub_sel & ((eE0[5:0] == 6'd31) | ((eE0[5:0] == 6'd30) & carry));

  // ---- flags and pack ----
  wire sub_zero = sub_sel & ~rnd & ~|x;
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_in | (ovr & ~zero_in));
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | sub_zero);

  wire [4:0] ey = (is_nan | is_inf) ? 5'd31 :
                  is_zero ? 5'd0 : e_field;
  wire [9:0] fy = is_nan ? 10'h200 :
                  (is_inf | is_zero) ? 10'd0 : frac;
  wire sy = is_nan ? 1'b0 : sgn;

  assign y = {sy, ey, fy};
endmodule

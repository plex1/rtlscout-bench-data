// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Structure: normalise the (at most one) subnormal significand before the
// multiply so the 11x11 product always lands in [2^20,2^22) -> 1-bit normalise.
// The subnormal result is derived from the normalised product with a small
// (10-bit) shifter driven by the result exponent.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  function [2:0] clz5;
    input [4:0] x;
    begin
      if      (x[4]) clz5 = 3'd0;
      else if (x[3]) clz5 = 3'd1;
      else if (x[2]) clz5 = 3'd2;
      else if (x[1]) clz5 = 3'd3;
      else           clz5 = 3'd4;
    end
  endfunction

  wire       sa = a[15];
  wire       sb = b[15];
  wire [4:0] ea = a[14:10];
  wire [4:0] eb = b[14:10];
  wire [9:0] fa = a[9:0];
  wire [9:0] fb = b[9:0];

  wire ea0 = (ea == 5'd0);
  wire eb0 = (eb == 5'd0);
  wire eaX = (ea == 5'd31);
  wire ebX = (eb == 5'd31);
  wire za  = ea0 & (fa == 10'd0);
  wire zb  = eb0 & (fb == 10'd0);
  wire na  = eaX & (fa != 10'd0);
  wire nb  = ebX & (fb != 10'd0);
  wire ia  = eaX & (fa == 10'd0);
  wire ib  = ebX & (fb == 10'd0);

  wire nan_in   = na | nb | (ia & zb) | (za & ib);
  wire inf_in   = ia | ib;
  wire zero_in  = za | zb;
  wire both_sub = ea0 & eb0;      // product is 0 or underflows to 0
  wire sgn      = sa ^ sb;

  // ---- normalise the (at most one) subnormal significand ----
  wire [9:0] fsub = ea0 ? fa : fb;
  wire [4:0] fh = fsub[9:5];
  wire [4:0] fl = fsub[4:0];
  wire       hnz = |fh;
  wire [2:0] czh = clz5(fh);
  wire [2:0] czl = clz5(fl);
  wire [3:0] clz = hnz ? {1'b0, czh} : (4'd5 + {1'b0, czl});
  wire [9:0] gn  = fsub << clz;

  wire [10:0] ma = ea0 ? {gn, 1'b0} : {1'b1, fa};
  wire [10:0] mb = eb0 ? {gn, 1'b0} : {1'b1, fb};
  wire [4:0]  lzab = (ea0 ? (5'd1 + {1'b0, clz}) : 5'd0)
                   + (eb0 ? (5'd1 + {1'b0, clz}) : 5'd0);

  // ---- significand product, leading one forced to bit 21 ----
  wire [21:0] pp  = ma * mb;
  wire        n1  = pp[21];
  wire [21:0] nrm = n1 ? pp : {pp[20:0], 1'b0};

  // ---- round the normal result ----
  wire [10:0] mpre  = nrm[21:11];
  wire        grdn  = nrm[10];
  wire        stkn  = |nrm[9:0];
  wire        rndn  = grdn & (stkn | mpre[0]);
  wire [11:0] msum  = {1'b0, mpre} + {11'd0, rndn};
  wire        carry = msum[11];
  wire [9:0]  fracn = msum[9:0];

  // ---- result exponent ----
  wire [6:0] exps = {2'd0, ea0 ? 5'd1 : ea} + {2'd0, eb0 ? 5'd1 : eb};
  wire [7:0] enc  = {1'b0, exps} - 8'd15 - {3'd0, lzab} + {7'd0, n1};
  wire [7:0] eE   = enc + {7'd0, carry};
  wire       sub_res = eE[7] | (eE == 8'd0);
  wire       ovr     = ~eE[7] & (eE > 8'd30);

  // ---- subnormal result path ----
  wire [7:0] negc = (~enc) + 8'd1;
  wire [4:0] dd   = enc[7] ? negc[4:0] : 5'd0;
  wire [9:0]  v   = nrm[21:12];
  wire [11:0] low = nrm[11:0];
  wire [9:0]  f10   = v >> dd;
  wire [9:0]  vsub  = v >> (dd - 5'd1);
  wire [9:0]  vleft = v << (5'd11 - dd);
  wire       grds = (dd == 5'd0) ? low[11] : vsub[0];
  wire       stks = (dd == 5'd0) ? |low[10:0] : (|low | |vleft);
  wire       rnds = grds & (stks | f10[0]);
  wire [10:0] fsum = {1'b0, f10} + {10'd0, rnds};
  wire        carry_s = fsum[10];
  wire [9:0]  fracs = fsum[9:0];

  // ---- pack ----
  wire is_nan   = nan_in;
  wire is_inf   = ~nan_in & (inf_in | (ovr & ~zero_in));
  wire is_sub   = ~nan_in & ~is_inf & sub_res;
  wire sub_zero = is_sub & ~carry_s & (fracs == 10'd0);
  wire is_zero  = ~nan_in & ~is_inf & (zero_in | both_sub | sub_zero);

  wire [4:0] ey = (is_nan | is_inf) ? 5'd31 :
                  is_zero ? 5'd0 :
                  is_sub  ? {4'd0, carry_s} : eE[4:0];
  wire [9:0] fy = is_nan ? 10'h200 :
                  (is_inf | is_zero) ? 10'd0 :
                  is_sub ? fracs : fracn;
  wire sy = is_nan ? 1'b0 : sgn;

  assign y = {sy, ey, fy};
endmodule

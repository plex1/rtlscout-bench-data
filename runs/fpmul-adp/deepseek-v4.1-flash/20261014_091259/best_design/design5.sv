// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction).
//
// Operand pre-normalisation: subnormal significands are shifted up so their
// leading one sits at bit 10, and their exponent is adjusted accordingly.  The
// 11x11 product of two such significands is then 22 bits with the leading one
// at bit 20 or 21, so a single small barrel shifter is enough:
//
//   result field = sec + carry - 15      (sec = pseudo exponent sum, biased
//                                         by 32 and saturating to 0 for a
//                                         zero operand so a signed zero falls
//                                         out of the arithmetic)
//   normal results : shift p right by 9+carry  (r4 = p[21:9+... ])
//   subnormal      : shift right so the mantissa lands on the 2^-24 grid,
//                    with the sticky bit gathered by the shifter stages.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a, b;
  output [15:0] y;

  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire za = ~|ea, zb = ~|eb;
  wire fza = za & ~|fa, fzb = zb & ~|fb;
  wire ia = (ea == 5'h1F), ib = (eb == 5'h1F);
  wire infa = ia & ~|fa, infb = ib & ~|fb;
  wire nana = ia & (|fa), nanb = ib & (|fb);
  wire sy    = sa ^ sb;
  wire isnan = nana | nanb | (infa & fzb) | (fza & infb);

  // left shift moving the leading one of a 10-bit significand to bit 10
  function [3:0] smt;
    input [9:0] f;
    if (|f[9:5]) begin
      if (|f[9:7]) begin
        if (f[9]) smt = 4'd1; else if (f[8]) smt = 4'd2; else smt = 4'd3;
      end else begin
        if (f[6]) smt = 4'd4; else smt = 4'd5;
      end
    end else begin
      if (|f[4:2]) begin
        if (f[4]) smt = 4'd6; else if (f[3]) smt = 4'd7; else smt = 4'd8;
      end else begin
        if (f[1]) smt = 4'd9; else smt = 4'd10;
      end
    end
  endfunction

  wire [3:0]  sha  = za ? smt(fa) : 4'd0;
  wire [3:0]  shb  = zb ? smt(fb) : 4'd0;
  wire [10:0] siga = {~za, fa} << sha;
  wire [10:0] sigb = {~zb, fb} << shb;

  // pseudo exponent sum: sa = ea+32 (normal), 33-sha (subnormal), 0 (zero)
  wire [6:0] Aa = fza ? 7'd0 : (za ? (7'd33 - {3'b0, sha}) : ({2'b0, ea} + 7'd32));
  wire [6:0] Ab = fzb ? 7'd0 : (zb ? (7'd33 - {3'b0, shb}) : ({2'b0, eb} + 7'd32));
  wire [7:0] Sec = {1'b0, Aa} + {1'b0, Ab};

  wire [21:0] p     = siga * sigb;
  wire        carry = p[21];

  wire [7:0]  S   = Sec + {7'b0, carry};
  wire        sub = (S <= 8'd79);          // result exponent <= 0
  wire [7:0]  shf = sub ? (8'd80 - S) : {7'b0, carry};
  wire        ge16 = |shf[7:4];

  // barrel shift of q = p[21:9]; sticky gathered stage by stage
  wire [12:0] q  = p[21:9];
  wire [12:0] r1 = shf[0] ? {1'b0, q[12:1]}   : q;   wire o1 = shf[0] &  q[0];
  wire [12:0] r2 = shf[1] ? {2'b0, r1[12:2]}  : r1;  wire o2 = shf[1] & (|r1[1:0]);
  wire [12:0] r3 = shf[2] ? {4'b0, r2[12:4]}  : r2;  wire o3 = shf[2] & (|r2[3:0]);
  wire [12:0] r4 = shf[3] ? {8'b0, r3[12:8]}  : r3;  wire o4 = shf[3] & (|r3[7:0]);
  wire [12:0] r5 = ge16 ? 13'd0 : r4;
  wire        o5 = ge16 & (|r4);

  wire [10:0] mant   = r5[11:1];
  wire        guard  = r5[0];
  wire        sticky = o1 | o2 | o3 | o4 | o5 | (|p[8:0]);
  wire        rup    = guard & (sticky | mant[0]);
  wire [11:0] M      = {1'b0, mant} + {11'b0, rup};

  wire [7:0]  fld  = sub ? {7'b0, M[10]} : (S - 8'd79 + {7'b0, M[11]});
  wire        isinf = (fld > 8'd30) | infa | infb | isnan;

  assign y = {sy & ~isnan,
              isinf  ? 5'h1F   : fld[4:0],
              isnan  ? 10'h200 : (isinf ? 10'd0 : M[9:0])};
endmodule

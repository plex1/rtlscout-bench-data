// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction).
//
// Subnormal significands are pre-normalised (shifted up so the leading one is
// at bit 10, exponent adjusted); a zero operand gets a pseudo exponent that
// makes the arithmetic itself produce a signed zero.  The 11x11 product p of
// two values in [1,2) is in [2^20, 2^22), so after removing the 9 low bits
// (q = p[21:9], whose bits [8:0] feed the sticky) only a small right shift is
// needed.  One round-to-nearest-even then feeds the packer.  The exponents are
// carried biased by +16, so sec = e_a + e_b_overhead + 32 and
//   result field (normal)  = sec + carry - 47 + rounding carry
//   subnormal              = sec + carry <= 47
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

  // left shift that moves the leading one of a 10-bit significand to bit 10
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

  // biased exponents (+16): normal = ea+16, subnormal = 17-sha, zero = 1
  wire [6:0] sea = fza ? 7'd1 : (za ? (7'd17 - {3'b0, sha}) : ({2'b0, ea} + 7'd16));
  wire [6:0] seb = fzb ? 7'd1 : (zb ? (7'd17 - {3'b0, shb}) : ({2'b0, eb} + 7'd16));
  wire [6:0] sec = sea + seb;

  wire [21:0] p     = siga * sigb;
  wire        carry = p[21];
  wire [7:0]  s     = {1'b0, sec} + {7'b0, carry};
  wire        sub   = ~(s > 8'd47);            // result exponent <= -15

  // right shift of q = p[21:9]; sticky gathered stage by stage
  wire [7:0]  dz  = 8'd48 - {1'b0, sec};
  wire [3:0]  sh  = sub ? ((dz > 8'd15) ? 4'd15 : dz[3:0]) : {3'b0, carry};
  wire [12:0] q   = p[21:9];
  wire [12:0] r1 = sh[0] ? {1'b0, q[12:1]}   : q;   wire o1 = sh[0] &  q[0];
  wire [12:0] r2 = sh[1] ? {2'b0, r1[12:2]}  : r1;  wire o2 = sh[1] & (|r1[1:0]);
  wire [12:0] r3 = sh[2] ? {4'b0, r2[12:4]}  : r2;  wire o3 = sh[2] & (|r2[3:0]);
  wire [12:0] r4 = sh[3] ? {8'b0, r3[12:8]}  : r3;  wire o4 = sh[3] & (|r3[7:0]);

  wire [10:0] mant   = r4[11:1];
  wire        guard  = r4[0];
  wire        sticky = o1 | o2 | o3 | o4 | (|p[8:0]);
  wire        rup    = guard & (sticky | mant[0]);
  wire [11:0] M      = {1'b0, mant} + {11'b0, rup};

  wire [7:0]  fld   = sub ? {7'b0, M[10]} : (s - 8'd47 + {7'b0, M[11]});
  wire        isinf = (fld > 8'd30) | infa | infb;

  assign y = {sy & ~isnan,
              isnan  ? 5'h1F   : (isinf ? 5'h1F : fld[4:0]),
              isnan  ? 10'h200 : (isinf ? 10'd0 : M[9:0])};
endmodule

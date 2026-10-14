// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction).
//
// Operand pre-normalisation (subnormals shifted so the leading significand bit
// sits at bit 10) -> 11x11 multiply -> small right-shift + sticky -> round to
// nearest even -> clamp/pack.  Zero, inf and overflow are handled by the
// exponent arithmetic itself so the output needs only one mux level.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a, b;
  output [15:0] y;

  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire za  = ~|ea;
  wire zb  = ~|eb;
  wire fza = za & ~|fa;
  wire fzb = zb & ~|fb;
  wire ia  = (ea == 5'h1F);
  wire ib  = (eb == 5'h1F);
  wire infa = ia & ~|fa;
  wire infb = ib & ~|fb;
  wire nana = ia & (|fa);
  wire nanb = ib & (|fb);
  wire sy     = sa ^ sb;
  wire isnan  = nana | nanb | (infa & fzb) | (fza & infb);
  wire anyinf = infa | infb;

  // left-shift amount that brings the leading one of fa to bit 10
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

  // Exponents, biased by +16 so everything stays unsigned:
  //   normal operand    -> ea + 16      (17 .. 46)
  //   subnormal operand -> 17 - sha     ( 7 .. 16)
  //   zero operand      -> 1            (forces a zero result through the
  //                                      exponent/round arithmetic alone)
  wire [6:0] sea = fza ? 7'd1 : (za ? (7'd17 - {3'b0, sha}) : {2'b0, ea} + 7'd16);
  wire [6:0] seb = fzb ? 7'd1 : (zb ? (7'd17 - {3'b0, shb}) : {2'b0, eb} + 7'd16);
  wire [6:0] sec = sea + seb;

  wire [21:0] p     = siga * sigb;
  wire        carry = p[21];

  // base exponent of a normal result = esum - 15 ; e = right shift of the
  // 13-bit window p[21:9] that maps the leading one into bit 10
  wire [6:0] base = (sec >= 7'd47) ? (sec - 7'd47) : 7'd0;
  // e = right shift of the 13-bit window q = p[21:9] such that the result
  // mantissa is formed from r4[11:1];  e = max(carry, 48 - sec)
  wire [6:0] tsub = 7'd48 - sec;
  wire [3:0] e    = (sec >= 7'd48) ? {3'b0, carry} :
                    (tsub > 7'd15) ? 4'd15 : tsub[3:0];

  wire [12:0] q  = p[21:9];
  wire [12:0] r1 = e[0] ? {1'b0,  q[12:1]} : q;
  wire        u1 = e[0] & q[0];
  wire [12:0] r2 = e[1] ? {2'b0,  r1[12:2]} : r1;
  wire        u2 = e[1] & (|r1[1:0]);
  wire [12:0] r3 = e[2] ? {4'b0,  r2[12:4]} : r2;
  wire        u3 = e[2] & (|r2[3:0]);
  wire [12:0] r4 = e[3] ? {8'b0,  r3[12:8]} : r3;
  wire        u4 = e[3] & (|r3[7:0]);

  wire [10:0] M      = r4[11:1];
  wire        guard  = r4[0];
  wire        sticky = u1 | u2 | u3 | u4 | (|p[8:0]);
  wire        rup    = guard & (sticky | M[0]);
  wire [10:0] F      = M[9:0] + {10'd0, rup};

  wire [6:0] expf = base + {6'b0, F[10]};
  wire ovf  = (expf > 7'sd30);
  wire spec = ovf | anyinf | isnan;

  assign y = {sy & ~isnan,
              spec  ? 5'h1F   : expf[4:0],
              isnan ? 10'h200 : (spec ? 10'd0 : F[9:0])};
endmodule

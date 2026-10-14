// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction).
//
// Operands are pre-normalised (subnormal significands are left shifted so the
// leading one sits at bit 10), then multiplied as 11x11 -> 22 bits.  Because
// the product of two normalised significands is always in [2^20, 2^22) the
// remaining normalisation is at most a 1-bit shift, so the only barrel shifter
// is the one that maps the product onto the subnormal grid.
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
  wire sy = sa ^ sb;
  wire isnan = nana | nanb | (infa & fzb) | (fza & infb);

  // left shift that moves the leading one of a 10 bit significand to bit 10
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

  // biased (+16) exponents; a zero operand gets 1 so the packer emits a zero
  wire [6:0] sea = fza ? 7'd1 : (za ? (7'd17 - {3'b0, sha}) : ({2'b0, ea} + 7'd16));
  wire [6:0] seb = fzb ? 7'd1 : (zb ? (7'd17 - {3'b0, shb}) : ({2'b0, eb} + 7'd16));
  wire [6:0] sec = sea + seb;

  wire [21:0] p     = siga * sigb;
  wire        carry = p[21];
  // product value = p * 2^(sec-82), unbiased result exponent = sec-62+carry
  wire [7:0]  s     = {1'b0, sec} + {7'b0, carry};
  wire        verge = ~(s > 8'd47);   // s<=47 -> result exponent <= 0 (subnormal)

  // ---- normal path: round p at bit (10+carry) ----
  wire [11:0] r   = carry ? p[21:10] : p[20:9];
  wire        stn = carry ? (|p[8:0]) : (|p[7:0]);
  wire [10:0] Mn  = r[11:1] + (r[0] & (stn | r[1]));   // in [1024,2048]
  wire [9:0]  fn  = Mn[9:0];
  wire [7:0]  En  = s - 8'd47 + {7'b0, Mn[11]};

  // ---- subnormal path: A = round(p / 2^(58-sec)) ----
  wire [7:0]  dd = 8'd48 - {1'b0, sec};      // shift of u = p[21:10]
  wire [5:0]  d1 = dd[5:0];
  wire [11:0] u  = p[21:10];
  wire        s0 = |p[9:0];
  wire [11:0] u1 = d1[0] ? {1'b0,  u[11:1]}  : u;   wire t1 = d1[0] & u[0];
  wire [11:0] u2 = d1[1] ? {2'b0,  u1[11:2]} : u1;  wire t2 = d1[1] & (|u1[1:0]);
  wire [11:0] u3 = d1[2] ? {4'b0,  u2[11:4]} : u2;  wire t3 = d1[2] & (|u2[3:0]);
  wire [11:0] u4 = d1[3] ? {8'b0,  u3[11:8]} : u3;  wire t4 = d1[3] & (|u3[7:0]);
  wire [11:0] u5 = d1[4] ? 12'b0 : u4;              wire t5 = d1[4] & (|u4);
  wire [11:0] u6 = d1[5] ? 12'b0 : u5;              wire t6 = d1[5] & (|u5);

  wire        stv = s0 | t1 | t2 | t3 | t4 | t5 | t6;
  wire [10:0] A   = u6[11:1] + (u6[0] & (stv | u6[1]));   // in [512,1024]
  wire [9:0]  fv  = A[9:0];
  wire        Ev  = A[10];

  // ---- pack ----
  wire [7:0]  E    = verge ? {7'b0, Ev} : En;
  wire [9:0]  F    = verge ? fv : fn;
  wire        isinf = (E > 8'd30) | infa | infb;
  wire [4:0]  expf  = isinf ? 5'h1F : E[4:0];
  wire [9:0]  fracf = isnan ? 10'h200 : (isinf ? 10'd0 : F);

  assign y = {sy & ~isnan, expf, fracf};
endmodule

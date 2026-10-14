// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction).
//
// Structure: operand pre-normalisation (subnormals shifted so the leading
// significand bit is at position 10) -> 11x11 multiply -> single small
// right-shift + sticky generation -> round-to-nearest-even -> pack.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a, b;
  output [15:0] y;

  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire za = ~|ea;
  wire zb = ~|eb;
  wire fza = za & ~|fa;
  wire fzb = zb & ~|fb;
  wire infa = (ea == 5'h1F) & ~|fa;
  wire infb = (eb == 5'h1F) & ~|fb;
  wire nana = (ea == 5'h1F) & (|fa);
  wire nanb = (eb == 5'h1F) & (|fb);
  wire sy   = sa ^ sb;

  // ---- leading-one position (0..9) ----
  function [3:0] leadpos;
    input [9:0] f;
    integer i;
    reg [3:0] p;
    begin
      p = 4'd0;
      for (i = 0; i < 10; i = i + 1) if (f[i]) p = i[3:0];
      leadpos = p;
    end
  endfunction

  // ---- normalise significands ----
  wire [3:0]  sha = za ? (4'd10 - leadpos(fa)) : 4'd0;
  wire [3:0]  shb = zb ? (4'd10 - leadpos(fb)) : 4'd0;
  wire [10:0] siga = {~za, fa} << sha;
  wire [10:0] sigb = {~zb, fb} << shb;

  wire signed [5:0] esa = za ? (6'sd1 - {2'b0, sha}) : {1'b0, ea};
  wire signed [5:0] esb = zb ? (6'sd1 - {2'b0, shb}) : {1'b0, eb};

  // ---- significand product ----
  wire [21:0] p = siga * sigb;

  // ---- exponent ----
  wire signed [6:0] esum  = esa + esb;
  wire              carry = p[21];
  wire signed [6:0] eres  = esum - 7'sd15 + $signed({6'b0, carry});
  wire [6:0]        base  = (eres > 0) ? eres : 7'd0;

  wire signed [6:0] cc = $signed({6'b0, carry});
  wire signed [6:0] t  = 7'sd16 - esum;
  wire signed [6:0] mx = (t > cc) ? t : cc;
  wire [3:0]        e  = (mx > 7'sd15) ? 4'd15 : mx[3:0];

  // ---- shift by e, collect sticky ----
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

  wire [10:0] F    = M[9:0] + {10'd0, rup};
  wire [9:0]  frac = F[9:0];
  wire [6:0]  expf = base + {6'd0, F[10]};
  wire        ovf  = (expf > 7'd30);

  assign y = (nana | nanb | (infa & fzb) | (fza & infb)) ? 16'h7E00 :
             (infa | infb) ? {sy, 5'h1F, 10'h0}                 :
             (fza  | fzb ) ? {sy, 5'h00, 10'h0}                 :
             ovf           ? {sy, 5'h1F, 10'h0}                 :
                             {sy, expf[4:0], frac};
endmodule

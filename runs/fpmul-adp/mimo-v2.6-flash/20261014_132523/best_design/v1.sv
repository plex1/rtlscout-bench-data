// fp_mul_e5f10 -- IEEE binary16 multiplier (1/5/10) with subnormals.
// Architecture: pre-normalise operand significands (LZC + left shift),
// 11x11 multiply, 0/1-bit normalise, round-to-nearest-even.
// Subnormal results share the same product with a single right shift.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- classify ----
  wire        sa = a[15],    sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  wire ea0 = (ea == 5'd0),  eb0 = (eb == 5'd0);
  wire eap = (ea == 5'd31), ebp = (eb == 5'd31);

  wire zero_a = ea0 & (fa == 10'd0), zero_b = eb0 & (fb == 10'd0);
  wire inf_a  = eap & (fa == 10'd0), inf_b  = ebp & (fb == 10'd0);
  wire nan_in = (eap & (fa != 10'd0)) | (ebp & (fb != 10'd0))
              | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ---- significands + pre-normalisation ----
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};

  wire [3:0] za = ma[10] ? 4'd0  : ma[9] ? 4'd1  : ma[8] ? 4'd2 : ma[7] ? 4'd3 :
                  ma[6]  ? 4'd4  : ma[5] ? 4'd5  : ma[4] ? 4'd6 : ma[3] ? 4'd7 :
                  ma[2]  ? 4'd8  : ma[1] ? 4'd9  : ma[0] ? 4'd10 : 4'd11;
  wire [3:0] zb = mb[10] ? 4'd0  : mb[9] ? 4'd1  : mb[8] ? 4'd2 : mb[7] ? 4'd3 :
                  mb[6]  ? 4'd4  : mb[5] ? 4'd5  : mb[4] ? 4'd6 : mb[3] ? 4'd7 :
                  mb[2]  ? 4'd8  : mb[1] ? 4'd9  : mb[0] ? 4'd10 : 4'd11;

  reg [10:0] man, mbn;
  integer k;
  always @* begin
    man = 11'd0;
    for (k = 0; k <= 10; k = k + 1)
      if (za == k[3:0]) man = ma << k;
  end
  always @* begin
    mbn = 11'd0;
    for (k = 0; k <= 10; k = k + 1)
      if (zb == k[3:0]) mbn = mb << k;
  end

  wire signed [6:0] ean  = (ea0 ? 7'sd1 : $signed({2'b00, ea})) - $signed({3'b000, za});
  wire signed [6:0] ebn  = (eb0 ? 7'sd1 : $signed({2'b00, eb})) - $signed({3'b000, zb});
  wire signed [6:0] esum = ean + ebn;

  // ---- product, 1-bit normalisation ----
  wire [21:0] prod = man * mbn;
  wire        s    = prod[21];
  wire [20:0] w    = s ? prod[21:1] : prod[20:0];   // w[20] = 1 for nonzero
  wire        tail = s & prod[0];                    // bit lost by the s-shift -> sticky

  // E = exponent of normalised product result (value = (w/2^20) * 2^(E-15))
  wire signed [6:0] eprod = esum + $signed({6'b000000, s}) - 7'sd15;
  wire signed [6:0] j     = -eprod;        // right shift for subnormal result
  wire [3:0]        jj    = j[3:0];
  wire juse = ~j[6] & (j <= 7'sd10);       // 0 <= j <= 10

  // ---- normal result: round w at 11-bit precision ----
  wire [10:0] rpre   = w[20:10];
  wire        gn     = w[9];
  wire        stn    = (|w[8:0]) | tail;
  wire        rup    = gn & (stn | rpre[0]);
  wire [11:0] rres   = {1'b0, rpre} + {11'd0, rup};
  wire        cr     = rres[11];                     // rounded up to 2^11
  wire signed [6:0] er = eprod + $signed({6'b000000, cr});
  wire ovf = ~er[6] & (er > 7'sd30);                 // er >= 31 -> overflow

  // ---- subnormal result: frac = round(w / 2^(11+j)) ----
  reg [20:0] z;
  always @* begin
    z = 21'd0;
    if (juse)
      for (k = 0; k <= 10; k = k + 1)
        if (j[3:0] == k[3:0]) z = w >> k;
  end
  wire [9:0]  t10  = z[20:11];
  wire        gs   = z[10];
  // sticky = |w[9+j:0] = |w[9:0] | (|(w[19:10] & thermo)), thermo[t] = (j >= t+1)
  wire [9:0]  th   = {jj >= 4'd10, jj >= 4'd9, jj >= 4'd8, jj >= 4'd7, jj >= 4'd6,
                      jj >= 4'd5,  jj >= 4'd4, jj >= 4'd3, jj >= 4'd2, jj >= 4'd1};
  wire        sts  = (|w[9:0]) | (|(w[19:10] & th)) | tail;
  wire        rs   = gs & (sts | t10[0]);
  wire [10:0] fres = {1'b0, t10} + {10'd0, rs};
  wire        cs   = fres[10];                       // rounded up to 2^10

  // ---- pack ----
  wire sub     = ~j[6];                              // E <= 0
  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (~sub & ovf & ~zero_in));
  wire is_zero = ~is_nan & ~is_inf & (zero_in | (sub & ~cs & (fres[9:0] == 10'd0)));
  wire is_sub  = ~is_nan & ~is_inf & ~is_zero & sub;

  assign y = is_nan  ? 16'h7E00 :
             is_inf  ? {sy, 5'd31, 10'd0} :
             is_zero ? {sy, 5'd0,  10'd0} :
             is_sub  ? {sy, 4'd0, cs, fres[9:0]} :
                       {sy, er[4:0], rres[9:0]};
endmodule

// v3: v2 with a leaner exponent/pack stage
module fp_mul_e5f10 (
  input  wire [15:0] a,
  input  wire [15:0] b,
  output wire [15:0] y
);
  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];

  wire ea0 = ~(|ea);
  wire eb0 = ~(|eb);
  wire ea1 = &ea;
  wire eb1 = &eb;
  wire fa0 = ~(|fa);
  wire fb0 = ~(|fb);
  wire nan_a = ea1 & ~fa0;
  wire nan_b = eb1 & ~fb0;
  wire inf_a = ea1 & fa0;
  wire inf_b = eb1 & fb0;
  wire zero_a = ea0 & fa0;
  wire zero_b = eb0 & fb0;
  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [21:0] P  = ma * mb;

  reg [4:0] lz;
  integer i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (P[i]) lz = 5'd21 - i[4:0];
  end
  wire [4:0] msb = 5'd21 - lz;

  wire [4:0] eae = ea0 ? 5'd1 : ea;
  wire [4:0] ebe = eb0 ? 5'd1 : eb;
  wire [5:0] S   = {1'b0, eae} + {1'b0, ebe};

  wire [6:0] Nsum   = {1'b0, S} + {2'b0, msb};
  wire       is_norm = (Nsum >= 7'd36);

  wire [4:0] k = is_norm ? (msb - 5'd10) : (5'd26 - S[4:0]);

  wire [24:0] PX = {3'b0, P};
  wire [4:0]  d  = (k == 5'd0) ? 5'd0 : (k - 5'd1);
  wire [24:0] Q  = PX >> d;
  wire [24:0] msk = (25'd1 << d) - 1;
  wire sticky = |(PX & msk);

  wire [10:0] M      = (k == 5'd0) ? P[10:0] : Q[11:1];
  wire        guard  = (k != 5'd0) & Q[0];
  wire        rnd    = guard & (sticky | M[0]);

  wire [11:0] Msum  = {1'b0, M} + {11'd0, rnd};
  wire        carry = Msum[11];

  // exponent field of a normal result: S + msb - 35 + carry  (signed, 7 bits)
  wire signed [6:0] En = {1'b0, S} + {2'b0, msb} - 7'sd35 + {6'd0, carry};
  wire       ovf = (En > 7'sd30);

  wire [4:0] e_fin = is_norm ? (zero_in ? 5'd0 : En[4:0]) : {4'd0, Msum[10]};
  wire [9:0] f_fin = Msum[9:0];

  wire sat = nan_in | inf_in | ovf;
  wire [4:0] e_y = sat ? 5'h1f : e_fin;
  wire [9:0] f_y = nan_in ? 10'h200 : (sat ? 10'd0 : f_fin);
  wire       s_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {s_y, e_y, f_y};
endmodule

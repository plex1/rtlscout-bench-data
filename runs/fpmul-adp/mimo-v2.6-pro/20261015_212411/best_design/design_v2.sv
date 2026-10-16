// fp_mul_e5f10 -- binary16 FP multiplier (1 sign, 5 exp, 10 frac), subnormals,
// round-to-nearest ties-to-even, canonical qNaN 16'h7E00.
//
// v2: all control arithmetic is done on the INPUT side (parallel to the
// multiply).  Post-multiply logic is only:
//   d      = ~prod[21-t]            (one bit of the product)
//   delta  = gt & d                 (1 gate)
//   the funnel shifter with EARLY control (S_pre), plus one late 1-bit stage
//   the rounding incrementer and 2 small output muxes
//
// Math:  lz = la + lb + d  (d in {0,1});  E_field = s - lz - 14
//        sh = max(11 - lz, 26 - s) = S_base - delta
//          S_base = gt ? (11-t) : (26-s),  delta = gt & d
//        normal result iff E_field >= 1;  same round formula for both paths:
//        round = guard & (sticky | lsb)

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ------------------------------------------------------------- unpack
  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'h1F);
  wire e_b_ones = (e_b == 5'h1F);
  wire zero_a = e_a_zero & (f_a == 10'd0);
  wire zero_b = e_b_zero & (f_b == 10'd0);
  wire nan_a  = e_a_ones & (f_a != 10'd0);
  wire nan_b  = e_b_ones & (f_b != 10'd0);
  wire inf_a  = e_a_ones & (f_a == 10'd0);
  wire inf_b  = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a = {~e_a_zero, f_a};
  wire [10:0] m_b = {~e_b_zero, f_b};
  wire [5:0] ea_eff = e_a_zero ? 6'd1 : {1'b0, e_a};
  wire [5:0] eb_eff = e_b_zero ? 6'd1 : {1'b0, e_b};

  wire [21:0] prod = m_a * m_b;

  // product of two subnormals is always < half of the min subnormal -> signed zero
  wire both_sub = (e_a_zero & (|f_a)) & (e_b_zero & (|f_b));

  // ------------------------------------------- operand CLZs (input side)
  function [3:0] clz11;
    input [10:0] v;
    begin
      if      (|v[10:8]) clz11 = v[10] ? 4'd0 : v[9] ? 4'd1 : 4'd2;
      else if (|v[7:5])  clz11 = v[7]  ? 4'd3 : v[6] ? 4'd4 : 4'd5;
      else if (|v[4:2])  clz11 = v[4]  ? 4'd6 : v[3] ? 4'd7 : 4'd8;
      else               clz11 = v[1]  ? 4'd9 : v[0] ? 4'd10 : 4'd11;
    end
  endfunction
  wire [3:0] la = clz11(m_a);
  wire [3:0] lb = clz11(m_b);

  // ------------------------------------------- input-side combined terms
  wire [4:0] t     = {1'b0, la} + {1'b0, lb};        // 0..22
  wire [5:0] s_sum = ea_eff + eb_eff;                // 2..60
  wire [5:0] pa    = ea_eff - {2'b0, la} + 6'd16;    // 6..47
  wire [5:0] pb    = eb_eff - {2'b0, lb} + 6'd16;
  wire [6:0] u0    = {1'b0, pa} + {1'b0, pb};        // = s - t + 32, 12..94

  wire gt = (u0 >= 7'd48);      // E_field >= 2 : normal for both d
  wire eq = (u0 == 7'd47);      // E_field == 1 : normal iff d == 0

  wire [4:0] Am1 = 5'd10 - t;          // (11 - t) - 1
  wire [4:0] Bm1 = 5'd25 - s_sum[4:0]; // (26 - s) - 1

  // ------------------------- one product bit: exact lz = t + d, d = ~prod[k]
  wire [4:0] k = (t >= 5'd21) ? 5'd0 : (5'd21 - t);
  wire [21:0] bp = prod;
  wire [11:0] v1;
  wire [7:0]  v2;
  wire [3:0]  v3;
  wire [1:0]  v4;
  wire        v5;
  genvar gi;
  generate
    for (gi = 0; gi < 11; gi = gi + 1) begin : g0
      assign v1[gi] = k[0] ? bp[2*gi+1] : bp[2*gi];
    end
    assign v1[11] = 1'b0;
    for (gi = 0; gi < 6; gi = gi + 1) begin : g1
      assign v2[gi] = k[1] ? v1[2*gi+1] : v1[2*gi];
    end
    assign v2[7:6] = 2'b00;
    for (gi = 0; gi < 3; gi = gi + 1) begin : g2
      assign v3[gi] = k[2] ? v2[2*gi+1] : v2[2*gi];
    end
    assign v3[3] = 1'b0;
    for (gi = 0; gi < 2; gi = gi + 1) begin : g3
      assign v4[gi] = k[3] ? v3[2*gi+1] : v3[2*gi];
    end
    assign v5 = k[4] ? v4[1] : v4[0];
  endgenerate
  wire d = ~v5;

  // ------------------------------------------------- shift control
  wire [4:0] S_pre = gt ? Am1 : Bm1;   // = S_base - 1
  wire       delta = gt & d;
  wire       L     = ~delta;           // late 1-bit extra shift

  // -------------------------- funnel shifter (early control) + sticky
  wire [22:0] w   = {prod, 1'b0};
  wire [22:0] qa0 = S_pre[0] ? {1'b0,  w[22:1]}   : w;
  wire [22:0] qa1 = S_pre[1] ? {2'b00, qa0[22:2]} : qa0;
  wire [22:0] qa2 = S_pre[2] ? {4'b0,  qa1[22:4]} : qa1;
  wire [22:0] qa3 = S_pre[3] ? {8'b0,  qa2[22:8]} : qa2;
  wire [22:0] qa4 = S_pre[4] ? {16'b0, qa3[22:16]}: qa3;

  wire sk1 = S_pre[1] & (|qa0[1:0]);
  wire sk2 = sk1 | (S_pre[2] & (|qa1[3:0]));
  wire sk3 = sk2 | (S_pre[3] & (|qa2[7:0]));
  wire sk4 = sk3 | (S_pre[4] & (|qa3[15:0]));

  wire [22:0] q  = L ? {1'b0, qa4[22:1]} : qa4;
  wire        sk = sk4 | (L & qa4[0]);

  // --------------------------------------------- shared round / significand
  wire [10:0] mant_pre = q[11:1];
  wire        guard    = q[0];
  wire        rnd      = guard & (sk | mant_pre[0]);
  wire [11:0] sum      = {1'b0, mant_pre} + {11'd0, rnd};
  wire        carry    = sum[11];

  // ------------------------------------- exponent candidates (input side)
  wire [6:0] E0v = u0 - 7'd46;
  wire [6:0] EM  = u0 - 7'd47;
  wire [6:0] EP  = u0 - 7'd45;
  wire ov0 = (E0v > 7'd30);
  wire ovM = (EM  > 7'd30);
  wire ovP = (EP  > 7'd30);

  wire selX = d ^ carry;
  wire [4:0] exp_norm = selX ? (carry ? EP[4:0] : EM[4:0]) : E0v[4:0];
  wire       ovf      = selX ? (carry ? ovP     : ovM    ) : ov0;
  wire       is_norm  = gt | (eq & ~d);
  wire [4:0] exp_sub  = sum[10] ? 5'd1 : 5'd0;
  wire [9:0] frac_y   = sum[9:0];

  // ------------------------------------------------------------- pack
  assign y = nan_in   ? 16'h7E00 :
             inf_in   ? {s_y, 5'h1F, 10'h0} :
             zero_in  ? {s_y, 15'h0} :
             both_sub ? {s_y, 15'h0} :
             is_norm  ? (ovf ? {s_y, 5'h1F, 10'h0} : {s_y, exp_norm, frac_y}) :
                        {s_y, exp_sub, frac_y};
endmodule

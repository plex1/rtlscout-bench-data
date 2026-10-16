// fp_mul_e5f10 -- binary16 FP multiplier (1 sign, 5 exp, 10 frac), subnormals,
// round-to-nearest ties-to-even, canonical qNaN 16'h7E00.
//
// v4 = v3 with normalized operands:
//   A' = m_a << la,  B' = m_b << lb   (input side, parallel to the multiply)
//   P' = A' * B' = P * 2^(la+lb)      (low la+lb bits of P' are zero)
//   => d = ~P'[21] is just the TOP product bit (no variable-index extraction!)
//   => window shift sh' = sh + t = max(11 - d, 58 - u0), no t dependence.
// Sticky unchanged: |P'[sh'-2:0] = |P[sh-2:0] since P'[t-1:0] = 0.

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

  // product of two subnormals always rounds to signed zero
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

  // ------------------------------- normalized significands (input side)
  wire [10:0] ma_s0 = la[0] ? {m_a[9:0],  1'b0}  : m_a;
  wire [10:0] ma_s1 = la[1] ? {ma_s0[8:0], 2'b0} : ma_s0;
  wire [10:0] ma_s2 = la[2] ? {ma_s1[6:0], 4'b0} : ma_s1;
  wire [10:0] A     = la[3] ? {ma_s2[2:0], 8'b0} : ma_s2;
  wire [10:0] mb_s0 = lb[0] ? {m_b[9:0],  1'b0}  : m_b;
  wire [10:0] mb_s1 = lb[1] ? {mb_s0[8:0], 2'b0} : mb_s0;
  wire [10:0] mb_s2 = lb[2] ? {mb_s1[6:0], 4'b0} : mb_s1;
  wire [10:0] B     = lb[3] ? {mb_s2[2:0], 8'b0} : mb_s2;

  wire [21:0] prod = A * B;          // = m_a*m_b << (la+lb)
  wire        d    = ~prod[21];      // lz(prod_unshifted) = la+lb+d

  // ------------------------------------------- input-side combined terms
  wire [5:0] pa = ea_eff - {2'b0, la} + 6'd16;    // 6..47
  wire [5:0] pb = eb_eff - {2'b0, lb} + 6'd16;
  wire [6:0] u0 = {1'b0, pa} + {1'b0, pb};        // = s - (la+lb) + 32

  wire gt = (u0 >= 7'd48);      // E_field >= 2 : normal for both d
  wire eq = (u0 == 7'd47);      // E_field == 1 : normal iff d == 0

  // ------------------------------------------------- shift control
  // sh' = max(11-d, 58-u0) = S_pre + L,  S_pre = S_base-1, L = ~delta
  wire [5:0] sp6  = gt ? 6'd10 : (6'd57 - u0[5:0]);
  wire [4:0] S_pre = (sp6 > 6'd31) ? 5'd31 : sp6[4:0];
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
  wire [4:0] S0 = ov0 ? 5'd31 : E0v[4:0];
  wire [4:0] SM = ovM ? 5'd31 : EM[4:0];
  wire [4:0] SP = ovP ? 5'd31 : EP[4:0];

  wire selX = d ^ carry;
  wire [4:0] exp_norm = selX ? (carry ? SP : SM) : S0;
  wire       ovf      = selX ? (carry ? ovP : ovM) : ov0;
  wire       is_norm  = gt | (eq & ~d);
  wire [4:0] exp_sub  = sum[10] ? 5'd1 : 5'd0;

  // ------------------------------------------------- parallel output pack
  wire z = zero_in | both_sub;
  wire [4:0] exp_data = is_norm ? exp_norm : exp_sub;
  wire [4:0] exp_y = (nan_in | inf_in) ? 5'd31 :
                     z                  ? 5'd0  : exp_data;
  wire [9:0] frac_y = nan_in           ? 10'h200 :
                      (inf_in | z | (is_norm & ovf)) ? 10'h0 : sum[9:0];
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

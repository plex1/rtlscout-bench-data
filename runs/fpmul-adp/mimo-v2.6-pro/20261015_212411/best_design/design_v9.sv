// fp_mul_e5f10 -- binary16 FP multiplier (1 sign, 5 exp, 10 frac), subnormals,
// round-to-nearest ties-to-even, canonical qNaN 16'h7E00.
//
// v9 = v5 with a cheaper fraction-zeroing cone: overflow-to-inf is detected
// as is_norm & (exp_norm == 31) (finite results never have exp field 31;
// candidates are saturated), instead of a separate mux cone over ov flags.

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

  // ---------------------- one product bit: exact lz = t + d, d = ~prod[21-t]
  wire [10:0] Mv = {1'b1,          (t >= 5'd1),  (t >= 5'd2),
                        (t >= 5'd3),  (t >= 5'd4),  (t >= 5'd5),
                        (t >= 5'd6),  (t >= 5'd7),  (t >= 5'd8),
                        (t >= 5'd9),  (t >= 5'd10)};
  wire d = ~|(prod[21:11] & Mv);

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

  // --------------------------------------------- round / significand
  wire [10:0] mp    = q[11:1];
  wire        guard = q[0];
  wire        rnd   = guard & (sk | mp[0]);      // carry into bit 0 only

  wire g01 = mp[0] & mp[1];
  wire g23 = mp[2] & mp[3];
  wire g45 = mp[4] & mp[5];
  wire g67 = mp[6] & mp[7];
  wire g89 = mp[8] & mp[9];
  wire A3  = g01 & g23;
  wire A7  = A3 & g45 & g67;
  wire A9  = A7 & g89;
  wire A10 = A9 & mp[10];
  wire A0  = mp[0];
  wire A1  = g01;
  wire A2  = g01 & mp[2];
  wire A4  = A3 & mp[4];
  wire A5  = A3 & g45;
  wire A6  = A5 & mp[6];
  wire A8  = A7 & mp[8];

  wire [9:0] frac;
  assign frac[0] = mp[0] ^ rnd;
  assign frac[1] = mp[1] ^ (guard & A0);
  assign frac[2] = mp[2] ^ (guard & A1);
  assign frac[3] = mp[3] ^ (guard & A2);
  assign frac[4] = mp[4] ^ (guard & A3);
  assign frac[5] = mp[5] ^ (guard & A4);
  assign frac[6] = mp[6] ^ (guard & A5);
  assign frac[7] = mp[7] ^ (guard & A6);
  assign frac[8] = mp[8] ^ (guard & A7);
  assign frac[9] = mp[9] ^ (guard & A8);

  wire carry   = guard & A10;                    // exp increment (normal)
  wire carry_s = mp[10] ^ (guard & A9);          // subnormal -> min normal

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
  wire       is_norm  = gt | (eq & ~d);
  wire [4:0] exp_sub  = carry_s ? 5'd1 : 5'd0;

  // ------------------------------------------------- parallel output pack
  wire z  = zero_in | both_sub;
  wire ovn = is_norm & (exp_norm == 5'd31);   // overflow -> inf
  wire fz = ~(nan_in | inf_in | z | ovn);     // frac passes through

  wire [4:0] exp_data = is_norm ? exp_norm : exp_sub;
  wire [4:0] exp_y = (nan_in | inf_in) ? 5'd31 :
                     z                  ? 5'd0  : exp_data;

  assign y[15] = ~nan_in & s_y;
  assign y[14] = exp_y[4];
  assign y[13] = exp_y[3];
  assign y[12] = exp_y[2];
  assign y[11] = exp_y[1];
  assign y[10] = exp_y[0];
  assign y[9]  = nan_in | (fz & frac[9]);
  assign y[8]  = fz & frac[8];
  assign y[7]  = fz & frac[7];
  assign y[6]  = fz & frac[6];
  assign y[5]  = fz & frac[5];
  assign y[4]  = fz & frac[4];
  assign y[3]  = fz & frac[3];
  assign y[2]  = fz & frac[2];
  assign y[1]  = fz & frac[1];
  assign y[0]  = fz & frac[0];
endmodule

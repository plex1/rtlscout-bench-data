// fp_mul_e5f10 -- binary16 FP multiplier (1 sign, 5 exp, 10 frac), subnormals,
// round-to-nearest ties-to-even, canonical qNaN 16'h7E00.
//
// v27 = v25 (narrow shifter, decoder Mv, flat sticky, balanced AND prefix,
// flat exp, AND-form frac gating) with the EM/EP subtractor exponent
// candidates (the form that consistently maps fast).

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
  wire t01 = t[1] & t[0];
  wire t0i = t[1] | t[0];
  wire lo2 = |t[2:0];
  wire [10:0] Mv = {1'b1,
                    |t[4:0],
                    |t[4:1],
                    |t[4:2] | t01,
                    |t[4:2],
                    |t[4:3] | (t[2] & t0i),
                    |t[4:3] | (t[2] & t[1]),
                    |t[4:3] | (t[2] & t01),
                    |t[4:3],
                    t[4] | (t[3] & lo2),
                    t[4] | (t[3] & (t[2] | t01))};
  wire d = ~|(prod[21:11] & Mv);

  // ------------------------------------------------- shift control
  wire [4:0] S_pre = gt ? Am1 : Bm1;   // = S_base - 1
  wire       delta = gt & d;
  wire       L     = ~delta;           // late 1-bit extra shift

  // ----------------------- narrow funnel shifter (early control) + sticky
  wire [27:0] wx  = {5'b0, prod, 1'b0};

  wire [27:0] v0 = S_pre[4] ? {16'b0, wx[27:16]} : wx;
  wire [19:0] v1 = S_pre[3] ? v0[27:8] : v0[19:0];
  wire [15:0] v2 = S_pre[2] ? v1[19:4] : v1[15:0];
  wire [13:0] v3 = S_pre[1] ? v2[15:2] : v2[13:0];
  wire [12:0] v4 = S_pre[0] ? v3[13:1] : v3[12:0];

  wire [11:0] q12 = L ? v4[12:1] : v4[11:0];

  wire sk = (S_pre[4] & (|wx[15:0])) | (S_pre[3] & (|v0[7:0])) |
            (S_pre[2] & (|v1[3:0]))  | (S_pre[1] & (|v2[1:0])) |
            (S_pre[0] & v3[0])       | (L & v4[0]);

  // --------------------------------------------- round / significand
  wire [10:0] mp    = q12[11:1];
  wire        guard = q12[0];
  wire        rnd   = guard & (sk | mp[0]);

  // balanced AND groups
  wire c01 = mp[0] & mp[1];
  wire c23 = mp[2] & mp[3];
  wire c45 = mp[4] & mp[5];
  wire c67 = mp[6] & mp[7];
  wire c89 = mp[8] & mp[9];

  wire A7  = c01 & c23 & c45 & c67;
  wire A3  = c01 & c23;
  wire A5  = A3 & c45;
  wire A8  = A7 & mp[8];
  wire A9  = A7 & c89;
  wire A10 = A9 & mp[10];
  wire A0  = mp[0];
  wire A1  = c01;
  wire A2  = c01 & mp[2];
  wire A4  = A3 & mp[4];
  wire A6  = A5 & mp[6];

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

  wire carry   = guard & A10;
  wire carry_s = mp[10] ^ (guard & A9);

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

  // ---------------------------- exponent output: weighted AND-OR (flat)
  wire sp  = nan_in | inf_in;
  wire z   = zero_in | both_sub;
  wire g_n = gt & ~sp & ~z;
  wire e_n = eq & ~sp & ~z;
  wire wSM = g_n & d & ~carry;
  wire wSP = (g_n | e_n) & ~d & carry;
  wire wS0 = ((g_n | e_n) & ~d & ~carry) | (g_n & d & carry);
  wire c_s = ~sp & ~z & ~gt & (~eq | d);

  assign y[14] = sp | (wSP & SP[4]) | (wSM & SM[4]) | (wS0 & S0[4]);
  assign y[13] = sp | (wSP & SP[3]) | (wSM & SM[3]) | (wS0 & S0[3]);
  assign y[12] = sp | (wSP & SP[2]) | (wSM & SM[2]) | (wS0 & S0[2]);
  assign y[11] = sp | (wSP & SP[1]) | (wSM & SM[1]) | (wS0 & S0[1]);
  assign y[10] = sp | (wSP & SP[0]) | (wSM & SM[0]) | (wS0 & S0[0]) | (c_s & carry_s);

  // ---------------------------------- fraction side: single AND gating
  wire selX = d ^ carry;
  wire ovf  = selX ? (carry ? ovP : ovM) : ov0;
  wire is_norm = gt | (eq & ~d);
  wire fz = ~(nan_in | inf_in | z | (is_norm & ovf));

  assign y[15] = ~nan_in & s_y;
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

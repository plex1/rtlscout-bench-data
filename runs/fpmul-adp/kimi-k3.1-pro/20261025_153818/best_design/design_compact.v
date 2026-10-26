// Binary16 multiplication with early leading-position detection.
module fp_mul_e5f10 (a, b, y);
  input [15:0] a;
  input [15:0] b;
  output [15:0] y;

  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0], f_b = b[9:0];
  wire ez_a = ~|e_a, ez_b = ~|e_b;
  wire eo_a = &e_a, eo_b = &e_b;
  wire fz_a = ~|f_a, fz_b = ~|f_b;
  wire zero_a = ez_a & fz_a, zero_b = ez_b & fz_b;
  wire inf_a = eo_a & fz_a, inf_b = eo_b & fz_b;
  wire nan_in = (eo_a & ~fz_a) | (eo_b & ~fz_b) |
                (inf_a & zero_b) | (inf_b & zero_a);
  wire zero_in = zero_a | zero_b;
  wire sign = a[15] ^ b[15];

  wire [10:0] m_a = {~ez_a, f_a};
  wire [10:0] m_b = {~ez_b, f_b};
  wire [21:0] prod = m_a * m_b;

  function [3:0] lead11;
    input [10:0] m;
    begin
      if      (m[10]) lead11 = 4'd10;
      else if (m[9])  lead11 = 4'd9;
      else if (m[8])  lead11 = 4'd8;
      else if (m[7])  lead11 = 4'd7;
      else if (m[6])  lead11 = 4'd6;
      else if (m[5])  lead11 = 4'd5;
      else if (m[4])  lead11 = 4'd4;
      else if (m[3])  lead11 = 4'd3;
      else if (m[2])  lead11 = 4'd2;
      else if (m[1])  lead11 = 4'd1;
      else            lead11 = 4'd0;
    end
  endfunction

  wire [4:0] base = {1'b0, lead11(m_a)} + {1'b0, lead11(m_b)};
  // Equivalent to effective_ea + effective_eb + base, with subnormal
  // effective exponent one, but without two exponent muxes.
  wire [7:0] pre_total = {3'b0, e_a} + {3'b0, e_b} +
                         {3'b0, base} + {7'b0, ez_a} + {7'b0, ez_b};
  wire [6:0] exp_sum = {2'b0, e_a} + {2'b0, e_b} +
                       {6'b0, ez_a} + {6'b0, ez_b};

  wire [22:0] aligned_n = {1'b0, prod} << (5'd21 - base);
  wire carry_n = aligned_n[22];
  wire [10:0] q_n = carry_n ? aligned_n[22:12] : aligned_n[21:11];
  wire g_n = carry_n ? aligned_n[11] : aligned_n[10];
  wire s_n = carry_n ? |aligned_n[10:0] : |aligned_n[9:0];

  wire [5:0] sub_shift = 6'd26 - exp_sum[5:0];
  wire [46:0] aligned_s = {13'b0, prod, 12'b0} << (6'd24 - sub_shift);
  wire [10:0] q_s = aligned_s[46:36];
  wire g_s = aligned_s[35];
  wire s_s = |aligned_s[34:0];

  wire normal = (pre_total > 8'd35) |
                ((pre_total == 8'd35) & carry_n);
  wire round_n = g_n & (s_n | q_n[0]);
  wire round_s = g_s & (s_s | q_s[0]);
  wire [11:0] rounded_n = {1'b0, q_n} + {11'b0, round_n};
  wire [11:0] rounded_s = {1'b0, q_s} + {11'b0, round_s};
  wire [11:0] rounded = normal ? rounded_n : rounded_s;
  wire round_carry_n = round_n & (&q_n);
  wire [7:0] final_total = pre_total + {7'b0, carry_n} + {7'b0, round_carry_n};
  wire overflow = final_total > 8'd65;
  wire [4:0] exp_n = final_total[4:0] - 5'd3;
  wire [4:0] exp_s = {4'b0, rounded_s[10]};

  assign y = nan_in ? 16'h7e00 :
             (inf_a | inf_b) ? {sign, 5'h1f, 10'b0} :
             zero_in ? {sign, 15'b0} :
             (normal & overflow) ? {sign, 5'h1f, 10'b0} :
             {sign, normal ? exp_n : exp_s, rounded[9:0]};
endmodule

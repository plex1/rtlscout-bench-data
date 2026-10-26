// IEEE binary16 multiplier with gradual underflow.
// Leading positions are found before multiplication; the product then needs
// only a left alignment and a one-bit normalization choice.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];
  wire ez_a = ~|e_a;
  wire ez_b = ~|e_b;
  wire eo_a = &e_a;
  wire eo_b = &e_b;
  wire fz_a = ~|f_a;
  wire fz_b = ~|f_b;
  wire zero_a = ez_a & fz_a;
  wire zero_b = ez_b & fz_b;
  wire inf_a = eo_a & fz_a;
  wire inf_b = eo_b & fz_b;
  wire nan_a = eo_a & ~fz_a;
  wire nan_b = eo_b & ~fz_b;
  wire nan_in = nan_a | nan_b | (inf_a & zero_b) | (inf_b & zero_a);
  wire inf_in = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sign = a[15] ^ b[15];

  wire [10:0] m_a = {~ez_a, f_a};
  wire [10:0] m_b = {~ez_b, f_b};
  wire [4:0] ea = ez_a ? 5'd1 : e_a;
  wire [4:0] eb = ez_b ? 5'd1 : e_b;
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

  wire [3:0] la = lead11(m_a);
  wire [3:0] lb = lead11(m_b);
  wire [4:0] base = {1'b0, la} + {1'b0, lb};
  wire [6:0] q_sum = {2'b0, ea} + {2'b0, eb};
  wire [7:0] pre_total = {1'b0, q_sum} + {3'b0, base};

  // Return {11 retained bits, guard, sticky} for prod >> sh.
  function [12:0] sub_bits;
    input [21:0] p;
    input [5:0] sh;
    reg [46:0] ext;
    begin
      ext = {13'b0, p, 12'b0} << (6'd24 - sh);
      sub_bits = {ext[46:36], ext[35], |ext[34:0]};
    end
  endfunction

  // The leading product bit is at base or base+1.  Left alignment retains
  // every product bit, so guard and sticky need no second variable shifter.
  wire [22:0] aligned_n = {1'b0, prod} << (5'd21 - base);
  wire carry_n = aligned_n[22];
  wire [10:0] q_n = carry_n ? aligned_n[22:12] : aligned_n[21:11];
  wire g_n = carry_n ? aligned_n[11] : aligned_n[10];
  wire s_n = carry_n ? |aligned_n[10:0] : |aligned_n[9:0];

  wire [5:0] sub_shift = 6'd26 - q_sum[5:0];
  wire [12:0] packed_s = sub_bits(prod, sub_shift);

  // pre_total is 35 + the result exponent before product normalization.
  wire normal = (pre_total > 8'd35) |
                ((pre_total == 8'd35) & carry_n);
`ifndef SHARED_ROUND
  wire round_n = g_n & (s_n | q_n[0]);
  wire round_s = packed_s[1] & (packed_s[0] | packed_s[2]);
  wire [11:0] rounded_n = {1'b0, q_n} + {11'b0, round_n};
  wire [11:0] rounded_s = {1'b0, packed_s[12:2]} + {11'b0, round_s};
  wire [11:0] rounded = normal ? rounded_n : rounded_s;
  wire [7:0] final_total = pre_total + {7'b0, carry_n} + {7'b0, rounded_n[11]};
  wire overflow = final_total > 8'd65;
  wire [4:0] exp_n = final_total[4:0] - 5'd3;
  wire [4:0] exp_s = {4'b0, rounded_s[10]};
`else
  wire [10:0] q = normal ? q_n : packed_s[12:2];
  wire g = normal ? g_n : packed_s[1];
  wire st = normal ? s_n : packed_s[0];
  wire round_up = g & (st | q[0]);
  wire [11:0] rounded = {1'b0, q} + {11'b0, round_up};

  wire [7:0] final_total = pre_total + {7'b0, carry_n} + {7'b0, rounded[11]};
  wire overflow = final_total > 8'd65;
  wire [4:0] exp_n = final_total[4:0] - 5'd3;
  wire [4:0] exp_s = {4'b0, rounded[10]};
`endif
  wire [4:0] exp_y = normal ? exp_n : exp_s;
  wire [9:0] frac_y = rounded[9:0];

  assign y = nan_in ? 16'h7e00 :
             inf_in ? {sign, 5'h1f, 10'b0} :
             zero_in ? {sign, 15'b0} :
             (normal & overflow) ? {sign, 5'h1f, 10'b0} :
             {sign, exp_y, frac_y};
endmodule

// Shared-alignment binary16 multiplier, including gradual underflow.
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

  wire [4:0] base = {1'b0, lead11(m_a)} + {1'b0, lead11(m_b)};
  function product_carry;
    input [21:0] p;
    input [4:0] n;
    begin
      case (n)
        5'd0:  product_carry = p[1];
        5'd1:  product_carry = p[2];
        5'd2:  product_carry = p[3];
        5'd3:  product_carry = p[4];
        5'd4:  product_carry = p[5];
        5'd5:  product_carry = p[6];
        5'd6:  product_carry = p[7];
        5'd7:  product_carry = p[8];
        5'd8:  product_carry = p[9];
        5'd9:  product_carry = p[10];
        5'd10: product_carry = p[11];
        5'd11: product_carry = p[12];
        5'd12: product_carry = p[13];
        5'd13: product_carry = p[14];
        5'd14: product_carry = p[15];
        5'd15: product_carry = p[16];
        5'd16: product_carry = p[17];
        5'd17: product_carry = p[18];
        5'd18: product_carry = p[19];
        5'd19: product_carry = p[20];
        default: product_carry = p[21];
      endcase
    end
  endfunction
  wire [6:0] q_sum = {2'b0, ea} + {2'b0, eb};
  wire [7:0] pre_total = {1'b0, q_sum} + {3'b0, base};

  // One lossless alignment supplies both paths.  For normal results the
  // leading one is at bit 35 or 36; for subnormals the binary point is fixed.
  wire [6:0] sub_amount = q_sum[6:0] - 7'd2;
  wire [6:0] norm_amount = 7'd33 - {2'b0, base};
  wire [6:0] shift = normal ? norm_amount : sub_amount;
  wire [46:0] aligned = {13'b0, prod, 12'b0} << shift;
  wire carry_n = aligned[46];
  wire [10:0] q_n = carry_n ? aligned[46:36] : aligned[45:35];
  wire g_n = carry_n ? aligned[35] : aligned[34];
  wire s_n = carry_n ? |aligned[34:0] : |aligned[33:0];
  wire [10:0] q_s = aligned[35:25];
  wire g_s = aligned[24];
  wire s_s = |aligned[23:0];
  wire normal = (pre_total > 8'd35) |
                ((pre_total == 8'd35) & product_carry(prod, base));

  wire round_n = g_n & (s_n | q_n[0]);
  wire round_s = g_s & (s_s | q_s[0]);
  wire [11:0] rounded_n = {1'b0, q_n} + {11'b0, round_n};
  wire [11:0] rounded_s = {1'b0, q_s} + {11'b0, round_s};
  wire [11:0] rounded = normal ? rounded_n : rounded_s;
  wire [7:0] final_total = pre_total + {7'b0, carry_n} + {7'b0, rounded_n[11]};
  wire overflow = final_total > 8'd65;
  wire [4:0] exp_n = final_total[4:0] - 5'd3;
  wire [4:0] exp_s = {4'b0, rounded_s[10]};

  assign y = nan_in ? 16'h7e00 :
             (inf_a | inf_b) ? {sign, 5'h1f, 10'b0} :
             zero_in ? {sign, 15'b0} :
             (normal & overflow) ? {sign, 5'h1f, 10'b0} :
             {sign, normal ? exp_n : exp_s, rounded[9:0]};
endmodule

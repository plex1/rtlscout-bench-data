// fp_mul_e5f10 -- binary16 floating-point multiplier (1 sign, 5 exponent, 10 fraction)
// Subnormals supported, round-to-nearest ties-to-even, canonical qNaN 16'h7E00.
//
// Optimized structure:
//   * operand CLZs are computed on the input side (parallel to the multiply);
//     only one product bit (d = ~prod[21-la-lb]) fixes the exact product CLZ:
//     lz = la + lb + d
//   * ONE funnel barrel shifter (with sticky accumulator) serves both the
//     normal path (shift 11-lz) and the subnormal path (shift 26-s):
//     sh = max(11-lz, 26-s)
//   * one shared rounding incrementer: round = guard & (sticky | lsb) is the
//     same formula for normal and subnormal results

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
  wire [6:0] s_sum  = {1'b0, ea_eff} + {1'b0, eb_eff};   // 2..60

  wire [21:0] prod = m_a * m_b;

  // ----------------------------------------- input-side CLZ of significands
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
  wire [4:0] t  = {1'b0, la} + {1'b0, lb};                    // 0..22
  wire [4:0] k  = (t >= 5'd21) ? 5'd0 : (5'd21 - t);          // index into prod

  // One product bit fixes the exact leading-zero count:  lz = t + d, d = ~prod[k]
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
  wire d = ~v5;                       // lz = t + d

  // -------------------------------------------------- exponent / shift control
  // u = s_sum - lz + 32;  normal result iff E_field = u-46 >= 1
  wire [6:0] u      = s_sum + 7'd32 - {2'b00, t} - {6'd0, d};
  wire       is_norm = (u >= 7'd47);
  wire [4:0] sh      = is_norm ? (5'd11 - t - {4'd0, d}) : (5'd26 - s_sum[4:0]);

  // --------------------------------------- single funnel shifter + sticky
  wire [22:0] w  = {prod, 1'b0};
  wire [22:0] q0 = sh[0] ? {1'b0,     w[22:1]}  : w;
  wire [22:0] q1 = sh[1] ? {2'b00,    q0[22:2]} : q0;
  wire [22:0] q2 = sh[2] ? {4'b0,     q1[22:4]} : q1;
  wire [22:0] q3 = sh[3] ? {8'b0,     q2[22:8]} : q2;
  wire [22:0] q4 = sh[4] ? {16'b0,    q3[22:16]} : q3;

  wire sk1 = sh[1] & (|q0[1:0]);
  wire sk2 = sk1 | (sh[2] & (|q1[3:0]));
  wire sk3 = sk2 | (sh[3] & (|q2[7:0]));
  wire sk4 = sk3 | (sh[4] & (|q3[15:0]));

  wire [10:0] mant_pre = q4[11:1];
  wire        guard    = q4[0];
  wire        sticky   = sk4;
  wire        rnd      = guard & (sticky | mant_pre[0]);
  wire [11:0] sum      = {1'b0, mant_pre} + {11'd0, rnd};

  // ------------------------------------------------------------- pack
  wire [6:0] e7        = u + {6'd0, sum[11]} - 7'd46;   // normal exponent field
  wire       ovf       = (e7 > 7'd30);
  wire [4:0] exp_norm  = e7[4:0];
  wire [4:0] exp_sub   = sum[10] ? 5'd1 : 5'd0;
  wire [9:0] frac_y    = sum[9:0];

  assign y = nan_in  ? 16'h7E00 :
             inf_in  ? {s_y, 5'h1F, 10'h0} :
             zero_in ? {s_y, 15'h0} :
             is_norm ? (ovf ? {s_y, 5'h1F, 10'h0} : {s_y, exp_norm, frac_y}) :
                       {s_y, exp_sub, frac_y};
endmodule

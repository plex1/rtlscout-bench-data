// fp_mul_e5f10 -- IEEE binary16 multiplier, unified single-shifter formulation.
module fp_mul_e5f10 (a, b, y);
  input  wire [15:0] a, b;
  output wire [15:0] y;

  wire        s_a = a[15], s_b = b[15];
  wire [4:0]  e_a = a[14:10], e_b = b[14:10];
  wire [9:0]  f_a = a[9:0],   f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'd31);
  wire e_b_ones = (e_b == 5'd31);
  wire f_a_zero = (f_a == 10'd0);
  wire f_b_zero = (f_b == 10'd0);

  wire zero_a = e_a_zero & f_a_zero;
  wire zero_b = e_b_zero & f_b_zero;
  wire nan_a  = e_a_ones & ~f_a_zero;
  wire nan_b  = e_b_ones & ~f_b_zero;
  wire inf_a  = e_a_ones & f_a_zero;
  wire inf_b  = e_b_ones & f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a = {~e_a_zero, f_a};
  wire [10:0] m_b = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  wire [21:0] p = m_a * m_b;

  // leading-zero count of the product
  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (p[i]) lz = 5'd21 - i[4:0];
  end

  // prefix OR of the product: por[k] = |p[k:0]
  wire [21:0] por;
  assign por[0] = p[0];
  genvar k;
  generate for (k = 1; k < 22; k = k + 1) begin : g_por
    assign por[k] = por[k-1] | p[k];
  end endgenerate
  wire [22:0] below = {por[21], por};   // below[22] = below[21] = |p[21:0]

  // unified right-shift amount:  t = max(11-lz, 26-ea-eb)
  wire signed [6:0] e0 = $signed({2'b0, ea_eff}) + $signed({2'b0, eb_eff})
                       - 7'sd14 - $signed({2'b0, lz});
  wire [5:0] ta = 6'd11 - {1'b0, lz};
  wire [5:0] tb = 6'd26 - {1'b0, ea_eff} - {1'b0, eb_eff};
  wire [5:0] t  = (e0 >= 7'sd1) ? ta : tb;      // >= 1 always
  wire [4:0] s1 = t[4:0] - 5'd1;

  wire [21:0] sh = p >> s1;                     // shift by t-1: sh[0]=guard
  wire        lsb    = sh[1];
  wire        guard  = sh[0];
  wire [4:0]  t2     = t[4:0] - 5'd2;
  wire        sticky = (t >= 6'd2) ? below[t2] : 1'b0;
  wire        round_up = guard & (sticky | lsb);

  wire [11:0] M = {1'b0, sh[11:1]} + {11'd0, round_up};

  wire signed [7:0] ef1 = $signed({3'b0, ea_eff}) + $signed({3'b0, eb_eff})
                        - 8'sd14 - $signed({3'b0, lz}) + {7'b0, M[11]};
  wire signed [7:0] ef  = (ef1 >= 8'sd1) ? ef1 : (M[10] ? 8'sd1 : 8'sd0);

  wire is_inf  = ~nan_in & (inf_in | (ef1 > 8'sd30));
  wire is_zero = ~nan_in & ~is_inf & (zero_in | ((ef == 8'sd0) & (M[9:0] == 10'd0)));

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_zero ? 5'd0 : (ef >= 8'sd1) ? ef[4:0] : 5'd0;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_zero) ? 10'd0 : M[9:0];
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

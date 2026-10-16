// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire ea_z = ~|e_a;
  wire eb_z = ~|e_b;
  wire ea_f =  &e_a;
  wire eb_f =  &e_b;
  wire fa_z = ~|f_a;
  wire fb_z = ~|f_b;

  wire zero_a  = ea_z & fa_z;
  wire zero_b  = eb_z & fb_z;
  wire zero_in = zero_a | zero_b;

  // NaN: an all-ones exponent with non-zero fraction, or inf * 0
  wire is_nan = (ea_f & (~fa_z | zero_b)) | (eb_f & (~fb_z | zero_a));
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a    = {|e_a, f_a};
  wire [10:0] m_b    = {|e_b, f_b};
  wire [4:0]  ea_eff = e_a | {4'd0, ea_z};
  wire [4:0]  eb_eff = e_b | {4'd0, eb_z};

  wire [21:0] prod = m_a * m_b;

  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = (k == 0) ? prod[0] : (below[k-1] | prod[k]);
    end
  endgenerate

  wire [31:0] rv = {10'b0, prod} | {{8{below[21]}}, below, 2'b00};

  wire [7:0] exp_sum = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire       sel     = (exp_sum < 8'd15 + {3'd0, lz});
  wire [7:0] e_lim   = sel ? 8'd0 : (exp_sum - 8'd14 - {3'd0, lz});
  wire       ovf_hi  = (e_lim > 8'd30);

  wire [3:0]  shn = 4'd11 - lz[3:0];
  wire [11:0] w_n = {prod, 1'b0} >> shn;
  wire        r_n = rv[shn];

  wire [4:0]  shs = 5'd26 - exp_sum[4:0];
  wire [11:0] w_s = {12'b0, prod, 1'b0} >> shs;
  wire        r_s = rv[shs];

  wire [11:0] w   = sel ? w_s : w_n;
  wire        rvb = sel ? r_s : r_n;

  wire [10:0] mant_pre = w[11:1];
  wire        round_up = w[0] & rvb;

  wire [9:0] frac  = mant_pre[9:0] + {9'd0, round_up};
  wire       carry = round_up & (&mant_pre[9:0]);

  // ea_f / eb_f already cover NaN and inf inputs
  wire big = ea_f | eb_f | (ovf_hi & ~zero_in);

  wire [4:0] exp_n5 = e_lim[4:0] + {4'd0, carry};
  wire [4:0] exp_y  = big ? 5'd31 : (zero_in ? 5'd0 : exp_n5);
  wire [9:0] frac_y = is_nan ? 10'h200 : (big ? 10'd0 : frac);
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

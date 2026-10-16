// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// The rounding adder is a single incrementer:  A = prod>>(sh-1), B = A + c
// with c = |prod[sh-2:0]| | prod[sh] precomputed per index.  frac and the
// rounding carry both come straight out of B.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam MAX_E = 31;

  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == MAX_E);
  wire e_b_ones = (e_b == MAX_E);
  wire zero_a   = e_a_zero & (f_a == 10'd0);
  wire zero_b   = e_b_zero & (f_b == 10'd0);
  wire nan_a    = e_a_ones & (f_a != 10'd0);
  wire nan_b    = e_b_ones & (f_b != 10'd0);
  wire inf_a    = e_a_ones & (f_a == 10'd0);
  wire inf_b    = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  wire [21:0] prod = m_a * m_b;

  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // below[j] = |prod[j:0]|
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = (k == 0) ? prod[0] : (below[k-1] | prod[k]);
    end
  endgenerate

  // rv[j] = prod[j] | |prod[j-2:0]|   (round-up input if guard is at j-1)
  wire [31:0] bv = {{8{below[21]}}, below, 2'b00};
  wire [31:0] rv = {10'b0, prod} | bv;

  wire [7:0] exp_sum = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire       sel     = (exp_sum < 8'd15 + {3'd0, lz});
  wire [7:0] e_lim   = sel ? 8'd0 : (exp_sum - 8'd14 - {3'd0, lz});
  wire       ovf_hi  = (e_lim > 8'd30);

  // ---- normal path : shift 0..11 ----
  wire [3:0]  shn = 4'd11 - lz[3:0];
  wire [11:0] w_n = {prod, 1'b0} >> shn;
  wire        c_n = rv[shn];

  // ---- subnormal path : shift 0..24 ----
  wire [4:0]  shs = 5'd26 - exp_sum[4:0];
  wire [11:0] w_s = {12'b0, prod, 1'b0} >> shs;
  wire        c_s = rv[shs];

  wire [11:0] w = sel ? w_s : w_n;
  wire        c = sel ? c_s : c_n;

  // B = A + c :  bit0 = guard, bit1 = mant_pre[0], bits 11:1 = mant_pre
  wire [11:0] bsum  = w + {11'd0, c};
  wire [9:0]  frac  = bsum[10:1];
  wire        carry = bsum[11] ^ w[11];

  wire is_nan = nan_in;
  wire big    = is_nan | inf_in | (ovf_hi & ~zero_in);

  wire [4:0] exp_n5 = e_lim[4:0] + {4'd0, carry};
  wire [4:0] exp_y  = big ? 5'd31 : (zero_in ? 5'd0 : exp_n5);
  wire [9:0] frac_y = is_nan ? 10'h200 : (big ? 10'd0 : frac);
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

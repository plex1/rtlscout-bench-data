// fp_mul_e5f10 -- IEEE binary16 multiplier (1+5+10) with subnormal support.
// Parallel narrow shifters; round-up bit is precomputed per product bit so the
// rounding stage is just an incrementer.

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
  // dr[j]    = prod[j] & (prod[j+1] | below[j-1])   == round-up bit if the
  //            guard bit sits at position j
  wire [21:0] below;
  wire [21:0] dr;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_pre
      assign below[k] = (k == 0) ? prod[0] : (below[k-1] | prod[k]);
      assign dr[k]    = prod[k] & (((k < 21) ? prod[k+1] : 1'b0) |
                                   ((k > 0)  ? below[k-1]  : 1'b0));
    end
  endgenerate

  wire [31:0] drext = {9'b0, dr, 1'b0};   // drext[j] = dr[j-1]

  wire [7:0] exp_sum   = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire       underflow = (exp_sum < 8'd15 + {3'd0, lz});
  wire [7:0] e_norm0   = exp_sum - 8'd14 - {3'd0, lz};

  // ---- normal path : shift 0..11 ----
  wire [3:0]  shn = 4'd11 - lz[3:0];
  wire [10:0] mp_n = prod >> shn;
  wire        ru_n = drext[shn];

  // ---- subnormal path : shift 0..24 ----
  wire [4:0]  shs = 5'd26 - exp_sum[4:0];
  wire [10:0] mp_s = prod >> shs;
  wire        ru_s = drext[shs];

  wire        sel  = underflow;
  wire [10:0] mant_pre = sel ? mp_s : mp_n;
  wire        round_up = sel ? ru_s : ru_n;

  wire [9:0] frac  = mant_pre[9:0] + {9'd0, round_up};
  wire       carry = round_up & (&mant_pre[9:0]);

  wire [7:0] exp_raw  = sel ? {7'd0, carry} : (e_norm0 + {7'd0, carry});
  wire       overflow = (exp_raw > 8'd30);

  wire is_nan  = nan_in;
  wire is_inf  = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_zero = ~is_nan & ~is_inf & (zero_in | (sel & ~carry & (frac == 10'd0)));

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_raw[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac;
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

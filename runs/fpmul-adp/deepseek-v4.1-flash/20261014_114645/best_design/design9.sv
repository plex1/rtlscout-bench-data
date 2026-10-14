// fp_mul_e5f10 -- IEEE binary16 multiplier, subnormal support.
// design7 structure + explicit balanced partial-product adder tree.
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam MAX_E  = 31;
  localparam PROD_W = 22;

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

  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  // ---- 11 partial products, balanced 4-level adder tree ----
  wire [22:0] S [0:10];
  genvar k;
  generate
    for (k = 0; k < 11; k = k + 1) begin : g_pp
      assign S[k] = ({12'b0, m_a} & {23{m_b[k]}}) << k;
    end
  endgenerate

  wire [22:0] A0 = S[0] + S[1];
  wire [22:0] A1 = S[2] + S[3];
  wire [22:0] A2 = S[4] + S[5];
  wire [22:0] A3 = S[6] + S[7];
  wire [22:0] A4 = S[8] + S[9];
  wire [22:0] B0 = A0 + A1;
  wire [22:0] B1 = A2 + A3;
  wire [22:0] B2 = A4 + S[10];
  wire [22:0] C0 = B0 + B1;
  wire [22:0] DD = C0 + B2;
  wire [PROD_W-1:0] prod = DD[21:0];

  // leading-zero count (0 when prod == 0)
  reg [4:0] lz;
  integer   i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < PROD_W; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  wire [PROD_W-1:0] below;
  generate
    for (k = 0; k < PROD_W; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  function bit_at;
    input [PROD_W-1:0] v;
    input integer      idx;
    begin
      bit_at = (idx >= 0 && idx < PROD_W) ? v[idx] : 1'b0;
    end
  endfunction

  // ---- normal path ----
  wire              need_right = (lz <= 5'd11);
  wire [3:0]        sr         = 4'd11 - lz[3:0];
  wire [3:0]        sl         = lz[3:0] - 4'd11;
  wire [PROD_W-1:0] shifted    = need_right ? (prod >> sr) : (prod << sl);
  wire [10:0]       mant_pre   = shifted[10:0];

  wire guard    = need_right ? bit_at(prod,  sr - 1) : 1'b0;
  wire sticky   = need_right ? bit_at(below, sr - 2) : 1'b0;
  wire round_up = guard & (sticky | mant_pre[0]);

  wire [11:0] mant_round = {1'b0, mant_pre} + {11'd0, round_up};
  wire        carry      = mant_round[11];
  wire [9:0]  frac_norm  = carry ? mant_round[10:1] : mant_round[9:0];

  // ---- exponent ----
  wire [7:0] exp_sum = {3'd0, ea_eff} + {3'd0, eb_eff};
  wire [7:0] pse     = exp_sum + {3'd0, (5'd21 - lz)};
  wire [7:0] psec    = pse + {7'd0, carry};
  wire       is_sub  = (pse <= 8'd35);
  wire       overflow= (psec >= 8'd66);
  wire [7:0] e_norm  = psec - 8'd35;
  wire [4:0] exp_norm= e_norm[4:0];

  // ---- subnormal path ----
  wire [7:0] sub_shift_c = 8'd26 - exp_sum;
  integer    total_shift;
  always @* total_shift = sub_shift_c[5:0];
  wire [PROD_W-1:0] sub_shifted = (sub_shift_c < 8'd22) ? (prod >> total_shift) : {PROD_W{1'b0}};
  wire [9:0]  frac_trunc = sub_shifted[9:0];
  wire        guard_s    = bit_at(prod,  total_shift - 1);
  wire        sticky_s   = bit_at(below, total_shift - 2);
  wire        round_s    = guard_s & (sticky_s | frac_trunc[0]);
  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_s};
  wire        carry_s    = frac_sum[10];
  wire [9:0]  frac_sub   = frac_sum[9:0];

  // ---- pack ----
  wire [4:0] exp_ns  = is_sub ? {4'd0, carry_s} : exp_norm;
  wire [9:0] frac_ns = is_sub ? frac_sub : frac_norm;

  assign y = nan_in  ? 16'h7E00
           : zero_in ? {s_y, 5'd0, 10'd0}
           : (inf_in | (overflow & ~is_sub)) ? {s_y, 5'd31, 10'd0}
           : {s_y, exp_ns, frac_ns};
endmodule

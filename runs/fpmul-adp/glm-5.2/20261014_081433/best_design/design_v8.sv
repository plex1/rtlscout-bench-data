// fp_mul_e5f10 -- behavioral, left-shift normal path, synthesis-friendly
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS = 15;

  wire        sa = a[15];
  wire        sb = b[15];
  wire [4:0]  ea = a[14:10];
  wire [4:0]  eb = b[14:10];
  wire [9:0]  fa = a[9:0];
  wire [9:0]  fb = b[9:0];

  // classify
  wire a_zero = (ea == 0) && (fa == 0);
  wire b_zero = (eb == 0) && (fb == 0);
  wire a_inf  = (ea == 31) && (fa == 0);
  wire b_inf  = (eb == 31) && (fb == 0);
  wire a_nan  = (ea == 31) && (fa != 0);
  wire b_nan  = (eb == 31) && (fb != 0);

  wire is_nan = a_nan | b_nan | (a_inf & b_zero) | (a_zero & b_inf);
  wire inf_in = a_inf | b_inf;
  wire zero_in = a_zero | b_zero;
  wire sy = sa ^ sb;

  // effective significands and exponents
  wire [10:0] ma = {(ea != 0), fa};
  wire [10:0] mb = {(eb != 0), fb};
  wire [5:0]  ea_eff = (ea == 0) ? 6'd1 : {1'b0, ea};
  wire [5:0]  eb_eff = (eb == 0) ? 6'd1 : {1'b0, eb};

  wire [21:0] prod = ma * mb;

  // leading zero count
  reg [4:0] lz;
  integer i;
  always @* begin
    lz = 5'd0;
    for (i = 0; i < 22; i = i + 1)
      if (prod[i]) lz = 5'd21 - i[4:0];
  end

  // below[k] = |prod[k:0]
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  function bit_at;
    input [21:0] v;
    input integer idx;
    begin
      bit_at = (idx >= 0 && idx < 22) ? v[idx] : 1'b0;
    end
  endfunction

  // Normal path: left-shift by lz, leading 1 at bit 21
  wire [21:0] norm = prod << lz;
  wire [10:0] mant_pre = norm[21:11];
  wire        guard_n  = norm[10];
  wire        sticky_n = |norm[9:0];
  wire        round_n  = guard_n & (sticky_n | mant_pre[0]);
  wire [11:0] mant_rnd = {1'b0, mant_pre} + {11'd0, round_n};
  wire        carry_n  = mant_rnd[11];
  wire [9:0]  frac_n   = carry_n ? mant_rnd[10:1] : mant_rnd[9:0];

  // Normal exponent
  wire [6:0] exp_sum = {1'b0, ea_eff} + {1'b0, eb_eff};
  wire [7:0] lhs     = {1'b0, exp_sum} + 8'd1 + {7'd0, carry_n};
  wire [7:0] lim_u   = 8'd15 + {3'd0, lz};
  wire [7:0] lim_o   = 8'd45 + {3'd0, lz};  // BIAS + MAX_FINITE_E = 15+30=45
  wire        under   = (lhs <= lim_u);
  wire        over    = (lhs >  lim_o);
  wire [7:0] e_norm   = lhs - lim_u;
  wire [4:0] exp_n    = e_norm[4:0];

  // Subnormal path
  integer           tsh;
  reg  [21:0]       sub_sh;
  always @* begin
    tsh = 26 - $signed({1'b0, exp_sum});
    sub_sh = (tsh >= 0 && tsh < 22) ? (prod >> tsh) : 22'd0;
  end
  wire [9:0]  frac_s0  = sub_sh[9:0];
  wire        guard_s  = bit_at(prod,  tsh - 1);
  wire        sticky_s = bit_at(below, tsh - 2);
  wire        round_s  = guard_s & (sticky_s | frac_s0[0]);
  wire [10:0] frac_ss  = {1'b0, frac_s0} + {10'd0, round_s};
  wire        carry_s  = frac_ss[10];
  wire [9:0]  frac_s   = frac_ss[9:0];

  // pack
  wire is_inf  = ~is_nan & (inf_in | (over & ~zero_in));
  wire is_sub  = ~is_nan & ~is_inf & under;
  wire sub_z   = is_sub & ~carry_s & (frac_s == 0);
  wire is_zero = ~is_nan & ~is_inf & (zero_in | sub_z);

  wire [4:0] exp_out  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : is_sub ? {4'd0, carry_s} : exp_n;
  wire [9:0] frac_out = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : is_sub ? frac_s : frac_n;

  assign y = {(is_nan ? 1'b0 : sy), exp_out, frac_out};
endmodule

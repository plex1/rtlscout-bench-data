// fp_mul_e5f10 -- IEEE binary16 multiplier (1 sign, 5 exponent, 10 fraction bits) with subnormal support.
//
// Structure:
//   - 11x11 significand product
//   - leading-one index of the product derived from the operands:
//       L = i + j + prod[i+j+1]   (i,j = leading-one indices of the significands)
//     so lz = 21 - L without a 22-bit priority encoder.
//   - normalised mantissa = window prod[L-10 .. L] obtained with a funnel shifter
//     (guard = prod[L-11], sticky = |prod[L-12:0]| from a prefix-OR table).
//   - subnormal results: window prod[ts .. ts+9] with ts = 26 - (ea_eff + eb_eff).
//   - exponent: e_norm = exp_sum + L + carry - 35 (one 7-bit adder, mod 128).

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- 1. unpack / classify ----
  wire       s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0],  f_b = b[9:0];

  wire e_a_zero = ~|e_a, e_b_zero = ~|e_b;
  wire e_a_ones =  &e_a, e_b_ones =  &e_b;
  wire f_a_zero = ~|f_a, f_b_zero = ~|f_b;
  wire zero_a = e_a_zero & f_a_zero, zero_b = e_b_zero & f_b_zero;
  wire nan_a  = e_a_ones & ~f_a_zero, nan_b = e_b_ones & ~f_b_zero;
  wire inf_a  = e_a_ones &  f_a_zero, inf_b = e_b_ones &  f_b_zero;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  // ---- 2. significand product ----
  wire        h_a = ~e_a_zero, h_b = ~e_b_zero;
  wire [10:0] m_a = {h_a, f_a}, m_b = {h_b, f_b};
  wire [21:0] prod = m_a * m_b;

  // ---- 3. leading-one index of each significand (0 when it is zero) ----
  function [3:0] pri10;
    input [9:0] f;
    begin
      pri10 = f[9] ? 4'd9 :
              f[8] ? 4'd8 :
              f[7] ? 4'd7 :
              f[6] ? 4'd6 :
              f[5] ? 4'd5 :
              f[4] ? 4'd4 :
              f[3] ? 4'd3 :
              f[2] ? 4'd2 :
              f[1] ? 4'd1 : 4'd0;
    end
  endfunction

  wire [3:0] i_idx = h_a ? 4'd10 : pri10(f_a);
  wire [3:0] j_idx = h_b ? 4'd10 : pri10(f_b);

  // leading one of prod sits at ij or ij+1 (prod < 2^(ij+2))
  wire [4:0] ij = i_idx + j_idx;
  wire       p1 = prod[ij + 5'd1];
  wire [4:0] L  = ij + {4'd0, p1};            // index of the product's leading one

  // ---- 4. normalisation window: mantissa = prod[L-10 .. L], guard = prod[L-11] ----
  wire [4:0] n5 = L - 5'd10;                  // (L-10) mod 32
  wire [3:0] n4 = n5[4] ? 4'd0 : n5[3:0];     // clamp; L<10 only for sub*sub -> result 0

  wire [22:0] U  = {prod, 1'b0};
  wire [22:0] u1 = n4[0] ? {1'b0,  U[22:1]}  : U;
  wire [22:0] u2 = n4[1] ? {2'b0,  u1[22:2]} : u1;
  wire [22:0] u4 = n4[2] ? {4'b0,  u2[22:4]} : u2;
  wire [22:0] u8 = n4[3] ? {8'b0,  u4[22:8]} : u4;

  wire [10:0] mant_pre = u8[11:1];            // prod[L-10 .. L]
  wire        guard    = u8[0];               // prod[L-11]

  // below[k] = |prod[k:0]|  (sticky lookup table)
  wire [21:0] below;
  genvar k;
  generate
    for (k = 0; k < 22; k = k + 1) begin : g_below
      assign below[k] = |prod[k:0];
    end
  endgenerate

  wire [3:0] sn_idx = n4 - 4'd2;
  wire       sticky = (n4 > 4'd1) ? below[sn_idx] : 1'b0;   // |prod[L-12:0]|

  // ---- 5. round the normal result (nearest, ties to even) ----
  wire round_up = guard & (sticky | mant_pre[0]);
  wire carry    = round_up & &mant_pre;       // mantissa overflowed to 2.0

  wire [10:0] cinc;
  generate
    for (k = 0; k < 11; k = k + 1) begin : g_inc
      if (k == 0) assign cinc[0] = round_up;
      else        assign cinc[k] = round_up & (&mant_pre[k-1:0]);
    end
  endgenerate
  wire [10:0] mant_post = carry ? {1'b1, 10'd0} : (mant_pre ^ cinc);

  // ---- 6. exponent bookkeeping ----
  wire [4:0] ea_eff = e_a | {4'd0, e_a_zero};         // subnormal -> effective exponent 1
  wire [4:0] eb_eff = e_b | {4'd0, e_b_zero};
  wire [5:0] exp_sum = ea_eff + eb_eff;

  // e_norm = exp_sum + L + carry - 35, computed mod 128 (range [-33, 47] fits)
  wire [6:0] e7 = {1'b0, exp_sum} + {2'b0, L} + {6'd0, carry} + 7'd93;
  wire underflow = e7[6] | (e7 == 7'd0);              // e_norm <= 0
  wire overflow  = ~e7[6] & (e7[5] | (e7[4] & &e7[3:0]));  // e_norm >= 31

  // ---- 7. subnormal-result path: frac = prod >> (26 - exp_sum), rounded once ----
  wire [7:0] ts8   = 8'd26 - {2'b0, exp_sum};         // 26 - exp_sum (signed)
  wire [7:0] ts1_8 = ts8 - 8'd1;
  wire [4:0] ts1   = ts1_8[7] ? 5'd0 : ts1_8[4:0];    // clamp to [0,23]

  wire [21:0] v1  = ts1[0] ? {1'b0,  prod[21:1]}  : prod;
  wire [21:0] v2  = ts1[1] ? {2'b0,  v1[21:2]}    : v1;
  wire [21:0] v4  = ts1[2] ? {4'b0,  v2[21:4]}    : v2;
  wire [21:0] v8  = ts1[3] ? {8'b0,  v4[21:8]}    : v4;
  wire [21:0] v16 = ts1[4] ? {16'b0, v8[21:16]}   : v8;

  wire [9:0] frac_trunc = v16[10:1];                  // prod[ts .. ts+9]
  wire       guard_s    = v16[0];                     // prod[ts-1]
  wire [4:0] ss_idx     = ts1 - 5'd1;
  wire       sticky_s   = (ts1 > 5'd0 && ts1 < 5'd23) ? below[ss_idx] : 1'b0;

  wire round_s = guard_s & (sticky_s | frac_trunc[0]);
  wire carry_s = round_s & &frac_trunc;               // rounds up to smallest normal

  wire [9:0] csub;
  generate
    for (k = 0; k < 10; k = k + 1) begin : g_sub
      if (k == 0) assign csub[0] = round_s;
      else        assign csub[k] = round_s & (&frac_trunc[k-1:0]);
    end
  endgenerate
  wire [9:0] frac_sub = carry_s ? 10'd0 : (frac_trunc ^ csub);
  wire [4:0] exp_sub  = {4'd0, carry_s};

  // ---- 8. pack ----
  wire is_nan      = nan_in;
  wire is_inf      = ~is_nan & (inf_in | (overflow & ~zero_in));
  wire is_sub      = ~is_nan & ~is_inf & underflow;
  wire sub_is_zero = is_sub & (exp_sub == 5'd0) & (frac_sub == 10'd0);
  wire is_zero     = ~is_nan & ~is_inf & (zero_in | sub_is_zero);

  wire [4:0] exp_y  = (is_nan | is_inf) ? 5'd31 : is_zero ? 5'd0 : is_sub ? exp_sub : e7[4:0];
  wire [9:0] frac_y = is_nan ? 10'h200 : (is_inf | is_zero) ? 10'd0 : is_sub ? frac_sub : mant_post[9:0];
  wire       sign_y = is_nan ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

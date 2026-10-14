// fp_mul_e5f10 -- IEEE binary16 multiplier (1/5/10), subnormal support.
// Structure:
//   p = m_a*m_b (22b).  With at least one normal operand p >= 2^10 so the
//   leading one sits in p[21:10]  ->  lz in 0..11.
//   result mantissa window = p >> S,  S = normal ? (11-lz) : (26-E)
//   (both branches: S in 0..24, 5 bits).  Rounding bits come from the same
//   shifted bundle, sticky from a prefix-OR (por) shifted alongside.
module fp_mul_e5f10 (a, b, y);
  input  wire [15:0] a, b;
  output wire [15:0] y;

  wire s_a = a[15], s_b = b[15];
  wire [4:0] e_a = a[14:10], e_b = b[14:10];
  wire [9:0] f_a = a[9:0],   f_b = b[9:0];

  wire eaz = ~|e_a, ebz = ~|e_b;
  wire eao =  &e_a, ebo =  &e_b;
  wire faz = ~|f_a, fbz = ~|f_b;
  wire z_a = eaz & faz, z_b = ebz & fbz;
  wire n_a = eao & ~faz, n_b = ebo & ~fbz;
  wire i_a = eao &  faz, i_b = ebo &  fbz;

  wire nan_in = n_a | n_b | (i_a & z_b) | (z_a & i_b);
  wire inf_in = i_a | i_b;
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a = {~eaz, f_a};
  wire [10:0] m_b = {~ebz, f_b};
  wire [5:0]  E   = (eaz ? 6'd1 : {1'b0, e_a}) + (ebz ? 6'd1 : {1'b0, e_b});

  wire [21:0] p = m_a * m_b;

  // ---- leading zero count of v = p[21:10] (v != 0 for all non-both-subnormal)
  function [3:0] lz4f;
    input [3:0] x;
    begin
      lz4f = x[3] ? 4'd0 : x[2] ? 4'd1 : x[1] ? 4'd2 : x[0] ? 4'd3 : 4'd4;
    end
  endfunction
  wire [11:0] v  = p[21:10];
  wire [3:0]  c0 = lz4f(v[11:8]);
  wire [3:0]  c1 = lz4f(v[7:4]);
  wire [3:0]  c2 = lz4f(v[3:0]);
  wire [4:0]  lz = (|v[11:8]) ? {1'b0, c0} :
                   (|v[7:4])  ? (5'd4 + {1'b0, c1}) :
                   (|v[3:0])  ? (5'd8 + {1'b0, c2}) : 5'd12;

  wire [5:0] lz6    = {1'b0, lz};
  wire       normal = (E >= (lz6 + 6'd15));
  wire [4:0] S      = normal ? (5'd11 - lz) : (5'd26 - E[4:0]);

  // ---- shifted bundle: sh = p<<1 >> S ; bits[11:1]=mantissa, [0]=guard, [1]=lsb
  wire [22:0] p2 = {p, 1'b0};
  wire [22:0] por2;
  assign por2[0] = p2[0];
  genvar g;
  generate for (g = 1; g < 23; g = g + 1) begin : g_por
    assign por2[g] = por2[g-1] | p2[g];
  end endgenerate
  wire [22:0] por2x = {por2[21:0], 1'b0};
  wire [22:0] sh  = p2    >> S;
  wire [22:0] sh2 = por2x >> S;

  wire [10:0] M_pre = sh[11:1];
  wire        guard = sh[0];
  wire        lsb   = sh[1];
  wire        sticky = sh2[0];
  wire        rup   = guard & (sticky | lsb);

  wire [11:0] Msum  = {1'b0, M_pre} + {11'd0, rup};
  wire        carry = Msum[11];

  wire signed [6:0] ef = $signed({1'b0, E}) - 7'sd14 - $signed({2'b0, lz}) + carry;
  wire [4:0] field = normal ? ef[4:0] : {4'd0, Msum[10]};
  wire [9:0] frac  = normal ? (carry ? Msum[10:1] : Msum[9:0]) : Msum[9:0];

  wire is_inf      = ~nan_in & (inf_in | (ef > 7'sd30));
  wire is_zero_res = ~nan_in & ~is_inf & (field == 5'd0) & (frac == 10'd0);

  wire zero_in     = z_a | z_b;
  wire is_zres     = ~nan_in & ~is_inf & (zero_in | ((field == 5'd0) & (frac == 10'd0)));
  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_zres ? 5'd0 : field;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_zres) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

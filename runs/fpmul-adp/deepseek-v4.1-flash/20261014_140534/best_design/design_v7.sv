// fp_mul_e5f10 -- IEEE binary16 multiplier (1/5/10), subnormal support.
// p = m_a*m_b; lz = leading zeros of p (0..22) computed from p[21:10].
// S = result shift = normal ? (11-lz) : (26-E)  (5 bits).
// sh = (p<<1) >> S : sh[11:1] mantissa, sh[0] guard, sh[1] lsb,
// sticky = OR of everything shifted out (via prefix-OR of (p<<1)).
// A zero result needs no special casing: p==0 gives field==0 and frac==0.
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
  wire [11:0] v = p[21:10];

  // 2-level leading-zero count of the 12-bit v (v==0 -> 12)
  function [2:0] lz6;
    input [5:0] x;
    begin
      lz6 = (|x[5:3]) ? (x[5] ? 3'd0 : x[4] ? 3'd1 : 3'd2)
                      : (|x[2:0]) ? (x[2] ? 3'd3 : x[1] ? 3'd4 : 3'd5) : 3'd6;
    end
  endfunction
  wire [2:0] c0 = lz6(v[11:6]);
  wire [2:0] c1 = lz6(v[5:0]);
  wire [4:0] lz = (|v[11:6]) ? {2'd0, c0} : (5'd6 + {2'd0, c1});

  wire [5:0] lz6e = {1'b0, lz};
  wire       normal = (E >= (lz6e + 6'd15)) & (|v);
  wire [4:0] S      = normal ? (5'd11 - lz) : (5'd26 - E[4:0]);

  wire [23:0] p2 = {p, 1'b0};
  wire [23:0] por2;
  assign por2[0] = p2[0];
  genvar g;
  generate for (g = 1; g < 24; g = g + 1) begin : g_por
    assign por2[g] = por2[g-1] | p2[g];
  end endgenerate
  wire [23:0] por2x = {por2[22:0], 1'b0};

  wire [23:0] sh  = p2    >> S[4:0];
  wire [23:0] sh2 = por2x >> S[4:0];

  wire [10:0] M_pre = sh[11:1];
  wire        guard = sh[0];
  wire        lsb   = sh[1];
  wire        sticky = sh2[0];
  wire        rup   = guard & (sticky | lsb);

  wire [11:0] Msum  = {1'b0, M_pre} + {11'd0, rup};
  wire        carry = Msum[11];
  wire signed [6:0] ef = $signed({1'b0, E}) - 7'sd14 - $signed({2'b0, lz}) + carry;

  wire [4:0] field  = normal ? ef[4:0] : {4'd0, Msum[10]};
  wire [9:0] frac   = normal ? (carry ? Msum[10:1] : Msum[9:0]) : Msum[9:0];

  wire is_inf = ~nan_in & (inf_in | (ef > 7'sd30));

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : field;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf ? 10'd0 : frac);
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

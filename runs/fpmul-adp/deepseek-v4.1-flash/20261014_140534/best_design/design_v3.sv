// fp_mul_e5f10 -- IEEE binary16 multiplier, single-shifter unified formulation.
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
  wire i_a = eao & faz,  i_b = ebo & fbz;

  wire nan_in = n_a | n_b | (i_a & z_b) | (z_a & i_b);
  wire inf_in = i_a | i_b;
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a = {~eaz, f_a};
  wire [10:0] m_b = {~ebz, f_b};
  wire [5:0]  ea  = eaz ? 6'd1 : {1'b0, e_a};
  wire [5:0]  eb  = ebz ? 6'd1 : {1'b0, e_b};
  wire [5:0]  E   = ea + eb;

  wire [21:0] p = m_a * m_b;

  // ---- balanced leading-zero count of p (0..22) -------------------------
  function [2:0] lz6;
    input [5:0] v;
    reg [2:0] hi, lo;
    begin
      hi = v[5] ? 3'd0 : v[4] ? 3'd1 : 3'd2;
      lo = v[2] ? 3'd3 : v[1] ? 3'd4 : 3'd5;
      lz6 = (|v[5:3]) ? hi : (|v[2:0]) ? lo : 3'd6;
    end
  endfunction
  function [2:0] lz4;
    input [3:0] v;
    begin
      lz4 = v[3] ? 3'd0 : v[2] ? 3'd1 : v[1] ? 3'd2 : v[0] ? 3'd3 : 3'd4;
    end
  endfunction

  wire [2:0] l3 = lz6(p[21:16]);
  wire [2:0] l2 = lz6(p[15:10]);
  wire [2:0] l1 = lz6(p[9:4]);
  wire [2:0] l0 = lz4(p[3:0]);
  wire nz3 = |p[21:16], nz2 = |p[15:10], nz1 = |p[9:4];
  wire [4:0] lz = nz3 ? {2'd0,l3} :
                  nz2 ? (5'd6  + {2'd0,l2}) :
                  nz1 ? (5'd12 + {2'd0,l1}) :
                        (5'd18 + {2'd0,l0}) ;

  // ---- prefix OR of p: bor[k] = |p[k:0] (k<=22) -------------------------
  wire [22:0] b0 = {1'b0, p};
  wire [22:0] b1 = b0 | (b0 << 1);
  wire [22:0] b2 = b1 | (b1 << 2);
  wire [22:0] b3 = b2 | (b2 << 4);
  wire [22:0] b4 = b3 | (b3 << 8);
  wire [22:0] bor = b4 | (b4 << 16);
  wire [8:0]  d2  = {4'd0, lz} * 2;
  wire [22:0] below = bor;

  // ---- exponent bookkeeping -------------------------------------------
  wire signed [6:0] ebase = $signed({1'b0, E}) - 7'sd14 - $signed({2'b0, lz});
  wire normal = (ebase >= 7'sd1);

  // ---- unified shift amount D = normal ? (11-lz) : (26-E) --------------
  wire [4:0] Dn = 5'd11 - lz;
  wire [4:0] Ds = 5'd26 - E[4:0];
  wire [4:0] D  = normal ? Dn : Ds;

  wire [4:0]  Dm1 = D - 5'd1;
  wire [21:0] sh  = p >> Dm1;
  wire        dz  = ~|D;
  wire [10:0] M_pre = dz ? p[10:0] : sh[11:1];
  wire        guard = ~dz & sh[0];
  wire        sticky = (D >= 5'd2) ? below[D-2] : 1'b0;
  wire        rup  = guard & (sticky | M_pre[0]);

  wire [11:0] Mfull = {1'b0, M_pre} + {11'd0, rup};
  wire        carry = Mfull[11];

  wire signed [6:0] ef = ebase + carry;
  wire [4:0] field = normal ? ef[4:0] : {4'd0, Mfull[10]};
  wire [9:0] frac  = normal ? (carry ? Mfull[10:1] : Mfull[9:0]) : Mfull[9:0];

  wire is_inf   = ~nan_in & (inf_in | (ef > 7'sd30));
  wire zero_res = ~nan_in & ~is_inf & (field == 5'd0) & (frac == 10'd0);

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : zero_res ? 5'd0 : field;
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | zero_res) ? 10'd0 : frac;
  wire       sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

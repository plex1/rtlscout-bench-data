// fp_mul_e5f10 -- IEEE binary16 multiplier with subnormal support.
// Strategy: compute the leading-one position of each input mantissa in parallel
// with the significand multiply, so the normalization shift amount is available
// almost immediately after the product.  Single rounding (round to nearest even).

module fp_mul_e5f10 (
  input  logic [15:0] a,
  input  logic [15:0] b,
  output logic [15:0] y
);

  wire        sgn = a[15] ^ b[15];
  wire [4:0]  ea  = a[14:10];
  wire [4:0]  eb  = b[14:10];
  wire [9:0]  fa  = a[9:0];
  wire [9:0]  fb  = b[9:0];

  wire eaz = ~|ea;
  wire ebz = ~|eb;
  wire faz = ~|fa;
  wire fbz = ~|fb;

  wire [10:0] ma = {~eaz, fa};
  wire [10:0] mb = {~ebz, fb};
  wire [21:0] prod = ma * mb;

  // ---------------------------------------------------------------- msb index
  function automatic [2:0] msb5 (input logic [4:0] v);
    logic h, l;
    h = v[4] | v[3];
    l = v[2] | v[1];
    msb5 = h ? (v[4] ? 3'd4 : 3'd3) : (l ? (v[2] ? 3'd2 : 3'd1) : 3'd0);
  endfunction

  function automatic [3:0] msb10 (input logic [9:0] v);
    msb10 = (|v[9:5]) ? ({1'b0, msb5(v[9:5])} + 4'd5) : {1'b0, msb5(v[4:0])};
  endfunction

  wire [3:0] la = eaz ? msb10(fa) : 4'd10;
  wire [3:0] lb = ebz ? msb10(fb) : 4'd10;
  wire [4:0] q  = la + lb;               // 0..20

  // ------------------------------------------------------------- exponents
  wire [5:0] t = (eaz ? 6'd1 : {1'b0, ea}) + (ebz ? 6'd1 : {1'b0, eb}); // 2..62
  wire [6:0] u = {1'b0, t} + {2'b0, q};  // 2..80
  wire [6:0] exp_base = u - 7'd35;       // 7-bit signed, -33..45

  wire normal_path = ~exp_base[6] & (|exp_base[5:0]);   // exp_base >= 1

  wire [6:0] A7 = {2'b0, q} - 7'd10;
  wire [6:0] B7 = 7'd26 - {1'b0, t};
  wire [4:0] s1 = normal_path ? A7[4:0] : B7[4:0];

  // ------------------------------------------------ stage-1 barrel shift
  // qx = {prod,1'b0} >> s1  (13 bits) : qx[0]=prod[s1-1], qx[12]=prod[s1+11]
  wire [22:0] pext = {prod, 1'b0};
  wire [22:0] sa = s1[0] ? (pext >> 1) : pext;
  wire [22:0] sb = s1[1] ? (sa   >> 2) : sa;
  wire [22:0] sc = s1[2] ? (sb   >> 4) : sb;
  wire [22:0] sd = s1[3] ? (sc   >> 8) : sc;
  wire [22:0] se = s1[4] ? (sd   >> 16) : sd;
  wire [12:0] qx = se[12:0];

  wire        ext   = normal_path & qx[12];
  wire [10:0] n_pre = ext ? qx[12:2] : qx[11:1];
  wire        guard = ext ? qx[1]   : qx[0];

  // ---------------------------------------------------------- sticky bit
  wire [22:0] bb = {1'b0, prod};
  wire [22:0] q1 = bb | {bb[21:0],  1'b0};
  wire [22:0] q2 = q1 | {q1[20:0],  2'b0};
  wire [22:0] q4 = q2 | {q2[18:0],  4'b0};
  wire [22:0] q8 = q4 | {q4[14:0],  8'b0};
  wire [22:0] q16 = q8 | {q8[6:0], 16'b0};
  wire [22:0] below = {q16[21:0], 1'b0};   // below[m] = |prod[m-1:0]

  wire [4:0] idx = (s1 == 5'd0) ? 5'd0 : (s1 - 5'd1);
  wire st0 = (s1 != 5'd0) & below[idx];
  wire sticky = ext ? (qx[0] | st0) : st0;

  wire rnd_up = guard & (sticky | n_pre[0]);

  // ------------------------------------------------------------- rounding
  wire [11:0] n_r = {1'b0, n_pre} + {11'd0, rnd_up};
  wire [9:0]  frac = n_r[11] ? n_r[10:1] : n_r[9:0];

  // ------------------------------------------------------------- exponent
  wire [6:0] exp_norm = exp_base + {6'd0, qx[12]} + {6'd0, n_r[11]};
  wire [6:0] exp_fin  = normal_path ? exp_norm : {4'b0, n_r[11], n_r[10]};
  wire       ovf      = exp_fin[6] | exp_fin[5] | (&exp_fin[4:0]);

  // ------------------------------------------------------------ specials
  wire eao = &ea;
  wire ebo = &eb;
  wire nan_in  = (eao & ~faz) | (ebo & ~fbz) | (eao & faz & ebz & fbz) | (ebo & fbz & eaz & faz);
  wire zero_in = (eaz & faz) | (ebz & fbz);
  wire inf_in  = ((eao & faz) | (ebo & fbz)) & ~nan_in;

  wire [15:0] y_val = (inf_in | ovf) ? {sgn, 5'h1F, 10'h0} : {sgn, exp_fin[4:0], frac};
  wire [15:0] y_fin = zero_in ? {sgn, 15'h0} : y_val;

  assign y = nan_in ? 16'h7E00 : y_fin;

endmodule

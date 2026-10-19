// fp_mul_e5f10 -- IEEE binary16 multiplier with subnormal support.
// v3: single barrel shifter with sticky injection, unified exponent path.

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

  // ---------------------------------------------------------- msb index
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

  wire [5:0] t   = (eaz ? 6'd1 : {1'b0, ea}) + (ebz ? 6'd1 : {1'b0, eb});
  wire [6:0] expb = {1'b0, t} + {2'b0, la} + {2'b0, lb} - 7'd35;

  wire is_norm = ~expb[6] & (|expb[5:0]);          // result exponent >= 1

  wire [5:0] A6 = {1'b0, la} + {1'b0, lb} - 6'd10; // q-10 : 0..10
  wire [5:0] B6 = 6'd26 - {1'b0, t};               // 26-t
  wire [4:0] s1 = is_norm ? A6[4:0] : B6[4:0];

  // ------------------------------------- barrel shift right with sticky
  wire [22:0] pext = {prod, 1'b0};
  wire [22:0] t1 = s1[0] ? (pext >> 1) : pext;
  wire        y1 = s1[0] & pext[0];
  wire [22:0] t2 = s1[1] ? (t1 >> 2) : t1;
  wire        y2 = y1 | (s1[1] & (|t1[1:0]));
  wire [22:0] t3 = s1[2] ? (t2 >> 4) : t2;
  wire        y3 = y2 | (s1[2] & (|t2[3:0]));
  wire [22:0] t4 = s1[3] ? (t3 >> 8) : t3;
  wire        y4 = y3 | (s1[3] & (|t3[7:0]));
  wire [22:0] t5 = s1[4] ? (t4 >> 16) : t4;
  wire        y5 = y4 | (s1[4] & (|t4[15:0]));

  wire [12:0] qx = t5[12:0];
  wire        ext = is_norm & qx[12];
  wire [10:0] n_pre = ext ? qx[12:2] : qx[11:1];
  wire        guard = ext ? qx[1]   : qx[0];
  wire        sticky = ext ? (qx[0] | y5) : y5;

  wire rnd = guard & (sticky | n_pre[0]);

  wire [11:0] n_r  = {1'b0, n_pre} + {11'd0, rnd};
  wire [9:0]  frac = n_r[11] ? n_r[10:1] : n_r[9:0];

  // ----------------------------------------------------------- exponent
  wire [6:0] c0  = expb[6] ? 7'd0 : expb;
  wire [6:0] en  = c0 + {6'd0, ext} + {6'd0, n_r[11]};
  wire [6:0] expf = is_norm ? en : {6'd0, n_r[10]};
  wire       ovf  = expf[6] | expf[5] | (&expf[4:0]);

  // ----------------------------------------------------------- specials
  wire eao = &ea;
  wire ebo = &eb;
  wire nan_in = (eao & ~faz) | (ebo & ~fbz)
              | (eao & faz & ebz & fbz) | (ebo & fbz & eaz & faz);
  wire zero_in = (eaz & faz) | (ebz & fbz);
  wire inf_in  = (eao & faz) | (ebo & fbz);

  wire sat = ovf | zero_in | inf_in;
  wire ones = ovf | inf_in;
  wire [4:0]  exp_o  = (nan_in | ones) ? 5'h1F : sat ? 5'd0 : expf[4:0];
  wire [9:0]  frac_o = nan_in ? 10'h200 : sat ? 10'd0 : frac;
  wire        sgn_o  = sgn & ~nan_in;

  assign y = {sgn_o, exp_o, frac_o};

endmodule

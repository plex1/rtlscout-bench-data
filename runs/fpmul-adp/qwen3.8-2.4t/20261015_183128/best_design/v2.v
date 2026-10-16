// fp_mul_e5f10 v2 -- input-side subnormal normalization.
// Subnormals are normalized before the multiply, so the 22-bit product has at
// most one leading zero (p[21]); no wide LZC or wide barrel shifter is needed.

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  // ---- unpack ----
  wire        sa = a[15], sb = b[15];
  wire [4:0]  ea = a[14:10], eb = b[14:10];
  wire [9:0]  fa = a[9:0],   fb = b[9:0];

  // ---- classify ----
  wire ez_a = ~|ea, ez_b = ~|eb;          // exponent field zero
  wire eo_a = &ea,  eo_b = &eb;           // exponent field all ones
  wire fz_a = ~|fa, fz_b = ~|fb;          // fraction zero
  wire zero_a = ez_a & fz_a,  zero_b = ez_b & fz_b;
  wire inf_a  = eo_a & fz_a,  inf_b  = eo_b & fz_b;
  wire nan_a  = eo_a & ~fz_a, nan_b  = eo_b & ~fz_b;
  wire sub_a  = ez_a & ~fz_a, sub_b  = ez_b & ~fz_b;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire sy      = sa ^ sb;

  // ---- leading-zero count of a 10-bit fraction ----
  function automatic [3:0] lzc10;
    input [9:0] x;
    begin
      lzc10 = x[9] ? 4'd0 : x[8] ? 4'd1 : x[7] ? 4'd2 : x[6] ? 4'd3 :
              x[5] ? 4'd4 : x[4] ? 4'd5 : x[3] ? 4'd6 : x[2] ? 4'd7 :
              x[1] ? 4'd8 : 4'd9;
    end
  endfunction

  wire [3:0] lz_a = lzc10(fa);
  wire [3:0] lz_b = lzc10(fb);

  // ---- normalized significands (leading 1 at bit 10) + effective biased exp ----
  wire [10:0] ma = sub_a ? ({fa, 1'b0} << lz_a) : {1'b1, fa};
  wire [10:0] mb = sub_b ? ({fb, 1'b0} << lz_b) : {1'b1, fb};

  wire signed [6:0] ea_e = sub_a ? -$signed({3'b000, lz_a}) : $signed({2'b00, ea});
  wire signed [6:0] eb_e = sub_b ? -$signed({3'b000, lz_b}) : $signed({2'b00, eb});

  // ---- multiply ----
  wire [21:0] p    = ma * mb;
  wire        p_hi = p[21];
  wire [20:0] mant = p_hi ? p[21:1] : p[20:0];   // leading 1 at bit 20

  wire signed [6:0] t = ea_e + eb_e + {6'd0, p_hi};

  wire        is_sub  = (t <= 7'sd15);   // E = t-15 <= 0 -> subnormal grid
  wire        ovf_pre = (t >= 7'sd46);
  wire        t_eq45  = (t == 7'sd45);

  // ---- subnormal result path: staged right shifter with sticky in parallel ----
  // sel = 15 - t  (= -E), clamped; result is zero once sel >= 11.
  wire [5:0] sel6 = 6'd15 - t[5:0];
  wire [4:0] sel  = (sel6 > 6'd15) ? 5'd15 : sel6[4:0];

  wire [10:0] e0 = mant[20:10];
  wire [10:0] e1 = sel[0] ? (e0 >> 1) : e0;
  wire        o1 = sel[0] & e0[0];
  wire [10:0] e2 = sel[1] ? (e1 >> 2) : e1;
  wire        o2 = sel[1] & (|e1[1:0]);
  wire [10:0] e4 = sel[2] ? (e2 >> 4) : e2;
  wire        o4 = sel[2] & (|e2[3:0]);
  wire [10:0] e8 = sel[3] ? (e4 >> 8) : e4;
  wire        o8 = sel[3] & (|e4[7:0]);

  // p[0] is dropped from 'mant' when p_hi is set; it must still feed the sticky
  wire p0 = p_hi & p[0];

  wire [9:0] frac_s  = e8[10:1];
  wire       guard_s = e8[0];
  wire       sticky_s = (|mant[9:0]) | o1 | o2 | o4 | o8 | p0;

  // ---- normal result path fields ----
  wire [9:0] frac_n  = mant[19:10];
  wire       guard_n = mant[9];
  wire       sticky_n = (|mant[8:0]) | p0;

  // ---- shared rounding (round to nearest, ties to even) ----
  wire [9:0] frac_pre = is_sub ? frac_s  : frac_n;
  wire       guard    = is_sub ? guard_s : guard_n;
  wire       sticky   = is_sub ? sticky_s : sticky_n;
  wire       round_up = guard & (sticky | frac_pre[0]);

  wire [10:0] frac11  = {1'b0, frac_pre} + {10'd0, round_up};
  wire        carry   = frac11[10];
  wire [9:0]  frac_o  = frac11[9:0];

  // ---- final exponent / flags ----
  wire [4:0] exp_norm = t[4:0] - 5'd15 + {4'd0, carry};
  wire [4:0] exp_o    = is_sub ? {4'd0, carry} : exp_norm;
  wire       ovf      = ovf_pre | (t_eq45 & carry);

  // ---- pack ----
  wire [15:0] core = ovf ? {sy, 5'h1F, 10'd0} : {sy, exp_o, frac_o};

  assign y = nan_in  ? 16'h7E00 :
             inf_in  ? {sy, 5'h1F, 10'd0} :
             zero_in ? {sy, 15'd0} :
             core;
endmodule

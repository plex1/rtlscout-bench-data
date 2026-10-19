// fp_mul_e5f10 -- IEEE binary16 multiplier (1s,5e,10f) with subnormal support.
//
//   value = ma*mb * 2^(esum-50),  p = ma*mb,  lz = lz(ma)+lz(mb)
//   Y = (2*p) >> s,  s = 10-lz (normal) or 26-esum (subnormal)
// The leading one of Y sits at bit 11 or 12 for a normal result (fixed up by the
// single select on Y[12]) and at bit 10 for a subnormal one; Y[0] is the round
// bit in both cases.  lz is computed from the raw significands, i.e. in parallel
// with the multiplier.
module fp_mul_e5f10 (input [15:0] a, input [15:0] b, output [15:0] y);

  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];

  wire ea_z = ~|ea, eb_z = ~|eb, ea_o = &ea, eb_o = &eb, fa_z = ~|fa, fb_z = ~|fb;
  wire inf_a  = ea_o &  fa_z;
  wire inf_b  = eb_o &  fb_z;
  wire zero_a = ea_z &  fa_z;
  wire zero_b = eb_z &  fb_z;
  wire nan_in = (ea_o & ~fa_z) | (eb_o & ~fb_z) | (inf_a & zero_b) | (zero_a & inf_b);

  // ---------------- leading-zero count (parallel with the multiplier) -------
  function [3:0] lz4(input [3:0] v);
    lz4 = v[3] ? 4'd0 : v[2] ? 4'd1 : v[1] ? 4'd2 : v[0] ? 4'd3 : 4'd4;
  endfunction

  function [3:0] lz11(input [10:0] v);
    lz11 = |v[10:7] ? lz4(v[10:7]) :
           |v[6:3]  ? (4'd4 + lz4(v[6:3])) :
           v[2] ? 4'd8 : v[1] ? 4'd9 : v[0] ? 4'd10 : 4'd11;
  endfunction

  wire [10:0] ma = {~ea_z, fa};
  wire [10:0] mb = {~eb_z, fb};
  wire [4:0]  lz = {1'b0, lz11(ma)} + {1'b0, lz11(mb)};

  wire [21:0] p = ma * mb;

  wire [5:0] esum = (ea_z ? 6'd1 : {1'b0, ea}) + (eb_z ? 6'd1 : {1'b0, eb});
  wire [6:0] W    = {1'b0, esum} - {2'd0, lz};
  wire       sub  = W[6] | ~|W[6:4];           // W <= 15 (or negative)

  wire [5:0] sn6 = 6'd10 - {1'b0, lz};
  wire [4:0] s_n = sn6[5] ? 5'd0 : sn6[4:0];   // clamp (zero operands only)
  wire [5:0] ss6 = 6'd26 - esum;               // 0..24 for subnormal results
  wire [4:0] s_s = ss6[5] ? 5'd31 : ss6[4:0];
  wire [4:0] s   = sub ? s_s : s_n;

  // ---------------- shifter with sticky accumulation ------------------------
  wire [22:0] X  = {p, 1'b0};
  wire [22:0] t1 = s[0] ? {1'b0,  X[22:1]}  : X;
  wire        q0 = s[0] & X[0];
  wire [22:0] t2 = s[1] ? {2'b0,  t1[22:2]} : t1;
  wire        q1 = s[1] & (|t1[1:0]);
  wire [22:0] t3 = s[2] ? {4'b0,  t2[22:4]} : t2;
  wire        q2 = s[2] & (|t2[3:0]);
  wire [22:0] t4 = s[3] ? {8'b0,  t3[22:8]} : t3;
  wire        q3 = s[3] & (|t3[7:0]);
  wire [22:0] t5 = s[4] ? {16'b0, t4[22:16]} : t4;
  wire        q4 = s[4] & (|t4[15:0]);

  wire [12:0] Y  = t5[12:0];
  wire        st = q0 | q1 | q2 | q3 | q4;

  // ---------------- round to nearest, ties to even -------------------------
  wire        e   = ~sub & Y[12];               // normal result needs one more shift
  wire [10:0] M   = e ? Y[12:2] : Y[11:1];
  wire        g   = e ? Y[1] : Y[0];
  wire        stk = st | (e & Y[0]);
  wire        rup = g & (stk | M[0]);
  wire [11:0] S   = {1'b0, M} + {11'd0, rup};

  // ---------------- exponent ------------------------------------------------
  wire [5:0] En = W[5:0] - 6'd15 + {5'd0, e} + {5'd0, S[11]};

  wire [6:0] exp_pre = sub ? {6'd0, S[10]} : {1'b0, En};
  wire       ovf = inf_a | inf_b | (~sub & (En[5] | (&En[4:0])));

  // ---------------- packing ------------------------------------------------
  wire [4:0] exp_y  = (nan_in | ovf) ? 5'd31 : (zero_a | zero_b) ? 5'd0 : exp_pre[4:0];
  wire [9:0] frac_y = nan_in ? 10'h200 : ovf ? 10'd0 : S[9:0];
  wire       sign_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {sign_y, exp_y, frac_y};
endmodule

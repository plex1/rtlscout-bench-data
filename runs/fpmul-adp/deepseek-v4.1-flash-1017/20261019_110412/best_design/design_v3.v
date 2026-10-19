// fp_mul_e5f10 -- IEEE binary16 multiplier, area/delay optimised.
// Significands (with implicit bit) are multiplied raw; the 22-bit product of two
// normalised significands has its leading one at bit 21 or 20, so post-multiply
// normalisation is a single-bit select.
module fp_mul_e5f10 (input [15:0] a, input [15:0] b, output [15:0] y);

  wire       sa = a[15], sb = b[15];
  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];

  wire ea_z = ~|ea, eb_z = ~|eb;
  wire ea_o = &ea,  eb_o = &eb;
  wire fa_z = ~|fa, fb_z = ~|fb;

  wire inf_a  = ea_o &  fa_z;
  wire inf_b  = eb_o &  fb_z;
  wire zero_a = ea_z &  fa_z;
  wire zero_b = eb_z &  fb_z;
  wire nan_a  = ea_o & ~fa_z;
  wire nan_b  = eb_o & ~fb_z;

  wire nan_in = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);

  // ---- significands (bit 10 = implicit one, cleared for subnormals) ----
  wire [10:0] ma = {~ea_z, fa};
  wire [10:0] mb = {~eb_z, fb};

  wire [21:0] p   = ma * mb;
  wire        ovf = p[21];

  wire [10:0] M  = ovf ? p[21:11] : p[20:10];
  wire        g0 = ovf ? p[10]    : p[9];
  wire        s0 = ovf ? |p[9:0]  : |p[8:0];

  // ---- exponent ------------------------------------------------------
  // biased exponent of the normalised result before rounding
  wire [6:0] ebr = {2'd0, ea} + {2'd0, eb} - 7'd15 + {6'd0, ovf}
                   + ({2'd0, ea_z} & 7'd1) + ({2'd0, eb_z} & 7'd1);
  wire       sub = ebr[6] | ~|ebr;          // ebr <= 0  -> subnormal result

  wire [5:0] s6 = 6'd1 - ebr[5:0];          // right shift for the subnormal grid
  wire [3:0] s4 = (|s6[5:4]) ? 4'd15 : s6[3:0];

  wire [10:0] Ms = M >> s4;
  wire [9:0]  frac_trunc = sub ? Ms[9:0] : M[9:0];

  wire guard  = sub ? M[s4-1] : g0;
  wire sticky = sub ? (|(M & ((11'd1 << (s4-1)) - 11'd1)) | g0 | s0) : s0;
  wire round_up = guard & (sticky | frac_trunc[0]);

  wire [10:0] frac_sum   = {1'b0, frac_trunc} + {10'd0, round_up};
  wire        carry      = frac_sum[10];
  wire [9:0]  frac_final = frac_sum[9:0];

  wire [6:0] exp_pre = (sub ? 7'd0 : ebr) + {6'd0, carry};

  wire is_inf  = ~nan_in & (inf_a | inf_b | (exp_pre > 7'd30));
  wire is_zero = ~nan_in & ~is_inf &
                 (zero_a | zero_b | ((exp_pre == 7'd0) & ~|frac_final));

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_pre[4:0];
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_final;
  wire       sign_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {sign_y, exp_y, frac_y};
endmodule

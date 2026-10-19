// fp_mul_e5f10 -- IEEE binary16 multiplier, restructured for area/delay.
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

  // inf*0 (either order) -> NaN ; any NaN input -> NaN
  wire nan_in = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);

  // ---- significands, normalised so that bit 10 is the leading one ----
  wire [10:0] ma = {~ea_z, fa};
  wire [10:0] mb = {~eb_z, fb};

  function [3:0] lz11(input [10:0] v);
    begin
      lz11 = v[10] ? 4'd0 : v[9] ? 4'd1 : v[8] ? 4'd2 : v[7] ? 4'd3 :
             v[6] ? 4'd4 : v[5] ? 4'd5 : v[4] ? 4'd6 : v[3] ? 4'd7 :
             v[2] ? 4'd8 : v[1] ? 4'd9 : v[0] ? 4'd10 : 4'd11;
    end
  endfunction

  wire [3:0]  lza = lz11(ma);
  wire [3:0]  lzb = lz11(mb);
  wire [10:0] na  = ma << lza;
  wire [10:0] nb  = mb << lzb;

  wire [21:0] p   = na * nb;      // leading one at bit 20 or 21
  wire        ovf = p[21];

  wire [10:0] M  = ovf ? p[21:11] : p[20:10];
  wire        g0 = ovf ? p[10]    : p[9];
  wire        s0 = ovf ? |p[9:0]  : |p[8:0];

  // ---- exponent ------------------------------------------------------
  wire [5:0] ea_e = ea_z ? 6'd1 : {1'b0, ea};
  wire [5:0] eb_e = eb_z ? 6'd1 : {1'b0, eb};
  // biased exponent of the normalised result, before rounding
  wire [6:0] ebr  = {1'b0, ea_e} + {1'b0, eb_e} - 7'd15 - {3'd0, lza} - {3'd0, lzb}
                    + {6'd0, ovf};
  wire       sub  = ebr[6] | ~|ebr;         // ebr <= 0  -> subnormal result

  wire [5:0] s6   = 6'd1 - ebr[5:0];        // right shift for the subnormal grid
  wire [3:0] s4   = (|s6[5:4]) ? 4'd15 : s6[3:0];

  wire [10:0] Ms  = M >> s4;
  wire [9:0]  frac_trunc = sub ? Ms[9:0] : M[9:0];

  wire guard  = sub ? M[s4-1] : g0;
  wire sticky = sub ? (|(M & ((11'd1 << (s4-1)) - 11'd1)) | g0 | s0) : s0;
  wire round_up = guard & (sticky | frac_trunc[0]);

  wire [10:0] frac_sum  = {1'b0, frac_trunc} + {10'd0, round_up};
  wire        carry     = frac_sum[10];
  wire [9:0]  frac_final= frac_sum[9:0];

  wire [6:0]  exp_pre   = (sub ? 7'd0 : ebr) + {6'd0, carry};

  wire is_inf  = ~nan_in & (inf_a | inf_b | (exp_pre > 7'd30));
  wire is_zero = ~nan_in & ~is_inf &
                 (zero_a | zero_b | ((exp_pre == 7'd0) & ~|frac_final));

  wire [4:0] exp_y  = (nan_in | is_inf) ? 5'd31 : is_zero ? 5'd0 : exp_pre[4:0];
  wire [9:0] frac_y = nan_in ? 10'h200 : (is_inf | is_zero) ? 10'd0 : frac_final;
  wire       sign_y = nan_in ? 1'b0 : (sa ^ sb);

  assign y = {sign_y, exp_y, frac_y};
endmodule

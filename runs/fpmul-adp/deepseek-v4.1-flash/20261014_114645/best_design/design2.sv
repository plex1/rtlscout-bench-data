// fp_mul_e5f10 -- IEEE binary16 multiplier, unified single-shift architecture.
//
// Core idea: P = m_a*m_b (22 bits).  Let p = position of P's leading one.
// Shifting P right by d = max(p-10, 26-esum) puts the significand on the
// normal grid (leading one at bit 10) or, when 26-esum dominates, directly on
// the subnormal grid.  The dominant term is exactly the "is subnormal result"
// predicate, so both cases share one shifter, one rounder and one pack stage.
module fp_mul_e5f10 (
  input  [15:0] a,
  input  [15:0] b,
  output [15:0] y
);
  wire        s_a = a[15], s_b = b[15];
  wire [4:0]  e_a = a[14:10], e_b = b[14:10];
  wire [9:0]  f_a = a[9:0],  f_b = b[9:0];

  wire e_az = ~|e_a, e_bz = ~|e_b;
  wire e_am = &e_a,  e_bm = &e_b;
  wire f_az = ~|f_a, f_bz = ~|f_b;

  wire zero_a = e_az & f_az;
  wire zero_b = e_bz & f_bz;
  wire inf_a  = e_am & f_az;
  wire inf_b  = e_bm & f_bz;
  wire nan_a  = e_am & ~f_az;
  wire nan_b  = e_bm & ~f_bz;

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a = {~e_az, f_a};
  wire [10:0] m_b = {~e_bz, f_b};
  wire [21:0] P   = m_a * m_b;

  // effective exponent sum (subnormal / zero uses 1)
  wire [5:0] esum = (e_az ? 6'd1 : {1'b0, e_a}) + (e_bz ? 6'd1 : {1'b0, e_b});

  // ---- leading-one position of P (0 if P==0): 5-level binary search ----
  reg  [4:0]  ppos;
  reg  [21:0] tt;
  always @* begin
    if (|P[21:11]) begin tt = P >> 11; ppos = 5'd11; end
    else           begin tt = P;        ppos = 5'd0;  end
    if (|tt[10:6]) begin tt = tt >> 5; ppos = ppos + 5'd5; end
    if (|tt[5:3])  begin tt = tt >> 3; ppos = ppos + 5'd3; end
    if (tt[2])     begin tt = tt >> 2; ppos = ppos + 5'd2; end
    if (tt[1])     ppos = ppos + 5'd1;
  end

  // ---- unified right shift d = max(p-10, 26-esum) ----
  wire signed [6:0] sh_n = {2'b00, ppos} - 7'sd10;
  wire signed [6:0] sh_s = 7'sd26 - {1'b0, esum};
  wire              is_sub = sh_n < sh_s;
  wire [4:0]        d      = is_sub ? sh_s[4:0] : sh_n[4:0];

  // ---- barrel shifter on P2={P,0}: R[k]=P[k+d-1], S=|P[d-2:0]| ----
  wire [22:0] P2 = {P, 1'b0};
  wire [22:0] q1 = d[0] ? (P2 >> 1)  : P2;
  wire [22:0] q2 = d[1] ? (q1 >> 2)  : q1;
  wire        k2 = d[1] & (|q1[1:0]);
  wire [22:0] q4 = d[2] ? (q2 >> 4)  : q2;
  wire        k4 = d[2] & (|q2[3:0]);
  wire [22:0] q8 = d[3] ? (q4 >> 8)  : q4;
  wire        k8 = d[3] & (|q4[7:0]);
  wire [22:0] R  = d[4] ? (q8 >> 16) : q8;
  wire        k16 = d[4] & (|q8[15:0]);
  wire        S  = k2 | k4 | k8 | k16;

  wire [10:0] mant    = R[11:1];
  wire        G       = R[0];
  wire        rup     = G & (S | mant[0]);
  wire [11:0] rounded = {1'b0, mant} + {11'd0, rup};
  wire [9:0]  frac    = rounded[9:0];

  // ---- exponent ----
  wire signed [7:0] E    = {3'b000, ppos} + {2'b00, esum} - 8'sd50;
  wire signed [7:0] yexp = E + 8'sd15 + {7'd0, rounded[11]};
  wire ovf    = ~is_sub & (yexp > 8'sd30);
  wire [4:0] e_res = is_sub ? {4'd0, rounded[10]} : (ovf ? 5'd31 : yexp[4:0]);

  // ---- pack ----
  wire is_nan = nan_in;
  wire is_inf = inf_in | ovf;
  wire [15:0] y_norm = {s_y, e_res, frac};
  wire [15:0] y_zero = {s_y, 5'd0,  10'd0};
  wire [15:0] y_inf  = {s_y, 5'd31, 10'd0};
  assign y = is_nan ? 16'h7E00 : is_inf ? y_inf : zero_in ? y_zero : y_norm;
endmodule

// fp_mul_e5f10 -- unified-shift, exponent derived from the shift amount d.
// E = esum + d - 40 + round_carry  holds for both the normal and the subnormal
// grid, so the normaliser position never enters the exponent path.
module fp_mul_e5f10 (
  input  [15:0] a,
  input  [15:0] b,
  output [15:0] y
);
  wire        s_a = a[15], s_b = b[15];
  wire [4:0]  e_a = a[14:10], e_b = b[14:10];
  wire [9:0]  f_a = a[9:0],  f_b = b[9:0];

  wire e_az = ~|e_a, e_bz = ~|e_b, e_am = &e_a, e_bm = &e_b;
  wire f_az = ~|f_a, f_bz = ~|f_b;

  wire zero_a = e_az & f_az, zero_b = e_bz & f_bz;
  wire inf_a  = e_am & f_az, inf_b  = e_bm & f_bz;
  wire nan_in = (e_am & ~f_az) | (e_bm & ~f_bz) | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y    = s_a ^ s_b;

  wire [10:0] m_a = {~e_az, f_a};
  wire [10:0] m_b = {~e_bz, f_b};
  wire [21:0] P   = m_a * m_b;

  // leading-one position of P (0 if P==0): 5-step binary search
  reg  [4:0]  p;
  reg  [21:0] t;
  always @* begin
    if (|P[21:11]) begin t = P >> 11; p = 5'd11; end
    else           begin t = P;       p = 5'd0;  end
    if (|t[10:6]) begin t = t >> 5; p = p + 5'd5; end
    if (|t[5:3])  begin t = t >> 3; p = p + 5'd3; end
    if (t[2])     begin t = t >> 2; p = p + 5'd2; end
    if (t[1])     p = p + 5'd1;
  end

  wire [5:0] esum = (e_az ? 6'd1 : {1'b0, e_a}) + (e_bz ? 6'd1 : {1'b0, e_b});
  wire [6:0] pse  = {2'b0, p} + {1'b0, esum};
  wire [4:0] dn   = p - 5'd10;
  wire [5:0] ds   = 6'd26 - esum;
  wire       is_sub = (pse < 7'd36);
  wire [4:0] d    = is_sub ? ds[4:0] : dn;

  // barrel shifter on P2={P,0}: R[k]=P[k+d-1], S=|P[d-2:0]|
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

  wire [10:0] mant = R[11:1];
  wire        G    = R[0];
  wire        rup  = G & (S | mant[0]);
  wire [11:0] I    = {1'b0, mant} + {11'd0, rup};
  wire        add  = is_sub ? I[10] : I[11];
  wire [9:0]  frac = (I[11] & ~is_sub) ? I[10:1] : I[9:0];

  // exponent: subnormal grid pins the biased exponent to 0; the normal case is
  // the biased E = p + esum + 15 - 50.  Rounding can add one more.
  wire signed [6:0] En = pse - 7'd35;
  wire signed [6:0] Ew = (is_sub ? 7'sd0 : En) + {6'd0, add};
  wire       ovf = ~is_sub & (Ew > 7'sd30);
  wire       oinf = inf_in | ovf;
  wire [4:0] e_res = Ew[4:0];

  assign y = nan_in ? 16'h7E00
           : oinf   ? {s_y, 5'd31, 10'd0}
           : zero_in? {s_y, 5'd0,  10'd0}
           : {s_y, e_res, frac};
endmodule

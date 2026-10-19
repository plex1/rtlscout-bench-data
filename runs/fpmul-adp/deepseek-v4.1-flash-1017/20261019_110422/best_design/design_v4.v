// fp_mul_e5f10 -- binary16 (e5f10) multiplier.
//
// value(a) = ma * 2^(ea_eff-25), ma = 11b (hidden bit 0 for subnormals)
// P  = ma*mb (22b)   ;  value = P * 2^(S-50),  S = ea_eff + eb_eff
// A  = lz(ma) + lz(mb)  -> Pp = P<<A has its MSB at bit 20 or 21 (ov = Pp[21])
// Normal   : frac = window[10:1], guard = window[0], E = S-15-A+ov+carry
// Subnormal: additional right shift m = 16+A-ov-S (>=1), E = carry
// One shared shift, one shared rounder.

module fp_mul_e5f10 (
  input  wire [15:0] a,
  input  wire [15:0] b,
  output wire [15:0] y
);

  // ---------------- leading-zero count of an 11-bit significand ----------
  function automatic [3:0] lz11 (input [10:0] v);
    begin
      if (|v[10:8]) begin
        if (v[10]) lz11 = 4'd0;
        else if (v[9]) lz11 = 4'd1;
        else lz11 = 4'd2;
      end else if (|v[7:5]) begin
        if (v[7]) lz11 = 4'd3;
        else if (v[6]) lz11 = 4'd4;
        else lz11 = 4'd5;
      end else if (|v[4:2]) begin
        if (v[4]) lz11 = 4'd6;
        else if (v[3]) lz11 = 4'd7;
        else lz11 = 4'd8;
      end else if (v[1]) lz11 = 4'd9;
      else if (v[0]) lz11 = 4'd10;
      else lz11 = 4'd11;
    end
  endfunction

  wire [4:0] ea = a[14:10], eb = b[14:10];
  wire [9:0] fa = a[9:0],   fb = b[9:0];
  wire ea0 = ~(|ea), eb0 = ~(|eb);
  wire ea1 =  (&ea), eb1 =  (&eb);
  wire fa0 = ~(|fa), fb0 = ~(|fb);

  // ---------------- specials (input only, early) ------------------------
  wire zero_a = ea0 & fa0;
  wire zero_b = eb0 & fb0;
  wire inf_a  = ea1 & fa0;
  wire inf_b  = eb1 & fb0;
  wire nan_in  = (ea1 & ~fa0) | (eb1 & ~fb0) | (inf_a & zero_b) | (inf_b & zero_a);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;

  // ---------------- significands, product, normalisation ----------------
  wire [10:0] ma = {~ea0, fa};
  wire [10:0] mb = {~eb0, fb};
  wire [3:0]  lza = lz11(ma);
  wire [3:0]  lzb = lz11(mb);
  wire [4:0]  A   = lza + lzb;
  wire [5:0]  S   = {1'b0, (ea0 ? 5'd1 : ea)} + {1'b0, (eb0 ? 5'd1 : eb)};

  wire [21:0] P  = ma * mb;
  wire [21:0] Pp = P << A;

  wire        ov = Pp[21];
  wire [12:0] Aw = ov ? {1'b0, Pp[21:10]} : Pp[20:9];
  wire        SL = (|Pp[8:0]) | (ov & Pp[9]);

  // extra subnormal right shift
  wire signed [6:0] md = 7'sd16 + {2'b0, A} - {6'd0, ov} - {1'b0, S};
  wire [4:0] m = zero_in ? 5'd16 :
                 md[6]   ? 5'd0  :
                 (|md[5:4]) ? 5'd16 : md[3:0];

  wire [12:0] B = Aw >> m;
  wire [9:0]  frac  = B[10:1];
  wire        guard = B[0];
  wire        mlsb  = B[1];

  // sticky: |Aw[m-1:0]
  wire [12:0] q1 = Aw | {Aw[11:0], 1'b0};
  wire [12:0] q2 = q1 | {q1[10:0], 2'b0};
  wire [12:0] q3 = q2 | {q2[8:0],  4'b0};
  wire [12:0] q4 = q3 | {q3[4:0],  8'b0};      // q4[i] = |Aw[i:0]
  wire [31:0] pt = {15'd0, q4, 1'b0};          // pt[i] = |Aw[i-1:0]
  wire        stk = pt[m] | SL;

  wire rnd = guard & (stk | mlsb);

  // ---------------- rounding (increment computed in parallel) -----------
  wire [10:0] fracp1 = {1'b0, frac} + 11'd1;
  wire        cout   = rnd & (&frac);
  wire [9:0]  frac_r = rnd ? fracp1[9:0] : frac;

  // ---------------- exponent --------------------------------------------
  wire signed [6:0] E0 = {1'b0, S} - 7'sd15 - {2'b0, A} + {6'd0, ov};
  wire signed [6:0] En = E0 + {6'd0, cout};
  wire ovf = (En > 7'sd30);
  wire [4:0] e_out = (m == 5'd0) ? En[4:0] : {4'd0, cout};

  // ---------------- pack ------------------------------------------------
  wire sat = nan_in | inf_in | ovf;
  wire [4:0] ey = sat ? 5'h1f : e_out;
  wire [9:0] fy = nan_in ? 10'h200 : (sat ? 10'd0 : frac_r);
  wire       sy = nan_in ? 1'b0 : (a[15] ^ b[15]);

  assign y = {sy, ey, fy};

endmodule

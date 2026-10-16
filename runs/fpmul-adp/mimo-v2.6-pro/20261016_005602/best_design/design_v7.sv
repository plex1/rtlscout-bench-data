// fp_mul_e5f10 -- IEEE binary16 (1+5+10) multiplier, subnormals supported, RNE.
// v7: raw-significand multiply (multiplier inputs are wiring), one unified
// window extractor serves both the normal and subnormal result paths.

module lzc11 (input [10:0] r, output [3:0] c);
  wire [15:0] F  = {r, 5'b0};
  wire        c3 = ~|F[15:8];
  wire [7:0]  F2 = c3 ? F[7:0] : F[15:8];
  wire        c2 = ~|F2[7:4];
  wire [3:0]  F3 = c2 ? F2[3:0] : F2[7:4];
  wire        c1 = ~|F3[3:2];
  wire [1:0]  F4 = c1 ? F3[1:0] : F3[3:2];
  wire        c0 = ~F4[1];
  assign c = {c3, c2, c1, c0};
endmodule

module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  wire        s_a = a[15];
  wire        s_b = b[15];
  wire [4:0]  e_a = a[14:10];
  wire [4:0]  e_b = b[14:10];
  wire [9:0]  f_a = a[9:0];
  wire [9:0]  f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == 5'd31);
  wire e_b_ones = (e_b == 5'd31);
  wire zero_a = e_a_zero & (f_a == 10'd0);
  wire zero_b = e_b_zero & (f_b == 10'd0);
  wire nan_a  = e_a_ones & (f_a != 10'd0);
  wire nan_b  = e_b_ones & (f_b != 10'd0);
  wire inf_a  = e_a_ones & (f_a == 10'd0);
  wire inf_b  = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire s_y     = s_a ^ s_b;

  // raw significands and their leading-zero counts (= pre-normalization shifts)
  wire [10:0] r_a = {~e_a_zero, f_a};
  wire [10:0] r_b = {~e_b_zero, f_b};
  wire [3:0]  lz_a, lz_b;
  lzc11 u_lza (.r(r_a), .c(lz_a));
  lzc11 u_lzb (.r(r_b), .c(lz_b));
  wire [5:0] k_a = zero_a ? 6'd21 : {2'd0, lz_a};
  wire [5:0] k_b = zero_b ? 6'd21 : {2'd0, lz_b};
  wire [5:0] k   = k_a + k_b;

  // effective exponent sum: subnormal operands contribute exponent 1
  wire [5:0] J = {1'b0, e_a_zero ? 5'd1 : e_a} + {1'b0, e_b_zero ? 5'd1 : e_b};

  wire [21:0] R = r_a * r_b;

  // norm = (R << k)[21]  (product's leading one at bit 21 after normalization)
  wire [33:0] Rsh  = {R, 12'b0};
  wire        norm = Rsh[33 - k[4:0]];

  // W0 = J - k + norm = pre-rounding biased exponent + 15
  wire signed [6:0] W0 = $signed({1'b0, J}) - $signed({1'b0, k}) + {6'd0, norm};

  // --- unified window extraction ------------------------------------------
  // u = window offset into R: 10-k+norm (normal result) or 26-J (subnormal result)
  // u2 = u + 9;  w12 = R[u-1 +: 12]  ->  G = w12[0], window = w12[11:1]
  wire sub_res = W0[6] | ~|W0[6:4];                    // W0 <= 15
  wire [5:0] u2 = sub_res ? (6'd35 - J) : (6'd19 - k + {5'd0, norm});

  wire [63:0] T3 = {32'b0, R, 10'b0};                  // T3[m] = R[m-10]
  wire [63:0] B;                                       // B[m] = |R[0 .. m-11]
  genvar m;
  generate
    for (m = 0; m < 64; m = m + 1) begin : gB
      if (m <= 10)       assign B[m] = 1'b0;
      else if (m >= 33)  assign B[m] = |R;
      else               assign B[m] = |R[m-11:0];
    end
  endgenerate

  wire [11:0] w12    = T3[u2 +: 12];
  wire [10:0] window = w12[11:1];
  wire        G      = w12[0];
  wire        S      = B[u2];

  // --- single rounding step (RNE) -----------------------------------------
  wire        round  = G & (S | window[0]);
  wire [11:0] sum12  = {1'b0, window} + {11'd0, round};
  wire        cy     = sum12[11];          // normal-path carry
  wire        fs_co  = sum12[10];          // subnormal carry (= exp LSB)
  wire [9:0]  frac_n = cy ? sum12[10:1] : sum12[9:0];
  wire [9:0]  frac_s = sum12[9:0];

  // --- pack ---------------------------------------------------------------
  wire ovf     = (W0 > 7'sd45) | ((W0 == 7'sd45) & cy);
  wire [4:0] ey5 = W0[4:0] + (cy ? 5'd18 : 5'd17);     // W0 + cy - 15

  wire spec_e = nan_in | inf_in | ovf;    // exponent field 31
  wire spec_f = inf_in | ovf;             // fraction field 0 (non-NaN specials)

  wire [4:0] exp_y  = spec_e ? 5'd31 :
                      sub_res ? {4'd0, fs_co} : ey5;
  wire [9:0] frac_y = nan_in  ? 10'h200 :
                      spec_f  ? 10'd0 :
                      sub_res ? frac_s : frac_n;
  wire sign_y = nan_in ? 1'b0 : s_y;

  assign y = {sign_y, exp_y, frac_y};
endmodule

// Optimized TPU MXU tile (logikbench tpu, cycle-accurate) - carry-save + Booth.
//
// Same structure/timing as the golden logikbench RTL (PE mesh, per-row input
// skew, eastward activation shift, southward psum accumulation, per-column
// output de-skew), with:
//   * 19-bit arithmetic: signed 8x8 products lie in [-16256, +16384] and a
//     column sums at most 8 of them, so every partial sum fits in 19 bits.
//   * Carry-save accumulation: each PE stores its partial sum as a (sum,
//     carry) pair; the full 19-bit CPA happens once, at the column output.
//   * Booth radix-4 recoding of the stationary weight: 4 partial products
//     (+1 correction) instead of 9, giving a shallow 4-level CSA tree.
//   * Reset only on the c_valid pipeline; datapath registers are naturally
//     flushed by the trace.

module csa3
  #(parameter W = 19)
   (input [W-1:0] x, y, z,
    output [W-1:0] s,
    output [W-1:0] c);

   assign s = x ^ y ^ z;
   assign c = ((x & y) | (x & z) | (y & z)) << 1;

endmodule

module tpu_pe
  #(parameter DW = 8,
    parameter W = 19)
   (
    input		  clk,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [W-1:0]	  cs_in,
    input [W-1:0]	  cc_in,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [W-1:0]	  cs_out,
    output reg [W-1:0]	  cc_out
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // activation sign-extended to accumulator width
   wire [W-1:0] ax = {{(W-DW){a_west[DW-1]}}, a_west};

   // Booth radix-4 recoding of w into 4 digits (k = 0, 2, 4, 6)
   // digit d_k = -2*b_{2k+1} + b_{2k} + b_{2k-1}
   wire b1_0 = w[1], b0_0 = w[0], bm_0 = 1'b0;
   wire b1_1 = w[3], b0_1 = w[2], bm_1 = w[1];
   wire b1_2 = w[5], b0_2 = w[4], bm_2 = w[3];
   wire b1_3 = w[7], b0_3 = w[6], bm_3 = w[5];

   wire neg1_0 = b1_0 & (b0_0 ^ bm_0);
   wire pos1_0 = ~b1_0 & (b0_0 ^ bm_0);
   wire neg2_0 = b1_0 & ~b0_0 & ~bm_0;
   wire pos2_0 = ~b1_0 & b0_0 & bm_0;
   wire zero_0 = ~(pos1_0 | neg1_0 | pos2_0 | neg2_0);

   wire neg1_1 = b1_1 & (b0_1 ^ bm_1);
   wire pos1_1 = ~b1_1 & (b0_1 ^ bm_1);
   wire neg2_1 = b1_1 & ~b0_1 & ~bm_1;
   wire pos2_1 = ~b1_1 & b0_1 & bm_1;
   wire zero_1 = ~(pos1_1 | neg1_1 | pos2_1 | neg2_1);

   wire neg1_2 = b1_2 & (b0_2 ^ bm_2);
   wire pos1_2 = ~b1_2 & (b0_2 ^ bm_2);
   wire neg2_2 = b1_2 & ~b0_2 & ~bm_2;
   wire pos2_2 = ~b1_2 & b0_2 & bm_2;
   wire zero_2 = ~(pos1_2 | neg1_2 | pos2_2 | neg2_2);

   wire neg1_3 = b1_3 & (b0_3 ^ bm_3);
   wire pos1_3 = ~b1_3 & (b0_3 ^ bm_3);
   wire neg2_3 = b1_3 & ~b0_3 & ~bm_3;
   wire pos2_3 = ~b1_3 & b0_3 & bm_3;
   wire zero_3 = ~(pos1_3 | neg1_3 | pos2_3 | neg2_3);

   wire dneg_0 = neg1_0 | neg2_0;
   wire dneg_1 = neg1_1 | neg2_1;
   wire dneg_2 = neg1_2 | neg2_2;
   wire dneg_3 = neg1_3 | neg2_3;

   wire dsel_0 = neg2_0 | pos2_0;   // |d| = 2
   wire dsel_1 = neg2_1 | pos2_1;
   wire dsel_2 = neg2_2 | pos2_2;
   wire dsel_3 = neg2_3 | pos2_3;

   wire [W-1:0] pp0 = zero_0 ? {W{1'b0}} :
                      (dneg_0 ? ~(dsel_0 ? (ax << 1) : ax)
                              : (dsel_0 ? (ax << 1) : ax));
   wire [W-1:0] pp1 = zero_1 ? {W{1'b0}} :
                      (dneg_1 ? ~(dsel_1 ? (ax << 3) : (ax << 2))
                              : (dsel_1 ? (ax << 3) : (ax << 2)));
   wire [W-1:0] pp2 = zero_2 ? {W{1'b0}} :
                      (dneg_2 ? ~(dsel_2 ? (ax << 5) : (ax << 4))
                              : (dsel_2 ? (ax << 5) : (ax << 4)));
   wire [W-1:0] pp3 = zero_3 ? {W{1'b0}} :
                      (dneg_3 ? ~(dsel_3 ? (ax << 7) : (ax << 6))
                              : (dsel_3 ? (ax << 7) : (ax << 6)));

   // +1 correction per negative digit (Ksum <= 4)
   wire [2:0] ksum = neg1_0 + neg2_0 + neg1_1 + neg2_1 +
                     neg1_2 + neg2_2 + neg1_3 + neg2_3;

   // 7 -> 2 compression tree (mod 2^W)
   wire [W-1:0] sA, cA, sB, cB, sC, cC, sD, cD;

   csa3 #(W) u1 (pp0, pp1, pp2, sA, cA);
   csa3 #(W) u2 (pp3, {{(W-3){1'b0}}, ksum}, cs_in, sB, cB);
   csa3 #(W) u3 (sA, sB, cc_in, sC, cC);
   csa3 #(W) u4 (sC, cA, cB, sD, cD);
   csa3 #(W) u5 (sD, cD, cC, cs_n, cc_n);

   wire [W-1:0] cs_n, cc_n;

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east <= a_west;
         cs_out <= cs_n;
         cc_out <= cc_n;
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,   // external result width
    parameter PACCW = 19)  // internal accumulator width
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i delayed i cycles; zero injected when a_valid low.
   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
		       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   // PE mesh
   wire [DW-1:0]    a_h  [0:N-1][0:N];
   wire [PACCW-1:0] cs_v [0:N][0:N-1];
   wire [PACCW-1:0] cc_v [0:N][0:N-1];
   wire [DW-1:0]    w_v  [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign cs_v[0][j] = {PACCW{1'b0}};
         assign cc_v[0][j] = {PACCW{1'b0}};
         assign w_v[0][j]  = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(PACCW)) u_pe
              (.clk        (clk),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .cs_in      (cs_v[i][j]),
               .cc_in      (cc_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .cs_out     (cs_v[i+1][j]),
               .cc_out     (cc_v[i+1][j]));
         end
      end
   endgenerate

   // Output: single CPA per column, then de-skew column j by N-1-j cycles
   // and sign-extend to 32 bits.
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = cs_v[N][j] + cc_v[N][j];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PACCW){cin[PACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PACCW){dr[D-1][PACCW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

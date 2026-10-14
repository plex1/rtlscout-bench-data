//#############################################################################
// Optimized weight-stationary systolic matrix-multiply tile (TPU MXU).
// Cycle-accurate re-implementation of the logikbench 'tpu' benchmark.
//
// Systolic data flow is IDENTICAL to the original logikbench tpu:
//   * weights shift south through per-column chains while w_valid is high
//     (bottom matrix row first); after 8 load cycles PE(i,j) holds B[i][j],
//   * activations hop west->east through one register per column,
//   * row i of A is delayed i cycles at the input (skew), column j of the
//     result is delayed N-1-j cycles at the output (de-skew),
//   * the whole array advances only while w_valid is low (stall semantics),
//   * result row m appears exactly 15 compute cycles after its A row was
//     accepted (c_valid = a_valid delayed 15 compute cycles).
//
// Micro-architecture optimizations:
//   1. Gated clocks (clk & enable) instead of enable feedback muxes.
//   2. No per-PE carry-propagate multiplier.  Each PE computes its product
//      as a Baugh-Wooley partial-product array (one AND/NAND level off the
//      activation register) and merges it into the redundant accumulator
//      with a 4-level carry-save reduction:
//          a*w = ARRAY - 2^15 + 2^8   (exact), so mod 2^19
//          a*w = ARRAY + (2^19 - 2^15) + 2^8.
//   3. Three-row redundant accumulator (sum, carry, carry2): the per-PE
//      tree reduces 11 rows (3 accumulator + 8 partial-product rows) to 3
//      rows in 4 carry-save levels -- one level less than a 2-row
//      accumulator, and the third row bypasses the final level.
//   4. Exact-width modular accumulation: every column result -- a sum of 8
//      signed 8x8 products -- lies in [-130048, 131072], inside the 19-bit
//      two's-complement range, so 19-bit registers reconstruct the exact
//      result:  V = S + C + D (mod 2^19), sign-extended once to 32 bits.
//   5. The output carry-propagate adder (3 rows -> 1) is pipelined across
//      the output de-skew stages: stage 1 registers the low 9 bits plus
//      both speculative high halves (carry-in 0/1), stage 2 selects, and
//      the remaining stages shift.  Column 7 (no de-skew delay) uses a
//      carry-save adder followed by a Kogge-Stone prefix adder straight to
//      the pins.
//   6. The eastward activation register of the last PE column is unused
//      and is removed.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,	 // async, active high
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;	// 15 compute cycles injection -> output
   localparam PW  = 19;		// redundant accumulator word width

   wire ld   = w_valid;
   wire en   = ~w_valid;
   wire gclk = clk & en;	// compute gated clock
   wire lclk = clk & ld;	// weight-load gated clock

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles; zeros
   // are injected while a_valid is low.
   //##########################################################################
   wire [N*DW-1:0] a_left;
   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge gclk or posedge rst) begin
               if (rst)
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               else
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= (k == 0) ? din : sr[k-1];
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // The systolic array with 3-row carry-save multiply-accumulate.
   //##########################################################################
   wire [DW-1:0] a_h [0:N-1][0:N];	// activation mesh (west -> east)
   wire [PW-1:0] s_v [0:N][0:N-1];	// accumulator row 1 (north -> south)
   wire [PW-1:0] c_v [0:N][0:N-1];	// accumulator row 2
   wire [PW-1:0] d_v [0:N][0:N-1];	// accumulator row 3
   wire [DW-1:0] w_v [0:N][0:N-1];	// weight shift mesh (top -> south)

   generate
      // west edge / top edge
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PW{1'b0}};
         assign c_v[0][j] = {PW{1'b0}};
         assign d_v[0][j] = {PW{1'b0}};
         assign w_v[0][j] = w_data[j*DW +: DW];
      end
      // PE mesh
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            // activation hop register (unused in the last column)
            if (j < N-1) begin : g_hop
               reg [DW-1:0] ae;
               always @(posedge gclk or posedge rst)
                 if (rst) ae <= {DW{1'b0}};
                 else     ae <= a_h[i][j];
               assign a_h[i][j+1] = ae;
            end

            tpu_pe #(.DW(DW), .PW(PW)) u_pe
              (.rst     (rst),
               .gclk    (gclk),
               .lclk    (lclk),
               .w_north (w_v[i][j]),
               .a_west  (a_h[i][j]),
               .s_in    (s_v[i][j]),
               .c_in    (c_v[i][j]),
               .d_in    (d_v[i][j]),
               .w_south (w_v[i+1][j]),
               .s_out   (s_v[i+1][j]),
               .c_out   (c_v[i+1][j]),
               .d_out   (d_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Output stage: carry-propagate conversion of the redundant triple
   // (S, C, D) -> V = S + C + D (mod 2^19), pipelined across the de-skew
   // chain so a result row leaves whole, 15 compute cycles after its A row.
   //
   //   column j: de-skew depth DSK = N-1-j.
   //   DSK >= 2 : stage 1 registers the speculative partial sums, stage 2
   //              selects, remaining DSK-2 stages shift.
   //   DSK == 1 : stage 1 registers the partial sums; the select feeds the
   //              output pins in the result cycle.
   //   DSK == 0 : carry-save adder + Kogge-Stone prefix adder to the pins.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
         localparam DSK = N - 1 - j;
         wire [PW-1:0] s_last = s_v[N][j];
         wire [PW-1:0] c_last = c_v[N][j];
         wire [PW-1:0] d_last = d_v[N][j];

         if (DSK == 0) begin : g_d0
            wire [PW-1:0] u  = s_last ^ c_last ^ d_last;
            wire [PW-1:0] mj = (s_last & c_last) | (s_last & d_last) |
                               (c_last & d_last);
            wire [PW-1:0] cc = {mj[PW-2:0], 1'b0};
            wire [PW-1:0] ksv;
            tpu_ks #(.W(PW)) u_ks (.x(u), .y(cc), .s(ksv));
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){ksv[PW-1]}}, ksv};
         end
         else begin : g_dn
            // low 9 bits: 3-input add (carry-out is two bits:
            // 3*511 = 1533 < 2^11)
            wire [10:0]	   lo_f  = {2'b00, s_last[8:0]} + {2'b00, c_last[8:0]} +
                                   {2'b00, d_last[8:0]};
            // high 10 bits, all three carry hypotheses (0, 1, 2)
            wire [9:0]	   hi0_f = {1'b0, s_last[PW-1:9]} +
                                   {1'b0, c_last[PW-1:9]} +
                                   {1'b0, d_last[PW-1:9]};
            wire [9:0]	   hi1_f = hi0_f + 1'b1;
            wire [9:0]	   hi2_f = hi0_f + 2'd2;

            reg [8:0]	 rlo;
            reg [1:0]	 rk;
            reg [9:0]	 rh0, rh1, rh2;
            reg [PW-1:0] dr [0:N-2];
            integer	 k;
            always @(posedge gclk or posedge rst) begin
               if (rst) begin
                  rlo <= {9{1'b0}};
                  rk  <= {2{1'b0}};
                  rh0 <= {10{1'b0}};
                  rh1 <= {10{1'b0}};
                  rh2 <= {10{1'b0}};
                  for (k = 0; k <= N-2; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else begin
                  rlo <= lo_f[8:0];
                  rk  <= lo_f[10:9];
                  rh0 <= hi0_f;
                  rh1 <= hi1_f;
                  rh2 <= hi2_f;
                  dr[0] <= {rk == 2'd2 ? rh2[9:0] :
			    rk == 2'd1 ? rh1[9:0] : rh0[9:0], rlo};
                  for (k = 1; k <= N-2; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            // DSK == 1: no room for the value register, select at the pins
            wire [PW-1:0] vout = (DSK == 1) ?
              {rk == 2'd2 ? rh2[9:0] : rk == 2'd1 ? rh1[9:0] : rh0[9:0],
	       rlo} : dr[DSK-2];
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-PW){vout[PW-1]}}, vout};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed 15 compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge gclk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Processing element: redundant multiply-accumulate cell.
//
//   (S,C,D) += a*w   (mod 2^PW), all non-clock logic carry-save.
//
// Weight load on lclk (column shift register); accumulate on gclk.
//#############################################################################
module tpu_pe
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		  rst,
    input		  gclk,
    input		  lclk,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PW-1:0]	  s_in,
    input [PW-1:0]	  c_in,
    input [PW-1:0]	  d_in,
    output [DW-1:0]	  w_south,
    output reg [PW-1:0]	  s_out,
    output reg [PW-1:0]	  c_out,
    output reg [PW-1:0]	  d_out
    );

   // stationary weight: shifts south on load edges
   reg [DW-1:0] w;
   always @(posedge lclk or posedge rst)
     if (rst) w <= {DW{1'b0}};
     else     w <= w_north;
   assign w_south = w;

   //########################################################################
   // Baugh-Wooley partial products of the signed 8x8 product, as PW-bit
   // rows.  With a_i/w_j the two's-complement bits (bit 7 = sign):
   //
   //   a*w = ARRAY - 2^15 + 2^8                       (exact identity)
   //   ARRAY = sum_{i<=6,j<=6} a_i&w_j << (i+j)
   //         + sum_{i<=6} !(a_i&w_7) << (i+7)
   //         + sum_{j<=6} !(a_7&w_j) << (j+7)
   //         + (a_7&w_7) << 14
   //
   // so mod 2^PW:  a*w = ARRAY + (2^PW - 2^15) + 2^8.
   // Row 0 carries the constant bits (8, 15, 16, 17, 18).
   //########################################################################
   reg [PW-1:0] V [0:7];
   integer	 ii, jj;
   always @* begin
      for (ii = 0; ii < 8; ii = ii + 1) begin
         V[ii] = {PW{1'b0}};
         for (jj = 0; jj < 7; jj = jj + 1)
           V[ii][ii+jj] = (ii < 7) ? (a_west[ii] & w[jj])
				   : ~(a_west[7] & w[jj]);
         V[ii][ii+7]  = (ii < 7) ? ~(a_west[ii] & w[7])
				 : (a_west[7] & w[7]);
      end
      V[0][8]  = 1'b1;   // +2^8
      V[0][15] = 1'b1;   // +2^PW - 2^15
      V[0][16] = 1'b1;
      V[0][17] = 1'b1;
      V[0][18] = 1'b1;
   end

   //########################################################################
   // Carry-save reduction of the 11 rows {s_in, c_in, d_in, V0..V7} down to
   // 3.  Every step is exact modulo 2^PW (the top carry is dropped), which
   // is valid because the true accumulated value fits in 19 bits.
   //########################################################################
   // pass 1: 3 CSAs -> 8 rows
   wire [PW-1:0] m1a = (V[0] & V[1]) | (V[0] & s_in) | (V[1] & s_in);
   wire [PW-1:0] s1a = V[0] ^ V[1] ^ s_in;
   wire [PW-1:0] k1a = {m1a[PW-2:0], 1'b0};
   wire [PW-1:0] m1b = (V[2] & V[3]) | (V[2] & c_in) | (V[3] & c_in);
   wire [PW-1:0] s1b = V[2] ^ V[3] ^ c_in;
   wire [PW-1:0] k1b = {m1b[PW-2:0], 1'b0};
   wire [PW-1:0] m1c = (V[4] & V[5]) | (V[4] & d_in) | (V[5] & d_in);
   wire [PW-1:0] s1c = V[4] ^ V[5] ^ d_in;
   wire [PW-1:0] k1c = {m1c[PW-2:0], 1'b0};
   // rows: s1a k1a s1b k1b s1c k1c V[6] V[7]

   // pass 2: 2 CSAs -> 6 rows
   wire [PW-1:0] m2a = (s1a & k1a) | (s1a & s1b) | (k1a & s1b);
   wire [PW-1:0] s2a = s1a ^ k1a ^ s1b;
   wire [PW-1:0] k2a = {m2a[PW-2:0], 1'b0};
   wire [PW-1:0] m2b = (k1b & s1c) | (k1b & k1c) | (s1c & k1c);
   wire [PW-1:0] s2b = k1b ^ s1c ^ k1c;
   wire [PW-1:0] k2b = {m2b[PW-2:0], 1'b0};
   // rows: s2a k2a s2b k2b V[6] V[7]

   // pass 3: 2 CSAs -> 4 rows
   wire [PW-1:0] m3a = (s2a & k2a) | (s2a & s2b) | (k2a & s2b);
   wire [PW-1:0] s3a = s2a ^ k2a ^ s2b;
   wire [PW-1:0] k3a = {m3a[PW-2:0], 1'b0};
   wire [PW-1:0] m3b = (k2b & V[6]) | (k2b & V[7]) | (V[6] & V[7]);
   wire [PW-1:0] s3b = k2b ^ V[6] ^ V[7];
   wire [PW-1:0] k3b = {m3b[PW-2:0], 1'b0};
   // rows: s3a k3a s3b k3b

   // pass 4: 1 CSA -> 3 rows
   wire [PW-1:0] m4 = (s3a & k3a) | (s3a & s3b) | (k3a & s3b);
   wire [PW-1:0] s4 = s3a ^ k3a ^ s3b;
   wire [PW-1:0] k4 = {m4[PW-2:0], 1'b0};
   // rows: s4 k4 k3b

   always @(posedge gclk or posedge rst) begin
      if (rst) begin
         s_out <= {PW{1'b0}};
         c_out <= {PW{1'b0}};
         d_out <= {PW{1'b0}};
      end
      else begin
         s_out <= s4;
         c_out <= k4;
         d_out <= k3b;
      end
   end

endmodule

//#############################################################################
// W-bit Kogge-Stone parallel-prefix adder (x + y mod 2^W).
//#############################################################################
module tpu_ks
  #(parameter W = 19)
   (input [W-1:0] x,
    input [W-1:0] y,
    output [W-1:0] s);

   wire [W-1:0] p0 = x ^ y;
   wire [W-1:0] g0 = x & y;

   // prefix levels: distance 1, 2, 4, 8, 16 (covers W <= 32)
   wire [W-1:0] p1, g1, p2, g2, p3, g3, p4, g4, p5, g5;
   genvar i;
   generate
      for (i = 0; i < W; i = i + 1) begin : lv1
         if (i == 0) begin : b0
            assign p1[0] = p0[0];
            assign g1[0] = g0[0];
         end
         else begin : bn
            assign p1[i] = p0[i] & p0[i-1];
            assign g1[i] = g0[i] | (p0[i] & g0[i-1]);
         end
      end
      for (i = 0; i < W; i = i + 1) begin : lv2
         if (i < 2) begin : pass
            assign p2[i] = p1[i];
            assign g2[i] = g1[i];
         end
         else begin : black
            assign p2[i] = p1[i] & p1[i-2];
            assign g2[i] = g1[i] | (p1[i] & g1[i-2]);
         end
      end
      for (i = 0; i < W; i = i + 1) begin : lv3
         if (i < 4) begin : pass
            assign p3[i] = p2[i];
            assign g3[i] = g2[i];
         end
         else begin : black
            assign p3[i] = p2[i] & p2[i-4];
            assign g3[i] = g2[i] | (p2[i] & g2[i-4]);
         end
      end
      for (i = 0; i < W; i = i + 1) begin : lv4
         if (i < 8) begin : pass
            assign p4[i] = p3[i];
            assign g4[i] = g3[i];
         end
         else begin : black
            assign p4[i] = p3[i] & p3[i-8];
            assign g4[i] = g3[i] | (p3[i] & g3[i-8]);
         end
      end
      for (i = 0; i < W; i = i + 1) begin : lv5
         if (i < 16) begin : pass
            assign p5[i] = p4[i];
            assign g5[i] = g4[i];
         end
         else begin : black
            assign p5[i] = p4[i] & p4[i-16];
            assign g5[i] = g4[i] | (p4[i] & g4[i-16]);
         end
      end
   endgenerate

   // carry into bit i is g5[i-1]; g5[i] spans [max(0,i-31) .. i] = [0..i]
   assign s[0] = p0[0];
   for (i = 1; i < W; i = i + 1) begin : sumgen
      assign s[i] = p0[i] ^ g5[i-1];
   end

endmodule

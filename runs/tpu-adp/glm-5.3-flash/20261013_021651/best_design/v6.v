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
//      as a Baugh-Wooley partial-product array and merges it directly into
//      the redundant (sum, carry) accumulator with a 5-level carry-save
//      reduction:
//          a*w = ARRAY - 2^15 + 2^8   (exact, ARRAY = 64 partial products
//               with complemented sign columns), so mod 2^19
//          a*w = ARRAY + (2^19 - 2^15) + 2^8
//      The per-PE datapath is thus AND gates + 5 carry-save adder levels,
//      with a single 19-bit carry-propagate adder per column at the output.
//   3. Exact-width modular accumulation: every column result -- a sum of 8
//      signed 8x8 products -- lies in [-130048, 131072], inside the 19-bit
//      two's-complement range, so 19-bit registers reconstruct the exact
//      result:  V = S + C (mod 2^19), sign-extended once to 32 bits.
//   4. The eastward activation register of the last PE column is unused
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
   // The systolic array with merged carry-save multiply-accumulate.
   //##########################################################################
   wire [DW-1:0] a_h [0:N-1][0:N];	// activation mesh (west -> east)
   wire [PW-1:0] s_v [0:N][0:N-1];	// accumulator sum mesh (north -> south)
   wire [PW-1:0] c_v [0:N][0:N-1];	// accumulator carry mesh
   wire [DW-1:0] w_v [0:N][0:N-1];	// weight shift mesh (top -> south)

   generate
      // west edge / top edge
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PW{1'b0}};
         assign c_v[0][j] = {PW{1'b0}};
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
               .w_south (w_v[i+1][j]),
               .s_out   (s_v[i+1][j]),
               .c_out   (c_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Column result conversion (S + C -> two's complement) + output de-skew:
   // column j is delayed N-1-j compute cycles so a result row leaves whole.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
         localparam D = N - 1 - j;
         wire [PW-1:0] val = s_v[N][j] + c_v[N][j];
         if (D == 0) begin : g_d0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){val[PW-1]}}, val};
         end
         else begin : g_dn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge gclk or posedge rst) begin
               if (rst)
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               else
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= (k == 0) ? val : dr[k-1];
            end
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
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
//   (S,C) += a*w   (mod 2^PW), all non-clock logic carry-save.
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
    output [DW-1:0]	  w_south,
    output reg [PW-1:0]	  s_out,
    output reg [PW-1:0]	  c_out
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
   // Carry-save reduction of the 10 rows {s_in, c_in, V0..V7} down to 2.
   // Every step is exact modulo 2^PW (the top carry is dropped), which is
   // valid because the true accumulated value fits in 19 bits.
   //########################################################################
   // pass 1: 3 CSAs -> 7 rows
   wire [PW-1:0] m1a = (V[0] & V[1]) | (V[0] & V[2]) | (V[1] & V[2]);
   wire [PW-1:0] s1a = V[0] ^ V[1] ^ V[2];
   wire [PW-1:0] k1a = {m1a[PW-2:0], 1'b0};
   wire [PW-1:0] m1b = (V[3] & V[4]) | (V[3] & V[5]) | (V[4] & V[5]);
   wire [PW-1:0] s1b = V[3] ^ V[4] ^ V[5];
   wire [PW-1:0] k1b = {m1b[PW-2:0], 1'b0};
   wire [PW-1:0] m1c = (V[6] & V[7]) | (V[6] & s_in) | (V[7] & s_in);
   wire [PW-1:0] s1c = V[6] ^ V[7] ^ s_in;
   wire [PW-1:0] k1c = {m1c[PW-2:0], 1'b0};
   // rows: s1a k1a s1b k1b s1c k1c c_in

   // pass 2: 2 CSAs -> 5 rows
   wire [PW-1:0] m2a = (s1a & k1a) | (s1a & s1b) | (k1a & s1b);
   wire [PW-1:0] s2a = s1a ^ k1a ^ s1b;
   wire [PW-1:0] k2a = {m2a[PW-2:0], 1'b0};
   wire [PW-1:0] m2b = (k1b & s1c) | (k1b & k1c) | (s1c & k1c);
   wire [PW-1:0] s2b = k1b ^ s1c ^ k1c;
   wire [PW-1:0] k2b = {m2b[PW-2:0], 1'b0};
   // rows: s2a k2a s2b k2b c_in

   // pass 3: 1 CSA -> 4 rows
   wire [PW-1:0] m3 = (s2a & k2a) | (s2a & s2b) | (k2a & s2b);
   wire [PW-1:0] s3 = s2a ^ k2a ^ s2b;
   wire [PW-1:0] k3 = {m3[PW-2:0], 1'b0};
   // rows: s3 k3 k2b c_in

   // pass 4: 1 CSA -> 3 rows
   wire [PW-1:0] m4 = (s3 & k3) | (s3 & k2b) | (k3 & k2b);
   wire [PW-1:0] s4 = s3 ^ k3 ^ k2b;
   wire [PW-1:0] k4 = {m4[PW-2:0], 1'b0};
   // rows: s4 k4 c_in

   // pass 5: 1 CSA -> 2 rows
   wire [PW-1:0] m5 = (s4 & k4) | (s4 & c_in) | (k4 & c_in);
   wire [PW-1:0] s5 = s4 ^ k4 ^ c_in;
   wire [PW-1:0] k5 = {m5[PW-2:0], 1'b0};

   always @(posedge gclk or posedge rst) begin
      if (rst) begin
         s_out <= {PW{1'b0}};
         c_out <= {PW{1'b0}};
      end
      else begin
         s_out <= s5;
         c_out <= k5;
      end
   end

endmodule

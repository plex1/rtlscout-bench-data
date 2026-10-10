//#############################################################################
// Weight-stationary 8x8 systolic matrix-multiply tile (TPU-style MXU).
//
// Merged-compressor-tree implementation, cycle-accurate with the logikbench
// `tpu` reference.
//
// Per PE, the product a*w (8x8 signed) is NEVER formed as a single number.
// Instead its Baugh-Wooley partial-product rows are merged directly with the
// incoming carry-save partial-sum pair (sum,carry) of the column chain and
// reduced by a wall of 3:2 compressors straight down to the new (sum,carry)
// register pair.  This removes both the multiplier's internal final
// carry-propagate adder and the separate accumulate adder from the critical
// path.
//
// Timing (in "en" cycles, i.e. cycles where w_valid is low):
//   * PE(i,j) multiplies/absorbs term A[m][i]*B[i][j] at en-index m+i+j,
//     using the activation byte at tap depth (i+j) and its live weight
//     register -- the exact operand schedule of the reference (so weight
//     reloads affect in-flight terms identically).
//   * result row m appears on c_data at en-index m + 2N-1.
//
// Chain width 20 bits (mod 2^20): 8 products of magnitude <= 2^14 sum to at
// most +/-2^17, so bits [18:0] of the mod-2^20 sum are exact.
//
// Baugh-Wooley bookkeeping: with a = a_u - 128*a7, w = w_u - 128*w7,
//   a*w = a_u*w_u - 128*a_u*w7 - 128*w_u*a7 + 2^14*a7*w7.
// Using ~(a_i & w7) etc. for the negated sign terms costs a constant
// (+2^7 - 2^14) per negated row, i.e. the tree sum equals the true product
// + 32512.  Over the 8 PEs of a column the offset is 260096, subtracted once
// in the column's final carry-propagate adder.
//#############################################################################

module tpu #(parameter N = 8, DW = 8, ACCW = 32) (
    input  wire              clk,
    input  wire              rst,
    input  wire              w_valid,
    input  wire [N*DW-1:0]   w_data,
    input  wire              a_valid,
    input  wire [N*DW-1:0]   a_data,
    output wire              c_valid,
    output wire [N*ACCW-1:0] c_data);

   localparam LAT = 2*N - 1;    // 15: accept -> result, in en-cycles
   localparam CW  = 19;         // carry-save chain width (mod 2^19)

   wire en = ~w_valid;
   wire ld =  w_valid;

   //#########################################################################
   // Weight shift chains (south), loaded from the top while ld is high.
   //#########################################################################
   reg [DW-1:0] wq [0:N*N-1];
   integer      a;
   always @(posedge clk) begin
      if (rst) begin
         for (a = 0; a < N*N; a = a + 1) wq[a] <= {DW{1'b0}};
      end else if (ld) begin
         for (a = 0; a < N; a = a + 1) wq[a] <= w_data[a*DW +: DW];
         for (a = N; a < N*N; a = a + 1) wq[a] <= wq[a-N];
      end
   end

   //#########################################################################
   // Activation tap lines: row i delays its byte by (i+j) en-cycles before
   // PE(i,j) multiplies it (the reference's a_west schedule).  Zero is
   // injected while a_valid is low.  Taps advance only on en.
   //#########################################################################
   wire [DW-1:0] taps [0:N*N-1];
   genvar        gi, gj;
   integer       gp;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_trow
         localparam DM = gi + N - 1;
         wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
         reg  [DW-1:0] sr  [1:DM];
         always @(posedge clk) begin
            if (rst) begin
               for (gp = 1; gp <= DM; gp = gp + 1) sr[gp] <= {DW{1'b0}};
            end else if (en) begin
               sr[1] <= din;
               for (gp = 2; gp <= DM; gp = gp + 1) sr[gp] <= sr[gp-1];
            end
         end
         for (gj = 0; gj < N; gj = gj + 1) begin : g_tap
            assign taps[gi*N + gj] = ((gi + gj) == 0) ? din : sr[gi+gj];
         end
      end
   endgenerate

   //#########################################################################
   // Per-PE merged compressor tree.
   //
   // BW rows (20-bit frames, row i = a-bit i):
   //   i in 0..6:  bits a_i&w_j  at 2^(i+j),  j = 0..6
   //               bit ~(a_i&w_7) at 2^(i+7)
   //   i = 7    :  bits ~(a_7&w_j) at 2^(j+7), j = 0..6
   //               bit  (a_7&w_7) at 2^14
   // Together with the incoming (sum,carry) pair: 10 rows -> 5 levels of 3:2
   // -> new (sum,carry).
   //#########################################################################
   reg [CW-1:0] sq [0:N*N-1];
   reg [CW-1:0] cq [0:N*N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_row
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pe
            localparam ID = gi*N + gj;
            localparam UP = ID - N;

            wire [DW-1:0] av = taps[ID];
            wire [DW-1:0] wv = wq[ID];

            // ---- Baugh-Wooley partial-product rows (20-bit frames) ----
            wire [CW-1:0] r [0:7];
            for (genvar ii = 0; ii < 8; ii = ii + 1) begin : g_bwr
               reg [CW-1:0] rr;
               integer jj;
               always @(*) begin
                  rr = {CW{1'b0}};
                  for (jj = 0; jj < 8; jj = jj + 1) begin
                     if (ii < 7) begin
                        if (jj < 7) rr[ii+jj] = av[ii] & wv[jj];
                        else        rr[ii+7]  = ~(av[ii] & wv[7]);
                     end else begin
                        if (jj < 7) rr[jj+7]  = ~(av[7] & wv[jj]);
                        else        rr[14]    = av[7] & wv[7];
                     end
                  end
               end
               assign r[ii] = rr;
            end

            // ---- incoming partial-sum pair ----
            wire [CW-1:0] sin = (gi == 0) ? {CW{1'b0}} : sq[UP];
            wire [CW-1:0] cin = (gi == 0) ? {CW{1'b0}} : cq[UP];

            // ---- 5-level 3:2 reduction of 10 rows -> 2 rows ----
            // L1: three parallel CSAs on rows 0..7 and the incoming sum
            wire [CW-1:0] a1s = r[0] ^ r[1] ^ r[2];
            wire [CW-1:0] a1c = ((r[0] & r[1]) | (r[0] & r[2]) | (r[1] & r[2])) << 1;
            wire [CW-1:0] a2s = r[3] ^ r[4] ^ r[5];
            wire [CW-1:0] a2c = ((r[3] & r[4]) | (r[3] & r[5]) | (r[4] & r[5])) << 1;
            wire [CW-1:0] a3s = r[6] ^ r[7] ^ sin;
            wire [CW-1:0] a3c = ((r[6] & r[7]) | (r[6] & sin) | (r[7] & sin)) << 1;
            // (cin deferred to the final stage)

            // L2: two CSAs
            wire [CW-1:0] b1s = a1s ^ a1c ^ a2s;
            wire [CW-1:0] b1c = ((a1s & a1c) | (a1s & a2s) | (a1c & a2s)) << 1;
            wire [CW-1:0] b2s = a2c ^ a3s ^ a3c;
            wire [CW-1:0] b2c = ((a2c & a3s) | (a2c & a3c) | (a3s & a3c)) << 1;
            // remaining: cin, b1s, b1c, b2s, b2c

            // L3
            wire [CW-1:0] d1s = b1s ^ b1c ^ b2s;
            wire [CW-1:0] d1c = ((b1s & b1c) | (b1s & b2s) | (b1c & b2s)) << 1;
            // remaining: cin, b2c, d1s, d1c

            // L4
            wire [CW-1:0] e1s = b2c ^ d1s ^ d1c;
            wire [CW-1:0] e1c = ((b2c & d1s) | (b2c & d1c) | (d1s & d1c)) << 1;
            // remaining: cin, e1s, e1c

            // L5 (final)
            wire [CW-1:0] ssum = cin ^ e1s ^ e1c;
            wire [CW-1:0] scry = ((cin & e1s) | (cin & e1c) | (e1s & e1c)) << 1;

            always @(posedge clk) begin
               if (rst) begin
                  sq[ID] <= {CW{1'b0}};
                  cq[ID] <= {CW{1'b0}};
               end else if (en) begin
                  sq[ID] <= ssum;
                  cq[ID] <= scry;
               end
            end
         end
      end
   endgenerate

   //#########################################################################
   // Per-column result: one carry-propagate add at the bottom of the chain
   // (including the -260096 BW offset for the column's 8 products), latched
   // by the first de-skew register.  Column j's word is then available at
   // en-index m+9+j, well before the m+15 result cycle for j <= 6; column 7
   // uses the combinational CPA directly at m+15.
   //#########################################################################
   genvar gj2;
   integer gq;
   generate
      for (gj2 = 0; gj2 < N; gj2 = gj2 + 1) begin : g_col
         localparam D  = N - 1 - gj2;
         localparam ID = (N-1)*N + gj2;
         wire [CW-1:0] val = sq[ID] + cq[ID] + 20'd788480; // -260096 mod 2^20
         wire [18:0]   res = val[18:0];
         if (D == 0) begin : g_d0
            assign c_data[gj2*ACCW +: ACCW] = {{(ACCW-19){res[18]}}, res};
         end else begin : g_dD
            reg [18:0] dr [1:D];
            always @(posedge clk) begin
               if (rst) begin
                  for (gq = 1; gq <= D; gq = gq + 1) dr[gq] <= 19'b0;
               end else if (en) begin
                  dr[1] <= res;
                  for (gq = 2; gq <= D; gq = gq + 1) dr[gq] <= dr[gq-1];
               end
            end
            assign c_data[gj2*ACCW +: ACCW] =
                 {{(ACCW-19){dr[D][18]}}, dr[D]};
         end
      end
   endgenerate

   //#########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT en-cycles.
   //#########################################################################
   reg [LAT-1:0] vpipe;
   integer       v;
   always @(posedge clk) begin
      if (rst) begin
         vpipe <= {LAT{1'b0}};
      end else if (en) begin
         vpipe[0] <= a_valid;
         for (v = 1; v < LAT; v = v + 1)
           vpipe[v] <= vpipe[v-1];
      end
   end
   assign c_valid = vpipe[LAT-1];

endmodule

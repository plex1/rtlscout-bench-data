//#############################################################################
// Weight-stationary 8x8 systolic matrix-multiply tile (TPU-style MXU).
//
// Carry-save retimed implementation, cycle-accurate with the logikbench `tpu`
// reference.
//
// Timing contract (in "en" cycles, i.e. cycles where w_valid is low):
//   * PE(i,j) multiplies term A[m][i]*B[i][j] at en-index m+i+j, using the
//     activation byte at tap depth (i+j) and its live weight register -- the
//     exact operand schedule of the reference (weight reloads therefore
//     affect in-flight terms identically, with no correction muxes).
//   * each product is absorbed into a per-column CARRY-SAVE accumulator
//     chain (sum/carry pair) through one 3:2 CSA stage; the only
//     carry-propagate adder per column sits at the bottom, timed to land
//     exactly inside the output de-skew budget.
//   * result row m appears on c_data at en-index m + 2N-1.
//
// All chain values fit in 18 signed bits (<= 8 products of magnitude 2^14),
// so a 20-bit mod-2^20 carry-save pair is lossless; results are
// sign-extended to the 32-bit port.
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
   localparam CW  = 20;         // carry-save chain width (mod 2^20)

   wire en = ~w_valid;   // array advances only while not loading weights
   wire ld =  w_valid;

   //#########################################################################
   // Weight shift chains (south), loaded from the top while ld is high.
   // Cycle k of a load carries B[N-1-k][*]; after N cycles PE(i,j) = B[i][j].
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
   // injected while a_valid is low.  Taps advance only on en, so they freeze
   // during weight loads exactly like the reference's a_east registers.
   //#########################################################################
   wire [DW-1:0] taps [0:N*N-1];
   genvar        gi, gj;
   integer       gp;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_trow
         localparam DM = gi + N - 1;      // deepest tap of this row
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
   // Per-PE carry-save accumulate: the (signed) product is absorbed into the
   // (sum,carry) chain with a single 3:2 CSA (mod 2^CW).  Row 0 starts a
   // fresh accumulation (zero north input, as in the reference).
   //#########################################################################
   reg [CW-1:0] sq [0:N*N-1];
   reg [CW-1:0] cq [0:N*N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_row
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pe
            localparam ID = gi*N + gj;
            localparam UP = ID - N;

            wire signed [2*DW-1:0] prod = $signed(taps[ID]) * $signed(wq[ID]);
            wire [CW-1:0] pz = {{(CW-2*DW){prod[2*DW-1]}}, prod};

            wire [CW-1:0] sin = (gi == 0) ? {CW{1'b0}} : sq[UP];
            wire [CW-1:0] cin = (gi == 0) ? {CW{1'b0}} : cq[UP];

            wire [CW-1:0] ssum = sin ^ cin ^ pz;
            wire [CW-1:0] cary = ((sin & cin) | (sin & pz) | (cin & pz)) << 1;

            always @(posedge clk) begin
               if (rst) begin
                  sq[ID] <= {CW{1'b0}};
                  cq[ID] <= {CW{1'b0}};
               end else if (en) begin
                  sq[ID] <= ssum;
                  cq[ID] <= cary;
               end
            end
         end
      end
   endgenerate

   //#########################################################################
   // Per-column result: one carry-propagate add at the bottom of the chain,
   // latched by the first de-skew register.  Column j's chain finishes at
   // en-index m+N-1+j; with de-skew depth (N-1-j) the word shows up at
   // en-index m + 2N-1, all columns together.
   //#########################################################################
   genvar gj2;
   integer gq;
   generate
      for (gj2 = 0; gj2 < N; gj2 = gj2 + 1) begin : g_col
         localparam D  = N - 1 - gj2;
         localparam ID = (N-1)*N + gj2;
         wire [CW-1:0] val = sq[ID] + cq[ID];   // mod 2^CW
         wire [18:0]   res = val[18:0];         // fits: |sum| <= 2^17
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

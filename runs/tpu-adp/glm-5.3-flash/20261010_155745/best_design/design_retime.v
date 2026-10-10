//#############################################################################
// Weight-stationary 8x8 systolic matrix-multiply tile (TPU-style MXU).
//
// Retimed implementation, cycle-accurate with the logikbench `tpu` reference.
//
// Timing contract (in "en" cycles, i.e. cycles where w_valid is low):
//   * term A[m][i]*B[i][j] is multiplied at en-index m+i+j using the weight
//     the PE holds at that cycle, so weight reloads affect in-flight terms
//     exactly as in the reference (including partially shifted weights).
//   * multiply and accumulate are split by a product pre-register that is
//     re-evaluated EVERY cycle (also during weight loads, against the weight
//     the PE is about to hold), so a mid-stream reload is picked up exactly
//     one cycle later, matching the reference.
//   * row 0 multiplies directly into its accumulator (reference row 0 has a
//     live input; PE(0,0) multiplies the current input byte).
//   * result row m appears on c_data at en-index m + 2N-1.
//
// Accumulators are 19 bits: 8 signed 8x8 products sum to at most +/-131072,
// exactly representable; outputs are sign-extended to the 32-bit port.
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

   localparam LAT   = 2*N - 1;    // 15: accept -> result, in en-cycles
   localparam PW    = 2*DW + 3;   // 19: accumulator width
   localparam PRODW = 2*DW;       // 16: raw product width

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
   // Activation tap lines.
   //
   // tap(i,j) presents, during en-index q, the activation byte that PE(i,j)
   // multiplies at en-index q+1 (rows >= 1, product pre-registered) or at
   // en-index q itself (row 0, multiplied straight into its accumulator):
   //      depth(row 0,   col j) = j        (PE(0,0): live input, depth 0)
   //      depth(row i>0, col j) = i + j - 1
   // depth 0 means the live (current-cycle) input byte.
   //#########################################################################
   wire [DW-1:0] taps [0:N*N-1];
   genvar        gi, gj;
   integer       gp;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_trow
         localparam DMAX = (gi == 0) ? N-1 : gi + N - 2;
         wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
         reg  [DW-1:0] sr  [1:DMAX];
         always @(posedge clk) begin
            if (rst) begin
               for (gp = 1; gp <= DMAX; gp = gp + 1)
                 sr[gp] <= {DW{1'b0}};
            end else if (en) begin
               sr[1] <= din;
               for (gp = 2; gp <= DMAX; gp = gp + 1)
                 sr[gp] <= sr[gp-1];
            end
         end
         for (gj = 0; gj < N; gj = gj + 1) begin : g_tap
            localparam D = (gi == 0) ? gj : (gi + gj - 1);
            wire [DW-1:0] tb = (D == 0) ? din : sr[D];
            assign taps[gi*N + gj] = tb;
         end
      end
   endgenerate

   //#########################################################################
   // Product pre-registers and the accumulator chain.
   //
   // prodq(ID) during en-index k = tap(k-1) * w(k) = A[k-i-j][i] * B[i][j],
   // i.e. exactly the term the reference adds at en-index k (row 0 adds its
   // product directly, reproducing the reference's live row-0 input).
   //
   // The product pre-register computes tap * wq (the weight register's
   // CURRENT output).  In the reference the product consumed at en-index k
   // used the weight captured at the end of en-index k-1; wq at en-index k-1
   // holds exactly that weight, so the pre-register matches the reference
   // even when a weight load interleaves with an in-flight computation.
   //#########################################################################
   reg [PRODW-1:0] prodq [0:N*N-1];
   reg [PW-1:0]    psumq [0:N*N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_row
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pe
            localparam ID = gi*N + gj;
            localparam UP = ID - N;

            wire signed [PRODW-1:0] pmul = $signed(taps[ID]) * $signed(wq[ID]);
            wire [PW-1:0] pex = {{(PW-PRODW){pmul[PRODW-1]}}, pmul};

            // product pre-register (re-evaluated every cycle, loads included)
            always @(posedge clk) begin
               if (rst) prodq[ID] <= {PRODW{1'b0}};
               else     prodq[ID] <= pmul;
            end

            // accumulator: row 0 takes the fresh product, deeper rows add the
            // pre-registered product to the partial sum from the row above.
            wire [PW-1:0] pin  = (gi == 0) ? {PW{1'b0}} : psumq[UP];
            wire [PW-1:0] pexd = {{(PW-PRODW){prodq[ID][PRODW-1]}}, prodq[ID]};
            wire [PW-1:0] pnew = (gi == 0) ? pex : (pin + pexd);

            always @(posedge clk) begin
               if (rst)     psumq[ID] <= {PW{1'b0}};
               else if (en) psumq[ID] <= pnew;
            end
         end
      end
   endgenerate

   //#########################################################################
   // Output de-skew: column j completed at en-index m+N-1+j; delay column j
   // by N-1-j en-cycles so a full result row leaves together at m + 2N-1.
   //#########################################################################
   genvar gj2;
   integer gq;
   generate
      for (gj2 = 0; gj2 < N; gj2 = gj2 + 1) begin : g_dsk
         localparam D  = N - 1 - gj2;
         localparam ID = (N-1)*N + gj2;
         if (D == 0) begin : g_d0
            assign c_data[gj2*ACCW +: ACCW] =
                 {{(ACCW-PW){psumq[ID][PW-1]}}, psumq[ID]};
         end else begin : g_dD
            reg [PW-1:0] dr [1:D];
            always @(posedge clk) begin
               if (rst) begin
                  for (gq = 1; gq <= D; gq = gq + 1)
                    dr[gq] <= {PW{1'b0}};
               end else if (en) begin
                  dr[1] <= psumq[ID];
                  for (gq = 2; gq <= D; gq = gq + 1)
                    dr[gq] <= dr[gq-1];
               end
            end
            assign c_data[gj2*ACCW +: ACCW] =
                 {{(ACCW-PW){dr[D][PW-1]}}, dr[D]};
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

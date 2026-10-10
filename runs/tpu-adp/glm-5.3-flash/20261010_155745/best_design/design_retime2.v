//#############################################################################
// Weight-stationary 8x8 systolic matrix-multiply tile (TPU-style MXU).
//
// Retimed implementation, cycle-accurate with the logikbench `tpu` reference.
//
// Timing contract (in "en" cycles, i.e. cycles where w_valid is low):
//   * term A[m][i]*B[i][j] is multiplied at en-index m+i+j using the weight
//     the PE holds at that cycle, exactly as in the reference (so weight
//     reloads affect in-flight terms identically, including partially
//     shifted weights during a reload).
//   * multiply and accumulate are split: the product is registered one cycle
//     before its accumulate cycle.  The pre-register re-evaluates EVERY
//     cycle.  During weight loads it samples one stage deeper in the
//     activation tap line (the byte frozen at the reference's multiply
//     point) and the post-shift weight, so a mid-stream reload is picked up
//     exactly as the reference does.
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

   // weight the PE will hold after the current cycle (tracks reloads)
   wire [DW-1:0] wnex [0:N*N-1];
   genvar        gv;
   generate
      for (gv = 0; gv < N*N; gv = gv + 1) begin : g_wnex
         assign wnex[gv] = ld ? ((gv < N) ? w_data[gv*DW +: DW] : wq[gv-N])
                              : wq[gv];
      end
   endgenerate

   //#########################################################################
   // Activation tap lines.
   //
   // taps(ID): during an en cycle, the activation byte that PE(i,j) multiplies
   //   (depth i+j for row 0 direct path convention below; see pmul mux).
   // tapl(ID): during a load cycle, the byte frozen at the multiply point
   //   (one stage deeper), used so an interrupted multiply resumes with the
   //   same activation byte the reference uses.
   //#########################################################################
   wire [DW-1:0] taps [0:N*N-1];
   wire [DW-1:0] tapl [0:N*N-1];
   genvar        gi, gj;
   integer       gp;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_trow
         wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
         if (gi == 0) begin : g_row0
            // row 0: depths 0..N-1 (depth 0 = live input); direct path
            reg [DW-1:0] sr [1:N-1];
            always @(posedge clk) begin
               if (rst) begin
                  for (gp = 1; gp <= N-1; gp = gp + 1) sr[gp] <= {DW{1'b0}};
               end else if (en) begin
                  sr[1] <= din;
                  for (gp = 2; gp <= N-1; gp = gp + 1) sr[gp] <= sr[gp-1];
               end
            end
            for (gj = 0; gj < N; gj = gj + 1) begin : g_t0
               assign taps[gi*N + gj] = (gj == 0) ? din : sr[gj];
               assign tapl[gi*N + gj] = taps[gi*N + gj]; // unused
            end
         end else begin : g_rown
            // rows >= 1: depths 1..gi+N-1 (ld path uses one extra stage)
            localparam DM = gi + N - 1;
            reg [DW-1:0] sr [1:DM];
            always @(posedge clk) begin
               if (rst) begin
                  for (gp = 1; gp <= DM; gp = gp + 1) sr[gp] <= {DW{1'b0}};
               end else if (en) begin
                  sr[1] <= din;
                  for (gp = 2; gp <= DM; gp = gp + 1) sr[gp] <= sr[gp-1];
               end
            end
            for (gj = 0; gj < N; gj = gj + 1) begin : g_tn
               localparam D = gi + gj - 1;
               assign taps[gi*N + gj] = (D == 0) ? din : sr[D];
               assign tapl[gi*N + gj] = sr[D+1];
            end
         end
      end
   endgenerate

   //#########################################################################
   // Product pre-registers and the accumulator chain.
   //#########################################################################
   reg [PRODW-1:0] prodq [0:N*N-1];
   reg [PW-1:0]    psumq [0:N*N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_row
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pe
            localparam ID = gi*N + gj;
            localparam UP = ID - N;

            // multiply operand selection:
            //  row 0  : live tap, current weight (direct into accumulator)
            //  rows>=1: en -> tap byte about to be multiplied, next weight
            //           ld -> frozen multiply-point byte, next weight
            wire [DW-1:0] am = (gi == 0) ? taps[ID]
                                         : (ld ? tapl[ID] : taps[ID]);
            wire [DW-1:0] wm = (gi == 0) ? wq[ID] : wnex[ID];
            wire signed [PRODW-1:0] pmul = $signed(am) * $signed(wm);
            wire [PW-1:0] pex = {{(PW-PRODW){pmul[PRODW-1]}}, pmul};

            // product pre-register (re-evaluated every cycle, loads included)
            always @(posedge clk) begin
               if (rst) prodq[ID] <= {PRODW{1'b0}};
               else     prodq[ID] <= pmul;
            end

            // accumulator chain (south)
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

//#############################################################################
// Weight-stationary 8x8 systolic matrix-multiply tile (TPU-style MXU).
//
// Retimed implementation, cycle-accurate with the logikbench `tpu` reference:
//
//   * Weights shift south through per-column chains on `ld` exactly like the
//     reference; after N load cycles PE(i,j) holds B[i][j].
//   * Activation bytes are delayed by per-byte tap lines so that every
//     product a*b is PRE-REGISTERED one cycle before its accumulate cycle.
//     This splits the multiply-accumulate critical path into a multiply
//     stage (tap reg -> 8x8 signed mul -> prod reg) and an accumulate stage
//     (prod reg + psum reg -> add -> psum reg).
//   * The product registers re-evaluate EVERY cycle (also during weight
//     load), so a mid-stream weight reload is automatically picked up
//     before the next compute cycle: they always hold
//     (pending activation) x (weight the PE is about to hold).
//   * Column accumulators chain south, one add per compute cycle. The
//     datapath latency is exactly LAT = 2N-1 en-cycles, so no output
//     de-skew / padding registers are needed at all.
//   * Accumulators are 19 bits: 8 signed 8x8 products sum to at most
//     +/-131072, exactly representable; results are sign-extended to the
//     32-bit port for free.
//#############################################################################

module tpu #(parameter N = 8, DW = 8, ACCW = 32) (
    input  wire            clk,
    input  wire            rst,
    input  wire            w_valid,
    input  wire [N*DW-1:0] w_data,
    input  wire            a_valid,
    input  wire [N*DW-1:0] a_data,
    output wire            c_valid,
    output wire [N*ACCW-1:0] c_data);

   localparam LAT    = 2*N - 1;   // 15: accept -> result, in en-cycles
   localparam DEL    = LAT - N - 1; // 6: base tap delay
   localparam PW     = 2*DW + 3;  // 19: accumulator width
   localparam PRODW  = 2*DW;      // 16: raw product width

   wire en = ~w_valid;   // array advances only while not loading weights
   wire ld =  w_valid;

   //#########################################################################
   // Per-byte activation tap lines.
   //
   // tapv[b] presents byte b of the activation row that row b of the PE grid
   // is about to multiply, one cycle ahead of its accumulate cycle:
   //   row 0 multiplies combinationally at accept+7  -> tap depth 7
   //   row i (i>=1) consumes a pre-registered product at accept+7+i, so its
   //   tap presents the byte at accept+6+i           -> tap depth 6+i
   // Taps advance only on en, so they hold their value across weight loads.
   //#########################################################################
   wire [DW-1:0] tapv [0:N-1];
   genvar        b, i, j;
   generate
      for (b = 0; b < N; b = b + 1) begin : g_tap
         localparam DB = (b == 0) ? (DEL + 1) : (DEL + b);
         wire [DW-1:0] din = a_valid ? a_data[b*DW +: DW] : {DW{1'b0}};
         reg  [DW-1:0] sr [0:DB-1];
         integer       kk;
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < DB; kk = kk + 1) sr[kk] <= {DW{1'b0}};
            end
            else if (en) begin
               sr[0] <= din;
               for (kk = 1; kk < DB; kk = kk + 1) sr[kk] <= sr[kk-1];
            end
         end
         assign tapv[b] = sr[DB-1];
      end
   endgenerate

   //#########################################################################
   // PE grid.  Flattened (i,j) indexing for interconnect wires.
   //#########################################################################
   wire [DW-1:0] wcur  [0:N*N-1];  // stationary weight register outputs
   wire [DW-1:0] wnext [0:N*N-1];  // weight the PE will hold after this cycle
   wire [PW-1:0] psumq [0:N*N-1];  // accumulator register outputs

   generate
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_pe
            localparam ID = i*N + j;
            localparam UP = ID - N;

            // ----- stationary weight (per-column shift chain on ld) -----
            reg  [DW-1:0] w;
            wire [DW-1:0] win = (i == 0) ? w_data[j*DW +: DW] : wcur[UP];
            always @(posedge clk) begin
               if (rst)     w <= {DW{1'b0}};
               else if (ld) w <= win;
            end
            assign wcur[ID]  = w;
            assign wnext[ID] = ld ? win : w;

            // ----- pre-registered product -----
            // Uses the weight the PE is about to hold, so the register is
            // always up to date, including across weight reloads.
            wire signed [PRODW-1:0] prodn =
                 $signed(tapv[i]) * $signed(wnext[ID]);
            wire [PW-1:0] prodse = {{(PW-PRODW){prodn[PRODW-1]}}, prodn};

            reg [PRODW-1:0] prodreg;   // unused (pruned) for row 0
            always @(posedge clk) begin
               if (rst) prodreg <= {PRODW{1'b0}};
               else     prodreg <= prodn;
            end

            // ----- accumulator -----
            // Row 0 has a zero north input: its accumulator is just the
            // (combinationally multiplied, registered) product.
            reg  [PW-1:0] psum;
            wire [PW-1:0] pin = (i == 0) ? {PW{1'b0}} : psumq[UP];
            wire signed [PW-1:0] padd = $signed(pin) + $signed(prodreg);
            wire [PW-1:0] pnext = (i == 0) ? prodse : padd;

            always @(posedge clk) begin
               if (rst)     psum <= {PW{1'b0}};
               else if (en) psum <= pnext;
            end
            assign psumq[ID] = psum;
         end
      end
   endgenerate

   //#########################################################################
   // Result row: bottom-row accumulators, sign-extended to the 32-bit port.
   //#########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
         assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-PW){psumq[(N-1)*N + j][PW-1]}}, psumq[(N-1)*N + j]};
      end
   endgenerate

   //#########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT en-cycles.
   //#########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)     vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

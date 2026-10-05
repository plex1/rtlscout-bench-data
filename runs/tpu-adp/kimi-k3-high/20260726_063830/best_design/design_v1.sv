//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) - optimized.
//
// Cycle-accurate re-implementation of the logikbench tpu. Key transform vs.
// the reference: the per-PE multiply and accumulate are split into two
// pipeline stages. Each PE registers prod_reg <= a_west * w (the deferred
// product) and accumulates the REGISTERED product into the partial sum on the
// next compute cycle. This halves the critical path (mul-only and add-only
// stages instead of mul+add). The extra cycle of latency is recovered by
// making the bottom row's accumulate combinational (no psum register on the
// last row), so c_bot is bit-identical to the reference on every cycle.
//
// Stall semantics: prod_reg is gated by en (~w_valid) exactly like the other
// state. On a mid-drain weight reload the deferred product keeps the OLD
// weight, matching the reference which commits a_west*w_old on the last
// compute cycle before the stall. Invariant: psum[i][j] at cycle t equals the
// reference psum_south[i][j] at the previous compute cycle, and prod_reg[i][j]
// holds the product the reference committed at that cycle; hence
//   c_bot[j] = psum[6][j] + prod_reg[7][j] == reference c_bot[j] always.
//
// Width reduction: |C[m][j]| <= 8*2^14 = 2^17, so 19-bit signed partial sums
// are exact. De-skew registers are also 19-bit; sign-extension to 32 bits
// happens only at the output pins.
//#############################################################################

module tpu
  #(parameter N = 8,     // array dimension (N x N)
    parameter DW = 8,    // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input                 clk,
    input                 rst,     // active high (asserted before streaming)
    input                 w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]      w_data,  // weight row, packed { col N-1 ... col 0 }
    input                 a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]      a_data,  // activation row, packed { col N-1 ... col 0 }
    output                c_valid, // result-row valid
    output [N*ACCW-1:0]   c_data   // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;  // injection-to-output latency (compute cycles)
   localparam PW  = 19;       // partial-sum width: 8*2^14 = 2^17 max magnitude

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles.
   //##########################################################################
   wire [DW-1:0] a_left [0:N-1];
   genvar        gi;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_skew
         wire [DW-1:0] din =
                       a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
         if (gi == 0) begin : g_skew0
            assign a_left[0] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:gi-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < gi; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < gi; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[gi] = sr[gi-1];
         end
      end
   endgenerate

   //##########################################################################
   // PE mesh state.
   //##########################################################################
   reg [DW-1:0]   w  [0:N-1][0:N-1];  // stationary weights
   reg [DW-1:0]   ae [0:N-1][0:N-2];  // activation regs (east out of last col unused)
   reg [2*DW-1:0] pr [0:N-1][0:N-1];  // deferred products
   reg [PW-1:0]   ps [0:N-2][0:N-1];  // partial sums (bottom row combinational)

   wire [DW-1:0]  aw [0:N-1][0:N-1];  // a_west per PE
   wire [PW-1:0]  pn [0:N-1][0:N-1];  // psum_north per PE
   wire [PW-1:0]  cbot [0:N-1];       // column results (narrow)

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            // west port
            if (j == 0) begin : g_w0
               assign aw[i][0] = a_left[i];
            end
            else begin : g_wn
               assign aw[i][j] = ae[i][j-1];
            end
            // north port
            if (i == 0) begin : g_n0
               assign pn[0][j] = {PW{1'b0}};
            end
            else begin : g_nn
               assign pn[i][j] = ps[i-1][j];
            end

            wire signed [2*DW-1:0] prod =
                 $signed(aw[i][j]) * $signed(w[i][j]);

            always @(posedge clk) begin
               if (rst) begin
                  w[i][j]  <= {DW{1'b0}};
                  pr[i][j] <= {2*DW{1'b0}};
               end
               else begin
                  if (ld)
                    w[i][j] <= (i == 0) ? w_data[j*DW +: DW] : w[i-1][j];
                  if (en)
                    pr[i][j] <= prod;
               end
            end

            // activation east register (skip unused last column)
            if (j < N-1) begin : g_ae
               always @(posedge clk) begin
                  if (rst)
                    ae[i][j] <= {DW{1'b0}};
                  else if (en)
                    ae[i][j] <= aw[i][j];
               end
            end

            // partial-sum register (skip bottom row: combinational there)
            if (i < N-1) begin : g_ps
               always @(posedge clk) begin
                  if (rst)
                    ps[i][j] <= {PW{1'b0}};
                  else if (en)
                    ps[i][j] <= pn[i][j] +
                                {{(PW-2*DW){pr[i][j][2*DW-1]}}, pr[i][j]};
               end
            end
         end
      end
   endgenerate

   // bottom row accumulate (combinational) recovers the reference c_bot
   generate
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign cbot[j] = pn[N-1][j] +
                          {{(PW-2*DW){pr[N-1][j][2*DW-1]}}, pr[N-1][j]};
      end
   endgenerate

   //##########################################################################
   // Output de-skew (narrow registers): column j delayed N-1-j cycles.
   // Sign-extend to ACCW only at the pins.
   //##########################################################################
   genvar gj;
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_deskew
         localparam D = N - 1 - gj;
         if (D == 0) begin : g_deskew0
            assign c_data[gj*ACCW +: ACCW] =
                   {{(ACCW-PW){cbot[gj][PW-1]}}, cbot[gj]};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cbot[gj];
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[gj*ACCW +: ACCW] =
                   {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
//
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Architecture notes (vs. classic PE mesh):
//  - Activations are broadcast combinationally across each row (no east-hop
//    registers).  Row i of the activation is delayed by (N-1+i) compute
//    cycles on the cheap 8-bit input side, which both time-aligns the
//    partial-sum chain down each column AND removes the need for any
//    output de-skew: the bottom of every column finishes simultaneously.
//  - Internal accumulators are 19 bits (max |sum of 8 signed 8x8 products|
//    = 8*128*128 = 131072); results are sign-extended to 32 bits only at
//    the output pins.
//  - Latency: A row accepted at compute cycle m -> PE(i,j) accumulates at
//    m + (N-1) + i, bottom row registers at end of m + 2N-2, result visible
//    at m + 2N-1 = 15, matching c_valid = a_valid delayed 15 compute cycles.
//  - The whole datapath advances only while w_valid is low (en); asserting
//    w_valid stalls in-flight computation and shifts new weights in.
//#############################################################################

module tpu #(parameter N = 8, DW = 8, ACCW = 32) (
    input  wire          clk,
    input  wire          rst,
    input  wire          w_valid,
    input  wire [N*DW-1:0] w_data,   // weight row, { col N-1 ... col 0 }
    input  wire          a_valid,
    input  wire [N*DW-1:0] a_data,   // activation row, { col N-1 ... col 0 }
    output wire          c_valid,
    output wire [N*ACCW-1:0] c_data
);

  localparam LAT = 2*N - 1;   // 15 compute cycles injection-to-output
  localparam PW  = 19;        // internal accumulator width (signed)

  wire en = ~w_valid;

  //##########################################################################
  // Activation distribution network.
  // Row i is delayed by (N-1+i) compute cycles; zeros are injected whenever
  // a_valid is low so idle/drain cycles add nothing.  This depth makes the
  // last psum-chain stage double as the final alignment register, which is
  // what makes the original's output de-skew network unnecessary.
  //##########################################################################
  wire [DW-1:0] ape [0:N-1];

  genvar gi, gj;
  generate
    for (gi = 0; gi < N; gi = gi + 1) begin : g_askew
      localparam D = (N-1) + gi;
      reg [DW-1:0] sr [0:D-1];
      integer k;
      always @(posedge clk) begin
        if (rst) begin
          for (k = 0; k < D; k = k + 1) sr[k] <= {DW{1'b0}};
        end else if (en) begin
          sr[0] <= a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
          for (k = 1; k < D; k = k + 1) sr[k] <= sr[k-1];
        end
      end
      assign ape[gi] = sr[D-1];
    end
  endgenerate

  //##########################################################################
  // Column machines: one chain of N accumulate stages per output column.
  // Stage i of column j adds A[m][i]*B[i][j] into the partial sum flowing
  // from the stage above; the top stage starts a fresh accumulation.
  //##########################################################################
  generate
    for (gj = 0; gj < N; gj = gj + 1) begin : g_col

      reg [DW-1:0]    w    [0:N-1];   // stationary weights (shift on load)
      reg [PW-1:0]    p    [0:N-1];   // partial sums
      wire [2*DW-1:0] prod [0:N-1];

      genvar i;
      for (i = 0; i < N; i = i + 1) begin : g_pr
        assign prod[i] = $signed(ape[i]) * $signed(w[i]);
      end

      integer s;
      always @(posedge clk) begin
        if (rst) begin
          for (s = 0; s < N; s = s + 1) begin
            w[s] <= {DW{1'b0}};
            p[s] <= {PW{1'b0}};
          end
        end else begin
          if (w_valid) begin
            w[0] <= w_data[gj*DW +: DW];
            for (s = 1; s < N; s = s + 1) w[s] <= w[s-1];
          end
          if (en) begin
            p[0] <= {{(PW-2*DW){prod[0][2*DW-1]}}, prod[0]};
            for (s = 1; s < N; s = s + 1)
              p[s] <= p[s-1] + {{(PW-2*DW){prod[s][2*DW-1]}}, prod[s]};
          end
        end
      end

      assign c_data[gj*ACCW +: ACCW] =
          {{(ACCW-PW){p[N-1][PW-1]}}, p[N-1]};
    end
  endgenerate

  //##########################################################################
  // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
  //##########################################################################
  reg [LAT-1:0] vpipe;
  always @(posedge clk) begin
    if (rst) vpipe <= {LAT{1'b0}};
    else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
  end
  assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Schedule (equivalent to the original PE mesh, which forms product (i,j)
// of activation row m at compute index m+i+j, sampling the stationary
// weight as of that same compute index):
//
//   Here every product (i,j) of row m is formed at compute index m+i+7
//   (independent of j), with:
//     - activation A[m][i] delayed by (7+i) en-gated stages (broadcast to
//       all columns combinationally), and
//     - the weight B[i][*] of column j read through a (7-j)-deep en-gated
//       delay line, so the weight value seen at compute index m+i+7 is
//       exactly the one the original mesh sampled at m+i+j.  Both the
//       activation skew and the weight delay lines advance only on compute
//       cycles, so they freeze together during weight-load stalls and the
//       compute-index alignment (hence the reload semantics) is exact.
//
//   Per column j a chain of 8 accumulate stages then runs:
//         p[0] <= prod(0,j);   p[i] <= p[i-1] + prod(i,j)
//   with stage i firing at m+i+7, so the finished row-m sum sits in p[7]
//   from compute index m+15 onward -- for EVERY column simultaneously.
//   No output de-skew network is needed at all.
//
//   Internal accumulators are tapered to the exact two's-complement width
//   needed for a sum of k signed 8x8 products (16/17/17/18/18/18/18/19);
//   results are sign-extended to 32 bits only at the output pins.
//#############################################################################

module tpu #(parameter N = 8, DW = 8, ACCW = 32) (
    input  wire            clk,
    input  wire            rst,
    input  wire            w_valid,
    input  wire [N*DW-1:0] w_data,
    input  wire            a_valid,
    input  wire [N*DW-1:0] a_data,
    output wire            c_valid,
    output wire [N*ACCW-1:0] c_data
);

  localparam LAT = 2*N-1;    // 15
  localparam PW  = 2*DW+3;   // 19
  wire en = ~w_valid;

  genvar gi, gj;

  // ---------------- activation distribution ----------------
  // row i delayed by (N-1+i) stages; zero injected while a_valid low
  wire [DW-1:0] a_h [0:N-1];

  generate
    for (gi = 0; gi < N; gi = gi + 1) begin : g_askew
      localparam D = (N-1) + gi;
      wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      reg [DW-1:0] sr [0:D-1];
      integer k;
      always @(posedge clk) begin
        if (rst) begin
          for (k = 0; k < D; k = k + 1) sr[k] <= {DW{1'b0}};
        end else if (en) begin
          sr[0] <= din;
          for (k = 1; k < D; k = k + 1) sr[k] <= sr[k-1];
        end
      end
      assign a_h[gi] = sr[D-1];
    end
  endgenerate

  // ---------------- per-column weights, delay line, psum chain ----------------
  generate
    for (gj = 0; gj < N; gj = gj + 1) begin : g_col

      localparam WD = (N-1) - gj;   // weight delay stages for this column

      // stationary weight row shift register (loads from the top)
      reg [DW-1:0] wj [0:N-1];
      integer s;
      always @(posedge clk) begin
        if (rst) begin
          for (s = 0; s < N; s = s + 1) wj[s] <= {DW{1'b0}};
        end else if (w_valid) begin
          wj[0] <= w_data[gj*DW +: DW];
          for (s = 1; s < N; s = s + 1) wj[s] <= wj[s-1];
        end
      end

      // en-gated weight delay line: tap(t) = wj value (WD compute cycles ago)
      wire [DW-1:0] wtap;
      if (WD == 0) begin : g_w0
        assign wtap = wj[N-1];
      end else begin : g_wn
        reg [DW-1:0] wdl [0:WD-1];
        integer d;
        always @(posedge clk) begin
          if (rst) begin
            for (d = 0; d < WD; d = d + 1) wdl[d] <= {DW{1'b0}};
          end else if (en) begin
            wdl[0] <= wj[N-1];
            for (d = 1; d < WD; d = d + 1) wdl[d] <= wdl[d-1];
          end
        end
        assign wtap = wdl[WD-1];
      end

      // products: all formed at the same compute index for this column
      wire [2*DW-1:0] prod [0:N-1];
      for (gi = 0; gi < N; gi = gi + 1) begin : g_pr
        assign prod[gi] = $signed(a_h[gi]) * $signed(wtap);
      end

      // psum chain (tapered widths)
      wire [PW-1:0] pchain [0:N-1];

      begin : g_st0
        reg [2*DW-1:0] p0;
        always @(posedge clk) begin
          if (rst) p0 <= {2*DW{1'b0}};
          else if (en) p0 <= prod[0];
        end
        assign pchain[0] = {{(PW-2*DW){p0[2*DW-1]}}, p0};
      end

      for (gi = 1; gi < N; gi = gi + 1) begin : g_st
        localparam WP = (gi >= 7) ? PW : (gi >= 3) ? 18 : 17;
        reg [WP-1:0] p;
        always @(posedge clk) begin
          if (rst) p <= {WP{1'b0}};
          else if (en) p <= pchain[gi-1][WP-1:0] +
                            {{(WP-2*DW){prod[gi][2*DW-1]}}, prod[gi]};
        end
        if (WP == PW) begin : g_pc0
          assign pchain[gi] = p;
        end else begin : g_pcs
          assign pchain[gi] = {{(PW-WP){p[WP-1]}}, p};
        end
      end

      assign c_data[gj*ACCW +: ACCW] =
          {{(ACCW-PW){pchain[N-1][PW-1]}}, pchain[N-1]};
    end
  endgenerate

  // ---------------- valid pipeline ----------------
  reg [LAT-1:0] vpipe;
  always @(posedge clk) begin
    if (rst) vpipe <= {LAT{1'b0}};
    else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
  end
  assign c_valid = vpipe[LAT-1];

endmodule

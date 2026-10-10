//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Structure (identical per-product schedule to the original PE mesh):
//   - Activation network: per row i an en-gated skew of i stages, then an
//     east chain of N-1 en-gated registers.  Stage (i,j) sees a_data[i]
//     delayed by i+j compute cycles (row 0 / column 0 taps the live input,
//     exactly like the original).  Zeros are injected while a_valid is low.
//   - Weights: per column j an N-deep shift register loaded from the top
//     while w_valid is high; after N load cycles wj[i] = B[i][j].
//   - Per column j, a chain of N accumulate stages:
//         p0 <= a_h[0][j]*w[0];   p[i] <= p[i-1] + a_h[i][j]*w[i]
//     so product (i,j) of activation row m is formed at compute cycle
//     m+i+j -- the same instant as in the original mesh, which is what makes
//     mid-stream weight reloads behave identically.
//   - Output de-skew: column j delayed N-1-j cycles so a full result row
//     leaves together, 2N-1 compute cycles after its A row was accepted.
//   - Internal accumulators are TAPERED to the exact two's-complement width
//     needed for a sum of k signed 8x8 products (16/17/17/18/18/18/18/19);
//     results are sign-extended to 32 bits only at the output pins.
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

  localparam LAT = 2*N-1;
  localparam PW  = 2*DW+3;   // 19: sum of N signed 8x8 products
  wire en = ~w_valid;

  genvar gi, gj, gk;

  // ---------------- activation distribution ----------------
  wire [DW-1:0] a_h [0:N-1][0:N-1];

  generate
    for (gi = 0; gi < N; gi = gi + 1) begin : g_arow
      wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      wire [DW-1:0] s0;
      if (gi == 0) begin : g_sk0
        assign s0 = din;
      end else begin : g_sk
        reg [DW-1:0] sr [0:gi-1];
        integer k;
        always @(posedge clk) begin
          if (rst) begin
            for (k = 0; k < gi; k = k + 1) sr[k] <= {DW{1'b0}};
          end else if (en) begin
            sr[0] <= din;
            for (k = 1; k < gi; k = k + 1) sr[k] <= sr[k-1];
          end
        end
        assign s0 = sr[gi-1];
      end
      assign a_h[gi][0] = s0;

      // east chain: one register per hop, column 1 .. N-1
      reg [DW-1:0] ech [0:N-2];
      integer e;
      always @(posedge clk) begin
        if (rst) begin
          for (e = 0; e < N-1; e = e + 1) ech[e] <= {DW{1'b0}};
        end else if (en) begin
          ech[0] <= s0;
          for (e = 1; e < N-1; e = e + 1) ech[e] <= ech[e-1];
        end
      end
      for (gk = 1; gk < N; gk = gk + 1) begin : g_ej
        assign a_h[gi][gk] = ech[gk-1];
      end
    end
  endgenerate

  // ---------------- per-column weight registers + psum chains ----------------
  generate
    for (gj = 0; gj < N; gj = gj + 1) begin : g_col

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

      wire [2*DW-1:0] prod [0:N-1];
      for (gi = 0; gi < N; gi = gi + 1) begin : g_pr
        assign prod[gi] = $signed(a_h[gi][gj]) * $signed(wj[gi]);
      end

      wire [PW-1:0] pchain [0:N-1];

      // stage 0: a single product fits in 2*DW bits
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
        wire [WP-1:0] px = {{(WP-2*DW){prod[gi][2*DW-1]}}, prod[gi]};
        always @(posedge clk) begin
          if (rst) p <= {WP{1'b0}};
          else if (en) p <= pchain[gi-1][WP-1:0] + px;
        end
        if (WP == PW) begin : g_pc0
          assign pchain[gi] = p;
        end else begin : g_pcs
          assign pchain[gi] = {{(PW-WP){p[WP-1]}}, p};
        end
      end

      // output de-skew: align all columns of one result row
      wire [PW-1:0] cout;
      if (N-1-gj == 0) begin : g_d0
        assign cout = pchain[N-1];
      end else begin : g_d
        reg [PW-1:0] dr [0:N-2-gj];
        integer d;
        always @(posedge clk) begin
          if (rst) begin
            for (d = 0; d < N-1-gj; d = d + 1) dr[d] <= {PW{1'b0}};
          end else if (en) begin
            dr[0] <= pchain[N-1];
            for (d = 1; d < N-1-gj; d = d + 1) dr[d] <= dr[d-1];
          end
        end
        assign cout = dr[N-2-gj];
      end
      assign c_data[gj*ACCW +: ACCW] = {{(ACCW-PW){cout[PW-1]}}, cout};
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

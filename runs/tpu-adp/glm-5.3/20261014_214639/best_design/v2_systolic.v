// Optimized weight-stationary systolic MXU, cycle-accurate to logikbench tpu.
//
// Key optimizations vs starting_point.v:
//  * partial sums narrowed to the minimum provable width per depth
//    (16,17,17,18,18,18,18,19 instead of 32) -> smaller adders + fewer FFs
//  * depth 0 needs no adder (its "partial sum" is just the product)
//  * output de-skew registers narrowed to 19 bits (result fits), sign
//    extended to 32 bits only at the pins
//  * activation delay network is one shared per-row shift chain with taps
//    (drops the dead east-exit register of the mesh)
//  * weight storage is the same per-column south-shifting register file
module tpu (
  input clk, rst, w_valid, a_valid,
  input [63:0] w_data, a_data,
  output c_valid,
  output [255:0] c_data
);
  wire en = ~w_valid;

  // ---------------- weights: per-column south-shifting registers --------
  reg [7:0] wq [0:7][0:7];
  integer ii, jj;
  always @(posedge clk) begin
    if (rst) begin
      for (jj = 0; jj < 8; jj = jj+1)
        for (ii = 0; ii < 8; ii = ii+1)
          wq[ii][jj] <= 8'h00;
    end else if (w_valid) begin
      for (jj = 0; jj < 8; jj = jj+1) begin
        wq[0][jj] <= w_data[jj*8 +: 8];
        for (ii = 1; ii < 8; ii = ii+1)
          wq[ii][jj] <= wq[ii-1][jj];
      end
    end
  end

  // ---------------- activation delay network ----------------------------
  // act[i][j] = injected row byte i delayed by (i+j) compute cycles
  wire [7:0] act [0:7][0:7];
  genvar gi, gj;
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_arow
      wire [7:0] din = a_valid ? a_data[gi*8 +: 8] : 8'h00;
      reg [7:0] c [0:13];
      integer q;
      always @(posedge clk) begin
        if (rst) begin
          for (q = 0; q < 14; q = q+1) c[q] <= 8'h00;
        end else if (en) begin
          c[0] <= din;
          for (q = 1; q < 14; q = q+1) c[q] <= c[q-1];
        end
      end
      for (gj = 0; gj < 8; gj = gj+1) begin : g_atap
        if (gi + gj == 0) begin : g_d0
          assign act[gi][gj] = din;
        end else begin : g_dn
          assign act[gi][gj] = c[gi+gj-1];
        end
      end
    end
  endgenerate

  // ---------------- products ---------------------------------------------
  wire signed [15:0] prod [0:7][0:7];
  generate
    for (gi = 0; gi < 8; gi = gi+1)
      for (gj = 0; gj < 8; gj = gj+1) begin : g_prod
        assign prod[gi][gj] = $signed(act[gi][gj]) * $signed(wq[gi][gj]);
      end
  endgenerate

  // ---------------- per-column accumulate chain --------------------------
  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_col
      reg signed [15:0] p0;
      reg signed [16:0] p1, p2;
      reg signed [17:0] p3, p4, p5, p6;
      reg signed [18:0] p7;

      // sign-extended products where the accumulator grows
      wire signed [16:0] e1 = prod[1][gj];
      wire signed [17:0] e3 = prod[3][gj];
      wire signed [18:0] e7 = prod[7][gj];

      always @(posedge clk) begin
        if (rst) begin
          p0 <= 16'sd0; p1 <= 17'sd0; p2 <= 17'sd0; p3 <= 18'sd0;
          p4 <= 18'sd0; p5 <= 18'sd0; p6 <= 18'sd0; p7 <= 19'sd0;
        end else if (en) begin
          p0 <= prod[0][gj];
          p1 <= p0 + e1;
          p2 <= p1 + prod[2][gj];
          p3 <= p2 + e3;
          p4 <= p3 + prod[4][gj];
          p5 <= p4 + prod[5][gj];
          p6 <= p5 + prod[6][gj];
          p7 <= p6 + e7;
        end
      end

      // de-skew: column j is delayed (7-j) compute cycles
      if (7 - gj == 0) begin : g_dd0
        assign c_data[gj*32 +: 32] = {{13{p7[18]}}, p7};
      end else begin : g_ddn
        reg [18:0] dr [0:7-gj-1];
        integer q;
        wire [18:0] last;
        always @(posedge clk) begin
          if (rst) begin
            for (q = 0; q < 7-gj; q = q+1) dr[q] <= 19'h0;
          end else if (en) begin
            dr[0] <= p7;
            for (q = 1; q < 7-gj; q = q+1) dr[q] <= dr[q-1];
          end
        end
        assign last = dr[7-gj-1];
        assign c_data[gj*32 +: 32] = {{13{last[18]}}, last};
      end
    end
  endgenerate

  // ---------------- valid pipeline ---------------------------------------
  reg [14:0] vpipe;
  always @(posedge clk) begin
    if (rst) vpipe <= 15'h0;
    else if (en) vpipe <= {vpipe[13:0], a_valid};
  end
  assign c_valid = vpipe[14];

endmodule

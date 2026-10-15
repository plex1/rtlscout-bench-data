// Optimized weight-stationary systolic MXU, cycle-accurate to logikbench tpu.
//
// vs v2: multiply and add are split into separate pipeline stages (products
// are registered in q[], then added one cycle later), so no cycle contains
// both a multiplier and an adder. The head of each column chain gets one
// extra pipeline register (p0) and the tail saves one deskew stage; column 7
// (zero deskew slack) uses a combinational final add at the output pin.
module tpu (
  input clk, rst, w_valid, a_valid,
  input [63:0] w_data, a_data,
  output c_valid,
  output [255:0] c_data
);
  wire en = ~w_valid;

  // ---------------- weights: per-column south-shifting registers ----------
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

  // ---------------- activation skew --------------------------------------
  // act[i][j] = injected row byte i delayed by (i+j) compute cycles
  wire [7:0] act [0:7][0:7];
  genvar gi, gj;
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_arow
      wire [7:0] din = a_valid ? a_data[gi*8 +: 8] : 8'h00;
      reg [7:0] c [0:gi+6];
      integer q;
      always @(posedge clk) begin
        if (rst) begin
          for (q = 0; q <= gi+6; q = q+1) c[q] <= 8'h00;
        end else if (en) begin
          c[0] <= din;
          for (q = 1; q <= gi+6; q = q+1) c[q] <= c[q-1];
        end
      end
      for (gj = 0; gj < 8; gj = gj+1) begin : g_tap
        if (gi + gj == 0) begin : g_t0
          assign act[gi][gj] = din;
        end else begin : g_tn
          assign act[gi][gj] = c[gi+gj-1];
        end
      end
    end
  endgenerate

  // ---------------- products (registered one cycle) -----------------------
  wire signed [15:0] prod [0:7][0:7];
  reg  signed [15:0] q    [0:7][0:7];
  generate
    for (gi = 0; gi < 8; gi = gi+1)
      for (gj = 0; gj < 8; gj = gj+1) begin : g_mul
        assign prod[gi][gj] = $signed(act[gi][gj]) * $signed(wq[gi][gj]);
      end
  endgenerate

  integer di, dj;
  always @(posedge clk) begin
    if (rst) begin
      for (dj = 0; dj < 8; dj = dj+1)
        for (di = 0; di < 8; di = di+1)
          q[di][dj] <= 16'sd0;
    end else if (en) begin
      for (dj = 0; dj < 8; dj = dj+1)
        for (di = 0; di < 8; di = di+1)
          q[di][dj] <= prod[di][dj];
    end
  end

  // ---------------- per-column accumulation chain -------------------------
  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_col
      reg signed [15:0] p0;
      reg signed [16:0] p1, p2;
      reg signed [17:0] p3, p4, p5, p6;

      always @(posedge clk) begin
        if (rst) begin
          p0 <= 16'sd0; p1 <= 17'sd0; p2 <= 17'sd0;
          p3 <= 18'sd0; p4 <= 18'sd0; p5 <= 18'sd0; p6 <= 18'sd0;
        end else if (en) begin
          p0 <= q[0][gj];
          p1 <= p0 + q[1][gj];
          p2 <= p1 + q[2][gj];
          p3 <= p2 + q[3][gj];
          p4 <= p3 + q[4][gj];
          p5 <= p4 + q[5][gj];
          p6 <= p5 + q[6][gj];
        end
      end

      if (gj < 7) begin : g_out_reg
        reg signed [18:0] p7;
        always @(posedge clk) begin
          if (rst) p7 <= 19'sd0;
          else if (en) p7 <= p6 + q[7][gj];
        end
        if (gj == 6) begin : g_d0
          assign c_data[gj*32 +: 32] = {{13{p7[18]}}, p7};
        end else begin : g_dn
          reg [18:0] dr [0:6-gj-1];
          integer k;
          always @(posedge clk) begin
            if (rst) begin
              for (k = 0; k < 6-gj; k = k+1) dr[k] <= 19'h0;
            end else if (en) begin
              dr[0] <= p7;
              for (k = 1; k < 6-gj; k = k+1) dr[k] <= dr[k-1];
            end
          end
          wire [18:0] last = dr[6-gj-1];
          assign c_data[gj*32 +: 32] = {{13{last[18]}}, last};
        end
      end else begin : g_out_comb
        wire signed [18:0] sum = p6 + q[7][gj];
        assign c_data[gj*32 +: 32] = {{13{sum[18]}}, sum};
      end
    end
  endgenerate

  // ---------------- valid pipeline ----------------------------------------
  reg [14:0] vpipe;
  always @(posedge clk) begin
    if (rst) vpipe <= 15'h0;
    else if (en) vpipe <= {vpipe[13:0], a_valid};
  end
  assign c_valid = vpipe[14];

endmodule

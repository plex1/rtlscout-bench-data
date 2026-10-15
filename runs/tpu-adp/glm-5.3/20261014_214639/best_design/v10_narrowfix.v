// v10: v9 with correct sign-extension in the narrowed accumulate stage
// (uses <<< 4 / direct assignment instead of concatenation).
module tpu (
  input clk, rst, w_valid, a_valid,
  input [63:0] w_data, a_data,
  output c_valid,
  output [255:0] c_data
);
  wire gclk = clk & ~w_valid;
  wire wclk = clk &  w_valid;
  wire en   = ~w_valid;

  reg [7:0] wq [0:7][0:7];
  integer ii, jj;
  always @(posedge wclk) begin
    for (jj = 0; jj < 8; jj = jj+1) begin
      wq[0][jj] <= w_data[jj*8 +: 8];
      for (ii = 1; ii < 8; ii = ii+1)
        wq[ii][jj] <= wq[ii-1][jj];
    end
  end

  wire [7:0] act [0:7][0:7];
  genvar gi, gj;
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_arow
      wire [7:0] din = a_valid ? a_data[gi*8 +: 8] : 8'h00;
      reg [7:0] c [0:gi+6];
      integer q;
      always @(posedge gclk) begin
        c[0] <= din;
        for (q = 1; q <= gi+6; q = q+1) c[q] <= c[q-1];
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

  wire signed [11:0] tlo [0:7][0:7];
  wire signed [11:0] thi [0:7][0:7];
  generate
    for (gi = 0; gi < 8; gi = gi+1)
      for (gj = 0; gj < 8; gj = gj+1) begin : g_mul
        wire signed [7:0] av = act[gi][gj];
        wire signed [7:0] wv = wq[gi][gj];
        assign tlo[gi][gj] = av * $signed({1'b0, wv[3:0]});
        assign thi[gi][gj] = av * $signed(wv[7:4]);
      end
  endgenerate

  reg signed [11:0] qlo [0:7][0:7];
  reg signed [11:0] qhi [0:7][0:7];
  integer di, dj;
  always @(posedge gclk) begin
    for (dj = 0; dj < 8; dj = dj+1)
      for (di = 0; di < 8; di = di+1) begin
        qlo[di][dj] <= tlo[di][dj];
        qhi[di][dj] <= thi[di][dj];
      end
  end

  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_col
      reg signed [15:0] p0;
      reg signed [16:0] p1, p2;
      reg signed [17:0] p3, p4, p5, p6;

      // y = (qhi <<< 4) and z = qlo, sign-extended to the target width
      wire signed [16:0] y1 = qhi[1][gj] <<< 4;
      wire signed [16:0] z1 = qlo[1][gj];
      wire signed [16:0] y2 = qhi[2][gj] <<< 4;
      wire signed [16:0] z2 = qlo[2][gj];
      wire signed [17:0] y3 = qhi[3][gj] <<< 4;
      wire signed [17:0] z3 = qlo[3][gj];
      wire signed [17:0] y4 = qhi[4][gj] <<< 4;
      wire signed [17:0] z4 = qlo[4][gj];
      wire signed [17:0] y5 = qhi[5][gj] <<< 4;
      wire signed [17:0] z5 = qlo[5][gj];
      wire signed [17:0] y6 = qhi[6][gj] <<< 4;
      wire signed [17:0] z6 = qlo[6][gj];
      wire signed [18:0] y7 = qhi[7][gj] <<< 4;
      wire signed [18:0] z7 = qlo[7][gj];

      wire signed [15:0] prod0 = (qhi[0][gj] <<< 4) + qlo[0][gj];

      wire [16:0] s1 = p0 + y1 + z1;
      wire [16:0] s2 = p1 + y2 + z2;
      wire [17:0] s3 = p2 + y3 + z3;
      wire [17:0] s4 = p3 + y4 + z4;
      wire [17:0] s5 = p4 + y5 + z5;
      wire [17:0] s6 = p5 + y6 + z6;

      always @(posedge gclk) begin
        p0 <= prod0;
        p1 <= s1;
        p2 <= s2;
        p3 <= s3;
        p4 <= s4;
        p5 <= s5;
        p6 <= s6;
      end

      wire signed [18:0] r7 = p6 + y7 + z7;

      if (gj < 7) begin : g_out_reg
        reg [18:0] p7;
        always @(posedge gclk) p7 <= r7;
        if (gj == 6) begin : g_d0
          assign c_data[gj*32 +: 32] = {{13{p7[18]}}, p7};
        end else begin : g_dn
          reg [18:0] dr [0:6-gj-1];
          integer k;
          always @(posedge gclk) begin
            dr[0] <= p7;
            for (k = 1; k < 6-gj; k = k+1) dr[k] <= dr[k-1];
          end
          wire [18:0] last = dr[6-gj-1];
          assign c_data[gj*32 +: 32] = {{13{last[18]}}, last};
        end
      end else begin : g_out_comb
        assign c_data[255-31 +: 32] = {{13{r7[18]}}, r7};
      end
    end
  endgenerate

  reg [14:0] vpipe;
  always @(posedge clk or posedge rst) begin
    if (rst) vpipe <= 15'h0;
    else if (en) vpipe <= {vpipe[13:0], a_valid};
  end
  assign c_valid = vpipe[14];

endmodule

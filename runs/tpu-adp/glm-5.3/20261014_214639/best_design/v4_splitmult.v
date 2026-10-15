// v4: split multiply (two 8x4 partial multiplies, registered), recombination
// folded into the accumulate step as an explicit CSA (3:2) + CPA.
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

  // ---------------- partial products (two 8x4 multiplies) ----------------
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
  always @(posedge clk) begin
    if (rst) begin
      for (dj = 0; dj < 8; dj = dj+1)
        for (di = 0; di < 8; di = di+1) begin
          qlo[di][dj] <= 12'sd0;
          qhi[di][dj] <= 12'sd0;
        end
    end else if (en) begin
      for (dj = 0; dj < 8; dj = dj+1)
        for (di = 0; di < 8; di = di+1) begin
          qlo[di][dj] <= tlo[di][dj];
          qhi[di][dj] <= thi[di][dj];
        end
    end
  end

  // ---------------- per-column accumulation -------------------------------
  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_col
      reg signed [18:0] p [0:6];

      // depth 0: p0 = (qhi<<<4) + qlo
      wire signed [15:0] prod0 = (qhi[0][gj] <<< 4) + qlo[0][gj];
      always @(posedge clk) begin
        if (rst) begin
          p[0] <= 19'sd0; p[1] <= 19'sd0; p[2] <= 19'sd0; p[3] <= 19'sd0;
          p[4] <= 19'sd0; p[5] <= 19'sd0; p[6] <= 19'sd0;
        end else if (en) begin
          p[0] <= {{3{prod0[15]}}, prod0};
        end
      end

      // depths 1..6 : p[i] = p[i-1] + (qhi<<<4) + qlo   (CSA + CPA)
      for (gi = 1; gi < 7; gi = gi+1) begin : g_acc
        wire [18:0] y = {{3{qhi[gi][gj][11]}}, qhi[gi][gj], 4'b0};
        wire [18:0] z = {{7{qlo[gi][gj][11]}}, qlo[gi][gj]};
        wire [18:0] s = p[gi-1] ^ y ^ z;
        wire [18:0] c = (p[gi-1] & y) | (z & (p[gi-1] | y));
        always @(posedge clk) begin
          if (rst) p[gi] <= 19'sd0;
          else if (en) p[gi] <= s + {c[17:0], 1'b0};
        end
      end

      // depth 7
      wire [18:0] y7 = {{3{qhi[7][gj][11]}}, qhi[7][gj], 4'b0};
      wire [18:0] z7 = {{7{qlo[7][gj][11]}}, qlo[7][gj]};
      wire [18:0] s7 = p[6] ^ y7 ^ z7;
      wire [18:0] c7 = (p[6] & y7) | (z7 & (p[6] | y7));
      wire [18:0] r7 = s7 + {c7[17:0], 1'b0};

      if (gj < 7) begin : g_out_reg
        reg [18:0] p7;
        always @(posedge clk) begin
          if (rst) p7 <= 19'sd0;
          else if (en) p7 <= r7;
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
        assign c_data[255-31 +: 32] = {{13{r7[18]}}, r7};
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

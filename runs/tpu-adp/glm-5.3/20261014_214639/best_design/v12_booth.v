// v12: Booth-coded weight-stationary systolic array.
//  * weights are stored as pre-decoded radix-4 Booth codes (12b) computed
//    during the load shift, so the per-cycle multiply is only a mux + XOR.
//  * mult cycle: 4 Booth partials + correction vector compressed by a 3-level
//    CSA tree into two rows (S,C); NO CPA in this cycle.
//  * acc cycle: p_i = p_{i-1} + S + 2C  (one CSA + CPA), narrowed widths.
//  * all arithmetic modulo 2^19 (result fits in 19 bits), so carry truncation
//    is exact; per-depth moduli are respected by reading low bits only.
module tpu (
  input clk, rst, w_valid, a_valid,
  input [63:0] w_data, a_data,
  output c_valid,
  output [255:0] c_data
);
  wire gclk = clk & ~w_valid;
  wire wclk = clk &  w_valid;
  wire en   = ~w_valid;

  // radix-4 Booth codes of a weight byte:
  //   [3g+0]=one_g, [3g+1]=two_g, [3g+2]=neg_g  (group g covers w[2g+1:2g-1])
  function [11:0] bcode;
    input [7:0] w;
    begin
      bcode[0]  = w[0];
      bcode[1]  = w[1] & ~w[0];
      bcode[2]  = w[1];
      bcode[3]  = w[2] ^ w[1];
      bcode[4]  = (w[3] & ~w[2] & ~w[1]) | (~w[3] & w[2] & w[1]);
      bcode[5]  = w[3];
      bcode[6]  = w[4] ^ w[3];
      bcode[7]  = (w[5] & ~w[4] & ~w[3]) | (~w[5] & w[4] & w[3]);
      bcode[8]  = w[5];
      bcode[9]  = w[6] ^ w[5];
      bcode[10] = (w[7] & ~w[6] & ~w[5]) | (~w[7] & w[6] & w[5]);
      bcode[11] = w[7];
    end
  endfunction

  reg [11:0] bc [0:7][0:7];
  integer ii, jj;
  always @(posedge wclk) begin
    for (jj = 0; jj < 8; jj = jj+1) begin
      bc[0][jj] <= bcode(w_data[jj*8 +: 8]);
      for (ii = 1; ii < 8; ii = ii+1)
        bc[ii][jj] <= bc[ii-1][jj];
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

  // ---------------- PE: Booth partials + CSA tree -> (S,C) ----------------
  reg [18:0] S [0:7][0:7];
  reg [17:0] C [0:7][0:7];
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_prow
      for (gj = 0; gj < 8; gj = gj+1) begin : g_pcol
        wire [7:0]  a  = act[gi][gj];
        wire [11:0] cd = bc[gi][gj];
        wire [15:0] sa  = {{8{a[7]}}, a};
        wire [15:0] sa2 = {sa[14:0], 1'b0};
        // magnitude mux (one/two are mutually exclusive), negation via XOR
        wire [15:0] m0 = (sa  & {16{cd[0]}}) | (sa2 & {16{cd[1]}});
        wire [15:0] m1 = (sa  & {16{cd[3]}}) | (sa2 & {16{cd[4]}});
        wire [15:0] m2 = (sa  & {16{cd[6]}}) | (sa2 & {16{cd[7]}});
        wire [15:0] m3 = (sa  & {16{cd[9]}}) | (sa2 & {16{cd[10]}});
        wire [15:0] p0 = m0 ^ {16{cd[2]}};
        wire [15:0] p1 = m1 ^ {16{cd[5]}};
        wire [15:0] p2 = m2 ^ {16{cd[8]}};
        wire [15:0] p3 = m3 ^ {16{cd[11]}};
        // 19-bit rows
        wire [18:0] e0 = {{3{p0[15]}}, p0};
        wire [18:0] e1 = {{3{p1[15]}}, p1, 2'b00};
        wire [18:0] e2 = {{3{p2[15]}}, p2, 4'b0};
        wire [18:0] e3 = {{3{p3[15]}}, p3, 6'b0};
        wire [18:0] cr = {11'b0, cd[11], 1'b0, cd[8], 1'b0, cd[5], 1'b0, cd[2]};
        // CSA tree 5 -> 2
        wire [18:0] x1 = e0 ^ e1 ^ cr;
        wire [18:0] y1 = (e0 & e1) | (cr & (e0 | e1));
        wire [18:0] b1 = {y1[17:0], 1'b0};
        wire [18:0] x2 = x1 ^ b1 ^ e2;
        wire [18:0] y2 = (x1 & b1) | (e2 & (x1 | b1));
        wire [18:0] b2 = {y2[17:0], 1'b0};
        wire [18:0] x3 = x2 ^ b2 ^ e3;
        wire [18:0] y3 = (x2 & b2) | (e3 & (x2 | b2));
        always @(posedge gclk) begin
          S[gi][gj] <= x3;
          C[gi][gj] <= y3[17:0];
        end
      end
    end
  endgenerate

  // ---------------- per-column accumulate chain ---------------------------
  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_acc
      reg [15:0] p0;
      reg [16:0] p1, p2;
      reg [17:0] p3, p4, p5, p6;

      wire [15:0] c2_0 = {C[0][gj][14:0], 1'b0};
      wire [16:0] c2_1 = {C[1][gj][15:0], 1'b0};
      wire [16:0] c2_2 = {C[2][gj][15:0], 1'b0};
      wire [17:0] c2_3 = {C[3][gj][16:0], 1'b0};
      wire [17:0] c2_4 = {C[4][gj][16:0], 1'b0};
      wire [17:0] c2_5 = {C[5][gj][16:0], 1'b0};
      wire [17:0] c2_6 = {C[6][gj][16:0], 1'b0};
      wire [18:0] c2_7 = {C[7][gj][17:0], 1'b0};

      always @(posedge gclk) begin
        p0 <= S[0][gj][15:0] + c2_0;
        p1 <= {p0[15], p0}    + S[1][gj][16:0] + c2_1;
        p2 <= {p1[16], p1}    + S[2][gj][16:0] + c2_2;
        p3 <= {p2[16], p2}    + S[3][gj][17:0] + c2_3;
        p4 <= {p3[17], p3}    + S[4][gj][17:0] + c2_4;
        p5 <= {p4[17], p4}    + S[5][gj][17:0] + c2_5;
        p6 <= {p5[17], p5}    + S[6][gj][17:0] + c2_6;
      end

      wire [18:0] r7 = {p6[17], p6} + S[7][gj][18:0] + c2_7;

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

  // ---------------- valid pipeline ----------------------------------------
  reg [14:0] vpipe;
  always @(posedge clk or posedge rst) begin
    if (rst) vpipe <= 15'h0;
    else if (en) vpipe <= {vpipe[13:0], a_valid};
  end
  assign c_valid = vpipe[14];

endmodule

// Sklansky parallel-prefix adder used by the ALU.
// Depth: clog2(W) prefix levels + carry select + sum XOR.
module alu_pfx_add #(
    parameter W = 32
) (
    input  logic [W-1:0] a_i,
    input  logic [W-1:0] b_i,
    input  logic         cin_i,
    output logic [W-1:0] sum_o,
    output logic         cout_o
);

  localparam L = $clog2(W);

  logic [L:0][W-1:0] G, P;

  assign G[0] = a_i & b_i;
  assign P[0] = a_i | b_i;

  genvar l, i;
  generate
    for (l = 0; l < L; l++) begin : g_lvl
      localparam int SB = 1 << l;
      localparam int BLK = 1 << (l + 1);
      for (i = 0; i < W; i++) begin : g_bit
        if ((i % BLK) >= SB) begin : g_up
          localparam int J = (i / BLK) * BLK + SB - 1;
          assign G[l+1][i] = G[l][i] | (P[l][i] & G[l][J]);
          assign P[l+1][i] = P[l][i] & P[l][J];
        end else begin : g_lo
          assign G[l+1][i] = G[l][i];
          assign P[l+1][i] = P[l][i];
        end
      end
    end
  endgenerate

  logic [W:0] c;
  assign c[0] = cin_i;
  generate
    for (i = 0; i < W; i++) begin : g_c
      assign c[i+1] = G[L][i] | (P[L][i] & cin_i);
    end
  endgenerate

  assign sum_o  = (a_i ^ b_i) ^ c[W-1:0];
  assign cout_o = c[W];

endmodule

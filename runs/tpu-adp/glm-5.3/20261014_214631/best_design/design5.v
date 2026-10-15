// Optimized weight-stationary systolic matmul tile (TPU MXU clone).
// Same cycle behavior as the logikbench tpu; structure: per-row activation
// skew/flow shift register, products at the original cycles, carry-save
// partial-sum chain down each column (3:2 compressors), one CPA per column,
// and a single-chain output de-skew (newest result enters at the bottom of
// the chain, oldest leaves at the top), which needs only 6 registers per
// column instead of 28.
module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		     clk,
    input		     rst,
    input		     w_valid,
    input [N*DW-1:0]	     w_data,
    input		     a_valid,
    input [N*DW-1:0]	     a_data,
    output		     c_valid,
    output [N*ACCW-1:0]     c_data
   );

   wire en = ~w_valid;

   //##########################################################################
   // Weight storage: per-column south shift while w_valid.
   //##########################################################################
   reg [DW-1:0] w [0:N-1][0:N-1];
   integer       wi, wj;
   always @(posedge clk) begin
      if (w_valid) begin
         for (wj = 0; wj < N; wj = wj + 1) begin
            w[0][wj] <= w_data[wj*DW +: DW];
            for (wi = 1; wi < N; wi = wi + 1)
              w[wi][wj] <= w[wi-1][wj];
         end
      end
   end

   //##########################################################################
   // Activation skew + east flow: sr[i][k] is delay stage k+1 of row i; the
   // multiplier tap of column j is stage i+j (combinational for i=j=0).
   //##########################################################################
   reg [DW-1:0] sr [0:N-1][0:N+5];
   integer       si, sk;
   always @(posedge clk) begin
      if (en) begin
         for (si = 0; si < N; si = si + 1) begin
            sr[si][0] <= a_valid ? a_data[si*DW +: DW] : {DW{1'b0}};
            for (sk = 1; sk <= si+6; sk = sk + 1)
              sr[si][sk] <= sr[si][sk-1];
         end
      end
   end

   wire [DW-1:0] tap [0:N-1][0:N-1];
   genvar	       gi, gj;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_taprow
         for (gj = 0; gj < N; gj = gj + 1) begin : g_tapcol
            if (gi == 0 && gj == 0) begin : g_c00
               assign tap[0][0] = a_valid ? a_data[0 +: DW] : {DW{1'b0}};
            end
            else begin : g_cnn
               assign tap[gi][gj] = sr[gi][gi+gj-1];
            end
         end
      end
   endgenerate

   //##########################################################################
   // Products (formed at the original cycles, with the current weights).
   //##########################################################################
   wire signed [2*DW-1:0] prd [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_prrow
         for (gj = 0; gj < N; gj = gj + 1) begin : g_prcol
            assign prd[gi][gj] = $signed(tap[gi][gj]) * $signed(w[gi][gj]);
         end
      end
   endgenerate

   //##########################################################################
   // Carry-save partial-sum chain down each column (fixed 19-bit modulus):
   //   value = S + C,  (S[i], C[i]) = 3:2(S[i-1], C[i-1], prd[i])
   //##########################################################################
   reg [18:0] S [0:N-1][0:N-1];
   reg [18:0] C [1:N-1][0:N-1];
   wire [18:0] colres [0:N-1];

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
         always @(posedge clk) begin
            if (en)
              S[0][gj] <= {{3{prd[0][gj][15]}}, prd[0][gj]};
         end
         for (gi = 1; gi < N; gi = gi + 1) begin : g_row
            wire [18:0] a = S[gi-1][gj];
            wire [18:0] b = (gi > 1) ? C[gi-1][gj] : {19{1'b0}};
            wire [18:0] c = {{3{prd[gi][gj][15]}}, prd[gi][gj]};
            wire [18:0] snew = a ^ b ^ c;
            wire [18:0] cnew = {((a & b) | (a & c) | (b & c)), 1'b0};
            always @(posedge clk) begin
               if (en) begin
                  S[gi][gj] <= snew;
                  C[gi][gj] <= cnew;
               end
            end
         end
         assign colres[gj] = S[N-1][gj] + C[N-1][gj];
      end
   endgenerate

   //##########################################################################
   // Output de-skew: column j needs N-1-j extra cycles. One chain per column,
   // newest result enters the bottom, the oldest leaves at the top after the
   // required number of cycles. Column 7 (no delay) is just registered.
   //##########################################################################
   reg [18:0] dq [0:N-2][0:N-2];   // dq[j][0] = newest, up the chain
   integer     qi, qj;
   always @(posedge clk) begin
      if (en) begin
         for (qj = 0; qj < N-1; qj = qj + 1) begin
            dq[qj][0] <= colres[qj];
            for (qi = 1; qi < N-1; qi = qi + 1)
              dq[qj][qi] <= dq[qj][qi-1];
         end
      end
   end

   genvar oj;
   generate
      for (oj = 0; oj < N; oj = oj + 1) begin : g_out
         if (oj == N-1) begin : g_last
            assign c_data[oj*ACCW +: ACCW] =
	      {{(ACCW-19){colres[N-1][18]}}, colres[N-1]};
         end
         else begin : g_oth
            assign c_data[oj*ACCW +: ACCW] =
	      {{(ACCW-19){dq[oj][N-2-oj][18]}}, dq[oj][N-2-oj]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed 2N-1 compute cycles.
   //##########################################################################
   localparam LAT = 2*N - 1;
   reg [LAT-1:0] vpipe;
   always @(posedge clk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

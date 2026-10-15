// Area/delay-optimized weight-stationary systolic matmul tile (TPU MXU clone).
//
// Cycle-exact clone of the logikbench tpu. Structure:
//  * per-row activation skew + east flow = one shift register chain per row,
//  * products formed at exactly the original cycles with the weight value
//    current at that cycle (so mid-drain weight reloads behave identically),
//  * the partial-sum chain down each column is kept in CARRY-SAVE form
//    (one 3:2 compressor per cell) so the multiplier output only drives a
//    compressor, not a carry-propagate adder; a single CPA at the bottom of
//    each column produces the final 19-bit result,
//  * width-tapered carry-save stages, 19-bit output de-skew,
//  * no datapath resets (self-flushing); only the valid pipe is reset.
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

   // the array advances only while w_valid is low
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
   // Activation skew + east flow. sr[i][k] is delay stage k+1 of row i; the
   // multiplier tap of column j is stage i+j (combinational input for i=j=0).
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
   // Carry-save partial-sum chain down each column.
   //   value of row i  = S[i] (signed, 16+i bits) + C[i] (unsigned, 17+i bits)
   //   row 0 : S[0] = product, C[0] = 0
   //   row i : (S[i], C[i]) = 3:2 compress( S[i-1], C[i-1], prd[i] )
   //##########################################################################
   // fixed width of the carry-save representation: the final value needs 19
   // bits signed, and the 3:2 identity S + C = value (mod 2^W) needs the SAME
   // modulus at every stage, so all stages use W = 19.
   reg [18:0] S [0:N-1][0:N-1];
   reg [18:0] C [1:N-1][0:N-1];
   wire [18:0] colres [0:N-1];

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
         // row 0: no carry-in, just the product
         always @(posedge clk) begin
            if (en)
              S[0][gj] <= {{3{prd[0][gj][15]}}, prd[0][gj]};
         end
         // rows 1..7: 3:2 compress (prev sum, prev carry, this product)
         for (gi = 1; gi < N; gi = gi + 1) begin : g_row
            wire [18:0] a = S[gi-1][gj];
            wire [18:0] b = (gi > 1) ? C[gi-1][gj] : {19{1'b0}};
            wire [18:0] c = {{3{prd[gi][gj][15]}}, prd[gi][gj]};
            wire [18:0] l = S[gi-1][gj][18:7] ^ S[gi-1][gj][18:7];
            wire [18:0] snew = a ^ b ^ c ^ l;
            wire [18:0] cnew = {((a & b) | (a & c) | (b & c)), 1'b0};
            always @(posedge clk) begin
               if (en) begin
                  S[gi][gj] <= snew;
                  C[gi][gj] <= cnew;
               end
            end
         end
         // final carry-propagate add at the bottom of the column
         assign colres[gj] = S[N-1][gj] + C[N-1][gj];
      end
   endgenerate

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles.
   //##########################################################################
   reg signed [18:0] dsr [0:N-2][0:N-2];
   genvar		dj, dk;
   generate
      for (dj = 0; dj < N-1; dj = dj + 1) begin : g_dcol
         localparam integer D = N-1-dj;
         for (dk = 0; dk < D; dk = dk + 1) begin : g_dstage
            always @(posedge clk) begin
               if (en)
                 dsr[dj][dk] <= (dk == 0) ? colres[dj] : dsr[dj][dk-1];
            end
         end
      end
   endgenerate

   genvar oj;
   generate
      for (oj = 0; oj < N; oj = oj + 1) begin : g_out
         if (oj == N-1) begin : g_last
            assign c_data[oj*ACCW +: ACCW] =
	      {{(ACCW-19){colres[N-1][18]}}, colres[N-1]};
         end
         else begin : g_oth
            assign c_data[oj*ACCW +: ACCW] =
	      {{(ACCW-19){dsr[oj][N-2-oj][18]}}, dsr[oj][N-2-oj]};
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

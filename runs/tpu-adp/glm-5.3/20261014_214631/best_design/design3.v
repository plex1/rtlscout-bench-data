// Area-optimized weight-stationary systolic matmul tile (TPU MXU clone).
//
// Cycle-exact clone of the logikbench tpu, restructured for area:
//  * per-row activation skew and the east flow are one shift register chain
//    (row i needs stages 1..i+7; the tap for column j is stage i+j), which is
//    the same register count as the PE a_east pipeline but without duplicate
//    edge logic,
//  * the partial-sum chain is width-tapered: after k products the value needs
//    at most 16..19 bits (|sum| <= k*16384), so the 32-bit accumulators of the
//    original shrink to 16/17/17/18/18/18/18/19 bits down each column,
//  * the output de-skew and c_data are 19-bit, sign-extended at the pins,
//  * no resets on the datapath (it is fully self-flushing: the psum chain is
//    seeded from a constant zero and the activation chain is flushed by the
//    zeros injected while a_valid is low); only the valid pipe is reset.
// Products are still formed at exactly the original cycles with the weight
// value current at that cycle, so mid-drain weight reloads behave identically.
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
   // Weight storage: per-column south shift while w_valid (original semantics,
   // including partial shifts).
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
   // Activation skew + east flow. sr[i][k] is delay stage k+1 of array row i;
   // the multiplier input of column j is stage i+j (combinational input for
   // i=j=0). Row i needs stages 1..i+7.
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
   // Products: formed at exactly the original cycles with the weight current
   // at that cycle.
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
   // Partial-sum chain, width-tapered down each column.
   //##########################################################################
   wire signed [18:0] psw [0:N-1][0:N-1];
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_pscol
         for (gi = 0; gi < N; gi = gi + 1) begin : g_psrow
            localparam integer W =
			       (gi==0) ? 16 :
			       (gi==1) ? 17 :
			       (gi==2) ? 17 :
			       (gi==3) ? 18 :
			       (gi==4) ? 18 :
			       (gi==5) ? 18 :
			       (gi==6) ? 18 : 19;
            reg signed [W-1:0] ps;
            if (gi == 0) begin : g_r0
               always @(posedge clk) begin
                  if (en) ps <= prd[0][gj];
               end
            end
            else begin : g_rn
               wire signed [W-1:0] prev = psw[gi-1][gj][W-1:0];
               wire signed [W-1:0] pext = {{(W-16){prd[gi][gj][15]}}, prd[gi][gj]};
               always @(posedge clk) begin
                  if (en) ps <= prev + pext;
               end
            end
            assign psw[gi][gj] = {{(19-W){ps[W-1]}}, ps};
         end
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
                 dsr[dj][dk] <= (dk == 0) ? psw[N-1][dj] : dsr[dj][dk-1];
            end
         end
      end
   endgenerate

   genvar oj;
   generate
      for (oj = 0; oj < N; oj = oj + 1) begin : g_out
         if (oj == N-1) begin : g_last
            assign c_data[oj*ACCW +: ACCW] =
	      {{(ACCW-19){psw[N-1][N-1][18]}}, psw[N-1][N-1]};
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

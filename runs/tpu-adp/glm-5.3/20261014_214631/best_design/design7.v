// Optimized TPU MXU clone (cycle-exact vs logikbench tpu).
//
// Structure: per-cell Baugh-Wooley partial-product arrays are reduced
// TOGETHER with the column's carry-save accumulator state in one balanced
// 5-level 3:2 tree (10 rows -> 2), so no carry-propagate adder and no
// separate multiplier adder is on the critical path.  The Baugh-Wooley array
// sums to (product + 32512); eight products are accumulated per column, so a
// constant is subtracted at the final (off-path) adder (mod 2^19).
module cs3
  (
   input [18:0]	a,
   input [18:0]	b,
   input [18:0]	c,
   output [18:0] s,
   output [18:0] cy
   );
   assign s	 = a ^ b ^ c;
   assign cy	 = {((a & b) | (a & c) | (b & c)), 1'b0};
endmodule

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
   // Activation skew + east flow.
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
   // Carry-save partial-sum chain with fused multipliers.
   //##########################################################################
   reg [18:0] S [0:N-1][0:N-1];
   reg [18:0] C [0:N-1][0:N-1];
   wire [18:0] colres [0:N-1];

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
         for (gi = 0; gi < N; gi = gi + 1) begin : g_row
            // Baugh-Wooley rows for this cell: sum = a*w + 32512
            wire [7:0] aa = tap[gi][gj];
            wire [7:0] ww = w[gi][gj];
            wire [7:0]  rb0 = {~(aa[7] & ww[0]), aa[6:0] & {7{ww[0]}}};
            wire [7:0]  rb1 = {~(aa[7] & ww[1]), aa[6:0] & {7{ww[1]}}};
            wire [7:0]  rb2 = {~(aa[7] & ww[2]), aa[6:0] & {7{ww[2]}}};
            wire [7:0]  rb3 = {~(aa[7] & ww[3]), aa[6:0] & {7{ww[3]}}};
            wire [7:0]  rb4 = {~(aa[7] & ww[4]), aa[6:0] & {7{ww[4]}}};
            wire [7:0]  rb5 = {~(aa[7] & ww[5]), aa[6:0] & {7{ww[5]}}};
            wire [7:0]  rb6 = {~(aa[7] & ww[6]), aa[6:0] & {7{ww[6]}}};
            wire [7:0]  rb7 = {(aa[7] & ww[7]), ~(aa[6:0] & {7{ww[7]}})};
            wire [18:0] r0  = {11'b0, rb0};
            wire [18:0] r1  = {11'b0, rb1} << 1;
            wire [18:0] r2  = {11'b0, rb2} << 2;
            wire [18:0] r3  = {11'b0, rb3} << 3;
            wire [18:0] r4  = {11'b0, rb4} << 4;
            wire [18:0] r5  = {11'b0, rb5} << 5;
            wire [18:0] r6  = {11'b0, rb6} << 6;
            wire [18:0] r7  = {11'b0, rb7} << 7;

            if (gi == 0) begin : g_first
               // 8 rows -> 2 in 4 levels
               wire [18:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5;
               wire [18:0] nS, nC;
               cs3 u1 (.a(r0), .b(r1), .c(r2), .s(s1), .cy(c1));
               cs3 u2 (.a(r3), .b(r4), .c(r5), .s(s2), .cy(c2));
               cs3 u3 (.a(r6), .b(r7), .c(s1), .s(s3), .cy(c3));
               cs3 u4 (.a(c1), .b(s2), .c(c2), .s(s4), .cy(c4));
               cs3 u5 (.a(s3), .b(c3), .c(s4), .s(s5), .cy(c5));
               cs3 u6 (.a(s5), .b(c5), .c(c4), .s(nS), .cy(nC));
               always @(posedge clk) begin
                  if (en) begin
                     S[0][gj] <= nS;
                     C[0][gj] <= nC;
                  end
               end
            end
            else begin : g_rest
               // (S,C) + 8 rows -> 2 in 5 levels; S,C enter early
               wire [18:0] Sp = S[gi-1][gj];
               wire [18:0] Cp = C[gi-1][gj];
               wire [18:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5;
               wire [18:0] s6, c6, s7, c7;
               wire [18:0] nS, nC;
               cs3 v1 (.a(r0), .b(r1), .c(r2), .s(s1), .cy(c1));
               cs3 v2 (.a(r3), .b(r4), .c(r5), .s(s2), .cy(c2));
               cs3 v3 (.a(r6), .b(r7), .c(Sp), .s(s3), .cy(c3));
               cs3 v4 (.a(s1), .b(c1), .c(s2), .s(s4), .cy(c4));
               cs3 v5 (.a(c2), .b(s3), .c(c3), .s(s5), .cy(c5));
               cs3 v6 (.a(s4), .b(c4), .c(s5), .s(s6), .cy(c6));
               cs3 v7 (.a(s6), .b(c6), .c(c5), .s(s7), .cy(c7));
               cs3 v8 (.a(s7), .b(c7), .c(Cp), .s(nS), .cy(nC));
               always @(posedge clk) begin
                  if (en) begin
                     S[gi][gj] <= nS;
                     C[gi][gj] <= nC;
                  end
               end
            end
         end
      end
   endgenerate

   //##########################################################################
   // Final adder (off the critical path) with the Baugh-Wooley constant:
   // each product is off by +32512, 8 products per column -> subtract 260096
   // (mod 2^19: add 2^19 - 260096 = 264192).
   //##########################################################################
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_fin
         assign colres[gj] = S[N-1][gj] + C[N-1][gj] + 19'd264192;
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
   // Valid pipeline.
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

// Optimized TPU MXU clone (cycle-exact vs logikbench tpu).
//
// Key optimization vs a plain PE mesh: each cell's multiplier is written as an
// explicit Baugh-Wooley partial-product array reduced by 3:2 compressors to a
// CARRY-SAVE pair (no internal carry-propagate adder), and the column
// partial-sum chain compresses (S, C, ms, mc) with two more 3:2 levels.  So
// the multiplier's ripple-carry adder is removed from the critical path.
// The Baugh-Wooley array sums to (product + 32512); eight products are
// accumulated per column, so a constant 8*32512 is subtracted at the final
// adder (mod 2^19).
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

module bw_mul
  (
   input [7:0]	 a,
   input [7:0]	 w,
   output [18:0] ms,
   output [18:0] mc
   );
   // Baugh-Wooley rows: sum of rows = a*w + 32512
   wire [7:0]  rb0 = {~(a[7] & w[0]), a[6:0] & {7{w[0]}}};
   wire [7:0]  rb1 = {~(a[7] & w[1]), a[6:0] & {7{w[1]}}};
   wire [7:0]  rb2 = {~(a[7] & w[2]), a[6:0] & {7{w[2]}}};
   wire [7:0]  rb3 = {~(a[7] & w[3]), a[6:0] & {7{w[3]}}};
   wire [7:0]  rb4 = {~(a[7] & w[4]), a[6:0] & {7{w[4]}}};
   wire [7:0]  rb5 = {~(a[7] & w[5]), a[6:0] & {7{w[5]}}};
   wire [7:0]  rb6 = {~(a[7] & w[6]), a[6:0] & {7{w[6]}}};
   wire [7:0]  rb7 = {(a[7] & w[7]), ~(a[6:0] & {7{w[7]}})};
   wire [18:0] r0  = {11'b0, rb0};
   wire [18:0] r1  = {11'b0, rb1} << 1;
   wire [18:0] r2  = {11'b0, rb2} << 2;
   wire [18:0] r3  = {11'b0, rb3} << 3;
   wire [18:0] r4  = {11'b0, rb4} << 4;
   wire [18:0] r5  = {11'b0, rb5} << 5;
   wire [18:0] r6  = {11'b0, rb6} << 6;
   wire [18:0] r7  = {11'b0, rb7} << 7;
   // balanced reduction: 8 rows -> 2 rows in 4 levels
   wire [18:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5;
   cs3 u1 (.a(r0), .b(r1), .c(r2), .s(s1), .cy(c1));
   cs3 u2 (.a(r3), .b(r4), .c(r5), .s(s2), .cy(c2));
   cs3 u3 (.a(r6), .b(r7), .c(s1), .s(s3), .cy(c3));
   cs3 u4 (.a(c1), .b(s2), .c(c2), .s(s4), .cy(c4));
   cs3 u5 (.a(s3), .b(c3), .c(s4), .s(s5), .cy(c5));
   cs3 u6 (.a(s5), .b(c5), .c(c4), .s(ms), .cy(mc));
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
   // Carry-save products (Baugh-Wooley, no internal CPA).
   //##########################################################################
   wire [18:0] ms [0:N-1][0:N-1];
   wire [18:0] mc [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_mrow
         for (gj = 0; gj < N; gj = gj + 1) begin : g_mcol
            bw_mul u_mul (.a(tap[gi][gj]), .w(w[gi][gj]),
                          .ms(ms[gi][gj]), .mc(mc[gi][gj]));
         end
      end
   endgenerate

   //##########################################################################
   // Carry-save partial-sum chain down each column (mod 2^19).
   //##########################################################################
   reg [18:0] S [0:N-1][0:N-1];
   reg [18:0] C [0:N-1][0:N-1];
   wire [18:0] colres [0:N-1];

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
         always @(posedge clk) begin
            if (en) begin
               S[0][gj] <= ms[0][gj];
               C[0][gj] <= mc[0][gj];
            end
         end
         for (gi = 1; gi < N; gi = gi + 1) begin : g_row
            wire [18:0] t1, u1;
            cs3 u_c1 (.a(S[gi-1][gj]), .b(C[gi-1][gj]), .c(ms[gi][gj]),
                      .s(t1), .cy(u1));
            wire [18:0] t2, u2;
            cs3 u_c2 (.a(t1), .b(u1), .c(mc[gi][gj]), .s(t2), .cy(u2));
            always @(posedge clk) begin
               if (en) begin
                  S[gi][gj] <= t2;
                  C[gi][gj] <= u2;
               end
            end
         end
         // final CPA (off the critical path) with the Baugh-Wooley constant:
         // 8 products per column sum, each off by +32512 -> subtract 260096.
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

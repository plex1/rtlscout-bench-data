// Area-optimized weight-stationary systolic matmul tile (TPU MXU clone).
//
// Same cycle behavior as the logikbench tpu but the internal structure is
// levelized: activations are captured directly in the PE register row above
// the weight they meet (no separate skew SRs), each column keeps its own
// accumulator, and the weight load path is a broadcast (no column shift
// chains). Only the small output de-skew (max 7 stages of 32b) and the valid
// pipeline are retained, together with a 3-cycle weight load enable counter.
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
   // Weight storage. Broadcast load: on the k-th w_valid cycle the bottom
   // matrix row B[7-k][*] is presented, so with a 3..0 down-counter the value
   // written into row r is w_data of cycle 7-r, exactly as if the rows had
   // been shifted south for 8 cycles.
   //##########################################################################
   reg [DW-1:0] w [0:N-1][0:N-1];

   reg [2:0] wcnt;   // 8 load cycles -> counts 7..0 during w_valid
   always @(posedge clk or posedge rst) begin
      if (rst) begin
         wcnt <= 3'd7;
      end
      else begin
         if (w_valid)
           wcnt <= (wcnt == 3'd0) ? 3'd7 : wcnt - 3'd1;
         else
           wcnt <= 3'd7;   // rearm so a later reload starts from row 7
      end
   end

   integer wi, wj;
   always @(posedge clk) begin
      if (w_valid) begin
         for (wj = 0; wj < N; wj = wj + 1)
                   w[wcnt][wj] <= w_data[wj*DW +: DW];
      end
   end
   // (weights need no reset: they are always fully loaded before compute)

   //##########################################################################
   // Activation registers. a[i] is the value entering column 0 of row i on
   // THIS cycle (i.e. it meets the weight of row i this cycle); the original
   // injects row i i cycles late, so this slot carries A[m][i]. The array
   // input itself is that value delayed one more cycle, which is exactly the
   // original PE a_east pipeline register.
   //##########################################################################
   reg [DW-1:0] a_in [0:N-1];      // value meeting row i's weight this cycle
   reg [DW-1:0] a_arr [0:N-1];     // arrival registers (original column 0)
   integer      ai;
   always @(posedge clk) begin
      if (en) begin
         for (ai = 0; ai < N; ai = ai + 1) begin
            a_in[ai]  <= a_valid ? a_data[ai*DW +: DW] : {DW{1'b0}};
            a_arr[ai] <= a_in[ai];
         end
      end
   end

   // activation entering column j of row i during compute (combinational)
   wire [DW-1:0] ain [0:N-1][0:N-1];
   genvar	       gi, gj;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ain
         assign ain[gi][0] = a_arr[gi];
         for (gj = 1; gj < N; gj = gj + 1) begin : g_ainc
            assign ain[gi][gj] = (gj == gi + 1) ? a_in[gi] : {DW{1'b0}};
         end
      end
   endgenerate

   //##########################################################################
   // Column accumulators: col[j] += sum_i A[m][i]*B[i][j] for the row being
   // consumed this cycle. Pairs add nothing, so only the diagonal terms feed
   // the adder tree. A zero row (a_valid low on that arrival) leaves the
   // accumulator unchanged, matching idle/drain cycles.
   //##########################################################################
   reg signed [ACCW-1:0] col [0:N-1];
   wire signed [ACCW-1:0] nxt [0:N-1];

   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_col
         wire signed [2*DW-1:0] pr [0:N-1];
         for (gi = 0; gi < N; gi = gi + 1) begin : g_pr
            assign pr[gi] = $signed(ain[gi][gj]) * $signed(w[gi][gj]);
         end
         wire signed [ACCW-1:0] psum =
                   $signed({{(ACCW-2*DW){pr[0][2*DW-1]}}, pr[0]})
                 + $signed({{(ACCW-2*DW){pr[1][2*DW-1]}}, pr[1]})
                 + $signed({{(ACCW-2*DW){pr[2][2*DW-1]}}, pr[2]})
                 + $signed({{(ACCW-2*DW){pr[3][2*DW-1]}}, pr[3]})
                 + $signed({{(ACCW-2*DW){pr[4][2*DW-1]}}, pr[4]})
                 + $signed({{(ACCW-2*DW){pr[5][2*DW-1]}}, pr[5]})
                 + $signed({{(ACCW-2*DW){pr[6][2*DW-1]}}, pr[6]})
                 + $signed({{(ACCW-2*DW){pr[7][2*DW-1]}}, pr[7]});
         assign nxt[gj] = col[gj] + psum;
      end
   endgenerate

   integer cj;
   always @(posedge clk) begin
      if (rst) begin
         for (cj = 0; cj < N; cj = cj + 1)
                   col[cj] <= {ACCW{1'b0}};
      end
      else if (en)
        for (cj = 0; cj < N; cj = cj + 1)
                   col[cj] <= nxt[cj];
   end

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles so all columns
   // of a result row leave together. Column N-1 (j=7) needs no delay, so
   // register only columns 0..6.
   //##########################################################################
   reg signed [ACCW-1:0] dr [0:N-2][0:N-2];   // dr[j][k], j=0..N-2
   genvar	       dj, dk;
   generate
      for (dj = 0; dj < N-1; dj = dj + 1) begin : g_deskew
         localparam D = N - 1 - dj;
         if (D == 1) begin : g_d1
            always @(posedge clk) begin
               if (rst) dr[dj][0] <= {ACCW{1'b0}};
               else if (en) dr[dj][0] <= col[dj];
            end
         end
         else begin : g_dn
            for (dk = 0; dk < D; dk = dk + 1) begin : g_dsr
               always @(posedge clk) begin
                  if (rst) dr[dj][dk] <= {ACCW{1'b0}};
                  else if (en) dr[dj][dk] <= (dk == 0) ? col[dj] : dr[dj][dk-1];
               end
            end
         end
      end
   endgenerate

   genvar oj;
   generate
      for (oj = 0; oj < N; oj = oj + 1) begin : g_out
         if (oj == N-1) begin : g_o7
            assign c_data[oj*ACCW +: ACCW] = col[N-1];
         end
         else begin : g_on
            assign c_data[oj*ACCW +: ACCW] = dr[oj][N-2-oj];
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
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

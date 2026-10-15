// Area-optimized TPU MXU clone.
//
// Instead of 64 PEs with 32-bit partial-sum chains + input skew + output
// de-skew, the activation row is simply delayed 15 cycles in an 8-bit-wide
// delay line (the narrowest place to absorb the mandatory latency) and the
// whole dot product per column is computed in parallel with an adder tree.
// Weights keep the exact shift-south load semantics of the original.
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
   // Weight storage: shift south on every w_valid cycle (identical semantics
   // to the original per-column shift register chain).
   //##########################################################################
   reg [DW-1:0] w [0:N-1][0:N-1];
   integer       wi, wj;
   always @(posedge clk or posedge rst) begin
      if (rst) begin
         for (wj = 0; wj < N; wj = wj + 1)
           for (wi = 0; wi < N; wi = wi + 1)
             w[wi][wj] <= {DW{1'b0}};
      end
      else if (w_valid) begin
         for (wj = 0; wj < N; wj = wj + 1) begin
            w[0][wj] <= w_data[wj*DW +: DW];
            for (wi = 1; wi < N; wi = wi + 1)
              w[wi][wj] <= w[wi-1][wj];
         end
      end
   end

   //##########################################################################
   // Activation delay: row m is captured at the end of cycle m (zero when
   // a_valid is low) and held for 15 cycles so that it is present at the
   // multiplier inputs during cycle m+15, exactly when its result must be
   // presented on c_data.
   //##########################################################################
   reg [DW-1:0] dly [0:N-1][0:2*N-2];
   integer       di, dk;
   always @(posedge clk) begin
      if (en) begin
         for (di = 0; di < N; di = di + 1) begin
            dly[di][0] <= a_valid ? a_data[di*DW +: DW] : {DW{1'b0}};
            for (dk = 1; dk < 2*N-1; dk = dk + 1)
              dly[di][dk] <= dly[di][dk-1];
         end
      end
   end

   wire [DW-1:0] e [0:N-1];
   genvar	       gi;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_e
         assign e[gi] = dly[gi][2*N-2];
      end
   endgenerate

   //##########################################################################
   // 64 fixed multipliers + one 8-input adder tree per column.
   //##########################################################################
   wire signed [2*DW-1:0] p [0:N-1][0:N-1];
   genvar		       pi, pj;
   generate
      for (pi = 0; pi < N; pi = pi + 1) begin : g_prow
         for (pj = 0; pj < N; pj = pj + 1) begin : g_pcol
            assign p[pi][pj] = $signed(e[pi]) * $signed(w[pi][pj]);
         end
      end
   endgenerate

   wire signed [19:0] t3 [0:N-1];
   generate
      for (pj = 0; pj < N; pj = pj + 1) begin : g_tree
         wire signed [19:0] t1 [0:3];
         wire signed [19:0] t2 [0:1];
         assign t1[0] = p[0][pj] + p[1][pj];
         assign t1[1] = p[2][pj] + p[3][pj];
         assign t1[2] = p[4][pj] + p[5][pj];
         assign t1[3] = p[6][pj] + p[7][pj];
         assign t2[0] = t1[0] + t1[1];
         assign t2[1] = t1[2] + t1[3];
         assign t3[pj] = t2[0] + t2[1];
         assign c_data[pj*ACCW +: ACCW] = {{(ACCW-20){t3[pj][19]}}, t3[pj]};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed 15 compute cycles.
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

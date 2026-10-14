//#############################################################################
// TPU MXU - Non-systolic pipelined design
// All 64 multipliers in parallel, pipelined adder tree, 15-cycle latency
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;  // 15

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Weight storage: 8x8 shift register array, loads from top (shift south)
   // After 8 load cycles, B_reg[i][j] = B[i][j]
   //##########################################################################
   reg [DW-1:0] B_reg [0:N-1][0:N-1];

   integer ii, jj;
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             B_reg[ii][jj] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (jj = 0; jj < N; jj = jj + 1) begin
            B_reg[0][jj] <= w_data[jj*DW +: DW];
            for (ii = 1; ii < N; ii = ii + 1)
              B_reg[ii][jj] <= B_reg[ii-1][jj];
         end
      end
   end

   //##########################################################################
   // Pipeline stage 0: capture activation row (inject zeros when !a_valid)
   //##########################################################################
   reg [DW-1:0] A_cap [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N; ii = ii + 1)
           A_cap[ii] <= {DW{1'b0}};
      end
      else if (en) begin
         for (ii = 0; ii < N; ii = ii + 1)
           A_cap[ii] <= a_valid ? a_data[ii*DW +: DW] : {DW{1'b0}};
      end
   end

   //##########################################################################
   // Pipeline stage 1: multiply (64 multipliers, 16-bit products)
   //##########################################################################
   reg [2*DW-1:0] prod [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             prod[ii][jj] <= {2*DW{1'b0}};
      end
      else if (en) begin
         for (ii = 0; ii < N; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             prod[ii][jj] <= $signed(A_cap[ii]) * $signed(B_reg[ii][jj]);
      end
   end

   //##########################################################################
   // Pipeline stage 2: add tree level 1 (pairs of products -> 17-bit)
   // 4 pairs per column, 8 columns = 32 sums
   //##########################################################################
   reg [2*DW:0] sum1 [0:N/2-1][0:N-1];  // 17-bit
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N/2; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             sum1[ii][jj] <= {2*DW+1{1'b0}};
      end
      else if (en) begin
         for (ii = 0; ii < N/2; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             sum1[ii][jj] <= $signed(prod[2*ii][jj]) + $signed(prod[2*ii+1][jj]);
      end
   end

   //##########################################################################
   // Pipeline stage 3: add tree level 2 (pairs of sum1 -> 18-bit)
   // 2 pairs per column, 8 columns = 16 sums
   //##########################################################################
   reg [2*DW+1:0] sum2 [0:N/4-1][0:N-1];  // 18-bit
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N/4; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             sum2[ii][jj] <= {2*DW+2{1'b0}};
      end
      else if (en) begin
         for (ii = 0; ii < N/4; ii = ii + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             sum2[ii][jj] <= $signed(sum1[2*ii][jj]) + $signed(sum1[2*ii+1][jj]);
      end
   end

   //##########################################################################
   // Pipeline stage 4: add tree level 3 (pair of sum2 -> 19-bit)
   // 1 pair per column, 8 columns = 8 sums
   //##########################################################################
   reg [2*DW+2:0] sum3 [0:N-1];  // 19-bit (actually 18-bit, but let's use 19 for safety)
   always @(posedge clk) begin
      if (rst) begin
         for (jj = 0; jj < N; jj = jj + 1)
           sum3[jj] <= {2*DW+3{1'b0}};
      end
      else if (en) begin
         for (jj = 0; jj < N; jj = jj + 1)
           sum3[jj] <= $signed(sum2[0][jj]) + $signed(sum2[1][jj]);
      end
   end

   //##########################################################################
   // Pipeline stages 5-14: delay (10 stages of 8 × 19-bit)
   //##########################################################################
   // We need 10 more stages to reach 15 total
   localparam PACCW = 2*DW + 3;  // 19
   reg [PACCW-1:0] delay_pipe [0:9][0:N-1];
   integer s;
   always @(posedge clk) begin
      if (rst) begin
         for (s = 0; s <= 9; s = s + 1)
           for (jj = 0; jj < N; jj = jj + 1)
             delay_pipe[s][jj] <= {PACCW{1'b0}};
      end
      else if (en) begin
         for (jj = 0; jj < N; jj = jj + 1) begin
            delay_pipe[0][jj] <= sum3[jj];
            for (s = 1; s <= 9; s = s + 1)
              delay_pipe[s][jj] <= delay_pipe[s-1][jj];
         end
      end
   end

   //##########################################################################
   // Output: sign-extend 19-bit to 32-bit
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
         wire [PACCW-1:0] result = delay_pipe[9][j];
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-PACCW){result[PACCW-1]}}, result};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT cycles
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

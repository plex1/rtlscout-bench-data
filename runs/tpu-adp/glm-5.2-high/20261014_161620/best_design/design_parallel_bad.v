//#############################################################################
// TPU MXU - Pipelined parallel multiply with CSA accumulation
// Stage 0: capture A row (with skew), Stage 1: multiply, 
// Stages 2-8: CSA accumulate, Stages 9-14: delay
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
   localparam PACCW = 19;

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Weight storage: 8x8 shift register, loads from top (shift south)
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
   // Input skew: row i of A is delayed i cycles
   //##########################################################################
   reg [DW-1:0] a_skew [0:N-1];
   reg [DW-1:0] skew_pipe [0:N-2][0:N-1];  // 7 stages, 8 rows each

   integer m, s;
   always @(posedge clk) begin
      if (rst) begin
         for (ii = 0; ii < N; ii = ii + 1)
           a_skew[ii] <= {DW{1'b0}};
         for (s = 0; s < N-1; s = s + 1)
           for (ii = 0; ii < N; ii = ii + 1)
             skew_pipe[s][ii] <= {DW{1'b0}};
      end
      else if (en) begin
         // Row 0: direct input
         a_skew[0] <= a_valid ? a_data[0*DW +: DW] : {DW{1'b0}};
         // Row i: delayed i cycles through skew_pipe
         skew_pipe[0][1] <= a_valid ? a_data[1*DW +: DW] : {DW{1'b0}};
         a_skew[1] <= skew_pipe[0][1];
         for (ii = 2; ii < N; ii = ii + 1) begin
            skew_pipe[0][ii] <= a_valid ? a_data[ii*DW +: DW] : {DW{1'b0}};
         end
         // Shift all skew pipes
         for (s = 1; s < N-1; s = s + 1) begin
            for (ii = s+1; ii < N; ii = ii + 1) begin
               skew_pipe[s][ii] <= skew_pipe[s-1][ii];
            end
         end
         // Tap outputs
         for (ii = 2; ii < N; ii = ii + 1) begin
            a_skew[ii] <= skew_pipe[ii-1][ii];
         end
      end
   end

   // This skew implementation is getting complicated. Let me go back to the
   // original generate-based skew, which is known to work.
   
   // Actually, let me just use the parallel multiply approach with the
   // original skew and CSA add tree.

endmodule

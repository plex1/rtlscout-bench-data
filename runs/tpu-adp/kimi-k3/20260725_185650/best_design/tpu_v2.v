//#############################################################################
// Optimized weight-stationary systolic matrix-multiply tile (8x8).
// Cycle-accurate drop-in for the logikbench tpu benchmark.
//
// Optimizations vs. the original RTL:
//   - Partial sums carried in PSW = 2*DW + clog2(N) = 19 bits (the exact
//     worst-case range of a sum of N signed 2*DW-bit products); sign-extended
//     to ACCW only at the output. Shrinks 64 psum regs, 28 deskew regs and
//     all 64 accumulators.
//   - Resets removed from all datapath/shift registers (only the valid pipe
//     keeps one). Any power-on garbage flushes out of the array and the
//     skew/de-skew chains well before c_valid can first rise (LAT cycles).
//   - East activation registers of the last column (unused) are omitted.
//#############################################################################

module tpu
  #(parameter N = 8,               // array dimension (N x N)
    parameter DW = 8,              // operand width (signed)
    parameter ACCW = 32)           // result width (signed)
   (
    input               clk,
    input               rst,       // synchronous, active high
    input               w_valid,   // weight-row valid (load phase)
    input [N*DW-1:0]    w_data,    // weight row, packed { col N-1 ... col 0 }
    input               a_valid,   // activation-row valid (compute phase)
    input [N*DW-1:0]    a_data,    // activation row, packed { col N-1 ... col 0 }
    output              c_valid,   // result-row valid
    output [N*ACCW-1:0] c_data     // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;             // injection-to-output latency
   localparam PSW = 2*DW + $clog2(N);    // exact partial-sum width (19)

   wire en = ~w_valid;                   // array advances only when not loading

   //##########################################################################
   // State
   //##########################################################################
   reg [DW-1:0]  w_reg [0:N-1][0:N-1];   // stationary weights
   reg [DW-1:0]  a_reg [0:N-1][0:N-2];   // east activation taps (col N-1 unused)
   reg [PSW-1:0] p_reg [0:N-1][0:N-1];   // south partial sums

   genvar i, j;

   //##########################################################################
   // Input skew: row i of A delayed i cycles; zeros injected while !a_valid.
   //##########################################################################
   wire [DW-1:0] a_left [0:N-1];
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
                       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // Weight load: per-column shift registers, load from the top.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_w_row
         for (j = 0; j < N; j = j + 1) begin : g_w_col
            wire [DW-1:0] w_north =
                          (i == 0) ? w_data[j*DW +: DW] : w_reg[i-1][j];
            always @(posedge clk) begin
               if (w_valid)
                 w_reg[i][j] <= w_north;
            end
         end
      end
   endgenerate

   //##########################################################################
   // The systolic mesh: activation east, partial sum south (19-bit).
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [DW-1:0] a_west =
                          (j == 0) ? a_left[i] : a_reg[i][j-1];
            wire signed [2*DW-1:0] prod =
                          $signed(a_west) * $signed(w_reg[i][j]);
            wire [PSW-1:0] psum_north =
                          (i == 0) ? {PSW{1'b0}} : p_reg[i-1][j];
            always @(posedge clk) begin
               if (en) begin
                  if (j < N-1)
                    a_reg[i][j] <= a_west;
                  p_reg[i][j] <= psum_north +
                                 {{(PSW-2*DW){prod[2*DW-1]}}, prod};
               end
            end
         end
      end
   endgenerate

   //##########################################################################
   // Output de-skew (19-bit) + sign extension to ACCW.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSW-1:0] cin = p_reg[N-1][j];
         wire [PSW-1:0] cout;
         if (D == 0) begin : g_deskew0
            assign cout = cin;
         end
         else begin : g_deskewn
            reg [PSW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign cout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSW){cout[PSW-1]}}, cout};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
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

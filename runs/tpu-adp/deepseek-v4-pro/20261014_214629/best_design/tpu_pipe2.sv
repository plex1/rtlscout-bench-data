module tpu
  #(parameter N = 8, DW = 8, ACCW = 32)
   (
    input  logic clk,
    input  logic rst,
    input  logic w_valid,
    input  logic [N*DW-1:0]   w_data,
    input  logic a_valid,
    input  logic [N*DW-1:0]   a_data,
    output logic c_valid,
    output logic [N*ACCW-1:0] c_data
    );

   localparam LAT = 2*N - 1;

   // Weight matrix register file, same load/reload/stall behavior.
   logic [DW-1:0] Breg [0:N-1][0:N-1];
   integer bi, bj;
   always_ff @(posedge clk) begin
      if (rst) begin
         for (bi = 0; bi < N; bi = bi + 1)
           for (bj = 0; bj < N; bj = bj + 1)
             Breg[bi][bj] <= '0;
      end else if (w_valid) begin
         for (bj = 0; bj < N; bj = bj + 1)
           Breg[0][bj] <= w_data[bj*DW +: DW];
         for (bi = 1; bi < N; bi = bi + 1)
           for (bj = 0; bj < N; bj = bj + 1)
             Breg[bi][bj] <= Breg[bi-1][bj];
      end
   end

   // Activation row delay line. Depth 16: stage 0 is the current input.
   logic [N*DW-1:0] adly [0:LAT];
   integer ai;
   always_ff @(posedge clk) begin
      if (rst) begin
         for (ai = 0; ai <= LAT; ai = ai + 1)
           adly[ai] <= '0;
      end else if (!w_valid) begin
         adly[0] <= a_valid ? a_data : {N*DW{1'b0}};
         for (ai = 1; ai <= LAT; ai = ai + 1)
           adly[ai] <= adly[ai-1];
      end
   end

   // products for the eight A rows currently in the pipeline.
   // The row i partial products are built from adly[i+1], i.e. with the
   // systolic arrival skew (activation row r is at stage r+1 after r en-cycles).
   wire signed [15:0] prod_i [0:N-1][0:N-1];
   wire signed [ACCW-1:0] partial_i [0:N-1][0:N-1];
   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            assign prod_i[i][j] = $signed(adly[i+1][i*DW +: DW])
                                * $signed(Breg[i][j]);
            assign partial_i[i][j] =
                {{(ACCW-2*DW){prod_i[i][j][15]}}, prod_i[i][j]};
         end
      end

      // Add partial products along the A-row dimension. This is intentionally
      // left as a simple expression so Yosys sees a wide carry-save tree.
      for (j = 0; j < N; j = j + 1) begin : g_out
         assign c_data[j*ACCW +: ACCW] =
             partial_i[0][j] + partial_i[1][j] + partial_i[2][j] +
             partial_i[3][j] + partial_i[4][j] + partial_i[5][j] +
             partial_i[6][j] + partial_i[7][j];
      end
   endgenerate

   // Result row r is valid LAT cycles after its input row, still in slot LAT.
   logic [LAT:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (!w_valid)
        vpipe <= {vpipe[LAT-1:0], a_valid};
   end
   assign c_valid = vpipe[LAT];
endmodule
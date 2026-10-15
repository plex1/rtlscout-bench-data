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

   // Weight matrix B[i][j], i = row, j = column. Load shifts south:
   // bottom matrix row first, so B[0] ends up in Breg[0] after N cycles.
   logic [DW-1:0] Breg [0:N-1][0:N-1];

   integer bi, bj;
   always_ff @(posedge clk) begin
      if (rst) begin
         for (bi = 0; bi < N; bi = bi + 1) begin
            for (bj = 0; bj < N; bj = bj + 1)
              Breg[bi][bj] <= '0;
         end
      end else if (w_valid) begin
         for (bj = 0; bj < N; bj = bj + 1)
           Breg[0][bj] <= w_data[bj*DW +: DW];
         for (bi = 1; bi < N; bi = bi + 1)
           for (bj = 0; bj < N; bj = bj + 1)
             Breg[bi][bj] <= Breg[bi-1][bj];
      end
   end

   // Activation row shift register: depth LAT (15). Shifts only while ~w_valid,
   // and zero is injected when a_valid is low (same as the systolic datapath).
   logic [N*DW-1:0] adly [0:LAT-1];
   integer ai;
   always_ff @(posedge clk) begin
      if (rst) begin
         for (ai = 0; ai < LAT; ai = ai + 1)
           adly[ai] <= '0;
      end else if (!w_valid) begin
         adly[0] <= a_valid ? a_data : {N*DW{1'b0}};
         for (ai = 1; ai < LAT; ai = ai + 1)
           adly[ai] <= adly[ai-1];
      end
   end

   // Compute C = A * B from the delayed activation row. Balanced addition tree.
   wire signed [15:0] prod [0:7][0:7];
   wire signed [ACCW-1:0] sum01 [0:3][0:7];
   wire signed [ACCW-1:0] sum23 [0:1][0:7];
   wire signed [ACCW-1:0] dot   [0:7];
   genvar j, i;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col_prod
         for (i = 0; i < N; i = i + 1) begin : g_prod_i
            assign prod[i][j] = $signed(adly[N-1+i][i*DW +: DW])
                              * $signed(Breg[i][j]);
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_col_s01
         for (i = 0; i < N/2; i = i + 1) begin : g_s01_i
            assign sum01[i][j] =
                $signed({{(ACCW-2*DW){prod[2*i][j][15]}}, prod[2*i][j]})
              + $signed({{(ACCW-2*DW){prod[2*i+1][j][15]}}, prod[2*i+1][j]});
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_col_s23
         for (i = 0; i < N/4; i = i + 1) begin : g_s23_i
            assign sum23[i][j] = sum01[2*i][j] + sum01[2*i+1][j];
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_col_out
         assign dot[j] = sum23[0][j] + sum23[1][j];
         assign c_data[j*ACCW +: ACCW] = dot[j];
      end
   endgenerate

   // Valid pipeline: a_valid delayed 15 en-cycles.
   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (!w_valid)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];
endmodule
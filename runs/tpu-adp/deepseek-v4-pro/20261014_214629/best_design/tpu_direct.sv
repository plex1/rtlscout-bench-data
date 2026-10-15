module tpu
  #(parameter N = 8, DW = 8, ACCW_TOP = 32, DOTW = 19)
   (
    input  logic clk,
    input  logic rst,
    input  logic w_valid,
    input  logic [N*DW-1:0]   w_data,
    input  logic a_valid,
    input  logic [N*DW-1:0]   a_data,
    output logic c_valid,
    output logic [N*ACCW_TOP-1:0] c_data
    );

   localparam LAT = 2*N - 1;

   //////////////////////////////////////////////////////////////////////////
   // Stationary weight matrix. Load/reload/stall behavior identical to the
   // systolic array.
   //////////////////////////////////////////////////////////////////////////
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

   //////////////////////////////////////////////////////////////////////////
   // Dot-product matrix: c_row[j] = sum_i a_row[i] * B[i][j].
   // Zeros are injected when a_valid is low, exactly as in the reference.
   //////////////////////////////////////////////////////////////////////////
   wire signed [DW-1:0] a_cur [0:N-1];
   wire signed [15:0]   prod  [0:N-1][0:N-1];
   wire signed [DOTW-1:0] psum01 [0:N/2-1][0:N-1];
   wire signed [DOTW-1:0] psum23 [0:N/4-1][0:N-1];
   wire signed [DOTW-1:0] dot    [0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_a
         assign a_cur[i] = a_valid ? $signed(a_data[i*DW +: DW])
                                   : $signed({DW{1'b0}});
      end

      for (i = 0; i < N; i = i + 1) begin : g_products
         for (j = 0; j < N; j = j + 1) begin : g_col
            assign prod[i][j] = a_cur[i] * $signed(Breg[i][j]);
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_l1
         for (i = 0; i < N/2; i = i + 1) begin : g_pairs
            assign psum01[i][j] =
                $signed({{(DOTW-2*DW){prod[2*i][j][15]}}, prod[2*i][j]})
              + $signed({{(DOTW-2*DW){prod[2*i+1][j][15]}}, prod[2*i+1][j]});
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_l2
         for (i = 0; i < N/4; i = i + 1) begin : g_pairs2
            assign psum23[i][j] = psum01[2*i][j] + psum01[2*i+1][j];
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_l3
         assign dot[j] = psum23[0][j] + psum23[1][j];
      end
   endgenerate

   //////////////////////////////////////////////////////////////////////////
   // Output de-skew / latency pipeline. Dot rows are captured and advanced
   // only on compute cycles (~w_valid), matching the systolic stall behavior.
   //////////////////////////////////////////////////////////////////////////
   logic [DOTW-1:0] pipe [0:LAT-1][0:N-1];
   integer aip, cj;
   always_ff @(posedge clk) begin
      if (rst) begin
         for (aip = 0; aip < LAT; aip = aip + 1)
           for (cj = 0; cj < N; cj = cj + 1)
             pipe[aip][cj] <= '0;
      end else if (!w_valid) begin
         for (cj = 0; cj < N; cj = cj + 1)
           pipe[0][cj] <= dot[cj];
         for (aip = 1; aip < LAT; aip = aip + 1)
           for (cj = 0; cj < N; cj = cj + 1)
             pipe[aip][cj] <= pipe[aip-1][cj];
      end
   end

   genvar oj;
   generate
      for (oj = 0; oj < N; oj = oj + 1) begin : g_out
         assign c_data[oj*ACCW_TOP +: ACCW_TOP] =
             {{(ACCW_TOP-DOTW){pipe[LAT-1][oj][DOTW-1]}}, pipe[LAT-1][oj]};
      end
   endgenerate

   //////////////////////////////////////////////////////////////////////////
   // Valid pipeline.
   //////////////////////////////////////////////////////////////////////////
   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (!w_valid)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];
endmodule
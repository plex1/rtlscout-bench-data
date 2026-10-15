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

   localparam LAT = 2*N - 1;   // 15
   localparam PIPE_W = LAT*N*2*DW;  // 15*8*16 = 1920 bits

   //////////////////////////////////////////////////////////////////////////
   // Static weight matrix. Load/reload/stall behavior identical to the
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
   // Activation row delay line, depth 15. Stage t is the input seen t compute
   // cycles ago (zero-padded when a_valid was low).
   //////////////////////////////////////////////////////////////////////////
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

   //////////////////////////////////////////////////////////////////////////
   // Deferred scalar products. Each delayed A-column byte is multiplied on
   // first use by its stationary B value, then the 16-bit product is pushed
   // through the remaining stages. The result row c_row[j] =
   // sum_i partial_i_j(t = N-1-i + (LAT-1)) = sum_i B[i][j]*A_row[i]
   // is valid together with vpipe[LAT-1].
   //////////////////////////////////////////////////////////////////////////
   logic [PIPE_W-1:0] pipe;

   always_ff @(posedge clk) begin
      if (rst) begin
         pipe <= '0;
      end else if (!w_valid) begin
         for (int s = 0; s < LAT; s = s + 1)
           for (int j = 0; j < N; j = j + 1)
             pipe[(s*N + j)*2*DW +: 2*DW] <=
                 (s < N-1) ? pipe[((s+1)*N + j)*2*DW +: 2*DW] :
                             ($signed(adly[s][j*DW +: DW]) * $signed(Breg[j][j]));
      end
   end

   generate
      for (genvar j = 0; j < N; j = j + 1) begin : g_out
         assign c_data[j*ACCW_TOP +: ACCW_TOP] =
             {{(ACCW_TOP-DOTW){pipe[(LAT-1)*N*2*DW + j*2*DW + 15]}}, 
              pipe[(LAT-1)*N*2*DW + j*2*DW +: 2*DW]};
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
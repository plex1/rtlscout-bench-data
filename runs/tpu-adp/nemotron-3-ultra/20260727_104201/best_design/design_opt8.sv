//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU)
// Optimized with 16-bit group carry-lookahead adder (2 groups for 32-bit)
//#############################################################################

module tpu_pe
  #(parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // accumulator / partial-sum width (signed)
   (
    input  logic              clk,
    input  logic              rst,         // synchronous, active high
    input  logic              ld,          // weight-load shift enable
    input  logic              en,          // compute enable
    input  logic [DW-1:0]     w_north,     // weight in from PE above
    input  logic [DW-1:0]     a_west,      // activation in from west
    input  logic [ACCW-1:0]   psum_north,  // partial sum in from above
    output logic [DW-1:0]     w_south,     // stored weight out to PE below
    output logic [DW-1:0]     a_east,      // activation out to east
    output logic [ACCW-1:0]   psum_south   // partial sum out to below
    );

   // stationary weight register
   logic [DW-1:0] w;
   assign w_south = w;

   // local product: signed multiply
   logic [2*DW-1:0] prod;
   assign prod = $signed(a_west) * $signed(w);

   // Sign-extended product
   logic [ACCW-1:0] prod_sext;
   assign prod_sext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // 16-bit group carry-lookahead adder (2 groups for 32-bit)
   logic [ACCW-1:0] sum;
   logic [1:0] group_carry;
   assign group_carry[0] = 1'b0;
   
   // Group 0 (bits 0-15)
   logic [15:0] g0, p0;
   assign g0 = psum_north[15:0] & prod_sext[15:0];
   assign p0 = psum_north[15:0] ^ prod_sext[15:0];
   
   logic [3:0] gg0_4, gp0_4;
   assign gg0_4[0] = g0[3] | (p0[3] & (g0[2] | (p0[2] & (g0[1] | (p0[1] & g0[0])))));
   assign gp0_4[0] = &p0[3:0];
   assign gg0_4[1] = g0[7] | (p0[7] & (g0[6] | (p0[6] & (g0[5] | (p0[5] & g0[4])))));
   assign gp0_4[1] = &p0[7:4];
   assign gg0_4[2] = g0[11] | (p0[11] & (g0[10] | (p0[10] & (g0[9] | (p0[9] & g0[8])))));
   assign gp0_4[2] = &p0[11:8];
   assign gg0_4[3] = g0[15] | (p0[15] & (g0[14] | (p0[14] & (g0[13] | (p0[13] & g0[12])))));
   assign gp0_4[3] = &p0[15:12];
   
   logic gg0, gp0;
   assign gg0 = gg0_4[3] | (gp0_4[3] & (gg0_4[2] | (gp0_4[2] | (gg0_4[1] | (gp0_4[1] & gg0_4[0])))));
   assign gp0 = &gp0_4;
   
   logic [15:0] c0;
   assign c0[0] = group_carry[0];
   assign c0[1] = g0[0] | (p0[0] & c0[0]);
   assign c0[2] = g0[1] | (p0[1] & c0[1]);
   assign c0[3] = g0[2] | (p0[2] & c0[2]);
   assign c0[4] = g0[3] | (p0[3] & c0[3]);
   
   logic c0_4;
   assign c0_4 = gg0_4[0] | (gp0_4[0] & c0[0]);
   assign c0[5] = g0[4] | (p0[4] & c0_4);
   assign c0[6] = g0[5] | (p0[5] & c0[5]);
   assign c0[7] = g0[6] | (p0[6] & c0[6]);
   assign c0[8] = g0[7] | (p0[7] & c0[7]);
   
   logic c0_8;
   assign c0_8 = gg0_4[1] | (gp0_4[1] & c0_4);
   assign c0[9]  = g0[8]  | (p0[8]  & c0_8);
   assign c0[10] = g0[9]  | (p0[9]  & c0[9]);
   assign c0[11] = g0[10] | (p0[10] & c0[10]);
   assign c0[12] = g0[11] | (p0[11] & c0[11]);
   
   logic c0_12;
   assign c0_12 = gg0_4[2] | (gp0_4[2] & c0_8);
   assign c0[13] = g0[12] | (p0[12] & c0_12);
   assign c0[14] = g0[13] | (p0[13] & c0[13]);
   assign c0[15] = g0[14] | (p0[14] & c0[14]);
   
   assign sum[15:0] = p0 ^ c0;
   assign group_carry[1] = gg0 | (gp0 & group_carry[0]);
   
   // Group 1 (bits 16-31)
   logic [15:0] g1, p1;
   assign g1 = psum_north[31:16] & prod_sext[31:16];
   assign p1 = psum_north[31:16] ^ prod_sext[31:16];
   
   logic [3:0] gg1_4, gp1_4;
   assign gg1_4[0] = g1[3] | (p1[3] & (g1[2] | (p1[2] & (g1[1] | (p1[1] & g1[0])))));
   assign gp1_4[0] = &p1[3:0];
   assign gg1_4[1] = g1[7] | (p1[7] & (g1[6] | (p1[6] & (g1[5] | (p1[5] & g1[4])))));
   assign gp1_4[1] = &p1[7:4];
   assign gg1_4[2] = g1[11] | (p1[11] & (g1[10] | (p1[10] & (g1[9] | (p1[9] & g1[8])))));
   assign gp1_4[2] = &p1[11:8];
   assign gg1_4[3] = g1[15] | (p1[15] & (g1[14] | (p1[14] & (g1[13] | (p1[13] & g1[12])))));
   assign gp1_4[3] = &p1[15:12];
   
   logic gg1, gp1;
   assign gg1 = gg1_4[3] | (gp1_4[3] & (gg1_4[2] | (gp1_4[2] & (gg1_4[1] | (gp1_4[1] & gg1_4[0])))));
   assign gp1 = &gp1_4;
   
   logic [15:0] c1;
   assign c1[0] = group_carry[1];
   assign c1[1] = g1[0] | (p1[0] & c1[0]);
   assign c1[2] = g1[1] | (p1[1] & c1[1]);
   assign c1[3] = g1[2] | (p1[2] & c1[2]);
   assign c1[4] = g1[3] | (p1[3] & c1[3]);
   
   logic c1_4;
   assign c1_4 = gg1_4[0] | (gp1_4[0] & c1[0]);
   assign c1[5] = g1[4] | (p1[4] & c1_4);
   assign c1[6] = g1[5] | (p1[5] & c1[5]);
   assign c1[7] = g1[6] | (p1[6] & c1[6]);
   assign c1[8] = g1[7] | (p1[7] & c1[7]);
   
   logic c1_8;
   assign c1_8 = gg1_4[1] | (gp1_4[1] & c1_4);
   assign c1[9]  = g1[8]  | (p1[8]  & c1_8);
   assign c1[10] = g1[9]  | (p1[9]  & c1[9]);
   assign c1[11] = g1[10] | (p1[10] & c1[10]);
   assign c1[12] = g1[11] | (p1[11] & c1[11]);
   
   logic c1_12;
   assign c1_12 = gg1_4[2] | (gp1_4[2] & c1_8);
   assign c1[13] = g1[12] | (p1[12] & c1_12);
   assign c1[14] = g1[13] | (p1[13] & c1[13]);
   assign c1[15] = g1[14] | (p1[14] & c1[14]);
   
   assign sum[31:16] = p1 ^ c1;

   always_ff @(posedge clk) begin
      if (rst) begin
         w          <= '0;
         a_east     <= '0;
         psum_south <= '0;
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= sum;
         end
      end
   end

endmodule


//#############################################################################
// N x N weight-stationary systolic array (the TPU MXU core)
//#############################################################################

module tpu_array
  #(parameter int N = 8,      // array dimension (N x N PEs)
    parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // accumulator width (signed)
   (
    input  logic              clk,
    input  logic              rst,        // synchronous, active high
    input  logic              ld,         // weight-load shift enable
    input  logic              en,         // compute enable
    input  logic [N*DW-1:0]   w_top,      // weight inputs, one per column (load)
    input  logic [N*DW-1:0]   a_left,     // activation inputs, one per row
    output logic [N*ACCW-1:0] c_bot       // column results, one per column
    );

   // mesh nets
   logic [DW-1:0]   a_h    [0:N-1][0:N];
   logic [ACCW-1:0] psum_v [0:N][0:N-1];
   logic [DW-1:0]   w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i++) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j++) begin : g_edge
         assign psum_v[0][j] = '0;
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
      end
      for (i = 0; i < N; i++) begin : g_pe_row
         for (j = 0; j < N; j++) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j]));
         end
      end
   endgenerate

endmodule


//#############################################################################
// Weight-stationary systolic matrix-multiply tile (a Google-TPU-style MXU)
//#############################################################################

module tpu
  #(parameter int N = 8,      // array dimension (N x N)
    parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // result / accumulator width (signed)
   (
    input  logic              clk,
    input  logic              rst,       // synchronous, active high
    input  logic              w_valid,   // weight-row valid (load phase)
    input  logic [N*DW-1:0]   w_data,    // weight row, packed { col N-1 ... col 0 }
    input  logic              a_valid,   // activation-row valid (compute phase)
    input  logic [N*DW-1:0]   a_data,    // activation row, packed { col N-1 ... col 0 }
    output logic              c_valid,   // result-row valid
    output logic [N*ACCW-1:0] c_data     // result row, packed { col N-1 ... col 0 }
    );

   localparam int LAT = 2*N - 1;
   logic ld, en;
   assign ld = w_valid;
   assign en = ~w_valid;

   //##########################################################################
   // Input skew
   //##########################################################################
   logic [N*DW-1:0] a_left;
   
   generate
      for (genvar i = 0; i < N; i++) begin : g_skew
         logic [DW-1:0] din;
         assign din = a_valid ? a_data[i*DW +: DW] : '0;
         
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            logic [DW-1:0] sr [0:i-1];
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (int k = 0; k < i; k++)
                    sr[k] <= '0;
               end
               else if (en) begin
                  sr[0] <= din;
                  for (int k = 1; k < i; k++)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // The systolic array
   //##########################################################################
   logic [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew
   //##########################################################################
   generate
      for (genvar j = 0; j < N; j++) begin : g_deskew
         localparam int D = N - 1 - j;
         logic [ACCW-1:0] cin;
         assign cin = c_bot[j*ACCW +: ACCW];
         
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin;
         end
         else begin : g_deskewn
            logic [ACCW-1:0] dr [0:D-1];
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (int k = 0; k < D; k++)
                    dr[k] <= '0;
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (int k = 1; k < D; k++)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr[D-1];
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline
   //##########################################################################
   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

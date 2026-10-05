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
   
   genvar grp;
   generate
      for (grp = 0; grp < 2; grp++) begin : cla_16bit_groups
         localparam int BASE = grp * 16;
         
         // Generate/propagate for 16 bits
         logic [15:0] g, p;
         assign g = psum_north[BASE+15:BASE] & prod_sext[BASE+15:BASE];
         assign p = psum_north[BASE+15:BASE] ^ prod_sext[BASE+15:BASE];
         
         // Group generate/propagate - computed hierarchically
         logic [3:0] gg_4, gp_4;
         genvar sg;
         generate
            for (sg = 0; sg < 4; sg++) begin : sub_groups
               localparam int SBASE = sg * 4;
               logic [3:0] sg_g, sg_p;
               assign sg_g = g[SBASE+3:SBASE];
               assign sg_p = p[SBASE+3:SBASE];
               assign gg_4[sg] = sg_g[3] | (sg_p[3] & (sg_g[2] | (sg_p[2] & (sg_g[1] | (sg_p[1] & sg_g[0])))));
               assign gp_4[sg] = &sg_p;
            end
         endgenerate
         
         logic gg, gp;
         assign gg = gg_4[3] | (gp_4[3] & (gg_4[2] | (gp_4[2] & (gg_4[1] | (gp_4[1] & gg_4[0])))));
         assign gp = &gp_4;
         
         // Carry into group
         logic cin;
         assign cin = group_carry[grp];
         
         // Carry out of group
         assign group_carry[grp+1] = gg | (gp & cin);
         
         // Internal carries - 4-bit sub-groups
         logic [15:0] c_int;
         assign c_int[0] = cin;
         
         // First 4-bit sub-group
         assign c_int[1] = g[0] | (p[0] & cin);
         assign c_int[2] = g[1] | (p[1] & c_int[1]);
         assign c_int[3] = g[2] | (p[2] & c_int[2]);
         assign c_int[4] = g[3] | (p[3] & c_int[3]);
         
         // Second 4-bit sub-group
         logic c4;
         assign c4 = gg_4[0] | (gp_4[0] & cin);
         assign c_int[5] = g[4] | (p[4] & c4);
         assign c_int[6] = g[5] | (p[5] & c_int[5]);
         assign c_int[7] = g[6] | (p[6] & c_int[6]);
         assign c_int[8] = g[7] | (p[7] & c_int[7]);
         
         // Third 4-bit sub-group
         logic c8;
         assign c8 = gg_4[1] | (gp_4[1] & c4);
         assign c_int[9]  = g[8]  | (p[8]  & c8);
         assign c_int[10] = g[9]  | (p[9]  & c_int[9]);
         assign c_int[11] = g[10] | (p[10] & c_int[10]);
         assign c_int[12] = g[11] | (p[11] & c_int[11]);
         
         // Fourth 4-bit sub-group
         logic c12;
         assign c12 = gg_4[2] | (gp_4[2] & c8);
         assign c_int[13] = g[12] | (p[12] & c12);
         assign c_int[14] = g[13] | (p[13] & c_int[13]);
         assign c_int[15] = g[14] | (p[14] & c_int[14]);
         
         // Sum bits
         assign sum[BASE+15:BASE] = p ^ c_int;
      end
   endgenerate

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

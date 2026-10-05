//#############################################################################
// Optimized TPU - using a hybrid adder (carry-lookahead at 8-bit boundaries)
//#############################################################################

module hybrid_adder_32 (
    input  [31:0] a,
    input  [31:0] b,
    output [31:0] sum
);
    // 4 x 8-bit carry-lookahead adders with ripple carry between them
    logic [3:0] carry;
    
    // Byte 0 (bits 0-7)
    logic [7:0] g0, p0, c0;
    assign g0 = a[7:0] & b[7:0];
    assign p0 = a[7:0] ^ b[7:0];
    assign c0[0] = 1'b0;
    assign c0[1] = g0[0];
    assign c0[2] = g0[1] | (p0[1] & g0[0]);
    assign c0[3] = g0[2] | (p0[2] & g0[1]) | (p0[2] & p0[1] & g0[0]);
    assign c0[4] = g0[3] | (p0[3] & g0[2]) | (p0[3] & p0[2] & g0[1]) | (p0[3] & p0[2] & p0[1] & g0[0]);
    assign c0[5] = g0[4] | (p0[4] & g0[3]) | (p0[4] & p0[3] & g0[2]) | (p0[4] & p0[3] & p0[2] & g0[1]) | (p0[4] & p0[3] & p0[2] & p0[1] & g0[0]);
    assign c0[6] = g0[5] | (p0[5] & g0[4]) | (p0[5] & p0[4] & g0[3]) | (p0[5] & p0[4] & p0[3] & g0[2]) | (p0[5] & p0[4] & p0[3] & p0[2] & g0[1]) | (p0[5] & p0[4] & p0[3] & p0[2] & p0[1] & g0[0]);
    assign c0[7] = g0[6] | (p0[6] & g0[5]) | (p0[6] & p0[5] & g0[4]) | (p0[6] & p0[5] & p0[4] & g0[3]) | (p0[6] & p0[5] & p0[4] & p0[3] & g0[2]) | (p0[6] & p0[5] & p0[4] & p0[3] & p0[2] & g0[1]) | (p0[6] & p0[5] & p0[4] & p0[3] & p0[2] & p0[1] & g0[0]);
    assign carry[0] = g0[7] | (p0[7] & g0[6]) | (p0[7] & p0[6] & g0[5]) | (p0[7] & p0[6] & p0[5] & g0[4]) | (p0[7] & p0[6] & p0[5] & p0[4] & g0[3]) | (p0[7] & p0[6] & p0[5] & p0[4] & p0[3] & g0[2]) | (p0[7] & p0[6] & p0[5] & p0[4] & p0[3] & p0[2] & g0[1]) | (p0[7] & p0[6] & p0[5] & p0[4] & p0[3] & p0[2] & p0[1] & g0[0]);
    assign sum[7:0] = p0 ^ c0;
    
    // Byte 1 (bits 8-15)
    logic [7:0] g1, p1, c1;
    assign g1 = a[15:8] & b[15:8];
    assign p1 = a[15:8] ^ b[15:8];
    assign c1[0] = carry[0];
    assign c1[1] = g1[0] | (p1[0] & carry[0]);
    assign c1[2] = g1[1] | (p1[1] & g1[0]) | (p1[1] & p1[0] & carry[0]);
    assign c1[3] = g1[2] | (p1[2] & g1[1]) | (p1[2] & p1[1] & g1[0]) | (p1[2] & p1[1] & p1[0] & carry[0]);
    assign c1[4] = g1[3] | (p1[3] & g1[2]) | (p1[3] & p1[2] & g1[1]) | (p1[3] & p1[2] & p1[1] & g1[0]) | (p1[3] & p1[2] & p1[1] & p1[0] & carry[0]);
    assign c1[5] = g1[4] | (p1[4] & g1[3]) | (p1[4] & p1[3] & g1[2]) | (p1[4] & p1[3] & p1[2] & g1[1]) | (p1[4] & p1[3] & p1[2] & p1[1] & g1[0]) | (p1[4] & p1[3] & p1[2] & p1[1] & p1[0] & carry[0]);
    assign c1[6] = g1[5] | (p1[5] & g1[4]) | (p1[5] & p1[4] & g1[3]) | (p1[5] & p1[4] & p1[3] & g1[2]) | (p1[5] & p1[4] & p1[3] & p1[2] & g1[1]) | (p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & g1[0]) | (p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & p1[0] & carry[0]);
    assign c1[7] = g1[6] | (p1[6] & g1[5]) | (p1[6] & p1[5] & g1[4]) | (p1[6] & p1[5] & p1[4] & g1[3]) | (p1[6] & p1[5] & p1[4] & p1[3] & g1[2]) | (p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & g1[1]) | (p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & g1[0]) | (p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & p1[0] & carry[0]);
    assign carry[1] = g1[7] | (p1[7] & g1[6]) | (p1[7] & p1[6] & g1[5]) | (p1[7] & p1[6] & p1[5] & g1[4]) | (p1[7] & p1[6] & p1[5] & p1[4] & g1[3]) | (p1[7] & p1[6] & p1[5] & p1[4] & p1[3] & g1[2]) | (p1[7] & p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & g1[1]) | (p1[7] & p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & g1[0]) | (p1[7] & p1[6] & p1[5] & p1[4] & p1[3] & p1[2] & p1[1] & p1[0] & carry[0]);
    assign sum[15:8] = p1 ^ c1;
    
    // Byte 2 (bits 16-23)
    logic [7:0] g2, p2, c2;
    assign g2 = a[23:16] & b[23:16];
    assign p2 = a[23:16] ^ b[23:16];
    assign c2[0] = carry[1];
    assign c2[1] = g2[0] | (p2[0] & carry[1]);
    assign c2[2] = g2[1] | (p2[1] & g2[0]) | (p2[1] & p2[0] & carry[1]);
    assign c2[3] = g2[2] | (p2[2] & g2[1]) | (p2[2] & p2[1] & g2[0]) | (p2[2] & p2[1] & p2[0] & carry[1]);
    assign c2[4] = g2[3] | (p2[3] & g2[2]) | (p2[3] & p2[2] & g2[1]) | (p2[3] & p2[2] & p2[1] & g2[0]) | (p2[3] & p2[2] & p2[1] & p2[0] & carry[1]);
    assign c2[5] = g2[4] | (p2[4] & g2[3]) | (p2[4] & p2[3] & g2[2]) | (p2[4] & p2[3] & p2[2] & g2[1]) | (p2[4] & p2[3] & p2[2] & p2[1] & g2[0]) | (p2[4] & p2[3] & p2[2] & p2[1] & p2[0] & carry[1]);
    assign c2[6] = g2[5] | (p2[5] & g2[4]) | (p2[5] & p2[4] & g2[3]) | (p2[5] & p2[4] & p2[3] & g2[2]) | (p2[5] & p2[4] & p2[3] & p2[2] & g2[1]) | (p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & g2[0]) | (p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & p2[0] & carry[1]);
    assign c2[7] = g2[6] | (p2[6] & g2[5]) | (p2[6] & p2[5] & g2[4]) | (p2[6] & p2[5] & p2[4] & g2[3]) | (p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & g2[1]) | (p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & g2[0]) | (p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & p2[0] & carry[1]);
    assign carry[2] = g2[7] | (p2[7] & g2[6]) | (p2[7] & p2[6] & g2[5]) | (p2[7] & p2[6] & p2[5] & g2[4]) | (p2[7] & p2[6] & p2[5] & p2[4] & g2[3]) | (p2[7] & p2[6] & p2[5] & p2[4] & p2[3] & g2[2]) | (p2[7] & p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & g2[1]) | (p2[7] & p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & g2[0]) | (p2[7] & p2[6] & p2[5] & p2[4] & p2[3] & p2[2] & p2[1] & p2[0] & carry[1]);
    assign sum[23:16] = p2 ^ c2;
    
    // Byte 3 (bits 24-31)
    logic [7:0] g3, p3, c3;
    assign g3 = a[31:24] & b[31:24];
    assign p3 = a[31:24] ^ b[31:24];
    assign c3[0] = carry[2];
    assign c3[1] = g3[0] | (p3[0] & carry[2]);
    assign c3[2] = g3[1] | (p3[1] & g3[0]) | (p3[1] & p3[0] & carry[2]);
    assign c3[3] = g3[2] | (p3[2] & g3[1]) | (p3[2] & p3[1] & g3[0]) | (p3[2] & p3[1] & p3[0] & carry[2]);
    assign c3[4] = g3[3] | (p3[3] & g3[2]) | (p3[3] & p3[2] & g3[1]) | (p3[3] & p3[2] & p3[1] & g3[0]) | (p3[3] & p3[2] & p3[1] & p3[0] & carry[2]);
    assign c3[5] = g3[4] | (p3[4] & g3[3]) | (p3[4] & p3[3] & g3[2]) | (p3[4] & p3[3] & p3[2] & g3[1]) | (p3[4] & p3[3] & p3[2] & p3[1] & g3[0]) | (p3[4] & p3[3] & p3[2] & p3[1] & p3[0] & carry[2]);
    assign c3[6] = g3[5] | (p3[5] & g3[4]) | (p3[5] & p3[4] & g3[3]) | (p3[5] & p3[4] & p3[3] & p3[2] & g3[1]) | (p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & g3[0]) | (p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & p3[0] & carry[2]);
    assign c3[7] = g3[6] | (p3[6] & g3[5]) | (p3[6] & p3[5] & g3[4]) | (p3[6] & p3[5] & p3[4] & p3[3] & g3[2]) | (p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & g3[1]) | (p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & g3[0]) | (p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & p3[0] & carry[2]);
    assign carry[3] = g3[7] | (p3[7] & g3[6]) | (p3[7] & p3[6] & g3[5]) | (p3[7] & p3[6] & p3[5] & g3[4]) | (p3[7] & p3[6] & p3[5] & p3[4] & g3[3]) | (p3[7] & p3[6] & p3[5] & p3[4] & p3[3] & g3[2]) | (p3[7] & p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & g3[1]) | (p3[7] & p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & g3[0]) | (p3[7] & p3[6] & p3[5] & p3[4] & p3[3] & p3[2] & p3[1] & p3[0] & carry[2]);
    assign sum[31:24] = p3 ^ c3;

endmodule


module tpu_pe
  #(parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // accumulator / partial-sum width (signed)
   (
    input                   clk,
    input                   rst,      // synchronous, active high
    input                   ld,       // weight-load shift enable
    input                   en,       // compute enable
    input [DW-1:0]          w_north,  // weight in from PE above
    input [DW-1:0]          a_west,   // activation in from west
    input [ACCW-1:0]        psum_north, // partial sum in from above
    output [DW-1:0]         w_south,  // stored weight out to PE below
    output reg [DW-1:0]     a_east,   // activation out to east
    output reg [ACCW-1:0]   psum_south // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [DW-1:0] a_signed = a_west;
   wire signed [DW-1:0] w_signed = w;
   wire signed [2*DW-1:0] prod = a_signed * w_signed;
   wire [ACCW-1:0] prod_sext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] sum;
   hybrid_adder_32 u_adder (
       .a   (psum_north),
       .b   (prod_sext),
       .sum (sum)
   );

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
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
// N x N weight-stationary systolic array (the TPU MXU core).
//#############################################################################

module tpu_array
  #(parameter N = 8,   // array dimension (N x N PEs)
    parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // accumulator width (signed)
   (
    input                   clk,
    input                   rst,  // synchronous, active high
    input                   ld,   // weight-load shift enable
    input                   en,   // compute enable
    input [N*DW-1:0]        w_top, // weight inputs, one per column (load)
    input [N*DW-1:0]        a_left, // activation inputs, one per row
    output [N*ACCW-1:0]     c_bot  // column results, one per column
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar          i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
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
// Weight-stationary systolic matrix-multiply tile (a Google-TPU-style MXU).
//#############################################################################

module tpu
  #(parameter N = 8,   // array dimension (N x N)
    parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // result / accumulator width (signed)
   (
    input                   clk,
    input                   rst,   // synchronous, active high
    input                   w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]        w_data,  // weight row, packed { col N-1 ... col 0 }
    input                   a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]        a_data,  // activation row, packed { col N-1 ... col 0 }
    output                  c_valid, // result-row valid
    output [N*ACCW-1:0]     c_data   // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;

   wire      ld = w_valid;
   wire      en = ~w_valid;

   // Input skew
   wire [N*DW-1:0] a_left;
   genvar          i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
               a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   // Systolic array
   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output de-skew
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr[D-1];
         end
      end
   endgenerate

   // Valid pipeline
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule
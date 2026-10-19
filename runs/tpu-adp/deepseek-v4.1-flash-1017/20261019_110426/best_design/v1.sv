//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU style), 8x8, 8-bit signed.
// Optimized: accumulator narrowed to 19 bits (exact: |sum| <= 8*16384 = 2^17,
// 19-bit signed covers +-262143), sign-extended to 32 only at the output.
//#############################################################################
module tpu_pe
  #(parameter DW = 8,
    parameter AW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [AW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [AW-1:0]	  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {AW{1'b0}};
      end
      else begin
         if (ld) w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north + {{(AW-2*DW){prod[2*DW-1]}}, prod};
         end
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19,
    parameter ACCW = 32)
   (
    input		  clk,
    input		  rst,
    input		  w_valid,
    input [N*DW-1:0]	  w_data,
    input		  a_valid,
    input [N*DW-1:0]	  a_data,
    output		  c_valid,
    output [N*ACCW-1:0]  c_data
    );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //--------------------------------------------------------------------
   // input skew: row i delayed i compute cycles, zeros injected when idle
   //--------------------------------------------------------------------
   wire [N*DW-1:0] a_left;
   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_s0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_sn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //--------------------------------------------------------------------
   // systolic array
   //--------------------------------------------------------------------
   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [AW-1:0]   psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]    = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
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

   //--------------------------------------------------------------------
   // output de-skew + sign extension to 32 bits
   //--------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [AW-1:0] cin = psum_v[N][j];
         if (D == 0) begin : g_d0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){cin[AW-1]}}, cin};
         end
         else begin : g_dn
            reg [AW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) dr[k] <= {AW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){dr[D-1][AW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //--------------------------------------------------------------------
   // valid pipeline
   //--------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)          vpipe <= {LAT{1'b0}};
      else if (en)      vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

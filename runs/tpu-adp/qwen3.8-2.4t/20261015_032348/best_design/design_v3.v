//#############################################################################
// Optimized weight-stationary 8x8 systolic matmul tile (cycle-exact with the
// logikbench golden tpu).
//
// Optimizations vs the original:
//  * Internal partial-sum width narrowed from 32 to 19 bits. The true sum of
//    eight signed 8x8 products lies in [-130048, +131072], which fits in 19
//    bits signed, so 19-bit accumulation is bit-exact vs 32-bit. Sign-extended
//    to 32 bits only at the output.
//  * Row-0 PEs omit the partial-sum adder (psum_north is the constant 0).
//#############################################################################

module tpu_pe
  #(parameter DW   = 8,  // operand width (signed)
    parameter PWW  = 19, // internal partial-sum width (signed)
    parameter ROW0 = 0)  // 1 for the top row: psum_north is always 0
   (
    input                  clk,
    input                  rst,
    input                  ld,        // weight-load shift enable
    input                  en,        // compute enable
    input [DW-1:0]         w_north,
    input [DW-1:0]         a_west,
    input [PWW-1:0]        psum_north,
    output [DW-1:0]        w_south,
    output reg [DW-1:0]    a_east,
    output reg [PWW-1:0]   psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PWW-1:0]  pext = {{(PWW-2*DW){prod[2*DW-1]}}, prod};

   always @(posedge clk) begin
      if (rst)
        w <= {DW{1'b0}};
      else if (ld)
        w <= w_north;
   end

   always @(posedge clk) begin
      if (rst)
        a_east <= {DW{1'b0}};
      else if (en)
        a_east <= a_west;
   end

   generate
      if (ROW0) begin : g_r0
         always @(posedge clk) begin
            if (rst)
              psum_south <= {PWW{1'b0}};
            else if (en)
              psum_south <= pext;
         end
      end
      else begin : g_acc
         always @(posedge clk) begin
            if (rst)
              psum_south <= {PWW{1'b0}};
            else if (en)
              psum_south <= psum_north + pext;
         end
      end
   endgenerate

endmodule

module tpu_array
  #(parameter N   = 8,
    parameter DW  = 8,
    parameter PWW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [N*DW-1:0]      w_top,
    input [N*DW-1:0]      a_left,
    output [N*PWW-1:0]    c_bot
    );

   wire [DW-1:0]  a_h    [0:N-1][0:N];
   wire [PWW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]  w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PWW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*PWW +: PWW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PWW(PWW), .ROW0(i == 0 ? 1 : 0)) u_pe
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

module tpu
  #(parameter N    = 8,
    parameter DW   = 8,
    parameter ACCW = 32)
   (
    input                 clk,
    input                 rst,
    input                 w_valid,
    input [N*DW-1:0]      w_data,
    input                 a_valid,
    input [N*DW-1:0]      a_data,
    output                c_valid,
    output [N*ACCW-1:0]   c_data
    );

   localparam LAT = 2*N - 1;
   localparam PWW = 2*DW + $clog2(N); // 19 bits: exact for N DWxDW products

   wire ld = w_valid;
   wire en = ~w_valid;

   //####################################################################
   // Input skew: row i delayed i compute cycles; zeros when !a_valid.
   //####################################################################
   wire [N*DW-1:0] a_left;
   genvar i, j;
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

   //####################################################################
   // Systolic array
   //####################################################################
   wire [N*PWW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PWW(PWW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //####################################################################
   // Output de-skew (PWW bits), then sign-extend to ACCW at the pins.
   //####################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PWW-1:0] cin = c_bot[j*PWW +: PWW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cin[PWW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PWW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PWW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PWW-1:0] cout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cout[PWW-1]}}, cout};
         end
      end
   endgenerate

   //####################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   //####################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

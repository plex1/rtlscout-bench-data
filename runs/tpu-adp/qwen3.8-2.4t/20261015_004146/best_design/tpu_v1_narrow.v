//#############################################################################
// tpu - weight-stationary systolic matmul tile, optimized variant v1.
//
// Same microarchitecture as the reference logikbench tpu (8x8 PE mesh,
// per-row input skew, per-column output de-skew, valid pipeline), but the
// internal partial-sum datapath is narrowed from 32 to IW=19 bits:
//
//   Each output is a sum of 8 signed 8x8 products. The extreme case is
//   8 * (-128 * -128) = +131072 = 2^17 (min is -130048), so a 19-bit signed
//   accumulator can hold every intermediate AND final value exactly. The
//   32-bit sign-extended result is formed only at the output ports (pure
//   wiring), so behavior is bit-exact versus the 32-bit reference.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,   // operand width (signed)
    parameter IW = 19)  // internal accumulator width (signed, exact)
   (
    input                 clk,
    input                 rst,
    input                 ld,        // weight-load shift enable
    input                 en,        // compute enable
    input [DW-1:0]        w_north,
    input [DW-1:0]        a_west,
    input [IW-1:0]        psum_north,
    output [DW-1:0]       w_south,
    output reg [DW-1:0]   a_east,
    output reg [IW-1:0]   psum_south
   );

   reg [DW-1:0] w;
   assign w_south = w;

   // signed 8x8 product, sign-extended to the accumulator width
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [IW-1:0]   pext = {{(IW-2*DW){prod[2*DW-1]}}, prod};

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {IW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north + pext;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [N*DW-1:0]      w_top,
    input [N*DW-1:0]      a_left,
    output [N*IW-1:0]     c_bot
   );

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [IW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {IW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*IW +: IW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .IW(IW)) u_pe
              (.clk(clk), .rst(rst), .ld(ld), .en(en),
               .w_north(w_v[i][j]),
               .a_west(a_h[i][j]),
               .psum_north(psum_v[i][j]),
               .w_south(w_v[i+1][j]),
               .a_east(a_h[i][j+1]),
               .psum_south(psum_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19,      // exact internal accumulator width
    parameter ACCW = 32)    // output width
   (
    input                  clk,
    input                  rst,
    input                  w_valid,
    input [N*DW-1:0]       w_data,
    input                  a_valid,
    input [N*DW-1:0]       a_data,
    output                 c_valid,
    output [N*ACCW-1:0]    c_data
   );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //#######################################################################
   // Input skew: row i delayed i compute cycles; zero when a_valid low.
   //#######################################################################
   wire [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
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

   //#######################################################################
   // Systolic array (narrow IW-bit partial sums).
   //#######################################################################
   wire [N*IW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .IW(IW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //#######################################################################
   // Output de-skew in the narrow domain, then sign-extend to 32 bits.
   //#######################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IW-1:0] cin = c_bot[j*IW +: IW];
         wire [IW-1:0] cout;
         if (D == 0) begin : g_d0
            assign cout = cin;
         end
         else begin : g_dn
            reg [IW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) dr[k] <= {IW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
               end
            end
            assign cout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] =
                {{(ACCW-IW){cout[IW-1]}}, cout};
      end
   endgenerate

   //#######################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   //#######################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

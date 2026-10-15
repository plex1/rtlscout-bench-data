// Variable-width carry-save accumulator.
//
// Carry-save psums are stored with the minimum signed word width needed at
// each vertical stage:
//   stage r terms: widths 1,16,17,17,18,18,18,18,19
// The final per-column carry-propagate add is only done once at the bottom.
// ---------------------------------------------------------------------------

module tpu_pe_csa
  #(parameter DW = 8,   // operand width (signed)
    parameter IW = 1,   // incoming carry-save word width
    parameter OW = 16)  // outgoing carry-save word width
   (
    input                    clk,
    input                    rst,
    input                    ld,
    input                    en,
    input [DW-1:0]          w_north,
    input [DW-1:0]          a_west,
    input [IW-1:0]          s_north,
    input [IW-1:0]          c_north,
    output [DW-1:0]         w_south,
    output reg [DW-1:0]     a_east,
    output reg [OW-1:0]     s_south,
    output reg [OW-1:0]     c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [OW-1:0]   p_ext = {{(OW-2*DW){prod[2*DW-1]}}, prod};

   wire [OW-1:0] s_ext = {{(OW-IW){s_north[IW-1]}}, s_north};
   wire [OW-1:0] c_ext = {{(OW-IW){c_north[IW-1]}}, c_north};

   wire [OW-1:0] s_next;
   wire [OW:0]   cy_full;
   assign cy_full[0] = 1'b0;

   genvar k;
   generate
      for (k = 0; k < OW; k = k + 1) begin : g_csa_bit
         assign s_next[k] = s_ext[k] ^ c_ext[k] ^ p_ext[k];
         assign cy_full[k+1] = (s_ext[k] & c_ext[k]) |
                               (s_ext[k] & p_ext[k])  |
                               (c_ext[k] & p_ext[k]);
      end
   endgenerate

   wire [OW-1:0] c_next = cy_full[OW-1:0];

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         s_south    <= {OW{1'b0}};
         c_south    <= {OW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east  <= a_west;
            s_south <= s_next;
            c_south <= c_next;
         end
      end
   end

endmodule

module tpu_array_csa
  #(parameter N = 8,
    parameter DW = 8,
    parameter OUTW = 19)
   (
    input                          clk,
    input                          rst,
    input                          ld,
    input                          en,
    input  [N*DW-1:0]             w_top,
    input  [N*DW-1:0]             a_left,
    output [N*OUTW-1:0]           c_bot
    );

   localparam integer W0 = 1;
   localparam integer W1 = 16;
   localparam integer W2 = 17;
   localparam integer W3 = 17;
   localparam integer W4 = 18;
   localparam integer W5 = 18;
   localparam integer W6 = 18;
   localparam integer W7 = 18;
   localparam integer W8 = 19;

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [DW-1:0] w_v [0:N][0:N-1];

   wire [N*W0-1:0] psum_s0, psum_c0;
   wire [N*W1-1:0] psum_s1, psum_c1;
   wire [N*W2-1:0] psum_s2, psum_c2;
   wire [N*W3-1:0] psum_s3, psum_c3;
   wire [N*W4-1:0] psum_s4, psum_c4;
   wire [N*W5-1:0] psum_s5, psum_c5;
   wire [N*W6-1:0] psum_s6, psum_c6;
   wire [N*W7-1:0] psum_s7, psum_c7;
   wire [N*W8-1:0] psum_s8, psum_c8;

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end

      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_s0[j*W0 +: W0] = {W0{1'b0}};
         assign psum_c0[j*W0 +: W0] = {W0{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];

         wire [W8-1:0] sum_final = psum_s8[j*W8 +: W8] +
                                   psum_c8[j*W8 +: W8];
         assign c_bot[j*OUTW +: OUTW] = sum_final[OUTW-1:0];
      end

      for (j = 0; j < N; j = j + 1) begin : g_row0
         tpu_pe_csa #(.DW(DW), .IW(W0), .OW(W1)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[0][j]), .a_west(a_h[0][j]),
            .s_north(psum_s0[j*W0 +: W0]),
            .c_north(psum_c0[j*W0 +: W0]),
            .w_south(w_v[1][j]), .a_east(a_h[0][j+1]),
            .s_south(psum_s1[j*W1 +: W1]),
            .c_south(psum_c1[j*W1 +: W1]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row1
         tpu_pe_csa #(.DW(DW), .IW(W1), .OW(W2)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[1][j]), .a_west(a_h[1][j]),
            .s_north(psum_s1[j*W1 +: W1]),
            .c_north(psum_c1[j*W1 +: W1]),
            .w_south(w_v[2][j]), .a_east(a_h[1][j+1]),
            .s_south(psum_s2[j*W2 +: W2]),
            .c_south(psum_c2[j*W2 +: W2]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row2
         tpu_pe_csa #(.DW(DW), .IW(W2), .OW(W3)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[2][j]), .a_west(a_h[2][j]),
            .s_north(psum_s2[j*W2 +: W2]),
            .c_north(psum_c2[j*W2 +: W2]),
            .w_south(w_v[3][j]), .a_east(a_h[2][j+1]),
            .s_south(psum_s3[j*W3 +: W3]),
            .c_south(psum_c3[j*W3 +: W3]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row3
         tpu_pe_csa #(.DW(DW), .IW(W3), .OW(W4)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[3][j]), .a_west(a_h[3][j]),
            .s_north(psum_s3[j*W3 +: W3]),
            .c_north(psum_c3[j*W3 +: W3]),
            .w_south(w_v[4][j]), .a_east(a_h[3][j+1]),
            .s_south(psum_s4[j*W4 +: W4]),
            .c_south(psum_c4[j*W4 +: W4]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row4
         tpu_pe_csa #(.DW(DW), .IW(W4), .OW(W5)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[4][j]), .a_west(a_h[4][j]),
            .s_north(psum_s4[j*W4 +: W4]),
            .c_north(psum_c4[j*W4 +: W4]),
            .w_south(w_v[5][j]), .a_east(a_h[4][j+1]),
            .s_south(psum_s5[j*W5 +: W5]),
            .c_south(psum_c5[j*W5 +: W5]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row5
         tpu_pe_csa #(.DW(DW), .IW(W5), .OW(W6)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[5][j]), .a_west(a_h[5][j]),
            .s_north(psum_s5[j*W5 +: W5]),
            .c_north(psum_c5[j*W5 +: W5]),
            .w_south(w_v[6][j]), .a_east(a_h[5][j+1]),
            .s_south(psum_s6[j*W6 +: W6]),
            .c_south(psum_c6[j*W6 +: W6]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row6
         tpu_pe_csa #(.DW(DW), .IW(W6), .OW(W7)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[6][j]), .a_west(a_h[6][j]),
            .s_north(psum_s6[j*W6 +: W6]),
            .c_north(psum_c6[j*W6 +: W6]),
            .w_south(w_v[7][j]), .a_east(a_h[6][j+1]),
            .s_south(psum_s7[j*W7 +: W7]),
            .c_south(psum_c7[j*W7 +: W7]));
      end
      for (j = 0; j < N; j = j + 1) begin : g_row7
         tpu_pe_csa #(.DW(DW), .IW(W7), .OW(W8)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[7][j]), .a_west(a_h[7][j]),
            .s_north(psum_s7[j*W7 +: W7]),
            .c_north(psum_c7[j*W7 +: W7]),
            .w_south(w_v[8][j]), .a_east(a_h[7][j+1]),
            .s_south(psum_s8[j*W8 +: W8]),
            .c_south(psum_c8[j*W8 +: W8]));
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19,
    parameter OUTW = 32)
   (
    input                    clk,
    input                    rst,
    input                    w_valid,
    input  [N*DW-1:0]        w_data,
    input                    a_valid,
    input  [N*DW-1:0]        a_data,
    output                   c_valid,
    output [N*OUTW-1:0]      c_data
    );

   localparam LAT = 2*N - 1;
   wire ld = w_valid;
   wire en = ~w_valid;

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

   wire [N*ACCW-1:0] c_bot;
   tpu_array_csa #(.N(N), .DW(DW), .OUTW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*OUTW +: OUTW] =
                  {{(OUTW-ACCW){cin[ACCW-1]}}, cin};
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
            assign c_data[j*OUTW +: OUTW] =
                  {{(OUTW-ACCW){dr[D-1][ACCW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule
// Variable-width accumulator variant.
// The signed 8x8 dot product bounds are:
//   term in [-16256, 16384], sum of r terms in [-16256*r, 16384*r].
// Minimal signed accumulator widths per vertical stage:
//   r: 0  1  2  3  4  5  6  7  8
//   w: 1 16 17 17 18 18 18 18 19
// This trims one bit from the early psum registers/adders while retaining
// exact two's-complement behavior and the original cycle timing.
// ---------------------------------------------------------------------------
module tpu_pe
  #(parameter DW = 8,      // operand width (signed)
    parameter PIW = 1,     // partial-sum input width (signed)
    parameter POW = 16)    // partial-sum output width (signed)
   (
    input                    clk,
    input                    rst,
    input                    ld,
    input                    en,
    input  [DW-1:0]          w_north,
    input  [DW-1:0]          a_west,
    input  [PIW-1:0]         psum_north,
    output [DW-1:0]          w_south,
    output reg [DW-1:0]      a_east,
    output reg [POW-1:0]     psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire signed [POW-1:0] psum_in_ext = {{(POW-PIW){psum_north[PIW-1]}}, psum_north};
   wire signed [POW-1:0] prod_ext   = {{(POW-2*DW){prod[2*DW-1]}}, prod};

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {POW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_in_ext + prod_ext;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8)
   (
    input                          clk,
    input                          rst,
    input                          ld,
    input                          en,
    input  [N*DW-1:0]             w_top,
    input  [N*DW-1:0]             a_left,
    output [N*19-1:0]             c_bot
    );

   // per-vertical-stage psum widths (0..N terms)
   localparam integer PW0 = 1;
   localparam integer PW1 = 16;
   localparam integer PW2 = 17;
   localparam integer PW3 = 17;
   localparam integer PW4 = 18;
   localparam integer PW5 = 18;
   localparam integer PW6 = 18;
   localparam integer PW7 = 18;
   localparam integer PW8 = 19;

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   // one flattened bus per psum stage; each column slice is PW[s] bits wide
   wire [N*PW0-1:0]             psum_s0;
   wire [N*PW1-1:0]             psum_s1;
   wire [N*PW2-1:0]             psum_s2;
   wire [N*PW3-1:0]             psum_s3;
   wire [N*PW4-1:0]             psum_s4;
   wire [N*PW5-1:0]             psum_s5;
   wire [N*PW6-1:0]             psum_s6;
   wire [N*PW7-1:0]             psum_s7;
   wire [N*PW8-1:0]             psum_s8;

   // this is used by generate instances; localparam expression is static
   // but cannot be directly used as a variable array, so we split it below.

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end

      for (j = 0; j < N; j = j + 1) begin : g_edge
         // top psum stage is zero (1 bit per column)
         assign psum_s0[j*PW0 +: PW0] = {PW0{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign c_bot[j*19 +: 19] = psum_s8[j*PW8 +: PW8];
      end

      // Row 0: PIW=1 -> POW=16
      for (j = 0; j < N; j = j + 1) begin : g_row0
         tpu_pe #(.DW(DW), .PIW(PW0), .POW(PW1)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[0][j]),
            .a_west     (a_h[0][j]),
            .psum_north (psum_s0[j*PW0 +: PW0]),
            .w_south    (w_v[1][j]),
            .a_east     (a_h[0][j+1]),
            .psum_south (psum_s1[j*PW1 +: PW1]));
      end

      // Row 1: PIW=16 -> POW=17
      for (j = 0; j < N; j = j + 1) begin : g_row1
         tpu_pe #(.DW(DW), .PIW(PW1), .POW(PW2)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[1][j]),
            .a_west     (a_h[1][j]),
            .psum_north (psum_s1[j*PW1 +: PW1]),
            .w_south    (w_v[2][j]),
            .a_east     (a_h[1][j+1]),
            .psum_south (psum_s2[j*PW2 +: PW2]));
      end

      // Row 2: PIW=17 -> POW=17
      for (j = 0; j < N; j = j + 1) begin : g_row2
         tpu_pe #(.DW(DW), .PIW(PW2), .POW(PW3)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[2][j]),
            .a_west     (a_h[2][j]),
            .psum_north (psum_s2[j*PW2 +: PW2]),
            .w_south    (w_v[3][j]),
            .a_east     (a_h[2][j+1]),
            .psum_south (psum_s3[j*PW3 +: PW3]));
      end

      // Row 3: PIW=17 -> POW=18
      for (j = 0; j < N; j = j + 1) begin : g_row3
         tpu_pe #(.DW(DW), .PIW(PW3), .POW(PW4)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[3][j]),
            .a_west     (a_h[3][j]),
            .psum_north (psum_s3[j*PW3 +: PW3]),
            .w_south    (w_v[4][j]),
            .a_east     (a_h[3][j+1]),
            .psum_south (psum_s4[j*PW4 +: PW4]));
      end

      // Row 4: PIW=18 -> POW=18
      for (j = 0; j < N; j = j + 1) begin : g_row4
         tpu_pe #(.DW(DW), .PIW(PW4), .POW(PW5)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[4][j]),
            .a_west     (a_h[4][j]),
            .psum_north (psum_s4[j*PW4 +: PW4]),
            .w_south    (w_v[5][j]),
            .a_east     (a_h[4][j+1]),
            .psum_south (psum_s5[j*PW5 +: PW5]));
      end

      // Row 5: PIW=18 -> POW=18
      for (j = 0; j < N; j = j + 1) begin : g_row5
         tpu_pe #(.DW(DW), .PIW(PW5), .POW(PW6)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[5][j]),
            .a_west     (a_h[5][j]),
            .psum_north (psum_s5[j*PW5 +: PW5]),
            .w_south    (w_v[6][j]),
            .a_east     (a_h[5][j+1]),
            .psum_south (psum_s6[j*PW6 +: PW6]));
      end

      // Row 6: PIW=18 -> POW=18
      for (j = 0; j < N; j = j + 1) begin : g_row6
         tpu_pe #(.DW(DW), .PIW(PW6), .POW(PW7)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[6][j]),
            .a_west     (a_h[6][j]),
            .psum_north (psum_s6[j*PW6 +: PW6]),
            .w_south    (w_v[7][j]),
            .a_east     (a_h[6][j+1]),
            .psum_south (psum_s7[j*PW7 +: PW7]));
      end

      // Row 7: PIW=18 -> POW=19
      for (j = 0; j < N; j = j + 1) begin : g_row7
         tpu_pe #(.DW(DW), .PIW(PW7), .POW(PW8)) u_pe
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[7][j]),
            .a_west     (a_h[7][j]),
            .psum_north (psum_s7[j*PW7 +: PW7]),
            .w_south    (w_v[8][j]),
            .a_east     (a_h[7][j+1]),
            .psum_south (psum_s8[j*PW8 +: PW8]));
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
   tpu_array #(.N(N), .DW(DW)) u_array
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
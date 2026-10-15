//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU).
// Booth variant v2:
//  * activation radix-4 Booth recoded once at the input ({t1,t2,neg} x 4);
//  * each PE: 4 term rows + correction row + psum pair = 7 rows, one balanced
//    3:2 CSA tree (4 levels); only the column-bottom CPAs propagate carries;
//  * psum pairs in OFFSET form: (s+c) mod 2^AW = partial sum + OFF, with
//    OFF = 16256*N so all values fit AW=18 unsigned bits; the offset is
//    removed by a constant bias in the bottom adder (result is 19-bit signed).
//  * term rows written as XORs of disjoint masked bases (flat mapping).
//#############################################################################

module tpu_pe
  #(parameter DW  = 8,
    parameter AW  = 18,
    parameter OFF = 130048)
   (
    input               clk,
    input               ld,
    input               en,
    input  [DW-1:0]     w_north,
    input  [11:0]       d_west,
    input  [AW-1:0]     s_north,
    input  [AW-1:0]     c_north,
    output [DW-1:0]     w_south,
    output reg [11:0]   d_east,
    output reg [AW-1:0] s_south,
    output reg [AW-1:0] c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [AW-1:0] ws = $signed(w);

   wire [AW-1:0] row [0:6];
   genvar k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_term
         assign row[k] = ((ws << (2*k))     & {AW{d_west[3*k]}})   ^
                         ((ws << (2*k + 1)) & {AW{d_west[3*k+1]}}) ^
                         {AW{d_west[3*k+2]}};
      end
   endgenerate
   wire [3:0] cnt4 = {3'b0, d_west[2]} + {3'b0, d_west[5]} +
                     {3'b0, d_west[8]} + {3'b0, d_west[11]};
   assign row[4] = {{(AW-3){1'b0}}, cnt4[2:0]};
   assign row[5] = s_north;
   assign row[6] = c_north;

   wire [AW-1:0] l1 [0:3];
   wire [AW-1:0] l2 [0:2];
   wire [AW-1:0] l3 [0:1];

   assign l1[0] = row[0] ^ row[1] ^ row[2];
   assign l1[1] = ((row[0] & row[1]) | (row[0] & row[2]) |
                   (row[1] & row[2])) << 1;
   assign l1[2] = row[3] ^ row[4] ^ row[5];
   assign l1[3] = ((row[3] & row[4]) | (row[3] & row[5]) |
                   (row[4] & row[5])) << 1;

   assign l2[0] = l1[0] ^ l1[1] ^ l1[2];
   assign l2[1] = ((l1[0] & l1[1]) | (l1[0] & l1[2]) |
                   (l1[1] & l1[2])) << 1;
   assign l2[2] = l1[3];

   assign l3[0] = l2[0] ^ l2[1] ^ l2[2];
   assign l3[1] = ((l2[0] & l2[1]) | (l2[0] & l2[2]) |
                   (l2[1] & l2[2])) << 1;

   wire [AW-1:0] s_next = l3[0] ^ l3[1] ^ row[6];
   wire [AW-1:0] c_next = ((l3[0] & l3[1]) | (l3[0] & row[6]) |
                           (l3[1] & row[6])) << 1;

   always @(posedge clk) begin
      if (ld) w <= w_north;
      if (en) begin
         d_east  <= d_west;
         s_south <= s_next;
         c_south <= c_next;
      end
   end

endmodule

module tpu_array
  #(parameter N   = 8,
    parameter DW  = 8,
    parameter AW  = 18,
    parameter OFF = 130048)
   (
    input               clk,
    input               ld,
    input               en,
    input  [N*DW-1:0]   w_top,
    input  [N*12-1:0]   d_left,
    output [N*(AW+1)-1:0] c_bot
    );

   wire [N*(N+1)*AW-1:0] s_v, c_v;
   wire [N*(N+1)*DW-1:0] w_v;
   wire [N*(N+1)*12-1:0] d_v;

   genvar i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_top
         assign s_v[j*AW +: AW] = OFF[AW-1:0];   // offset-form zero
         assign c_v[j*AW +: AW] = {AW{1'b0}};
      end
      for (j = 0; j < N; j = j + 1) begin : g_wtop
         assign w_v[j*DW +: DW] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_dwest
         assign d_v[(i*(N+1))*12 +: 12] = d_left[i*12 +: 12];
      end
      // bottom: remove the offset bias (19-bit signed result)
      for (j = 0; j < N; j = j + 1) begin : g_bot
         localparam BIAS = (1 << (AW+1)) - OFF;
         wire [AW-1:0] sv = s_v[(N*N+j)*AW +: AW];
         wire [AW-1:0] cv = c_v[(N*N+j)*AW +: AW];
         assign c_bot[j*(AW+1) +: (AW+1)] =
                {1'b0, sv + cv} + BIAS[AW:0];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe #(.DW(DW), .AW(AW), .OFF(OFF)) u_pe
              (.clk     (clk),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[(i*N+j)*DW +: DW]),
               .d_west  (d_v[(i*(N+1)+j)*12 +: 12]),
               .s_north (s_v[(i*N+j)*AW +: AW]),
               .c_north (c_v[(i*N+j)*AW +: AW]),
               .w_south (w_v[((i+1)*N+j)*DW +: DW]),
               .d_east  (d_v[(i*(N+1)+j+1)*12 +: 12]),
               .s_south (s_v[((i+1)*N+j)*AW +: AW]),
               .c_south (c_v[((i+1)*N+j)*AW +: AW]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N    = 8,
    parameter DW   = 8,
    parameter ACCW = 32)
   (
    input               clk,
    input               rst,
    input               w_valid,
    input  [N*DW-1:0]   w_data,
    input               a_valid,
    input  [N*DW-1:0]   a_data,
    output              c_valid,
    output [N*ACCW-1:0] c_data
    );

   localparam LAT = 2*N - 1;
   localparam PPMAX = (1 << (DW-1)) * (1 << (DW-1));
   localparam PPMIN = (1 << (DW-1)) * ((1 << (DW-1)) - 1);
   localparam OFF   = PPMIN * N;
   localparam AW    = $clog2(OFF + PPMAX*N + 1);        // 18
   localparam SUMW  = AW + 1;                            // 19-bit signed

   wire ld = w_valid;
   wire en = ~w_valid;

   function [11:0] booth (input [DW-1:0] a);
      integer k;
      reg y, y1, y2;
      begin
         for (k = 0; k < 4; k = k + 1) begin
            y  = a[2*k];
            y1 = (k == 0) ? 1'b0 : a[2*k-1];
            y2 = a[2*k+1];
            booth[3*k]   = y ^ y1;
            booth[3*k+1] = (y1 & y & ~y2) | (~y1 & ~y & y2);
            booth[3*k+2] = y2 & ~(y1 & y);
         end
      end
   endfunction

   wire [N*12-1:0] d_left;
   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [11:0] din = a_valid ? booth(a_data[i*DW +: DW]) : 12'b0;
         if (i == 0) begin : g_skew0
            assign d_left[i*12 +: 12] = din;
         end
         else begin : g_skewn
            reg [11:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign d_left[i*12 +: 12] = sr[i-1];
         end
      end
   endgenerate

   wire [N*SUMW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(AW), .OFF(OFF)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .d_left(d_left), .c_bot(c_bot));

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [SUMW-1:0] cin = c_bot[j*SUMW +: SUMW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-SUMW){cin[SUMW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [SUMW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-SUMW){dr[D-1][SUMW-1]}}, dr[D-1]};
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

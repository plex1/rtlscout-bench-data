//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU).
// Best design: one-hot signed-digit recode + carry-save psum pairs.
//  * a recoded once at the input: a = sum_{k<=2} d_k 4^k + d3*64,
//    d_k in {0,1,-1,-2}, d3 in {-2..2}; digits are ONE-HOT (p1,n1,n2)/(p1,p2,
//    n1,n2) so masked bases combine with OR (== ADD, masks exclusive) and no
//    per-digit base mux is needed;
//  * each PE stores w and -w (9-bit, computed at load);
//  * 4 term rows + psum pair = 6 rows, one balanced 3:2 CSA tree (3 levels);
//  * psum pairs store (partial sum + OFF) mod 2^AW, AW = 18 unsigned bits
//    (OFF = 16256*N); offset removed by a constant subtract on the output;
//  * carry row stored shifted right one (its LSB is always 0).
//#############################################################################

module tpu_pe
  #(parameter DW  = 8,
    parameter AW  = 18)
   (
    input               clk,
    input               ld,
    input               en,
    input  [DW-1:0]     w_north,
    input  [12:0]       d_west,
    input  [AW-1:0]     s_north,
    input  [AW-1:0]     c_north,
    output [DW-1:0]     w_south,
    output reg [12:0]   d_east,
    output reg [AW-1:0] s_south,
    output [AW-1:0]     c_south
    );

   reg [DW-1:0] w;
   reg [DW:0]   negw;
   reg [AW-2:0] c_hi;                  // c_south >> 1 (LSB is always 0)
   assign w_south = w;
   assign c_south = {c_hi, 1'b0};

   wire signed [AW-1:0] ws  = $signed(w);
   wire signed [AW-1:0] wns = $signed(negw);

   wire [AW-1:0] row [0:5];
   genvar k;
   generate
      for (k = 0; k < 3; k = k + 1) begin : g_term
         assign row[k] = ((ws         & {AW{d_west[3*k]}})   |
                          (wns        & {AW{d_west[3*k+1]}}) |
                          ((wns << 1) & {AW{d_west[3*k+2]}})) << (2*k);
      end
   endgenerate
   assign row[3] = ((ws         & {AW{d_west[9]}})   |
                    ((ws  << 1) & {AW{d_west[10]}})  |
                    (wns        & {AW{d_west[11]}})  |
                    ((wns << 1) & {AW{d_west[12]}})) << 6;
   assign row[4] = s_north;
   assign row[5] = c_north;

   // balanced 3:2 CSA tree: 6 -> 4 -> 3 -> 2
   wire [AW-1:0] l1 [0:3];
   wire [AW-1:0] l2 [0:1];

   assign l1[0] = row[0] ^ row[1] ^ row[2];
   assign l1[1] = ((row[0] & row[1]) | (row[0] & row[2]) |
                   (row[1] & row[2])) << 1;
   assign l1[2] = row[3] ^ row[4] ^ row[5];
   assign l1[3] = ((row[3] & row[4]) | (row[3] & row[5]) |
                   (row[4] & row[5])) << 1;

   assign l2[0] = l1[0] ^ l1[1] ^ l1[2];
   assign l2[1] = ((l1[0] & l1[1]) | (l1[0] & l1[2]) |
                   (l1[1] & l1[2])) << 1;

   wire [AW-1:0] s_next = l2[0] ^ l2[1] ^ l1[3];
   wire [AW-1:0] c_next = ((l2[0] & l2[1]) | (l2[0] & l1[3]) |
                           (l2[1] & l1[3])) << 1;

   always @(posedge clk) begin
      if (ld) begin
         w    <= w_north;
         negw <= -$signed({w_north[DW-1], w_north});
      end
      if (en) begin
         d_east  <= d_west;
         s_south <= s_next;
         c_hi    <= c_next[AW-1:1];
      end
   end

endmodule

module tpu_array
  #(parameter N   = 8,
    parameter DW  = 8,
    parameter AW  = 18,
    parameter OFF = 130048)
   (
    input                 clk,
    input                 ld,
    input                 en,
    input  [N*DW-1:0]     w_top,
    input  [N*13-1:0]     d_left,
    output [N*AW-1:0]     c_bot
    );

   wire [N*(N+1)*AW-1:0] s_v, c_v;
   wire [N*(N+1)*DW-1:0] w_v;
   wire [N*(N+1)*13-1:0] d_v;

   genvar i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_top
         assign s_v[j*AW +: AW] = OFF[AW-1:0];
         assign c_v[j*AW +: AW] = {AW{1'b0}};
      end
      for (j = 0; j < N; j = j + 1) begin : g_wtop
         assign w_v[j*DW +: DW] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_dwest
         assign d_v[(i*(N+1))*13 +: 13] = d_left[i*13 +: 13];
      end
      for (j = 0; j < N; j = j + 1) begin : g_bot
         wire [AW-1:0] sv = s_v[(N*N+j)*AW +: AW];
         wire [AW-1:0] cv = c_v[(N*N+j)*AW +: AW];
         assign c_bot[j*AW +: AW] = sv + cv;
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
              (.clk     (clk),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[(i*N+j)*DW +: DW]),
               .d_west  (d_v[(i*(N+1)+j)*13 +: 13]),
               .s_north (s_v[(i*N+j)*AW +: AW]),
               .c_north (c_v[(i*N+j)*AW +: AW]),
               .w_south (w_v[((i+1)*N+j)*DW +: DW]),
               .d_east  (d_v[(i*(N+1)+j+1)*13 +: 13]),
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

   localparam LAT   = 2*N - 1;
   localparam PPMAX = (1 << (DW-1)) * (1 << (DW-1));
   localparam PPMIN = (1 << (DW-1)) * ((1 << (DW-1)) - 1);
   localparam OFF   = PPMIN * N;
   localparam AW    = $clog2(OFF + PPMAX*N + 1);

   wire ld = w_valid;
   wire en = ~w_valid;

   function [12:0] recode (input [DW-1:0] a);
      reg signed [9:0] A0, A1, A2, A3, S;
      begin
         A0 = {{(10-DW){a[DW-1]}}, a};
         S  = A0 >>> 2;
         A1 = S + {9'b0, A0[1]};
         S  = A1 >>> 2;
         A2 = S + {9'b0, A1[1]};
         S  = A2 >>> 2;
         A3 = S + {9'b0, A2[1]};
         // d_k = A[0] - 2*A[1] in {0,1,-1,-2} -> one-hot (p1,n1,n2)
         recode[0] =  A0[0] & ~A0[1];
         recode[1] =  A0[0] &  A0[1];
         recode[2] = ~A0[0] &  A0[1];
         recode[3] =  A1[0] & ~A1[1];
         recode[4] =  A1[0] &  A1[1];
         recode[5] = ~A1[0] &  A1[1];
         recode[6] =  A2[0] & ~A2[1];
         recode[7] =  A2[0] &  A2[1];
         recode[8] = ~A2[0] &  A2[1];
         // d3 = A3 in {-2..2} -> one-hot (p1,p2,n1,n2)
         recode[9]  = (A3 ==  1);
         recode[10] = (A3 ==  2);
         recode[11] = (A3 == -1);
         recode[12] = (A3 == -2);
      end
   endfunction

   wire [N*13-1:0] d_left;
   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [12:0] din = a_valid ? recode(a_data[i*DW +: DW]) : 13'b0;
         if (i == 0) begin : g_skew0
            assign d_left[i*13 +: 13] = din;
         end
         else begin : g_skewn
            reg [12:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign d_left[i*13 +: 13] = sr[i-1];
         end
      end
   endgenerate

   wire [N*AW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(AW), .OFF(OFF)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .d_left(d_left), .c_bot(c_bot));

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [AW-1:0] cin = c_bot[j*AW +: AW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){1'b0}}, cin} - OFF;
         end
         else begin : g_deskewn
            reg [AW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-AW){1'b0}}, dr[D-1]} - OFF;
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

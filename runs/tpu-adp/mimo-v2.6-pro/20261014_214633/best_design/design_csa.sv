//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU).
// Carry-save variant: the vertical partial sums are kept in carry-save form
// ({s,c} pairs), so each PE performs one CSA level (no carry chain) per cycle;
// the carry-propagate add happens once at the bottom of each column.
//
// Invariant: (s + c) mod 2^AW = the partial sum (two's complement, AW bits).
// AW = 19 bits for DW=8, N=8 (sum of k<=8 signed products fits signed 19-bit).
// The psum pair is born as {0,0} at the top of each column, so garbage init
// values never reach a checked output.
//#############################################################################

module tpu_pe
  #(parameter DW  = 8,
    parameter AW  = 19)
   (
    input               clk,
    input               ld,
    input               en,
    input  [DW-1:0]     w_north,
    input  [DW-1:0]     a_west,
    input  [AW-1:0]     s_north,
    input  [AW-1:0]     c_north,
    output [DW-1:0]     w_south,
    output reg [DW-1:0] a_east,
    output reg [AW-1:0] s_south,
    output reg [AW-1:0] c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // one CSA level: {s,c} + product -> {s',c'}, c' always even
   wire [AW-1:0] p = {{(AW-2*DW){prod[2*DW-1]}}, prod};
   wire [AW-1:0] s_next = s_north ^ c_north ^ p;
   wire [AW-1:0] c_next = ((s_north & c_north) | (s_north & p) |
                           (c_north & p)) << 1;

   always @(posedge clk) begin
      if (ld) w <= w_north;
      if (en) begin
         a_east <= a_west;
         s_south <= s_next;
         c_south <= c_next;
      end
   end

endmodule

//#############################################################################
// N x N PE mesh. Slot packing: slot(b,j) at [ (b*N+j)*SLOT +: SLOT ]
//   - s/c psum bands b = 0..N (band 0 = {0,0}, band N = results)
//   - weight bands b = 0..N (band 0 = w_top, band b+1 = w_south of row b)
//   - activation slots i*(N+1)+j (j=0 = a_left[i], j+1 = a_east of PE(i,j))
//#############################################################################
module tpu_array
  #(parameter N    = 8,
    parameter DW   = 8,
    parameter AW   = 19)
   (
    input               clk,
    input               ld,
    input               en,
    input  [N*DW-1:0]   w_top,
    input  [N*DW-1:0]   a_left,
    output [N*AW-1:0]   c_bot
    );

   wire [N*(N+1)*AW-1:0] s_v, c_v;
   wire [N*(N+1)*DW-1:0] w_v, a_v;

   genvar i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_top
         assign s_v[j*AW +: AW] = {AW{1'b0}};
         assign c_v[j*AW +: AW] = {AW{1'b0}};
      end
      for (j = 0; j < N; j = j + 1) begin : g_wtop
         assign w_v[j*DW +: DW] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_awest
         assign a_v[(i*(N+1))*DW +: DW] = a_left[i*DW +: DW];
      end
      // bottom: one carry-propagate add per column
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign c_bot[j*AW +: AW] =
                s_v[(N*N+j)*AW +: AW] + c_v[(N*N+j)*AW +: AW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
              (.clk     (clk),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[(i*N+j)*DW +: DW]),
               .a_west  (a_v[(i*(N+1)+j)*DW +: DW]),
               .s_north (s_v[(i*N+j)*AW +: AW]),
               .c_north (c_v[(i*N+j)*AW +: AW]),
               .w_south (w_v[((i+1)*N+j)*DW +: DW]),
               .a_east  (a_v[(i*(N+1)+j+1)*DW +: DW]),
               .s_south (s_v[((i+1)*N+j)*AW +: AW]),
               .c_south (c_v[((i+1)*N+j)*AW +: AW]));
         end
      end
   endgenerate

endmodule

//#############################################################################
// Top: input skew -> systolic array -> output de-skew -> valid pipeline.
//#############################################################################
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
   localparam AW  = $clog2((((1<<(DW-1))*(1<<(DW-1)))*N) + 1) + 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i, j;
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
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   wire [N*AW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(AW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [AW-1:0] cin = c_bot[j*AW +: AW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){cin[AW-1]}}, cin};
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
                   {{(ACCW-AW){dr[D-1][AW-1]}}, dr[D-1]};
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

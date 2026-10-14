//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU)
// Cycle-accurate re-architecture of the logikbench tpu benchmark.
//
// Topology is the canonical one (input skew i, horizontal hop j, output
// de-skew 7-j): PE(i,j) captures its product for result-row m on compute
// edge T+i+j.  This per-(i,j) capture schedule is REQUIRED for exact
// equivalence -- a mid-flight weight reload must split old/new weights along
// the anti-diagonal i+j = const, which only the hop-based dataflow produces.
//
// Optimization vs the original: the accumulators only need to hold sums of
// k <= N signed 8x8 products, which lie in [-16256k, 16384k]; for k = 8 that
// is [-130048, +131072], exactly fitting 19 signed bits.  Partial sums and
// the output de-skew chain therefore use PW = 19 bits instead of 32, and
// c_data is formed by sign extension (checked only while c_valid, where the
// value is a true result row).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,    // operand width (signed)
    parameter PW = 19)   // partial-sum width (signed)
   (
    input                  clk,
    input                  rst,       // async, active high
    input                  ld,        // weight-load shift enable
    input                  en,        // compute enable
    input  [DW-1:0]        w_north,   // weight in from PE above
    input  [DW-1:0]        a_west,    // activation in from west
    input  [PW-1:0]        psum_north,// partial sum in from above
    output [DW-1:0]        w_south,   // stored weight out to PE below
    output reg [DW-1:0]    a_east,    // activation out to east
    output reg [PW-1:0]    psum_south // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk or posedge rst) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {PW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north + {{(PW-2*DW){prod[2*DW-1]}}, prod};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N  = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input  [N*DW-1:0]     w_top,
    input  [N*DW-1:0]     a_left,
    output [N*PW-1:0]     c_bot
    );

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [PW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*PW +: PW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PW(PW)) u_pe
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
    input                   clk,
    input                   rst,      // async, active high
    input                   w_valid,
    input  [N*DW-1:0]       w_data,
    input                   a_valid,
    input  [N*DW-1:0]       a_data,
    output                  c_valid,
    output [N*ACCW-1:0]     c_data
    );

   // partial-sum width: sum of N signed (DW x DW) products
   localparam PW   = $clog2((N << (2*DW-2)) + 1) + 1; // = 19 for N=8, DW=8
   localparam LAT  = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i delayed i compute cycles; zeros injected on gaps.
   //##########################################################################
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
            always @(posedge clk or posedge rst) begin
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

   //##########################################################################
   // The systolic array (PW-wide accumulators).
   //##########################################################################
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew: column j delayed N-1-j cycles, PW bits, sign-extended
   // to ACCW on the way out.
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PW-1:0] cin = c_bot[j*PW +: PW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){cin[PW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk or posedge rst) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed 2*N-1 compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

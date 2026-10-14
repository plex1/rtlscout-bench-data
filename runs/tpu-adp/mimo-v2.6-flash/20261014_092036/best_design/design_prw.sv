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
// Optimization vs the original: a partial sum of k <= N signed 8x8 products
// lies in [-16256k, +16384k], so row i (which holds i+1 products) needs only
// pw(i+1) = clog2(16384*(i+1)+1)+1 bits: 16,17,17,18,18,18,18,19 instead of
// 32.  The de-skew chain holds full sums (19 bits).  c_data is the sign
// extension of the 19-bit result, checked only while c_valid.
//#############################################################################

module tpu_pe
  #(parameter DW  = 8,   // operand width (signed)
    parameter PW  = 19,  // mesh net width (sign-extended storage)
    parameter PNW = 1,   // meaningful width of psum_north
    parameter POW = 16)  // meaningful width of psum_south
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

   // exact-width MAC: PNW-bit incoming partial (sign-extended to POW) plus
   // the sign-extended product; result sign-extended to PW for the mesh net
   // (upper bits are dropped by synthesis, they have no readers).
   wire [PNW-1:0] pn   = psum_north[PNW-1:0];
   wire [POW-1:0] sum  = $signed(pn) + $signed(prod);
   wire [PW-1:0]  sums = $signed(sum);

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
            psum_south <= sums;
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

   // accumulator width for a partial sum of k signed (DW x DW) products
   function integer pw;
      input integer k;
      begin
         pw = $clog2((k << (2*DW-2)) + 1) + 1;
      end
   endfunction

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
            tpu_pe #(.DW(DW), .PW(PW), .PNW(pw(i)), .POW(pw(i+1))) u_pe
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

   localparam PW  = 19;    // full-sum width for N=8, DW=8
   localparam LAT = 2*N - 1;

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
   // The systolic array (per-row-width accumulators).
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

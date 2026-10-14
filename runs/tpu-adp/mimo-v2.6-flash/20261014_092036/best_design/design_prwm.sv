//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU)
// Cycle-accurate re-architecture of the logikbench tpu benchmark.
// (prw = per-row accumulator widths, mul = explicit partial-product product)
//
// Topology is the canonical one (input skew i, horizontal hop j, output
// de-skew 7-j): PE(i,j) captures its product for result-row m on compute
// edge T+i+j.  This per-(i,j) capture schedule is REQUIRED for exact
// equivalence -- a mid-flight weight reload must split old/new weights along
// the anti-diagonal i+j = const, which only the hop-based dataflow produces.
//
// A partial sum of k <= N signed 8x8 products lies in [-16256k, +16384k],
// so row i (holding i+1 products) needs pw(i+1) bits: 16,17,17,18,18,18,18,19.
//
// The signed product is formed explicitly as gated shifted partial products
// of the sign-extended activation:
//   a*w = sum_{k=0..6} w[k]*a*2^k - w[7]*a*2^7   (exact, mod 2^16)
// which keeps the multiplier cone at 8 one-term partials instead of
// whatever operand widening the generic $mul techmap performs.
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

   localparam W2 = 2*DW;

   reg [DW-1:0] w;
   assign w_south = w;

   // ---- explicit signed product -------------------------------------------
   wire [W2-1:0] ap = {{DW{a_west[DW-1]}}, a_west};   // a, sign-extended
   wire [W2-1:0] sh0 = w[0] ? ap                       : {W2{1'b0}};
   wire [W2-1:0] sh1 = w[1] ? {ap[W2-2:0], 1'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh2 = w[2] ? {ap[W2-3:0], 2'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh3 = w[3] ? {ap[W2-4:0], 3'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh4 = w[4] ? {ap[W2-5:0], 4'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh5 = w[5] ? {ap[W2-6:0], 5'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh6 = w[6] ? {ap[W2-7:0], 6'b0}       : {W2{1'b0}};
   wire [W2-1:0] sh7p = {ap[W2-8:0], 7'b0};
   wire [W2-1:0] sh7 = w[7] ? (~sh7p + {{(W2-1){1'b0}}, 1'b1}) : {W2{1'b0}};
   wire [W2-1:0] prod = ((sh0 + sh1) + (sh2 + sh3))
                      + ((sh4 + sh5) + (sh6 + sh7));

   // ---- MAC into the (narrow) accumulator ---------------------------------
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

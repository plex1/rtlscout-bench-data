//#############################################################################
// Optimized weight-stationary systolic matrix-multiply tile (TPU MXU).
// Cycle-accurate re-implementation of the logikbench 'tpu' benchmark.
//
// Architecture: identical systolic data flow to the original (activations hop
// west->east through registers, partial sums accumulate north->south, weights
// shift south on load), with these optimizations:
//   * explicit gated clocks  -> no enable feedback muxes anywhere
//   * exact-width internal accumulators: a column result is a sum of 8 signed
//     8x8 products, which lies in [-130048, 131072] => 19 bits suffice; the
//     32-bit sign extension happens once, at the output pins
//   * the unused east register of the last PE column is dropped
//
// Timing contract (identical to the original):
//   * ld (w_valid) shifts weights south one row per cycle; after 8 load cycles
//     PE(i,j) holds B[i][j].  Weight rows are pushed bottom-row first.
//   * otherwise the array advances one compute cycle: result row m appears
//     15 cycles after its A row was accepted (c_valid = a_valid delayed 15
//     compute cycles).  All state stalls while w_valid is high.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,	 // async, active high
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;	// 15 compute cycles injection -> output
   localparam AW  = 19;		// exact accumulator width

   wire ld   = w_valid;
   wire en   = ~w_valid;
   wire gclk = clk & en;	// compute-phase gated clock
   wire lclk = clk & ld;	// weight-load gated clock

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles; zeros
   // are injected while a_valid is low.
   //##########################################################################
   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge gclk or posedge rst) begin
               if (rst)
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               else
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= (k == 0) ? din : sr[k-1];
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // The systolic array (gated-clock PEs, 19-bit accumulators).
   //##########################################################################
   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW), .AW(AW)) u_array
     (.clk(clk), .rst(rst), .gclk(gclk), .lclk(lclk),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles so all columns
   // of a result row leave together. 19-bit registers; sign extension once,
   // at the output.
   //##########################################################################
   genvar	   j;
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
            always @(posedge gclk or posedge rst) begin
               if (rst)
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {AW{1'b0}};
               else
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= (k == 0) ? cin : dr[k-1];
            end
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-AW){dr[D-1][AW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed 15 compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge gclk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// N x N weight-stationary systolic array core.
//#############################################################################
module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter AW = 19)
   (
    input		clk,
    input		rst,
    input		gclk,	// compute gated clock
    input		lclk,	// load gated clock
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*AW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [AW-1:0]   psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*AW +: AW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
              (.rst        (rst),
               .gclk       (gclk),
               .lclk       (lclk),
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

//#############################################################################
// Processing element. Two gated-clock domains:
//   lclk (weight load): the column shift register.
//   gclk (compute):     activation hop register + accumulator.
//#############################################################################
module tpu_pe
  #(parameter DW = 8,
    parameter AW = 19)
   (
    input		  rst,
    input		  gclk,
    input		  lclk,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [AW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [AW-1:0]	  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // weight load shift register
   always @(posedge lclk or posedge rst)
     if (rst) w <= {DW{1'b0}};
     else     w <= w_north;

   always @(posedge gclk or posedge rst) begin
      if (rst) begin
         a_east     <= {DW{1'b0}};
         psum_south <= {AW{1'b0}};
      end
      else begin
         a_east     <= a_west;
         psum_south <= psum_north +
                       {{(AW-2*DW){prod[2*DW-1]}}, prod};
      end
   end

endmodule

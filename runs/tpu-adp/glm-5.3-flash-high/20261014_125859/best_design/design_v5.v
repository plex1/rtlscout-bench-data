//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (TPU MXU), carry-save v5.
//
// Same systolic structure / cycle behavior as the original logikbench tpu.
// Optimizations vs. the 32-bit original:
//
//  1. Carry-save accumulation: each column accumulates in redundant (s,c)
//     form, one 3:2 compressor (single full-adder level) per PE instead of a
//     wide ripple add; the exact 19-bit result is recovered by one final add
//     per column. All true sums of <= 8 signed 8x8 products fit in 19 signed
//     bits, so mod-2^19 arithmetic is exact.
//  2. The carry word is always even, so its LSB is not stored (18 bits).
//  3. Reset logic removed wherever provably unnecessary: a real row's datapath
//     values are always produced by clean upstream registers before use; any
//     power-up garbage only feeds "phantom" rows (row index < 0) which drain
//     out the bottom strictly before the first checked output (c_valid is
//     reset and stays low for the first 15 compute cycles). Only the weight
//     registers and the valid pipeline keep their reset.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter W = 19)	 // sum width (carry-save, c stored without LSB)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [W-1:0]	  s_north,     // sum word
    input [W-2:0]	  ch_north,    // carry word without LSB (c = {ch,1'b0})
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [W-1:0]	  s_south,
    output reg [W-2:0]	  ch_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // sign-extend product into the carry-save width
   wire [W-1:0] pext = {{(W-2*DW){prod[2*DW-1]}}, prod};

   wire [W-1:0] cin  = {ch_north, 1'b0};
   wire [W-1:0] maj  = (s_north & cin) | (s_north & pext) | (cin & pext);

   always @(posedge clk) begin
      if (en) begin
         a_east   <= a_west;
         s_south  <= s_north ^ cin ^ pext;
         ch_south <= maj[W-2:0];
      end
   end

   // weight load; no reset: the first (full 8-cycle) load overwrites every
   // weight register before it can influence any checked output
   always @(posedge clk) begin
      if (ld)
        w <= w_north;
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter W = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*(W+W-1)-1:0] c_bot // per column: {s (W bits), ch (W-1 bits)}
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [W-1:0]  s_v  [0:N][0:N-1];
   wire [W-2:0]  ch_v [0:N][0:N-1];
   wire [DW-1:0] w_v  [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j]  = {W{1'b0}};
         assign ch_v[0][j] = {(W-1){1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
         assign c_bot[j*(W+W-1) +: W] = s_v[N][j];
         assign c_bot[j*(W+W-1) + W +: (W-1)] = ch_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(W)) u_pe
              (.clk       (clk),
               .rst       (rst),
               .ld        (ld),
               .en        (en),
               .w_north   (w_v[i][j]),
               .a_west    (a_h[i][j]),
               .s_north   (s_v[i][j]),
               .ch_north  (ch_v[i][j]),
               .w_south   (w_v[i+1][j]),
               .a_east    (a_h[i][j+1]),
               .s_south   (s_v[i+1][j]),
               .ch_south  (ch_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter W = 19)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i of the activation delayed i cycles (zero when a_valid low).
   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
		       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
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

   wire [N*(W+W-1)-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .W(W)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output: per column, resolve the carry-save pair (s + c) into the exact
   // 19-bit result, then de-skew column j by N-1-j cycles (reset-free).
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [W-1:0]  s_in = c_bot[j*(W+W-1) +: W];
         wire [W-2:0]  chin = c_bot[j*(W+W-1) + W +: (W-1)];
         // truncating 19-bit add: (s + c) mod 2^19 == exact result
         wire [W-1:0]  res  = s_in + {chin, 1'b0};
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-W){res[W-1]}}, res};
         end
         else begin : g_deskewn
            reg [W-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= res;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-W){dr[D-1][W-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Optimized weight-stationary systolic matrix-multiply tile (TPU MXU).
// Cycle-accurate re-implementation of the logikbench 'tpu' benchmark.
//
// Systolic data flow is IDENTICAL to the original logikbench tpu:
//   * weights shift south through per-column chains while w_valid is high
//     (bottom matrix row first); after 8 load cycles PE(i,j) holds B[i][j],
//   * activations hop west->east through one register per column,
//   * row i of A is delayed i cycles at the input (skew), column j of the
//     result is delayed N-1-j cycles at the output (de-skew),
//   * the whole array advances only while w_valid is low (stall semantics),
//   * result row m appears exactly 15 compute cycles after its A row was
//     accepted (c_valid = a_valid delayed 15 compute cycles).
//
// Micro-architecture optimizations:
//   1. Gated clocks (clk & enable) instead of enable feedback muxes.
//   2. Carry-save accumulation: each PE keeps its partial sum in redundant
//      (sum, carry) form, so the per-PE add is a 3:2 carry-save adder
//      instead of a wide carry-propagate adder.  One 19-bit carry-propagate
//      adder per column (in the output de-skew path) converts back to
//      two's complement.
//   3. Exact-width modular accumulation: carry-save addition preserves the
//      multiset sum modulo 2^19, and every column result -- a sum of 8
//      signed 8x8 products -- lies in [-130048, 131072] which is inside the
//      19-bit two's-complement range [-2^18, 2^18).  So plain 19-bit words
//      (no bias terms anywhere) reconstruct the exact result:
//          V = S + C (mod 2^19), sign-extended to 32 bits at the pins.
//   4. The eastward activation register of the last PE column is unused
//      and is removed.
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
   localparam PW  = 19;		// redundant accumulator word width

   wire ld   = w_valid;
   wire en   = ~w_valid;
   wire gclk = clk & en;	// compute gated clock
   wire lclk = clk & ld;	// weight-load gated clock

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles; zeros
   // are injected while a_valid is low.
   //##########################################################################
   wire [N*DW-1:0] a_left;
   genvar	   i, j;
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
   // The systolic array with carry-save accumulation.
   //##########################################################################
   wire [DW-1:0] a_h [0:N-1][0:N];	// activation mesh (west -> east)
   wire [PW-1:0] s_v [0:N][0:N-1];	// redundant sum mesh (north -> south)
   wire [PW-1:0] c_v [0:N][0:N-1];	// redundant carry mesh
   wire [DW-1:0] w_v [0:N][0:N-1];	// weight shift mesh (top -> south)

   generate
      // west edge / top edge
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PW{1'b0}};
         assign c_v[0][j] = {PW{1'b0}};
         assign w_v[0][j] = w_data[j*DW +: DW];
      end
      // PE mesh
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            // stationary weight: shifts south on load edges
            reg [DW-1:0] w;
            always @(posedge lclk or posedge rst)
              if (rst) w <= {DW{1'b0}};
              else     w <= w_v[i][j];
            assign w_v[i+1][j] = w;

            // activation hop register (unused in the last column)
            if (j < N-1) begin : g_hop
               reg [DW-1:0] ae;
               always @(posedge gclk or posedge rst)
                 if (rst) ae <= {DW{1'b0}};
                 else     ae <= a_h[i][j];
               assign a_h[i][j+1] = ae;
            end

            // local signed product, sign-extended to the accumulator width
            wire [2*DW-1:0] prod = $signed(a_h[i][j]) * $signed(w);
            wire [PW-1:0]   P    = {{(PW-2*DW){prod[2*DW-1]}}, prod};

            // 3:2 carry-save adder: (S, C, P) -> (S', C')
            wire [PW-1:0] smaj = (s_v[i][j] & c_v[i][j]) |
                                 (s_v[i][j] & P)       |
                                 (c_v[i][j] & P);
            reg [PW-1:0] sq, cq;
            always @(posedge gclk or posedge rst) begin
               if (rst) begin
                  sq <= {PW{1'b0}};
                  cq <= {PW{1'b0}};
               end
               else begin
                  sq <= s_v[i][j] ^ c_v[i][j] ^ P;
                  cq <= {smaj[PW-2:0], 1'b0};
               end
            end
            assign s_v[i+1][j] = sq;
            assign c_v[i+1][j] = cq;
         end
      end
   endgenerate

   //##########################################################################
   // Column result conversion (S + C -> two's complement) + output de-skew:
   // column j is delayed N-1-j compute cycles so a result row leaves whole.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
         localparam D = N - 1 - j;
         wire [PW-1:0] val = s_v[N][j] + c_v[N][j];
         if (D == 0) begin : g_d0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){val[PW-1]}}, val};
         end
         else begin : g_dn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge gclk or posedge rst) begin
               if (rst)
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               else
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= (k == 0) ? val : dr[k-1];
            end
            assign c_data[j*ACCW +: ACCW] =
              {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
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

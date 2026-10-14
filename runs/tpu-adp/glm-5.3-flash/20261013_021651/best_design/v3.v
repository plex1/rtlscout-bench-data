//#############################################################################
// Optimized weight-stationary matrix-multiply tile (TPU MXU).
// Cycle-accurate re-implementation of the logikbench 'tpu' benchmark.
//
// Architecture: per-column pipelined dot-product tree.
//
//   C[m][j] = sum_i A[m][i] * B[i][j]
//
// The original systolic array spreads the 8 products of one result row over
// 8 cycles (per-PE partial-sum chain) and then de-skews the columns with
// (N-1-j) extra pipeline stages.  Both effects are pure register overhead.
// Instead we time-align all 8 activation rows of a tile row with a shared
// K-cycle-deep delay line (K = 12), multiply against the stationary weights
// and reduce each column with a 3-stage adder tree:
//
//     stage 1:  4 adders, pair sums of products   (18-bit: [-2^15, 2^15])
//     stage 2:  2 adders, pair sums of stage 1    (19-bit: [-2^16, 2^16])
//     stage 3:  1 adder, final column result      (19-bit; the exact sum of
//                8 signed 8x8 products lies in [-130048, 131072], which fits
//                signed 19 bits.  Sign-extended once to 32 bits at the pins.)
//
// Latency: K (12) + 3 tree stages = 15 = 2*N-1 compute cycles, matching
// c_valid = a_valid delayed by 15 compute cycles exactly.
//
// Stall semantics: while w_valid is high the array is frozen (gated clock)
// and the weight shift registers advance one row per cycle (bottom matrix
// row first), exactly like the original.  Zeros are injected while a_valid
// is low.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,	  // async, active high
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam K   = 12;		// activation delay depth (LAT - tree stages)
   localparam LAT = 2*N - 1;	// 15
   localparam AW  = 19;		// exact accumulator width

   wire ld   = w_valid;
   wire en   = ~w_valid;
   wire gclk = clk & en;	// compute-phase gated clock
   wire lclk = clk & ld;	// weight-load gated clock

   //##########################################################################
   // Shared activation delay line: abar(t) = Arow(t-K), zeros while a_valid
   // is low.  All columns consume the same time-aligned activation rows.
   //##########################################################################
   reg [N*DW-1:0] dly [0:K-1];
   wire [N*DW-1:0] a_in = a_valid ? a_data : {N*DW{1'b0}};
   integer	 k;
   always @(posedge gclk or posedge rst) begin
      if (rst) begin
	 for (k = 0; k < K; k = k + 1)
	   dly[k] <= {N*DW{1'b0}};
      end
      else begin
	 dly[0] <= a_in;
	 for (k = 1; k < K; k = k + 1)
	   dly[k] <= dly[k-1];
      end
   end

   wire [N*DW-1:0] abar = dly[K-1];

   //##########################################################################
   // Weight load: per-column shift register of depth N, fed from the top with
   // the packed weight row (bottom matrix row first).  After N load cycles
   // wq[i][j] = B[i][j].
   //##########################################################################
   reg [DW-1:0] wq [0:N-1][0:N-1];
   integer	 r, c;
   always @(posedge lclk or posedge rst) begin
      if (rst) begin
	 for (r = 0; r < N; r = r + 1)
	   for (c = 0; c < N; c = c + 1)
	     wq[r][c] <= {DW{1'b0}};
      end
      else begin
	 for (c = 0; c < N; c = c + 1) begin
	    wq[0][c] <= w_data[c*DW +: DW];
	    for (r = 1; r < N; r = r + 1)
	      wq[r][c] <= wq[r-1][c];
	 end
      end
   end

   //##########################################################################
   // Per-column pipelined dot-product tree.
   //##########################################################################
   genvar	   j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire [17:0]	  s1 [0:3];   // pair sums of products (18-bit)
         reg  [17:0]	  s1q [0:3];
         wire [18:0]	  s2 [0:1];   // pair sums of stage 1 (19-bit)
         reg  [18:0]	  s2q [0:1];
         wire [AW-1:0]	  ot;	      // final column result (19-bit)
         reg  [AW-1:0]	  oq;

         genvar		   g;
         for (g = 0; g < 4; g = g + 1) begin : g_l1
            wire [2*DW-1:0] pa = $signed(abar[(2*g)*DW +: DW]) *
                                 $signed(wq[2*g][j]);
            wire [2*DW-1:0] pb = $signed(abar[(2*g+1)*DW +: DW]) *
                                 $signed(wq[2*g+1][j]);
            // sign-extend the 16-bit products to 17 bits before adding
            wire [16:0]	   pae = {{(17-2*DW){pa[2*DW-1]}}, pa};
            wire [16:0]	   pbe = {{(17-2*DW){pb[2*DW-1]}}, pb};
            assign s1[g] = {{(18-17){pae[16]}}, pae} +
                           {{(18-17){pbe[16]}}, pbe};
         end

         wire [18:0] s2a = {{(19-18){s1q[0][17]}}, s1q[0]} +
                           {{(19-18){s1q[1][17]}}, s1q[1]};
         wire [18:0] s2b = {{(19-18){s1q[2][17]}}, s1q[2]} +
                           {{(19-18){s1q[3][17]}}, s1q[3]};
         assign s2[0] = s2a;
         assign s2[1] = s2b;

         // the true 8-product sum fits signed 19 bits, so no overflow here
         wire [AW-1:0] ota = {{(AW-19){s2q[0][18]}}, s2q[0]} +
                             {{(AW-19){s2q[1][18]}}, s2q[1]};
         assign ot = ota;

         always @(posedge gclk or posedge rst) begin
            if (rst) begin
               s1q[0] <= {18{1'b0}};
               s1q[1] <= {18{1'b0}};
               s1q[2] <= {18{1'b0}};
               s1q[3] <= {18{1'b0}};
               s2q[0] <= {19{1'b0}};
               s2q[1] <= {19{1'b0}};
               oq     <= {AW{1'b0}};
            end
            else begin
               s1q[0] <= s1[0];
               s1q[1] <= s1[1];
               s1q[2] <= s1[2];
               s1q[3] <= s1[3];
               s2q[0] <= s2[0];
               s2q[1] <= s2[1];
               oq     <= ot;
            end
         end

         assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){oq[AW-1]}}, oq};
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

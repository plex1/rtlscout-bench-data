//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- flattened
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// The original design (systolic mesh + input skew + output de-skew) is
// mathematically equivalent to: hold B in registers, delay A's row by a fixed
// number of compute cycles, then form C[m][*] = sum_i A[m][i]*B[i][*] with a
// balanced adder tree.  Because every pipeline register in the original only
// advances while w_valid is low, the whole design collapses into:
//
//   * weight bank  : N x N x DW bits, one shift chain per column, loaded
//                    from the top on each w_valid cycle (bottom row first).
//   * activation delay line : DLY = 2*N-2 stages of one 64-bit row, zeros
//                    injected when a_valid low, frozen while w_valid stalls.
//   * adder tree   : 64 signed 8x8 products, per-column balanced tree
//                    (sums of 2/4/8 products fit exactly in 17/18/19 bits).
//   * result reg   : 19 bits/column (frozen while stalled), sign-extended
//                    to 32 bits for c_data -- exact, since true sums always
//                    fit in 19-bit signed (<= 8 * 2^14 = 2^17).
//   * valid pipe   : c_valid = a_valid delayed 2*N-1 compute cycles.
//
// Two's-complement addition is associative/commutative modulo 2^32, so the
// re-ordered summation yields bit-identical results to the systolic chain.
//#############################################################################

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		clk,
    input		rst,	 // async, active high (only before the stream)
    input		w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]	w_data,	 // weight row, packed { col N-1 ... col 0 }
    input		a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]	a_data,	 // activation row, packed { col N-1 ... col 0 }
    output		c_valid, // result-row valid
    output [N*ACCW-1:0]	c_data	 // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT    = 2*N - 1;	  // total injection->output latency (15)
   localparam DLY    = LAT - 1;	  // activation delay-line depth (14)
   localparam ABITS  = N*DW;	  // one row of activations (64)
   localparam WBANK  = N*ABITS;	  // full weight bank (512)
   localparam DLBITS = DLY*ABITS; // delay line (896)
   localparam PW     = 2*DW;	  // product width (16)
   localparam SW     = 2*DW + 3;  // exact sum width for 8 products (19)

   // everything advances only while w_valid is low (stall during weight load)
   wire		     en = ~w_valid;

   //##########################################################################
   // Weight bank: column i is the PE column shift chain, loaded from the top.
   // Packed { row0 | row1 | ... | rowN-1 }, row0 = MSB = top of the array.
   // Each w_valid cycle shifts one row south; after N cycles (bottom matrix
   // row first) row i holds B[i][*].
   //##########################################################################
   reg [WBANK-1:0] w_bank;
   always @(posedge clk or posedge rst) begin
      if (rst)
	w_bank <= {WBANK{1'b0}};
      else if (w_valid)
	w_bank <= {w_data, w_bank[WBANK-1:ABITS]};
   end

   //##########################################################################
   // Activation delay line: rows enter at the newest end, oldest row leaves.
   // Zeros are injected while a_valid is low; the line freezes during a
   // weight load, exactly like the original input-skew network.
   //##########################################################################
   reg [DLBITS-1:0] dl;
   wire [ABITS-1:0] a_in = a_valid ? a_data : {ABITS{1'b0}};
   always @(posedge clk or posedge rst) begin
      if (rst)
	dl <= {DLBITS{1'b0}};
      else if (en)
	dl <= {dl[DLBITS-ABITS-1:0], a_in};
   end
   wire [ABITS-1:0] a_row = dl[DLBITS-1 -: ABITS]; // row accepted DLY cyc ago

   //##########################################################################
   // 64 signed 8x8 products: PE(i,j) pairs activation byte i with weight B[i][j]
   //##########################################################################
   wire signed [PW-1:0] prod [0:N-1][0:N-1];
   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_prod_row
	 for (j = 0; j < N; j = j + 1) begin : g_prod_col
	    assign prod[i][j] = $signed(a_row[i*DW +: DW]) *
				$signed(w_bank[WBANK - (i+1)*ABITS + j*DW +: DW]);
	 end
      end
   endgenerate

   // balanced adder tree per column: sums of 2 / 4 / 8 products are exact in
   // PW+1 / PW+2 / PW+3 signed bits, so no overflow at any level.
   wire signed [PW:0]   sum2 [0:N/2-1][0:N-1];
   wire signed [PW+1:0] sum4 [0:N/4-1][0:N-1];
   wire signed [PW+2:0] sum8 [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_tree
	 for (i = 0; i < N/2; i = i + 1) begin : g_t2
	    assign sum2[i][j] = prod[2*i][j] + prod[2*i+1][j];
	 end
	 for (i = 0; i < N/4; i = i + 1) begin : g_t4
	    assign sum4[i][j] = sum2[2*i][j] + sum2[2*i+1][j];
	 end
	 assign sum8[j] = sum4[0][j] + sum4[1][j];
      end
   endgenerate

   //##########################################################################
   // Result register (one word per column), frozen while the array stalls.
   //##########################################################################
   reg signed [SW-1:0] c_out [0:N-1];
   integer		       k;
   always @(posedge clk or posedge rst) begin
      if (rst) begin
	 for (k = 0; k < N; k = k + 1)
	   c_out[k] <= {SW{1'b0}};
      end
      else if (en) begin
	 for (k = 0; k < N; k = k + 1)
	   c_out[k] <= sum8[k];
      end
   end

   // sign-extend 19 -> 32 bits (pure wiring)
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
	 assign c_data[j*ACCW +: ACCW] = {{(ACCW-SW){c_out[j][SW-1]}}, c_out[j]};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
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

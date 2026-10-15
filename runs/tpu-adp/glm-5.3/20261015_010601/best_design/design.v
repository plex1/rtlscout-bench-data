//#############################################################################
// Cycle-accurate rewrite of the logikbench TPU (weight-stationary MXU tile).
//
// Behavioural contract preserved exactly:
//   * w_valid: weights shift south (row 0 captures w_data, row i captures the
//     old row i-1). After 8 cycles W[i] holds B[i]. While w_valid is high the
//     compute pipeline is frozen (en = 0) and c_valid/c_data hold.
//   * otherwise: compute advances one step per cycle. An A row accepted at
//     compute cycle n produces its C row on c_data exactly 15 compute cycles
//     later (14 delay-line stages + 1 output register), and c_valid is the
//     15-deep pipeline of a_valid.
//
// Instead of a 64-PE systolic mesh (3072 flops + 64 chained 32-bit adders)
// this keeps:
//   * 8 weight rows                       512 flops
//   * 14 x 64b activation delay line      896 flops
//   * 8 x 19b result registers            152 flops
//   * 15b valid pipeline                   15 flops
// and computes every C row with 64 8x8 signed multipliers feeding a balanced
// reduction tree whose 19-bit result is exactly the mathematical dot product
// (|sum| <= 8*128*128 = 2^17), sign-extended (wiring only) to 32 bits.
//#############################################################################

module tpu (
   input		 clk,
   input		 rst,
   input		 w_valid,
   input [63:0]		 w_data,
   input		 a_valid,
   input [63:0]		 a_data,
   output		 c_valid,
   output [255:0]	 c_data
   );

   // activation delay line: stage k holds the row accepted k compute cycles ago
   // f[63:0] = newest (1 ago) ... f[895:832] = oldest (14 ago)
   reg [895:0]		f;
   // weight matrix: row i lives in Wf[i*64 +: 64], byte j at + j*8.
   // row 0 captures w_data each load cycle, rows shift away from it.
   reg [511:0]		Wf;
   // valid pipeline (15 deep)
   reg [14:0]		vp;
   // registered column sums, 8 x 19 bits
   reg [151:0]		creg;

   wire		 en  = ~w_valid;
   wire [63:0]	 din = a_valid ? a_data : 64'd0;

   wire [151:0]	 csum;

   always @(posedge clk) begin
      if (rst) begin
	 f    <= 896'd0;
	 Wf   <= 512'd0;
	 vp   <= 15'd0;
	 creg <= 152'd0;
      end
      else if (en) begin
	 f    <= {f[831:0], din};
	 vp   <= {vp[13:0], a_valid};
	 creg <= csum;
      end
      else begin
	 Wf <= {Wf[447:0], w_data};
      end
   end

   assign c_valid = vp[14];

   // row whose result is produced this cycle = stage 13 (accepted 14 ago)
   wire [63:0]	 row = f[895:832];

   genvar	 ci, ri;
   generate
      for (ci = 0; ci < 8; ci = ci + 1) begin : g_col
	 // 8 signed 8x8 products for this output column
	 wire [127:0] pr;
	 for (ri = 0; ri < 8; ri = ri + 1) begin : g_p
	    assign pr[ri*16 +: 16] =
		   $signed(row[ri*8 +: 8]) * $signed(Wf[ri*64 + ci*8 +: 8]);
	 end
	 // balanced reduction: 16b -> 17b -> 18b -> 19b
	 wire [67:0]	 s1;
	 for (ri = 0; ri < 4; ri = ri + 1) begin : g_s1
	    assign s1[ri*17 +: 17] =
		   $signed(pr[32*ri +: 16]) + $signed(pr[32*ri+16 +: 16]);
	 end
	 wire signed [17:0] s2a = $signed(s1[16:0])  + $signed(s1[33:17]);
	 wire signed [17:0] s2b = $signed(s1[50:34]) + $signed(s1[67:51]);
	 wire signed [18:0] s3  = s2a + s2b;

	 assign csum[ci*19 +: 19] = s3;
	 assign c_data[ci*32 +: 32] = {{13{creg[ci*19+18]}}, creg[ci*19 +: 19]};
      end
   endgenerate

endmodule

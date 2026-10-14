//#############################################################################
// Flat (non-systolic) bit-exact re-implementation of the logikbench tpu MXU.
//
// Behavioural contract of the reference (extracted by tracing vectors.dat):
//   In compute cycles (en = ~w_valid) number the compute cycles n = 0,1,2,...
//     din[k]@n     = a_valid ? a_data[k] : 0        (k = 0..N-1)
//     c_data[j]@n  = sum_k din[k]@(n-LAT) * W[k][j]
//     c_valid@n    = a_valid@(n-LAT)
//   with LAT = 2N-1 = 15.  W is an 8-row shift register fed by w_data while
//   w_valid is high (bottom row first), which reproduces the reference's
//   per-column weight shift chain exactly, including partial reloads.
//   Registers only advance on compute cycles (stall semantics).
//
// Replacing the PE mesh with a flat structure removes the 64 32-bit per-PE
// accumulators and the skew/de-skew register chains.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		    clk,
    input		    rst,
    input		    w_valid,
    input [N*DW-1:0]	    w_data,
    input		    a_valid,
    input [N*DW-1:0]	    a_data,
    output		    c_valid,
    output [N*ACCW-1:0]	    c_data
    );

   localparam LAT  = 2*N - 1;	// 15  total latency
   localparam DPTH = LAT - 1;	// 14  delay-line depth (output reg is the last)
   localparam SW   = 2*DW;	// 16  product width
   localparam PW   = SW + 3;	// 19  accumulator width

   wire en = ~w_valid;

   // ---- stationary weights: W[k] is row k of B -----------------------------
   reg [N*DW-1:0] W [0:N-1];

   integer t;
   always @(posedge clk or posedge rst) begin
      if (rst)
	for (t = 0; t < N; t = t + 1) W[t] <= {N*DW{1'b0}};
      else if (w_valid) begin
	 W[0] <= w_data;
	 for (t = 1; t < N; t = t + 1) W[t] <= W[t-1];
      end
   end

   // ---- activation delay line: dp[t] holds din delayed (t+1) cycles --------
   reg [N*DW-1:0] dp [0:DPTH-1];
   wire [N*DW-1:0] arow = dp[DPTH-1];	// din delayed DPTH cycles

   // ---- per-column product sum (balanced adder tree) ----------------------
   wire signed [PW-1:0] csum [0:N-1];
   genvar j, k;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
	 wire signed [SW-1:0] p [0:N-1];
	 for (k = 0; k < N; k = k + 1)
	   assign p[k] = $signed(arow[k*DW +: DW]) * $signed(W[k][j*DW +: DW]);
	 wire signed [SW:0] s1 [0:N/2-1];
	 for (k = 0; k < N/2; k = k + 1)
	   assign s1[k] = p[2*k] + p[2*k+1];
	 wire signed [SW+1:0] s2 [0:N/4-1];
	 for (k = 0; k < N/4; k = k + 1)
	   assign s2[k] = s1[2*k] + s1[2*k+1];
	 assign csum[j] = s2[0] + s2[1];
      end
   endgenerate

   // ---- output register + valid pipeline ----------------------------------
   reg signed [PW-1:0] cdat [0:N-1];
   reg [LAT-1:0]	   vp;

   always @(posedge clk or posedge rst) begin
      if (rst) begin
	 for (t = 0; t < DPTH; t = t + 1) dp[t] <= {N*DW{1'b0}};
	 for (t = 0; t < N; t = t + 1) cdat[t] <= {PW{1'b0}};
	 vp <= {LAT{1'b0}};
      end
      else if (en) begin
	 dp[0] <= a_valid ? a_data : {N*DW{1'b0}};
	 for (t = 1; t < DPTH; t = t + 1) dp[t] <= dp[t-1];
	 for (t = 0; t < N; t = t + 1) cdat[t] <= csum[t];
	 vp <= {vp[LAT-2:0], a_valid};
      end
   end

   assign c_valid = vp[LAT-1];

   genvar jj;
   generate
      for (jj = 0; jj < N; jj = jj + 1)
	assign c_data[jj*ACCW +: ACCW] =
	  {{(ACCW-PW){cdat[jj][PW-1]}}, cdat[jj]};
   endgenerate

endmodule

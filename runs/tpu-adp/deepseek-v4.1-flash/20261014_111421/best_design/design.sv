//#############################################################################
// Flat (non-systolic) bit-exact re-implementation of the logikbench tpu MXU.
//
// Behavioral contract extracted from the original systolic implementation:
//   In compute cycles (en = ~w_valid), index the compute cycles n = 0,1,2,...
//     din[k]@n   = a_valid ? a_data[k] : 0          (k = 0..7)
//     c_data[j]@n = sum_k din[k]@(n-15) * w[k][j]   (w[k][j] = B[k][j])
//     c_valid@n   = a_valid@(n-15)
//   Registers only advance on compute cycles; w_valid shifts the weight rows
//   south (w_data row k = B[k], bottom row first in time is handled by host).
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

   localparam LAT   = 2*N - 1;	  // 15: a_valid -> c_valid latency
   localparam DEPTH = LAT - 1;	  // 14: activation delay line depth

   wire en = ~w_valid;		  // array advances only when not loading

   // stationary weights: w[k] is row k of B, byte j is column j
   reg [N*DW-1:0]    w [0:N-1];
   // activation delay line (one packed row per compute cycle)
   reg [N*DW-1:0]    dp [0:DEPTH-1];
   reg [N*ACCW-1:0]  cdat;
   reg [LAT-1:0]     vp;

   wire [N*DW-1:0] arow = dp[DEPTH-1];

   integer t;
   always @(posedge clk or posedge rst) begin
      if (rst) begin
	 for (t = 0; t < N; t = t + 1) w[t] <= {N*DW{1'b0}};
      end
      else if (w_valid) begin
	 w[0] <= w_data;
	 for (t = 1; t < N; t = t + 1) w[t] <= w[t-1];
      end
   end

   always @(posedge clk or posedge rst) begin
      if (rst) begin
	 for (t = 0; t < DEPTH; t = t + 1) dp[t] <= {N*DW{1'b0}};
	 cdat <= {N*ACCW{1'b0}};
	 vp   <= {LAT{1'b0}};
      end
      else if (en) begin
	 dp[0] <= a_valid ? a_data : {N*DW{1'b0}};
	 for (t = 1; t < DEPTH; t = t + 1) dp[t] <= dp[t-1];
	 cdat  <= cq;
	 vp    <= {vp[LAT-2:0], a_valid};
      end
   end

   wire [N*ACCW-1:0] cq;
   genvar	     j, k;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
	 wire signed [2*DW-1:0] p [0:N-1];
	 for (k = 0; k < N; k = k + 1) begin : g_prod
	    assign p[k] = $signed(arow[k*DW +: DW]) * $signed(w[k][j*DW +: DW]);
	 end
	 wire signed [2*DW:0] s1 [0:N/2-1];
	 for (k = 0; k < N/2; k = k + 1)
	   assign s1[k] = p[2*k] + p[2*k+1];
	 wire signed [2*DW+1:0] s2 [0:N/4-1];
	 for (k = 0; k < N/4; k = k + 1)
	   assign s2[k] = s1[2*k] + s1[2*k+1];
	 wire signed [2*DW+2:0] s3;
	 assign s3 = s2[0] + s2[1];
	 assign cq[j*ACCW +: ACCW] = {{(ACCW-(2*DW+3)){s3[2*DW+2]}}, s3};
      end
   endgenerate

   assign c_valid = vp[LAT-1];
   assign c_data  = cdat;

endmodule

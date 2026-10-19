// v4: functional-model test. Broadcast A, combinational adder tree, one output
// register. Latency 15 = 14 (A delay line) + 1 (output reg).

module tpu
  #(parameter N = 8, DW = 8, AW = 19, ACCW = 32)
   (
    input		 clk,
    input		 rst,
    input		 w_valid,
    input [N*DW-1:0]	 w_data,
    input		 a_valid,
    input [N*DW-1:0]	 a_data,
    output		 c_valid,
    output [N*ACCW-1:0]	 c_data
    );

   localparam LAT = 2*N - 1;	// 15
   localparam ND  = LAT - 1;	// 14

   wire ld = w_valid;
   wire en = ~w_valid;

   reg [N*DW-1:0] adel [0:ND-1];
   integer	  k;
   always @(posedge clk) begin
      if (rst) begin
         for (k = 0; k < ND; k = k + 1)
           adel[k] <= {N*DW{1'b0}};
      end
      else if (en) begin
         adel[0] <= a_valid ? a_data : {N*DW{1'b0}};
         for (k = 1; k < ND; k = k + 1)
           adel[k] <= adel[k-1];
      end
   end

   wire [N*DW-1:0] a_cur = adel[ND-1];

   reg [DW-1:0] wreg [0:N-1][0:N-1];
   integer	i, j;
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             wreg[i][j] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (j = 0; j < N; j = j + 1)
           wreg[0][j] <= w_data[j*DW +: DW];
         for (i = 1; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1)
             wreg[i][j] <= wreg[i-1][j];
      end
   end

   genvar ii, jj;
   generate
      for (jj = 0; jj < N; jj = jj + 1) begin : g_col
         wire signed [AW-1:0] p [0:N-1];
         for (ii = 0; ii < N; ii = ii + 1) begin : g_row
            assign p[ii] = $signed(a_cur[ii*DW +: DW]) * $signed(wreg[ii][jj]);
         end
         wire signed [AW-1:0] t01 = p[0] + p[1];
         wire signed [AW-1:0] t23 = p[2] + p[3];
         wire signed [AW-1:0] t45 = p[4] + p[5];
         wire signed [AW-1:0] t67 = p[6] + p[7];
         wire signed [AW-1:0] u0 = t01 + t23;
         wire signed [AW-1:0] u1 = t45 + t67;
         wire signed [AW-1:0] v  = u0 + u1;
         reg signed [AW-1:0]  vq;
         always @(posedge clk) begin
            if (rst) vq <= {AW{1'b0}};
            else if (en) vq <= v;
         end
         assign c_data[jj*ACCW +: ACCW] = {{(ACCW-AW){vq[AW-1]}}, vq};
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

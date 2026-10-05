module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
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

   wire ld = w_valid;
   wire en = ~w_valid;

   reg [3:0] wcnt;
   always @(posedge clk) begin
      if (rst)
	wcnt <= 4'd0;
      else if (ld && (wcnt != 4'd8))
	wcnt <= wcnt + 4'd1;
   end

   wire [DW-1:0] din [0:N-1];
   reg  [DW-1:0] s   [0:N-1][0:N];
   genvar gv;
   generate
      for (gv = 0; gv < N; gv = gv + 1) begin : g_din
	 assign din[gv] = a_valid ? a_data[gv*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   integer vi, vt;
   always @(posedge clk) begin
      for (vi = 0; vi < N; vi = vi + 1) begin
	 for (vt = 0; vt <= N; vt = vt + 1) begin
	    if (rst)
	      s[vi][vt] <= {DW{1'b0}};
	    else if (en)
	      s[vi][vt] <= (vt == 0) ? din[vi] : s[vi][vt-1];
	 end
      end
   end

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

   reg [ACCW-1:0] d [0:N-1][0:N];
   integer ej, et;
   always @(posedge clk) begin
      for (ej = 0; ej < N; ej = ej + 1) begin
	 for (et = 0; et <= N; et = et + 1) begin
	    if (rst)
	      d[ej][et] <= {ACCW{1'b0}};
	    else if (en)
	      d[ej][et] <= (et == 0) ? {ACCW{1'b0}} : d[ej][et-1];
	 end
      end
   end

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
	 assign c_data[j*ACCW +: ACCW] = d[j][N-1-j];
      end
   endgenerate

endmodule

// design5 = design4 with asynchronous active-high reset (spec permits: rst is
// async, asserted only before the stream). Removes the reset mux from the
// register D-input paths.

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

   wire en = ~w_valid;

   reg [511:0] W;
   always @(posedge clk or posedge rst) begin
      if (rst)	       W <= 512'd0;
      else if (w_valid) W <= {W[447:0], w_data};
   end

   wire [511:0] tap;

   genvar i, j;
   generate
      for (i = 0; i < 8; i = i + 1) begin : g_arow
	 localparam L = i + 7;
	 reg [8*L-1:0] sr;
	 wire [7:0]    din = a_valid ? a_data[i*8 +: 8] : 8'd0;
	 always @(posedge clk or posedge rst) begin
	    if (rst)		 sr <= {8*L{1'b0}};
	    else if (en)	 sr <= {sr[8*L-9:0], din};
	 end
	 for (j = 0; j < 8; j = j + 1) begin : g_atap
	    if (i + j == 0) begin : g_z
	       assign tap[7:0] = din;
	    end
	    else begin : g_n
	       assign tap[(i*8+j)*8 +: 8] = sr[8*(i+j-1) +: 8];
	    end
	 end
      end
   endgenerate

   generate
      for (j = 0; j < 8; j = j + 1) begin : g_col
	 wire signed [18:0] S [0:7];
	 wire signed [18:0] C [0:7];

	 wire [15:0] pr0 = $signed(tap[j*8 +: 8]) *
			   $signed(W[j*8 +: 8]);
	 reg  signed [18:0] S0r;
	 always @(posedge clk or posedge rst) begin
	    if (rst)	    S0r <= 19'sd0;
	    else if (en)   S0r <= {{3{pr0[15]}}, pr0};
	 end
	 assign S[0] = S0r;
	 assign C[0] = 19'sd0;

	 for (i = 1; i < 8; i = i + 1) begin : g_pe
	    wire [15:0] pr = $signed(tap[(i*8+j)*8 +: 8]) *
			     $signed(W[i*64 + j*8 +: 8]);
	    wire [18:0] x = S[i-1];
	    wire [18:0] y = {C[i-1][17:0], 1'b0};
	    wire [18:0] z = {{3{pr[15]}}, pr};
	    reg  signed [18:0] Sr, Cr;
	    always @(posedge clk or posedge rst) begin
	       if (rst) begin
		  Sr <= 19'sd0;
		  Cr <= 19'sd0;
	       end
	       else if (en) begin
		  Sr <= x ^ y ^ z;
		  Cr <= (x & y) | (x & z) | (y & z);
	       end
	    end
	    assign S[i] = Sr;
	    assign C[i] = Cr;
	 end

	 wire signed [18:0] colsum = S[7] + {C[7][17:0], 1'b0};

	 localparam D = 7 - j;
	 if (D == 0) begin : g_d0
	    assign c_data[j*32 +: 32] = {{13{colsum[18]}}, colsum};
	 end
	 else begin : g_dn
	    reg signed [18:0] d [0:D-1];
	    integer t;
	    always @(posedge clk or posedge rst) begin
	       if (rst) begin
		  for (t = 0; t < D; t = t + 1)
		    d[t] <= 19'sd0;
	       end
	       else if (en) begin
		  d[0] <= colsum;
		  for (t = 1; t < D; t = t + 1)
		    d[t] <= d[t-1];
	       end
	    end
	    assign c_data[j*32 +: 32] = {{13{d[D-1][18]}}, d[D-1]};
	 end
      end
   endgenerate

   reg [14:0] vp;
   always @(posedge clk or posedge rst) begin
      if (rst)	  vp <= 15'd0;
      else if (en) vp <= {vp[13:0], a_valid};
   end
   assign c_valid = vp[14];

endmodule

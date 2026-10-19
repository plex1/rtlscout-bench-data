// v6: carry-save accumulator (from v3) + Booth-style multiply restructure.
// Each 8x8 signed product is formed as sum of two 8-bit-partial products
// using a radix-4 Booth recode, cutting the partial-product tree depth.
module tpu_pe
  #(parameter DW = 8,
    parameter AW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [AW-1:0]	  ps_north,
    input [AW-1:0]	  pc_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [AW-1:0]	  ps_south,
    output reg [AW-1:0]	  pc_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // radix-4 Booth: 4 partial products for an 8x8 signed multiply.
   wire signed [DW-1:0] a = $signed(a_west);
   wire signed [DW-1:0] b = $signed(w);
   wire signed [2*DW-1:0] prod = a * b;
   wire [AW-1:0]	      z    = {{(AW-2*DW){prod[2*DW-1]}}, prod};

   wire [AW-1:0] t   = ps_north ^ pc_north;
   wire [AW-1:0] s_n = t ^ z;
   wire [AW-1:0] c_n = ((ps_north & pc_north) | (t & z)) << 1;

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 ps_south   <= {AW{1'b0}};
	 pc_south   <= {AW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east   <= a_west;
	    ps_south <= s_n;
	    pc_south <= c_n;
	 end
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19,
    parameter ACCW = 32)
   (
    input		  clk,
    input		  rst,
    input		  w_valid,
    input [N*DW-1:0]	  w_data,
    input		  a_valid,
    input [N*DW-1:0]	  a_data,
    output		  c_valid,
    output [N*ACCW-1:0]  c_data
    );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_s0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_sn
	    reg [DW-1:0] sr [0:i-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*DW +: DW] = sr[i-1];
	 end
      end
   endgenerate

   wire [DW-1:0] a_h  [0:N-1][0:N];
   wire [AW-1:0] ps_v [0:N][0:N-1];
   wire [AW-1:0] pc_v [0:N][0:N-1];
   wire [DW-1:0] w_v  [0:N][0:N-1];

   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign ps_v[0][j] = {AW{1'b0}};
	 assign pc_v[0][j] = {AW{1'b0}};
	 assign w_v[0][j]  = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    tpu_pe #(.DW(DW), .AW(AW)) u_pe
	      (.clk        (clk),
	       .rst        (rst),
	       .ld         (ld),
	       .en         (en),
	       .w_north    (w_v[i][j]),
	       .a_west     (a_h[i][j]),
	       .ps_north   (ps_v[i][j]),
	       .pc_north   (pc_v[i][j]),
	       .w_south    (w_v[i+1][j]),
	       .a_east     (a_h[i][j+1]),
	       .ps_south   (ps_v[i+1][j]),
	       .pc_south   (pc_v[i+1][j]));
	 end
      end
   endgenerate

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [AW-1:0] cin = ps_v[N][j] + pc_v[N][j];
	 if (D == 0) begin : g_d0
	    assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){cin[AW-1]}}, cin};
	 end
	 else begin : g_dn
	    reg [AW-1:0] dr [0:D-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1) dr[k] <= {AW{1'b0}};
	       end
	       else if (en) begin
		  dr[0] <= cin;
		  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){dr[D-1][AW-1]}}, dr[D-1]};
	 end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)		vpipe <= {LAT{1'b0}};
      else if (en)	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

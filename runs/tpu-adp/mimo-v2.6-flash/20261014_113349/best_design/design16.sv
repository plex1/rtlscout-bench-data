//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) - variant 16.
// Behavior identical to design9. Tree-C, constant row paired with the
// negated bit-7 row in the first level-0 compressor:
//   t0 = (nx, one, r0), t1 = (r1,r2,r3), t2 = (r4,r5,r6).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [DW-1:0]	 w_north,
    input [DW-1:0]	 a_west,
    input [2*PW-1:0]	 psum_north,
    output [DW-1:0]	 w_south,
    output reg [DW-1:0]	 a_east,
    output reg [2*PW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   function [2*PW-1:0] csa3;
      input [PW-1:0] x, y, z;
      reg [PW-1:0]   sum, cry;
      begin
	 sum = x ^ y ^ z;
	 cry = ((x & y) | (y & z) | (x & z)) << 1;
	 csa3 = {cry, sum};
      end
   endfunction

   wire [PW-1:0] as = {{(PW-DW){a_west[DW-1]}}, a_west};
   wire [PW-1:0] r0 = w[0] ? as              : {PW{1'b0}};
   wire [PW-1:0] r1 = w[1] ? (as << 1)       : {PW{1'b0}};
   wire [PW-1:0] r2 = w[2] ? (as << 2)       : {PW{1'b0}};
   wire [PW-1:0] r3 = w[3] ? (as << 3)       : {PW{1'b0}};
   wire [PW-1:0] r4 = w[4] ? (as << 4)       : {PW{1'b0}};
   wire [PW-1:0] r5 = w[5] ? (as << 5)       : {PW{1'b0}};
   wire [PW-1:0] r6 = w[6] ? (as << 6)       : {PW{1'b0}};
   wire [PW-1:0] ex = w[7] ? (as << 7)       : {PW{1'b0}};
   wire [PW-1:0] nx = ~ex;
   wire [PW-1:0] one = {{(PW-1){1'b0}}, 1'b1};

   wire [PW-1:0] ns = psum_north[PW-1:0];
   wire [PW-1:0] nc = psum_north[2*PW-1:PW];

   // level 0
   wire [2*PW-1:0] t0 = csa3(nx, one, r0);
   wire [2*PW-1:0] t1 = csa3(r1, r2, r3);
   wire [2*PW-1:0] t2 = csa3(r4, r5, r6);
   // level 1
   wire [2*PW-1:0] u0 = csa3(t0[PW-1:0], t0[2*PW-1:PW], t1[PW-1:0]);
   wire [2*PW-1:0] u1 = csa3(t1[2*PW-1:PW], t2[PW-1:0], t2[2*PW-1:PW]);
   // level 2: psum pair split
   wire [2*PW-1:0] v0 = csa3(u0[PW-1:0], u0[2*PW-1:PW], ns);
   wire [2*PW-1:0] v1 = csa3(u1[PW-1:0], u1[2*PW-1:PW], nc);
   // level 3
   wire [2*PW-1:0] f0 = csa3(v0[PW-1:0], v0[2*PW-1:PW], v1[PW-1:0]);
   // level 4
   wire [2*PW-1:0] f1 = csa3(f0[PW-1:0], f0[2*PW-1:PW], v1[2*PW-1:PW]);

   wire [PW-1:0] ps = f1[PW-1:0];
   wire [PW-1:0] pc = f1[2*PW-1:PW];

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 psum_south <= {(2*PW){1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    psum_south <= {pc, ps};
	 end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [N*DW-1:0]	 w_top,
    input [N*DW-1:0]	 a_left,
    output [N*2*PW-1:0]	 c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [2*PW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign psum_v[0][j] = {(2*PW){1'b0}};
	 assign w_v[0][j]    = w_top[j*DW +: DW];
	 assign c_bot[j*2*PW +: 2*PW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	    tpu_pe #(.DW(DW), .PW(PW)) u_pe
	      (.clk        (clk),
	       .rst        (rst),
	       .ld         (ld),
	       .en         (en),
	       .w_north    (w_v[i][j]),
	       .a_west     (a_h[i][j]),
	       .psum_north (psum_v[i][j]),
	       .w_south    (w_v[i+1][j]),
	       .a_east     (a_h[i][j+1]),
	       .psum_south (psum_v[i+1][j]));
	 end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter PW = 19)
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

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din =
			a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < i; k = k + 1)
		    sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*DW +: DW] = sr[i-1];
	 end
      end
   endgenerate

   wire [N*2*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [2*PW-1:0] cin = c_bot[j*2*PW +: 2*PW];
	 wire [PW-1:0]   bin = cin[PW-1:0] + cin[2*PW-1:PW];
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-PW){bin[PW-1]}}, bin};
	 end
	 else begin : g_deskewn
	    reg [PW-1:0] dr [0:D-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1)
		    dr[k] <= {PW{1'b0}};
	       end
	       else if (en) begin
		  dr[0] <= bin;
		  for (k = 1; k < D; k = k + 1)
		    dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
	 end
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

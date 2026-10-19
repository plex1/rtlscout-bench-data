//#############################################################################
// v10: reference systolic dataflow + radix-4 Booth multiplier whose partial
//      products are compressed straight into the carry-save accumulator.
//
// Per PE the exact reference value is  psum_s(t+1) = psum_n(t) + a(t)*w(t).
// The accumulator is kept as a redundant (sum,carry) pair, so no carry-
// propagate adder is on the critical path; the 8x8 signed product is formed
// as 4 Booth digits (values 0,+-a,+-2a) and folded into the same 3:2
// compressor tree.  Critical path = Booth digit select + 4 compressor levels
// instead of a ~19-level array multiplier.
//#############################################################################
module tpu_csa
  #(parameter W = 19)
   (input [W-1:0] x,
    input [W-1:0] y,
    input [W-1:0] z,
    output [W-1:0] s,
    output [W-1:0] c);

   assign s = x ^ y ^ z;
   assign c = ((x & y) | (x & z) | (y & z)) << 1;

endmodule

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

   wire signed [DW-1:0] aa = a_west;
   wire signed [DW-1:0] bb = w;

   // sign-extended multiplicand and its double
   wire signed [AW-1:0] ae  = {{(AW-DW){aa[DW-1]}}, aa};
   wire signed [AW-1:0] a2e = ae <<< 1;

   // radix-4 Booth digits of bb:  digit k = (bb[2k+1],bb[2k],bb[2k-1])
   wire [3:0] bm1, bm2, bneg;
   wire [AW-1:0] pp0, pp1, pp2, pp3;

   genvar k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_booth
	 wire lo  = (k == 0) ? 1'b0 : bb[2*k-1];
	 wire mid = bb[2*k];
	 wire hi  = bb[2*k+1];
	 wire signed [AW-1:0] xk;
	 assign bm1[k]  = lo ^ mid;
	 assign bm2[k]  = (hi & ~mid & ~lo) | (~hi & mid & lo);
	 assign bneg[k] = hi;
	 assign xk = bm2[k] ? a2e : (bm1[k] ? ae : {AW{1'b0}});
	 // negative partial products are kept in one's-complement form; the
	 // +1 corrections are collected in 'corr' (bits 0,2,4,6 only)
	 if (k == 0) assign pp0 = (xk ^ {AW{bneg[k]}}) <<< 0;
	 if (k == 1) assign pp1 = (xk ^ {AW{bneg[k]}}) <<< 2;
	 if (k == 2) assign pp2 = (xk ^ {AW{bneg[k]}}) <<< 4;
	 if (k == 3) assign pp3 = (xk ^ {AW{bneg[k]}}) <<< 6;
      end
   endgenerate

   wire [AW-1:0] corr = {{(AW-7){1'b0}}, bneg[3], 1'b0, bneg[2], 1'b0,
			 bneg[1], 1'b0, bneg[0]};

   // 3:2 compressor tree:  (ps,pc) + pp0..pp3 + corr
   wire [AW-1:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5;

   tpu_csa #(.W(AW)) u_c1 (.x(ps_north), .y(pc_north), .z(pp0),
			   .s(s1), .c(c1));
   tpu_csa #(.W(AW)) u_c2 (.x(pp1), .y(pp2), .z(pp3), .s(s2), .c(c2));
   tpu_csa #(.W(AW)) u_c3 (.x(s1), .y(c1), .z(s2), .s(s3), .c(c3));
   tpu_csa #(.W(AW)) u_c4 (.x(s3), .y(c3), .z(c2), .s(s4), .c(c4));
   tpu_csa #(.W(AW)) u_c5 (.x(s4), .y(c4), .z(corr), .s(s5), .c(c5));

   always @(posedge clk) begin
      if (rst) begin
	 w        <= {DW{1'b0}};
	 a_east   <= {DW{1'b0}};
	 ps_south <= {AW{1'b0}};
	 pc_south <= {AW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east   <= a_west;
	    ps_south <= s5;
	    pc_south <= c5;
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

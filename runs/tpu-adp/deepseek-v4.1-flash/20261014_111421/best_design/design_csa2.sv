//#############################################################################
// TPU MXU variant "csa2": two-stage multiply + carry-save accumulation.
//
// Per PE (rows 0..N-2):
//   stage 1: register Pa = a_west[7:4] * w  and  Pb = a_west[3:0] * w
//   stage 2: carry-save compress {s_north, c_north, Pa<<4, Pb}
//            -> (s_south, c_south)   (two 3:2 compressors = 2 FA levels)
// Bottom row (N-1): stage 1 as above, then a combinational carry-save
//   compress + carry-propagate add yields c_bot directly.
//
// The carry-save invariant is  psum = s + c  (mod 2^PW), and the true partial
// sums always fit in signed PW bits, so the final carry-propagate add is exact.
// Timing: activation at PE(i,j) at t_m+i+j; (s,c) of row i+1 registered one
// cycle later per row; bottom resolve at t_m+j+8; deskew (N-1-j) -> t_m+15.
//#############################################################################

module tpu_pe_cs
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    s_north,
    input [PW-1:0]	    c_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output reg [PW-1:0]	    s_south,
    output reg [PW-1:0]	    c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // ---- stage 1: split the activation into two halves --------------------
   wire signed [DW/2-1:0]  ah = a_west[DW-1:DW/2];
   wire signed [DW/2:0]    al = {1'b0, a_west[DW/2-1:0]};
   wire signed [DW/2+DW-1:0]     pa = ah * $signed(w);	  // 12b signed
   wire signed [DW/2+1+DW-1:0]   pb = al * $signed(w);	  // 13b signed

   reg signed [DW/2+1+DW-1:0] pa_r;
   reg signed [DW/2+1+DW-1:0] pb_r;

   localparam CW = DW/2 + 1 + DW;	// 13

   // ---- stage 2: carry-save compress -------------------------------------
   wire [PW-1:0] A = {{(PW-CW){pa_r[CW-1]}}, pa_r} << 4;
   wire [PW-1:0] B = {{(PW-CW){pb_r[CW-1]}}, pb_r};

   wire [PW-1:0] e1 = s_north ^ c_north;
   wire [PW-1:0] s1 = e1 ^ A;
   wire [PW-1:0] c1 = ((s_north & c_north) | (e1 & A)) << 1;
   wire [PW-1:0] e2 = s1 ^ c1;
   wire [PW-1:0] s2 = e2 ^ B;
   wire [PW-1:0] c2 = ((s1 & c1) | (e2 & B)) << 1;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
	 pa_r   <= {CW{1'b0}};
	 pb_r   <= {CW{1'b0}};
	 s_south<= {PW{1'b0}};
	 c_south<= {PW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east  <= a_west;
	    pa_r    <= pa;
	    pb_r    <= pb;
	    s_south <= s2;
	    c_south <= c2;
	 end
      end
   end

endmodule

module tpu_pe_csbot
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    s_north,
    input [PW-1:0]	    c_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output [PW-1:0]	    c_bot
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [DW/2-1:0]  ah = a_west[DW-1:DW/2];
   wire signed [DW/2:0]    al = {1'b0, a_west[DW/2-1:0]};
   wire signed [DW/2+DW-1:0]     pa = ah * $signed(w);
   wire signed [DW/2+1+DW-1:0]   pb = al * $signed(w);

   reg signed [DW/2+1+DW-1:0] pa_r;
   reg signed [DW/2+1+DW-1:0] pb_r;

   localparam CW = DW/2 + 1 + DW;	// 13

   wire [PW-1:0] A = {{(PW-CW){pa_r[CW-1]}}, pa_r} << 4;
   wire [PW-1:0] B = {{(PW-CW){pb_r[CW-1]}}, pb_r};

   wire [PW-1:0] e1 = s_north ^ c_north;
   wire [PW-1:0] s1 = e1 ^ A;
   wire [PW-1:0] c1 = ((s_north & c_north) | (e1 & A)) << 1;
   wire [PW-1:0] e2 = s1 ^ c1;
   wire [PW-1:0] s2 = e2 ^ B;
   wire [PW-1:0] c2 = ((s1 & c1) | (e2 & B)) << 1;

   assign c_bot = s2 + c2;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
	 pa_r   <= {CW{1'b0}};
	 pb_r   <= {CW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east <= a_west;
	    pa_r   <= pa;
	    pb_r   <= pb;
	 end
      end
   end

endmodule

module tpu_array_cs
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PW-1:0]	c_bot
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [PW-1:0] s_v [0:N][0:N-1];
   wire [PW-1:0] c_v [0:N][0:N-1];
   wire [DW-1:0] w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PW{1'b0}};
         assign c_v[0][j] = {PW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            if (i == N-1) begin : g_bot
               tpu_pe_csbot #(.DW(DW), .PW(PW)) u_pe
                 (.clk(clk), .rst(rst), .ld(ld), .en(en),
                  .w_north(w_v[i][j]), .a_west(a_h[i][j]),
                  .s_north(s_v[i][j]), .c_north(c_v[i][j]),
                  .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
                  .c_bot(c_bot[j*PW +: PW]));
            end
            else begin : g_mid
               tpu_pe_cs #(.DW(DW), .PW(PW)) u_pe
                 (.clk(clk), .rst(rst), .ld(ld), .en(en),
                  .w_north(w_v[i][j]), .a_west(a_h[i][j]),
                  .s_north(s_v[i][j]), .c_north(c_v[i][j]),
                  .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
                  .s_south(s_v[i+1][j]), .c_south(c_v[i+1][j]));
            end
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
    input		    clk,
    input		    rst,
    input		    w_valid,
    input [N*DW-1:0]	    w_data,
    input		    a_valid,
    input [N*DW-1:0]	    a_data,
    output		    c_valid,
    output [N*ACCW-1:0]	    c_data
    );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
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

   wire [N*PW-1:0] c_bot;
   tpu_array_cs #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PW-1:0] cin = c_bot[j*PW +: PW];
         wire [PW-1:0] dout;
         if (D == 0) begin : g_deskew0
            assign dout = cin;
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
               end
            end
            assign dout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){dout[PW-1]}}, dout};
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

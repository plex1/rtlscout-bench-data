//#############################################################################
// TPU MXU variant "split2": exact split-multiply systolic array.
//
// The 8-bit signed activation is written as  a = ah*16 + al  where
// ah = a[7:4] (signed 4-bit) and al = a[3:0] (unsigned 4-bit = {0,a[3:0]}).
// Then  a*w = (ah*w)*16 + al*w, and the column accumulation splits into two
// independent exact accumulators:
//      C = 16 * (sum_i ah_i*w_i) + (sum_i al_i*w_i)
// Both partial chains use a 4x8 multiplier (about half the depth of an 8x8),
// which shortens the PE critical path.  The combine (shift + add) happens once
// per column at the bottom of the array, and the deskew is unchanged (the
// combined value fits in 19 bits exactly).
//
// accumulators:  |acc_hi| <= 8*8*128 = 8192, |acc_lo| <= 8*15*128 = 15360
//                -> 15 signed bits each.
//#############################################################################

module tpu_pe2
  #(parameter DW = 8,
    parameter PW = 15)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    hi_north,
    input [PW-1:0]	    lo_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output reg [PW-1:0]	    hi_south,
    output reg [PW-1:0]	    lo_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [3:0]	ah = a_west[7:4];
   wire signed [4:0]	al = {1'b0, a_west[3:0]};

   wire signed [11:0] ph = ah * $signed(w);	 // 12b
   wire signed [12:0] pl = al * $signed(w);	 // 13b

   wire signed [PW-1:0] the = {{(PW-12){ph[11]}}, ph};
   wire signed [PW-1:0] tle = {{(PW-13){pl[12]}}, pl};

   always @(posedge clk) begin
      if (rst) begin
	 w        <= {DW{1'b0}};
	 a_east   <= {DW{1'b0}};
	 hi_south <= {PW{1'b0}};
	 lo_south <= {PW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east   <= a_west;
	    hi_south <= hi_north + the;
	    lo_south <= lo_north + tle;
	 end
      end
   end

endmodule

module tpu_pe2_bot
  #(parameter DW = 8,
    parameter PW = 15)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    hi_north,
    input [PW-1:0]	    lo_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output [PW-1:0]	    hi_south,
    output [PW-1:0]	    lo_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [3:0]	ah = a_west[7:4];
   wire signed [4:0]	al = {1'b0, a_west[3:0]};
   wire signed [11:0] ph = ah * $signed(w);
   wire signed [12:0] pl = al * $signed(w);
   wire signed [PW-1:0] the = {{(PW-12){ph[11]}}, ph};
   wire signed [PW-1:0] tle = {{(PW-13){pl[12]}}, pl};

   assign hi_south = hi_north + the;
   assign lo_south = lo_north + tle;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) a_east <= a_west;
      end
   end

endmodule

module tpu_array2
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 15)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PW-1:0]	c_bot
    );

   wire [DW-1:0] a_h  [0:N-1][0:N];
   wire [PW-1:0] hi_v [0:N][0:N-1];
   wire [PW-1:0] lo_v [0:N][0:N-1];
   wire [DW-1:0] w_v  [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1)
	assign a_h[i][0] = a_left[i*DW +: DW];
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign hi_v[0][j] = {PW{1'b0}};
	 assign lo_v[0][j] = {PW{1'b0}};
	 assign w_v[0][j]  = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
	 for (j = 0; j < N; j = j + 1) begin : g_pe_col
	       tpu_pe2 #(.DW(DW), .PW(PW)) u_pe
		 (.clk(clk), .rst(rst), .ld(ld), .en(en),
		  .w_north(w_v[i][j]), .a_west(a_h[i][j]),
		  .hi_north(hi_v[i][j]), .lo_north(lo_v[i][j]),
		  .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
		  .hi_south(hi_v[i+1][j]), .lo_south(lo_v[i+1][j]));
	 end
      end
      for (j = 0; j < N; j = j + 1) begin : g_cout
	 wire signed [PW-1:0] h = hi_v[N][j];
	 wire signed [PW-1:0] l = lo_v[N][j];
	 wire signed [PW+3:0] comb = (h <<< 4) + l;
	 assign c_bot[j*PW +: PW] = comb[PW-1:0];
      end   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter PW = 15)
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
   tpu_array2 #(.N(N), .DW(DW), .PW(PW)) u_array
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

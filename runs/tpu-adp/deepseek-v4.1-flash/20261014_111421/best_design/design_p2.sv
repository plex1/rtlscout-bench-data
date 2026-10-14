//#############################################################################
// TPU MXU variant "p2": multiply and accumulate split into separate register
// stages, with the latency restored by making the bottom row's final add
// combinational.
//
//   PE(i,j) for i < 7:
//       preg       <= a_west * w
//       psum_south <= psum_north + preg      (registered add)
//   PE(7,j) (bottom row):
//       preg       <= a_west * w
//       c_bot      =  psum_north + preg      (combinational add)
//
// Contributions: activation for row i reaches PE(i,j) at cycle p+i+j; it is
// registered into preg (1 cycle), enters the psum register of row i (1 cycle)
// and falls 7-i rows -> psum_v[7] at p+j+8.  The bottom-row combinational add
// adds row 7's contribution at the same cycle, so c_bot is valid at p+j+8.
// The deskew (7-j) then lines every column up at p+15 = the required latency.
// Critical path per PE is now either the multiplier or the adder, not both.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    psum_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output reg [PW-1:0]	    psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // split multiply: a = ah*16 + al  =>  a*w = (ah*w)*16 + al*w
   wire signed [2*DW-1:0] pcomb = $signed(a_west) * $signed(w);
   reg  signed [2*DW-1:0] preg;
   wire [PW-1:0] pb = {{(PW-2*DW){preg[2*DW-1]}}, preg};

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 preg       <= {2*DW{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    preg       <= pcomb;
	    psum_south <= psum_north + pb;
	 end
      end
   end

endmodule

module tpu_pe_bot
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		    clk,
    input		    rst,
    input		    ld,
    input		    en,
    input [DW-1:0]	    w_north,
    input [DW-1:0]	    a_west,
    input [PW-1:0]	    psum_north,
    output [DW-1:0]	    w_south,
    output reg [DW-1:0]	    a_east,
    output [PW-1:0]	    psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [2*DW-1:0] pcomb = $signed(a_west) * $signed(w);
   reg  signed [2*DW-1:0] preg;
   wire [PW-1:0] pb = {{(PW-2*DW){preg[2*DW-1]}}, preg};

   assign psum_south = psum_north + pb;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
	 preg   <= {2*DW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east <= a_west;
	    preg   <= pcomb;
	 end
      end
   end

endmodule

module tpu_array
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

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [PW-1:0]   psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            if (i == N-1) begin : g_bot
               tpu_pe_bot #(.DW(DW), .PW(PW)) u_pe
                 (.clk(clk), .rst(rst), .ld(ld), .en(en),
                  .w_north(w_v[i][j]), .a_west(a_h[i][j]),
                  .psum_north(psum_v[i][j]),
                  .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
                  .psum_south(psum_v[i+1][j]));
            end
            else begin : g_mid
               tpu_pe #(.DW(DW), .PW(PW)) u_pe
                 (.clk(clk), .rst(rst), .ld(ld), .en(en),
                  .w_north(w_v[i][j]), .a_west(a_h[i][j]),
                  .psum_north(psum_v[i][j]),
                  .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
                  .psum_south(psum_v[i+1][j]));
            end
         end
      end
      for (j = 0; j < N; j = j + 1) begin : g_cout
         assign c_bot[j*PW +: PW] = psum_v[N][j];
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
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
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

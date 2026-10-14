//#############################################################################
// TPU MXU variant "p3": like p2 (mult and accumulate in separate register
// stages, bottom row combinational), but the 8x8 product is split into two
// smaller multipliers so neither the 8x8 multiply nor the final add alone
// dominates the PE critical path:
//     a = ah*16 + al        =>        a*w = (ah*w)*16 + al*w
//   stage 1: pH <= ah*w (4x8),  pL <= al*w (5x8)
//   stage 2: psum_south <= psum_north + (pH<<4) + pL
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

   localparam CW = DW/2 + DW;	// 12: width of ah*w (4x8 signed -> 12b)
   localparam CL = DW/2 + 1 + DW;	// 13: width of al*w (5x8 signed -> 13b)

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [DW/2-1:0] ah = a_west[DW-1:DW/2];
   wire signed [DW/2:0]	  al = {1'b0, a_west[DW/2-1:0]};
   wire signed [CW-1:0]	  pH = ah * $signed(w);
   wire signed [CL-1:0]	  pL = al * $signed(w);

   reg signed [CW-1:0] pH_r;
   reg signed [CL-1:0] pL_r;

   wire [PW-1:0] hi = {{(PW-CW){pH_r[CW-1]}}, pH_r};
   wire [PW-1:0] lo = {{(PW-CL){pL_r[CL-1]}}, pL_r};

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 pH_r       <= {CW{1'b0}};
	 pL_r       <= {CL{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    pH_r       <= pH;
	    pL_r       <= pL;
	    psum_south <= psum_north + (hi <<< (DW/2)) + lo;
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

   localparam CW = DW/2 + DW;
   localparam CL = DW/2 + 1 + DW;

   reg [DW-1:0] w;
   assign w_south = w;

   wire signed [DW/2-1:0] ah = a_west[DW-1:DW/2];
   wire signed [DW/2:0]	  al = {1'b0, a_west[DW/2-1:0]};
   wire signed [CW-1:0]	  pH = ah * $signed(w);
   wire signed [CL-1:0]	  pL = al * $signed(w);

   reg signed [CW-1:0] pH_r;
   reg signed [CL-1:0] pL_r;

   wire [PW-1:0] hi = {{(PW-CW){pH_r[CW-1]}}, pH_r};
   wire [PW-1:0] lo = {{(PW-CL){pL_r[CL-1]}}, pL_r};

   assign psum_south = psum_north + (hi <<< (DW/2)) + lo;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
	 pH_r   <= {CW{1'b0}};
	 pL_r   <= {CL{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east <= a_west;
	    pH_r <= pH;
	    pL_r <= pL;
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
      for (i = 0; i < N; i = i + 1)
	assign a_h[i][0] = a_left[i*DW +: DW];
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign psum_v[0][j] = {PW{1'b0}};
	 assign w_v[0][j]    = w_top[j*DW +: DW];
	 assign c_bot[j*PW +: PW] = psum_v[N][j];
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

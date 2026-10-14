//#############################################################################
// TPU MXU variant "booth": pipelined PE (mult and accumulate in separate
// register stages, bottom row combinational to restore 15-cycle latency) with
// an explicit radix-4 Booth multiplier (4 partial products instead of 8).
//#############################################################################

module booth_mul #(parameter DW = 8, parameter PW = 19)
   (input  signed [DW-1:0] x,	// multiplicand (activation, signed)
    input	  [DW-1:0] y,	// multiplier (weight)
    output signed [PW-1:0] p);

   wire signed [PW-1:0] xe = {{(PW-DW){x[DW-1]}}, x};
   wire signed [PW-1:0] x2 = xe <<< 1;

   wire signed [PW-1:0] pp [0:DW/2-1];
   genvar i;
   generate
      for (i = 0; i < DW/2; i = i + 1) begin : g_digit
	 wire a = y[2*i+1];
	 wire b = y[2*i];
	 wire c = (2*i == 0) ? 1'b0 : y[2*i-1];
	 wire mg2 = a ? (~b & ~c) : (b & c);
	 wire mg1 = b ^ c;
	 wire signed [PW-1:0] sel0 = mg2 ? x2 : (mg1 ? xe : {PW{1'b0}});
	 wire signed [PW-1:0] sel  = a ? (-sel0) : sel0;
	 assign pp[i] = sel <<< (2*i);
      end
   endgenerate

   wire signed [PW-1:0] s01 = pp[0] + pp[1];
   wire signed [PW-1:0] s23 = pp[2] + pp[3];
   assign p = s01 + s23;

endmodule

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

   wire signed [PW-1:0] prod;
   booth_mul #(.DW(DW), .PW(PW)) u_mul (.x(a_west), .y(w), .p(prod));

   reg signed [PW-1:0] preg;

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 preg       <= {PW{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    preg       <= prod;
	    psum_south <= psum_north + preg;
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

   wire signed [PW-1:0] prod;
   booth_mul #(.DW(DW), .PW(PW)) u_mul (.x(a_west), .y(w), .p(prod));

   reg signed [PW-1:0] preg;

   assign psum_south = psum_north + preg;

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 a_east <= {DW{1'b0}};
	 preg   <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east <= a_west;
	    preg   <= prod;
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

// Systolic TPU with per-stage narrowed accumulators.
// The psum register after row i holds a sum of (i+1) products:
//   [-2^15, 2^15-1] * (i+1) -> widths 16,17,17,18,18,18,18,19.
// Incoming psum at row i is sign-extended (zero-padded for row 0) to the
// local width and $signed operands keep the arithmetic correct.

module tpu_pe
  #(parameter DW = 8,
    parameter PSW = 19)      // width of this PE's psum_south register
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PSW-1:0]	  psum_north, // sign-extended to match PSW
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PSW-1:0]  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
	 w	  <= {DW{1'b0}};
	 a_east	  <= {DW{1'b0}};
	 psum_south <= {PSW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) psum_south <= $signed(psum_north) + $signed(prod);
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*19-1:0]	c_bot
    );

   function integer psw;   // width after accumulating k products
      input integer k;
      begin
	 if (k <= 1)    psw = 16;
	 else if (k <= 3) psw = 17;
	 else if (k <= 7) psw = 18;
	 else             psw = 19;
      end
   endfunction

   // psum nets, max width 19; PE(i,j) uses [psw(i+1)-1:0]
   wire [18:0] psum_v [0:N][0:N-1];

   wire [DW-1:0]   a_h [0:N-1][0:N];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign c_bot[j*19 +: 19] = psum_v[N-1][j][18:0];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            if (i == 0) begin : g_top
	       // top row: incoming psum is zero (fresh accumulation)
               wire [psw(1)-1:0] zero = {psw(1){1'b0}};
               tpu_pe #(.DW(DW), .PSW(psw(1))) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[0][j]),
               .a_west     (a_h[0][j]),
               .psum_north (zero),
               .w_south    (w_v[1][j]),
               .a_east     (a_h[0][j+1]),
               .psum_south (psum_v[0][j][psw(1)-1:0]));
	    end
	 end
      end
   endgenerate

endmodule

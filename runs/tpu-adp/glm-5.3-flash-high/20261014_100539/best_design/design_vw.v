// Systolic TPU with per-stage narrowed accumulators.
// The psum register leaving row i accumulates (i+1) products, each in
// [-2^15, 2^15-1], so widths 16,17,17,18,18,18,18,19 suffice for rows 0..7.
// PE(i,j) computes with only psw(i+1) bits of accumulator logic; the wider
// interconnect nets are sign-filled so the bottom row presents a full 19-bit
// result. The tpu top sign-extends 19 -> 32 at the output pins.

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
    input [PSW-1:0]	  psum_north, // zero/sign-extended by the array wiring
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
	 if (en) begin
	    a_east     <= a_west;
	    psum_south <= $signed(psum_north) + $signed(prod);
	 end
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

   function integer psw;   // accumulator width after (k) products
      input integer k;
      begin
	 if (k <= 1)      psw = 16;
	 else if (k <= 3) psw = 17;
	 else if (k <= 7) psw = 18;
	 else             psw = 19;
      end
   endfunction

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [18:0]   psum_v [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = 19'd0;	     // fresh accumulation at the top
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*19 +: 19] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            // incoming psum, cut to this row's width (row 0: zero)
            wire [psw(i+1)-1:0] pin;
            if (i == 0) begin : g_zin
               assign pin = {(psw(1)){1'b0}};
            end
            else begin : g_sin
               assign pin = psum_v[i][j][psw(i+1)-1:0];
            end
            tpu_pe #(.DW(DW), .PSW(psw(i+1))) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (pin),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j][psw(i+1)-1:0]));
         end
      end
      // sign-fill the unused high bits of every psum net
      for (i = 1; i <= N; i = i + 1) begin : g_fill
         for (j = 0; j < N; j = j + 1) begin : g_fillj
            if (psw(i) < 19) begin : g_fillw
               assign psum_v[i][j][18:psw(i)] =
		      {(19-psw(i)){psum_v[i][j][psw(i)-1]}};
            end
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter IACCW = 19)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
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

   wire [N*IACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IACCW-1:0] cin = c_bot[j*IACCW +: IACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-IACCW){cin[IACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [IACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {IACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-IACCW){dr[D-1][IACCW-1]}}, dr[D-1]};
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

//#############################################################################
// TPU MXU with variable-width accumulators per row
//#############################################################################

module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north + prod_ext;
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
    output [23*N-1:0]	c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [DW-1:0]   w_v    [0:N][0:N-1];
   wire [15:0]  psum_0 [0:N-1];
   wire [16:0]  psum_1 [0:N-1];
   wire [17:0]  psum_2 [0:N-1];
   wire [18:0]  psum_3 [0:N-1];
   wire [19:0]  psum_4 [0:N-1];
   wire [20:0]  psum_5 [0:N-1];
   wire [21:0]  psum_6 [0:N-1];
   wire [22:0]  psum_7 [0:N-1];

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_west
         assign a_h[0][0] = a_left[0*DW +: DW];
         assign a_h[1][0] = a_left[1*DW +: DW];
         assign a_h[2][0] = a_left[2*DW +: DW];
         assign a_h[3][0] = a_left[3*DW +: DW];
         assign a_h[4][0] = a_left[4*DW +: DW];
         assign a_h[5][0] = a_left[5*DW +: DW];
         assign a_h[6][0] = a_left[6*DW +: DW];
         assign a_h[7][0] = a_left[7*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign c_bot[j*23 +: 23] = psum_7[j];
      end
      // Row 0 PEs (ACCW=16)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r0
         tpu_pe #(.DW(DW), .ACCW(16)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[0][j]), .a_west(a_h[0][j]),
            .psum_north(16'b0),
            .w_south(w_v[1][j]), .a_east(a_h[0][j+1]),
            .psum_south(psum_0[j]));
      end
      // Row 1 PEs (ACCW=17)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r1
         tpu_pe #(.DW(DW), .ACCW(17)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[1][j]), .a_west(a_h[1][j]),
            .psum_north({psum_0[j][15], psum_0[j]}),
            .w_south(w_v[2][j]), .a_east(a_h[1][j+1]),
            .psum_south(psum_1[j]));
      end
      // Row 2 PEs (ACCW=18)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r2
         tpu_pe #(.DW(DW), .ACCW(18)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[2][j]), .a_west(a_h[2][j]),
            .psum_north({psum_1[j][16], psum_1[j]}),
            .w_south(w_v[3][j]), .a_east(a_h[2][j+1]),
            .psum_south(psum_2[j]));
      end
      // Row 3 PEs (ACCW=19)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r3
         tpu_pe #(.DW(DW), .ACCW(19)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[3][j]), .a_west(a_h[3][j]),
            .psum_north({psum_2[j][17], psum_2[j]}),
            .w_south(w_v[4][j]), .a_east(a_h[3][j+1]),
            .psum_south(psum_3[j]));
      end
      // Row 4 PEs (ACCW=20)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r4
         tpu_pe #(.DW(DW), .ACCW(20)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[4][j]), .a_west(a_h[4][j]),
            .psum_north({psum_3[j][18], psum_3[j]}),
            .w_south(w_v[5][j]), .a_east(a_h[4][j+1]),
            .psum_south(psum_4[j]));
      end
      // Row 5 PEs (ACCW=21)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r5
         tpu_pe #(.DW(DW), .ACCW(21)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[5][j]), .a_west(a_h[5][j]),
            .psum_north({psum_4[j][19], psum_4[j]}),
            .w_south(w_v[6][j]), .a_east(a_h[5][j+1]),
            .psum_south(psum_5[j]));
      end
      // Row 6 PEs (ACCW=22)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r6
         tpu_pe #(.DW(DW), .ACCW(22)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[6][j]), .a_west(a_h[6][j]),
            .psum_north({psum_5[j][20], psum_5[j]}),
            .w_south(w_v[7][j]), .a_east(a_h[6][j+1]),
            .psum_south(psum_6[j]));
      end
      // Row 7 PEs (ACCW=23)
      for (j = 0; j < N; j = j + 1) begin : g_pe_r7
         tpu_pe #(.DW(DW), .ACCW(23)) u_pe
           (.clk(clk), .rst(rst), .ld(ld), .en(en),
            .w_north(w_v[7][j]), .a_west(a_h[7][j]),
            .psum_north({psum_6[j][21], psum_6[j]}),
            .w_south(w_v[8][j]), .a_east(a_h[7][j+1]),
            .psum_south(psum_7[j]));
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
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
   localparam PACCW = 23;

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

   wire [N*PACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = c_bot[j*PACCW +: PACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PACCW){cin[PACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PACCW-1:0] dout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PACCW){dout[PACCW-1]}}, dout};
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

// Carry-save systolic array; the final carry-propagate addition is moved
// after the per-column output de-skew. The array emits both s and c streams,
// de-skew stores them separately (redundant columns stay time-aligned), then
// a combinational final add reconstructs each 32-bit result.

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
    input [ACCW-1:0]	  s_north,
    input [ACCW-1:0]	  c_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] s_south,
    output reg [ACCW-1:0] c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] s_ab = s_north ^ c_north;
   wire [ACCW-1:0] sum_w = s_ab ^ prod_ext;
   wire [ACCW-1:0] carry_w = {(s_north & c_north) | (s_ab & prod_ext), 1'b0};

   always @(posedge clk) begin
      if (rst) begin
         w      <= {DW{1'b0}};
         a_east <= {DW{1'b0}};
         s_south <= {ACCW{1'b0}};
         c_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east <= a_west;
            s_south <= sum_w;
            c_south <= carry_w[ACCW-1:0];
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*ACCW-1:0]	s_bot,
    output [N*ACCW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h [0:N-1][0:N];
   wire [ACCW-1:0] s_v [0:N][0:N-1];
   wire [ACCW-1:0] c_v [0:N][0:N-1];
   wire [DW-1:0]   w_v [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {ACCW{1'b0}};
         assign c_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign s_bot[j*ACCW +: ACCW] = s_v[N][j];
         assign c_bot[j*ACCW +: ACCW] = c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .s_north    (s_v[i][j]),
               .c_north    (c_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .s_south    (s_v[i+1][j]),
               .c_south    (c_v[i+1][j]));
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

   wire [N*IACCW-1:0] s_bot, c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(IACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .s_bot(s_bot), .c_bot(c_bot));

   wire [N*IACCW-1:0] s_delayed, c_delayed;
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IACCW-1:0] sin = s_bot[j*IACCW +: IACCW];
         wire [IACCW-1:0] cin = c_bot[j*IACCW +: IACCW];
         if (D == 0) begin : g_deskew0
            assign s_delayed[j*IACCW +: IACCW] = sin;
            assign c_delayed[j*IACCW +: IACCW] = cin;
         end
         else begin : g_deskewn
            reg [IACCW-1:0] sr [0:D-1];
            reg [IACCW-1:0] cr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) begin
                     sr[k] <= {IACCW{1'b0}};
                     cr[k] <= {IACCW{1'b0}};
                  end
               end
               else if (en) begin
                  sr[0] <= sin;
                  cr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) begin
                     sr[k] <= sr[k-1];
                     cr[k] <= cr[k-1];
                  end
               end
            end
            assign s_delayed[j*IACCW +: IACCW] = sr[D-1];
            assign c_delayed[j*IACCW +: IACCW] = cr[D-1];
         end
      end
   endgenerate

   // Final carry-propagate addition after de-skew, then sign-extension.
   genvar j2;
   generate
      for (j2 = 0; j2 < N; j2 = j2 + 1) begin : g_final
         wire [IACCW-1:0] s_d = s_delayed[j2*IACCW +: IACCW];
         wire [IACCW-1:0] c_d = c_delayed[j2*IACCW +: IACCW];
         wire [IACCW-1:0] sum_d = s_d + c_d;
         assign c_data[j2*ACCW +: ACCW] =
              {{(ACCW-IACCW){sum_d[IACCW-1]}}, sum_d};
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
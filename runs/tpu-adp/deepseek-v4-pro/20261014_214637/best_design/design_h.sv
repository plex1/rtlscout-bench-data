// Carry-save systolic variant with per-row adaptive accumulator width.
// Row i sums at most i+1 products, so a minimal signed accumulator width
// can be used at each stage: row widths are 16,17,17,18,18,18,18,19 for N=8.
// PEs sign-extend incoming partials to the local output width.

module tpu_pe
  #(parameter DW = 8,
    parameter INW = 16,
    parameter OUTW = 17)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [INW-1:0]	  s_north,
    input [INW-1:0]	  c_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [OUTW-1:0] s_south,
    output reg [OUTW-1:0] c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [OUTW-1:0] s_ext = {{(OUTW-INW){s_north[INW-1]}}, s_north};
   wire [OUTW-1:0] c_ext = {{(OUTW-INW){c_north[INW-1]}}, c_north};

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [OUTW-1:0] prod_ext = {{(OUTW-2*DW){prod[2*DW-1]}}, prod};

   wire [OUTW-1:0] s_ab = s_ext ^ c_ext;
   wire [OUTW-1:0] sum_w = s_ab ^ prod_ext;
   wire [OUTW:0]   carry_w = {(s_ext & c_ext) | (s_ab & prod_ext), 1'b0};

   always @(posedge clk) begin
      if (rst) begin
         w      <= {DW{1'b0}};
         a_east <= {DW{1'b0}};
         s_south <= {OUTW{1'b0}};
         c_south <= {OUTW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east <= a_west;
            s_south <= sum_w;
            c_south <= carry_w[OUTW-1:0];
         end
      end
   end

endmodule

// Minimum signed accumulator width that represents a sum of k signed
// 8x8 products. Max positive = 16384*k; bound 2^(W-1)-1 >= 16384*k.
function integer min_accw(input integer k);
   integer w;
   begin
      w = 16;
      while (w < 32) begin
         if (((1 << (w - 1)) - 1) >= 16384 * k)
           return w;
         w = w + 1;
      end
      min_accw = w;
   end
endfunction

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
         assign c_bot[j*ACCW +: ACCW] = s_v[N][j] + c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         localparam integer K    = i + 1;
         localparam integer WOUT = min_accw(K);
         localparam integer WIN  = (i == 0) ? WOUT : min_accw(i);
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .INW(WIN), .OUTW(WOUT)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .s_north    (s_v[i][j][WIN-1:0]),
               .c_north    (c_v[i][j][WIN-1:0]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .s_south    (s_v[i+1][j][WOUT-1:0]),
               .c_south    (c_v[i+1][j][WOUT-1:0]));
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
   tpu_array #(.N(N), .DW(DW), .ACCW(IACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   wire [N*ACCW-1:0] c_before_ext;
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IACCW-1:0] cin = c_bot[j*IACCW +: IACCW];
         if (D == 0) begin : g_deskew0
            assign c_before_ext[j*ACCW +: ACCW] =
                 {{(ACCW-IACCW){cin[IACCW-1]}}, cin};
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
            assign c_before_ext[j*ACCW +: ACCW] =
                 {{(ACCW-IACCW){dr[D-1][IACCW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   assign c_data = c_before_ext;

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule
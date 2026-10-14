// Optimized version 11: split multiply into two 4x8 partial products
// a * w = (a_lo + a_hi * 16) * w = a_lo * w + (a_hi * w) << 4
// where a_lo = a[3:0], a_hi = a[7:4] (both unsigned for the multiply,
// then correct for sign)
// Actually, let's keep $signed but let Yosys decide. Instead try:
// Pre-register the weight sign extension and just do unsigned multiply
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

   // Split a into low 4 and high 4 bits
   wire [3:0] a_lo = a_west[3:0];
   wire [3:0] a_hi = a_west[7:4];
   
   // Partial products: a_lo * w (unsigned 4x8 = 12 bits) and a_hi * w (4x8 = 12 bits)
   wire [11:0] pp_lo = a_lo * w;       // a_lo * w (unsigned)
   wire [11:0] pp_hi = a_hi * w;       // a_hi * w (unsigned)
   
   // Combine: full unsigned product = pp_lo + (pp_hi << 4) (16 bits)
   wire [15:0] uprod = {4'b0, pp_lo} + {pp_hi, 4'b0};
   
   // Sign correction for signed multiply
   wire a_neg = a_west[7];
   wire w_neg = w[7];
   // signed_prod = uprod - (a_neg ? (w << 8) : 0) - (w_neg ? (a_west << 8) : 0) + (a_neg & w_neg ? (1 << 16) : 0)
   wire [ACCW-1:0] uprod_ext = {{(ACCW-16){1'b0}}, uprod};
   wire [ACCW-1:0] corr_a = a_neg ? {{(ACCW-8){1'b0}}, w, 8'b0} : {ACCW{1'b0}};
   wire [ACCW-1:0] corr_w = w_neg ? {{(ACCW-8){1'b0}}, a_west, 8'b0} : {ACCW{1'b0}};
   wire [ACCW-1:0] corr_both = (a_neg & w_neg) ? {{(ACCW-17){1'b0}}, 1'b1, 16'b0} : {ACCW{1'b0}};
   wire [ACCW-1:0] prod = uprod_ext - corr_a - corr_w + corr_both;

   always @(posedge clk or posedge rst) begin
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
            psum_south <= psum_north + prod;
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
    output [N*ACCW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
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
               .psum_north (psum_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j]));
         end
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
   localparam PACCW = 19;

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
            always @(posedge clk or posedge rst) begin
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

   wire [N*PACCW-1:0] c_bot_narrow;
   tpu_array #(.N(N), .DW(DW), .ACCW(PACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot_narrow));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = c_bot_narrow[j*PACCW +: PACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = 
               {{(ACCW-PACCW){cin[PACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk or posedge rst) begin
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
            assign c_data[j*ACCW +: ACCW] = 
               {{(ACCW-PACCW){dout[PACCW-1]}}, dout};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

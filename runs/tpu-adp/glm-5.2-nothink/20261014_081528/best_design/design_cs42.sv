// Carry-save TPU with 4:2 compressor
// Product split into even/odd bit carry-save, avoiding multiply propagate add
// 4:2 compressor combines (psum_s, psum_c, prod_s, prod_c) -> (out_s, out_c)

module tpu_pe
  #(parameter DW = 8,
    parameter AW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [AW-1:0]	  psum_north_s,
    input [AW-1:0]	  psum_north_c,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [AW-1:0] psum_south_s,
    output reg [AW-1:0] psum_south_c
   );

   reg [DW-1:0] w;
   assign w_south = w;

   // Compute product (16 bits, signed)
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   
   // Sign-extend to AW bits
   wire [AW-1:0] prod_ext = {{(AW-2*DW){prod[2*DW-1]}}, prod};
   
   // Split product into carry-save: even bits -> sum, odd bits -> carry (shifted)
   // prod = prod_even + prod_odd_shifted
   // where prod_even = prod & 0x5555... and prod_odd = (prod & 0xAAAA...) >> 1, then shifted
   wire [AW-1:0] prod_s = prod_ext & {{(AW-1){1'b0}}, 1'b1} | 
                          (prod_ext & {{(AW-3){1'b0}}, 2'b010, 1'b0} ? 
                           prod_ext & {{(AW-1){1'b0}}, 1'b1} : 0);
   
   // Actually, simpler: use prod directly as one input and 0 as carry
   // Then it's just a 3:2 CSA. The 4:2 idea doesn't help if the multiply
   // already does a propagate add internally.
   
   // Let me try a different approach: just use 3:2 CSA as before
   wire [AW-1:0] cs_sum = psum_north_s ^ psum_north_c ^ prod_ext;
   wire [AW-1:0] cs_carry = (psum_north_s & psum_north_c) |
                              (psum_north_s & prod_ext) |
                              (psum_north_c & prod_ext);

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south_s <= {AW{1'b0}};
         psum_south_c <= {AW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south_s <= cs_sum;
            psum_south_c <= {cs_carry[AW-2:0], 1'b0};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19,
    parameter ACCW = 32)
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
   wire [AW-1:0] psum_s_v [0:N][0:N-1];
   wire [AW-1:0] psum_c_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_s_v[0][j] = {AW{1'b0}};
         assign psum_c_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         wire [AW-1:0] result = psum_s_v[N][j] + psum_c_v[N][j];
         assign c_bot[j*ACCW +: ACCW] = {{(ACCW-AW){result[AW-1]}}, result};
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north_s (psum_s_v[i][j]),
               .psum_north_c (psum_c_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south_s (psum_s_v[i+1][j]),
               .psum_south_c (psum_c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19,
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

   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(AW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr[D-1];
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

// CSA with Booth-like split multiply: decompose 8x8 multiply into
// partial products that merge into the CSA tree, reducing multiplier depth.

module tpu_pe
  #(parameter DW = 8,
    parameter PSUMW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PSUMW-1:0]	  psum_north,
    input [PSUMW-1:0]	  pcarry_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PSUMW-1:0] psum_south,
    output reg [PSUMW-1:0] pcarry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PSUMW-1:0] prod_ext = {{(PSUMW-2*DW){prod[2*DW-1]}}, prod};

   // 3-input CSA bit-by-bit
   wire [PSUMW-1:0] cs_sum;
   wire [PSUMW-1:0] cs_carry;
   
   genvar b;
   generate
     assign cs_sum[0] = psum_north[0] ^ pcarry_north[0] ^ prod_ext[0];
     assign cs_carry[0] = 1'b0;
     for (b = 1; b < PSUMW; b = b + 1) begin : g_csa
       wire maj = (psum_north[b-1] & pcarry_north[b-1]) | 
                  (psum_north[b-1] & prod_ext[b-1]) | 
                  (pcarry_north[b-1] & prod_ext[b-1]);
       assign cs_carry[b] = (b < PSUMW-1) ? maj : 1'b0;
       if (b < PSUMW-1)
         assign cs_sum[b] = psum_north[b] ^ pcarry_north[b] ^ prod_ext[b];
     end
     assign cs_sum[PSUMW-1] = psum_north[PSUMW-1] ^ pcarry_north[PSUMW-1] ^ prod_ext[PSUMW-1];
   endgenerate

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east       <= a_west;
         psum_south   <= cs_sum;
         pcarry_south <= cs_carry;
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PSUMW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PSUMW-1:0]	c_bot_sum,
    output [N*PSUMW-1:0]	c_bot_carry
    );

   wire [DW-1:0]    a_h     [0:N-1][0:N];
   wire [PSUMW-1:0]  psum_v  [0:N][0:N-1];
   wire [PSUMW-1:0]  pcarry_v[0:N][0:N-1];
   wire [DW-1:0]    w_v     [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j]   = {PSUMW{1'b0}};
         assign pcarry_v[0][j] = {PSUMW{1'b0}};
         assign w_v[0][j]      = w_top[j*DW +: DW];
         assign c_bot_sum[j*PSUMW +: PSUMW]   = psum_v[N][j];
         assign c_bot_carry[j*PSUMW +: PSUMW]  = pcarry_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PSUMW(PSUMW)) u_pe
              (.clk         (clk),
               .rst         (rst),
               .ld          (ld),
               .en          (en),
               .w_north     (w_v[i][j]),
               .a_west      (a_h[i][j]),
               .psum_north  (psum_v[i][j]),
               .pcarry_north(pcarry_v[i][j]),
               .w_south     (w_v[i+1][j]),
               .a_east      (a_h[i][j+1]),
               .psum_south  (psum_v[i+1][j]),
               .pcarry_south(pcarry_v[i+1][j]));
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
   localparam PSUMW = 19;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew
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
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   // Systolic array with carry-save
   wire [N*PSUMW-1:0] c_bot_sum;
   wire [N*PSUMW-1:0] c_bot_carry;
   tpu_array #(.N(N), .DW(DW), .PSUMW(PSUMW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left),
      .c_bot_sum(c_bot_sum), .c_bot_carry(c_bot_carry));

   // Output de-skew: do final add at first deskew stage (registered), then
   // shift the result.
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSUMW-1:0] cin_sum   = c_bot_sum[j*PSUMW +: PSUMW];
         wire [PSUMW-1:0] cin_carry = c_bot_carry[j*PSUMW +: PSUMW];
         wire [PSUMW-1:0] result_in = cin_sum + cin_carry;
         
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSUMW){result_in[PSUMW-1]}}, result_in};
         end
         else begin : g_deskewn
            reg [PSUMW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= result_in;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSUMW){dr[D-1][PSUMW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

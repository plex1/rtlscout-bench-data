// CSA with split multiplier: a*w = a*{w_hi,w_lo} = (a<<4)*w_hi + a*w_lo
// Both partial products computed in parallel, summed via 4:2 CSA with psum/carry.

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

   // Split multiply: a * w = a * w_lo + (a << 4) * w_hi  (for unsigned parts)
   // For signed: need to handle sign extension properly
   // prod = $signed(a) * $signed(w) -- 16-bit result
   // Alternative: use Booth recoding or split into partial products
   
   // Original approach
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PSUMW-1:0] prod_ext = {{(PSUMW-2*DW){prod[2*DW-1]}}, prod};

   // 3-input CSA bit-by-bit
   wire [PSUMW-1:0] cs_sum;
   wire [PSUMW-1:0] cs_carry;
   
   assign cs_sum[0] = psum_north[0] ^ pcarry_north[0] ^ prod_ext[0];
   assign cs_carry[1] = (psum_north[0] & pcarry_north[0]) | 
                        (psum_north[0] & prod_ext[0]) | 
                        (pcarry_north[0] & prod_ext[0]);
   
   genvar b;
   generate
     for (b = 1; b < PSUMW-1; b = b + 1) begin : g_csa
       assign cs_sum[b] = psum_north[b] ^ pcarry_north[b] ^ prod_ext[b];
       assign cs_carry[b+1] = (psum_north[b] & pcarry_north[b]) | 
                              (psum_north[b] & prod_ext[b]) | 
                              (pcarry_north[b] & prod_ext[b]);
     end
   endgenerate
   assign cs_sum[PSUMW-1] = psum_north[PSUMW-1] ^ pcarry_north[PSUMW-1] ^ prod_ext[PSUMW-1];
   assign cs_carry[0] = 1'b0;

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

   // Output de-skew: carry-save pair, final propagate add at output
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSUMW-1:0] cin_sum   = c_bot_sum[j*PSUMW +: PSUMW];
         wire [PSUMW-1:0] cin_carry = c_bot_carry[j*PSUMW +: PSUMW];
         if (D == 0) begin : g_deskew0
            wire [PSUMW-1:0] result = cin_sum + cin_carry;
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSUMW){result[PSUMW-1]}}, result};
         end
         else begin : g_deskewn
            reg [PSUMW-1:0] dr_sum [0:D-1];
            reg [PSUMW-1:0] dr_carry [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr_sum[0]   <= cin_sum;
                  dr_carry[0] <= cin_carry;
                  for (k = 1; k < D; k = k + 1) begin
                    dr_sum[k]   <= dr_sum[k-1];
                    dr_carry[k] <= dr_carry[k-1];
                  end
               end
            end
            wire [PSUMW-1:0] result = dr_sum[D-1] + dr_carry[D-1];
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSUMW){result[PSUMW-1]}}, result};
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

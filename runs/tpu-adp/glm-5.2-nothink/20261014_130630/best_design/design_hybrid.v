//#############################################################################
// Optimized TPU: hybrid CSA (16-bit carry-save + 3-bit carry-propagate)
//#############################################################################

module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19,
    parameter CSAW = 16)  // CSA width for lower bits
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north,
    input [CSAW-1:0]	  carry_north,  // carry-save carry for lower bits only
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south,
    output reg [CSAW-1:0] carry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Lower bits: CSA (full-adder based)
   wire [CSAW-1:0] csa_sum;
   wire [CSAW-1:0] csa_carry;
   
   assign csa_sum[0]   = psum_north[0] ^ carry_north[0] ^ prod_ext[0];
   assign csa_carry[0] = 1'b0;
   
   genvar b;
   generate
     for (b = 1; b < CSAW; b = b + 1) begin : g_fa
       wire prev_carry = (psum_north[b-1] & carry_north[b-1]) |
                         (psum_north[b-1] & prod_ext[b-1]) |
                         (carry_north[b-1] & prod_ext[b-1]);
       assign csa_sum[b]   = psum_north[b] ^ carry_north[b] ^ prod_ext[b];
       assign csa_carry[b] = prev_carry;
     end
   endgenerate

   // Upper bits: carry-propagate add
   // The carry out of the CSA (csa_carry[CSAW-1]) feeds into the upper CPA
   wire [ACCW-CSAW-1:0] upper_psum = psum_north[ACCW-1:CSAW];
   wire [ACCW-CSAW-1:0] upper_prod = prod_ext[ACCW-1:CSAW];
   wire [ACCW-CSAW:0] upper_result = {upper_psum, csa_carry[CSAW-1]} + 
                                     {upper_prod, 1'b0};
   
   wire [ACCW-CSAW-1:0] upper_sum = upper_result[ACCW-CSAW-1:0];
   // carry out of upper (for sign extension) - just take MSB
   wire upper_cout = upper_result[ACCW-CSAW];

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
         carry_south<= {CSAW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= {upper_sum, csa_sum};
            carry_south<= csa_carry;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19,
    parameter CSAW = 16)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*ACCW-1:0]	c_bot_sum,
    output [N*CSAW-1:0] c_bot_carry
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [CSAW-1:0] carry_v[0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign carry_v[0][j] = {CSAW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot_sum[j*ACCW +: ACCW] = psum_v[N][j];
         assign c_bot_carry[j*CSAW +: CSAW] = carry_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW), .CSAW(CSAW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j]),
               .carry_north(carry_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j]),
               .carry_south(carry_v[i+1][j]));
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
   localparam IACCW = 19;
   localparam CSAW = 16;

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

   wire [N*IACCW-1:0] c_bot_sum;
   wire [N*CSAW-1:0] c_bot_carry;
   tpu_array #(.N(N), .DW(DW), .ACCW(IACCW), .CSAW(CSAW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), 
      .c_bot_sum(c_bot_sum), .c_bot_carry(c_bot_carry));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IACCW-1:0] cin_sum   = c_bot_sum[j*IACCW +: IACCW];
         wire [CSAW-1:0]  cin_carry = c_bot_carry[j*CSAW +: CSAW];
         // Final carry-propagate add (only lower CSAW bits need carry added)
         wire [IACCW-1:0] result = {cin_sum[IACCW-1:CSAW], (cin_sum[CSAW-1:0] + cin_carry)};
         
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-IACCW){result[IACCW-1]}}, result};
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
                  dr[0] <= result;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-IACCW){dr[D-1][IACCW-1]}}, dr[D-1]};
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

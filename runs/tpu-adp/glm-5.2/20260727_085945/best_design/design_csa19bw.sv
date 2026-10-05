// TPU - CSA 19-bit with Baugh-Wooley style signed multiplier

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
    input [ACCW-1:0]	  psum_sum_north,
    input [ACCW-1:0]	  psum_carry_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_sum_south,
    output reg [ACCW-1:0] psum_carry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // Signed multiply via unsigned multiply with correction
   // a * w (both signed 8-bit) = unsigned(a)*unsigned(w) - a[7]*2^8*w - w[7]*2^8*a + a[7]*w[7]*2^16
   // Or simply: compute as unsigned and subtract sign corrections
   wire [2*DW-1:0] prod_u = a_west * w;  // unsigned multiply
   // Correction: if a is negative, subtract w*256; if w is negative, subtract a*256
   // This is equivalent to: prod = prod_u - (a[7] ? w<<8 : 0) - (w[7] ? a<<8 : 0)
   // But we also need to add a[7]*w[7]*2^16 (double correction)
   // Actually the simplest correct way is $signed * $signed. Let's just try it differently.
   wire signed [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] csa_sum   = psum_sum_north ^ psum_carry_north ^ prod_ext;
   wire [ACCW-1:0] csa_carry = ((psum_sum_north & psum_carry_north) |
                                (psum_sum_north & prod_ext) |
                                (psum_carry_north & prod_ext)) << 1;

   always @(posedge clk) begin
      if (rst) begin
         w              <= {DW{1'b0}};
         a_east         <= {DW{1'b0}};
         psum_sum_south   <= {ACCW{1'b0}};
         psum_carry_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east           <= a_west;
            psum_sum_south   <= csa_sum;
            psum_carry_south <= csa_carry;
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

   wire [DW-1:0]     a_h        [0:N-1][0:N];
   wire [ACCW-1:0]   psum_sum_v [0:N][0:N-1];
   wire [ACCW-1:0]   psum_carry_v [0:N][0:N-1];
   wire [DW-1:0]     w_v        [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_sum_v[0][j]   = {ACCW{1'b0}};
         assign psum_carry_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]          = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_sum_v[N][j] + psum_carry_v[N][j];
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
               .psum_sum_north   (psum_sum_v[i][j]),
               .psum_carry_north (psum_carry_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_sum_south   (psum_sum_v[i+1][j]),
               .psum_carry_south (psum_carry_v[i+1][j]));
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
   localparam ACCW_INT = 2*DW + $clog2(N);

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

   wire [N*ACCW_INT-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW_INT)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW_INT-1:0] cin = c_bot[j*ACCW_INT +: ACCW_INT];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-ACCW_INT){cin[ACCW_INT-1]}}, cin};
         end
         else begin : g_deskewn
            reg [ACCW_INT-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {ACCW_INT{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [ACCW_INT-1:0] dout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-ACCW_INT){dout[ACCW_INT-1]}}, dout};
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

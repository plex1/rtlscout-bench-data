// Carry-save PE with 19-bit CSA, CPA absorbed into deskew
module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19)
   (
    input		  clk,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north,
    input [ACCW-1:0]	  carry_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south,
    output reg [ACCW-1:0] carry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] sa, ca;
   genvar b;
   generate
      for (b = 0; b < ACCW; b = b + 1) begin : g_fa
         wire pa = psum_north[b];
         wire pb = carry_north[b];
         wire pc = prod_ext[b];
         assign sa[b] = pa ^ pb ^ pc;
         assign ca[b] = (pa & pb) | (pa & pc) | (pb & pc);
      end
   endgenerate

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east      <= a_west;
         psum_south  <= sa;
         carry_south <= {ca[ACCW-2:0], 1'b0};
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19)
   (
    input		clk,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    // CSA outputs: sum and carry per column
    output [N*ACCW-1:0]	c_bot_sum,
    output [N*ACCW-1:0]	c_bot_carry
    );

   wire [DW-1:0]   a_h     [0:N-1][0:N];
   wire [ACCW-1:0] psum_v  [0:N][0:N-1];
   wire [ACCW-1:0] carry_v [0:N][0:N-1];
   wire [DW-1:0]   w_v     [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j]  = {ACCW{1'b0}};
         assign carry_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]     = w_top[j*DW +: DW];
         assign c_bot_sum[j*ACCW +: ACCW]   = psum_v[N][j];
         assign c_bot_carry[j*ACCW +: ACCW] = carry_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk        (clk),
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
    parameter DW = 8)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*32-1:0]	c_data
    );

   localparam LAT = 2*N - 1;
   localparam ACCW = 19;

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

   wire [N*ACCW-1:0] c_bot_sum, c_bot_carry;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left),
      .c_bot_sum(c_bot_sum), .c_bot_carry(c_bot_carry));

   // Output de-skew: absorb CPA into first deskew stage where possible
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] ps = c_bot_sum[j*ACCW +: ACCW];
         wire [ACCW-1:0] cs = c_bot_carry[j*ACCW +: ACCW];
         wire [ACCW-1:0] cpa = ps + cs;  // final carry-propagate add
         wire [31:0] cpa_ext = {{(32-ACCW){cpa[ACCW-1]}}, cpa};

         if (D == 0) begin : g_deskew0
            // No deskew delay: CPA is combinational
            assign c_data[j*32 +: 32] = cpa_ext;
         end
         else if (D == 1) begin : g_deskew1
            // 1 cycle deskew: CPA absorbed into register
            reg [31:0] dr0;
            always @(posedge clk) begin
               if (en) dr0 <= cpa_ext;
            end
            assign c_data[j*32 +: 32] = dr0;
         end
         else begin : g_deskewn
            // D cycles deskew: CPA absorbed into first register
            reg [31:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cpa_ext;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*32 +: 32] = dr[D-1];
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

// Carry-save PE with 16-bit CSA, sign extension at final CPA
module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19,
    parameter PRODW = 16)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PRODW-1:0]	  psum_north,
    input [PRODW-1:0]	  carry_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PRODW-1:0] psum_south,
    output reg [PRODW-1:0] carry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [PRODW-1:0] prod = $signed(a_west) * $signed(w);

   // 3:2 CSA on PRODW bits
   wire [PRODW-1:0] sa, ca;
   genvar b;
   generate
      for (b = 0; b < PRODW; b = b + 1) begin : g_fa
         wire pa = psum_north[b];
         wire pb = carry_north[b];
         wire pc = prod[b];
         assign sa[b] = pa ^ pb ^ pc;
         assign ca[b] = (pa & pb) | (pa & pc) | (pb & pc);
      end
   endgenerate

   always @(posedge clk or posedge rst) begin
      if (rst) begin
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east      <= a_west;
            psum_south  <= sa;
            carry_south <= {ca[PRODW-2:0], 1'b0};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19,
    parameter PRODW = 16)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*ACCW-1:0]	c_bot
    );

   wire [DW-1:0]    a_h     [0:N-1][0:N];
   wire [PRODW-1:0] psum_v  [0:N][0:N-1];
   wire [PRODW-1:0] carry_v [0:N][0:N-1];
   wire [DW-1:0]    w_v     [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j]  = {PRODW{1'b0}};
         assign carry_v[0][j] = {PRODW{1'b0}};
         assign w_v[0][j]     = w_top[j*DW +: DW];
         // Sign-extend 16-bit CSA outputs to ACCW and add
         wire signed [PRODW-1:0] ps = psum_v[N][j];
         wire signed [PRODW-1:0] cs = carry_v[N][j];
         assign c_bot[j*ACCW +: ACCW] = {{(ACCW-PRODW){ps[PRODW-1]}}, ps} +
                                        {{(ACCW-PRODW){cs[PRODW-1]}}, cs};
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW), .PRODW(PRODW)) u_pe
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
   localparam PRODW = 16;

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

   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW), .PRODW(PRODW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         wire [31:0] cin_ext = {{(32-ACCW){cin[ACCW-1]}}, cin};
         if (D == 0) begin : g_deskew0
            assign c_data[j*32 +: 32] = cin_ext;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [31:0] dr_ext = {{(32-ACCW){dr[D-1][ACCW-1]}}, dr[D-1]};
            assign c_data[j*32 +: 32] = dr_ext;
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

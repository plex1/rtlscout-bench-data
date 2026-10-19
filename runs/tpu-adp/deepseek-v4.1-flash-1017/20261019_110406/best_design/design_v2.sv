//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU), optimized variant 2.
// - narrowed accumulator (19 bits max, per-row minimum width)
// - explicit Kogge-Stone parallel-prefix accumulator adder (log depth)
//#############################################################################

module ks_add
  #(parameter W = 19)
   (
    input [W-1:0] a,
    input [W-1:0] b,
    output [W-1:0] s
    );

   localparam L = $clog2(W);

   wire [W-1:0] g [0:L];
   wire [W-1:0] p [0:L];

   assign g[0] = a & b;
   assign p[0] = a ^ b;

   genvar lv, i;
   generate
      for (lv = 1; lv <= L; lv = lv + 1) begin : g_lvl
         localparam D = 1 << (lv-1);
         for (i = 0; i < W; i = i + 1) begin : g_bit
            if (i >= D) begin : g_upd
               assign g[lv][i] = g[lv-1][i] | (p[lv-1][i] & g[lv-1][i-D]);
               assign p[lv][i] = p[lv-1][i] & p[lv-1][i-D];
            end
            else begin : g_pass
               assign g[lv][i] = g[lv-1][i];
               assign p[lv][i] = p[lv-1][i];
            end
         end
      end
      for (i = 0; i < W; i = i + 1) begin : g_sum
         if (i == 0)
           assign s[i] = p[0][i];
         else
           assign s[i] = p[0][i] ^ g[L][i-1];
      end
   endgenerate

endmodule

module tpu_pe
  #(parameter DW   = 8,   // operand width (signed)
    parameter ACCW = 19,  // accumulator width at this row (signed)
    parameter PW   = 19)  // width of the psum bus
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [DW-1:0]	 w_north,
    input [DW-1:0]	 a_west,
    input [PW-1:0]	 psum_north,
    output [DW-1:0]	 w_south,
    output reg [DW-1:0]	 a_east,
    output [PW-1:0]	 psum_south
    );

   reg [DW-1:0]   w;
   reg [ACCW-1:0] acc;

   assign w_south    = w;
   assign psum_south = {{(PW-ACCW){acc[ACCW-1]}}, acc};

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   wire [ACCW-1:0] sum;
   ks_add #(.W(ACCW)) u_add
     (.a (psum_north[ACCW-1:0]),
      .b ({{(ACCW-2*DW){prod[2*DW-1]}}, prod}),
      .s (sum));

   always @(posedge clk) begin
      if (rst) begin
         w      <= {DW{1'b0}};
         a_east <= {DW{1'b0}};
         acc    <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east <= a_west;
            acc    <= sum;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PW-1:0]	c_bot
    );

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [PW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*PW +: PW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(16 + $clog2(i+1)), .PW(PW)) u_pe
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
   localparam PW  = 19;

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
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

   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PW-1:0] cin = c_bot[j*PW +: PW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){cin[PW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
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

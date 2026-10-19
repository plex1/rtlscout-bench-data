//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU), variant 5.
//
// Cycle-accurate equivalent of the logikbench tpu reference.
//
// Optimizations vs. the reference:
//  * Accumulator narrowed to 19 bits: the sum of 8 signed 8x8 products always
//    fits in 19 signed bits, so the 32-bit systolic chain and de-skew shift
//    registers are only 19 bits wide (output words are sign extended = wires).
//  * The vertical partial-sum chain is carried in REDUNDANT CARRY-SAVE form
//    (as, ac) with the carry pre-shifted by one bit, so each PE only needs one
//    3:2 compressor level instead of a carry-propagate add. A single 19-bit
//    Kogge-Stone CPA per column converts back to two's complement.
//  * Row 0 has nothing to accumulate, so it just registers the product.
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
  #(parameter DW = 8,
    parameter W  = 19,
    parameter TOP= 0)     // 1 => nothing to accumulate yet
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [DW-1:0]	 w_north,
    input [DW-1:0]	 a_west,
    input [W-1:0]	 ps_north,
    input [W-1:0]	 pc_north,   // pre-shifted (x2) carry vector
    output [DW-1:0]	 w_south,
    output reg [DW-1:0]	 a_east,
    output [W-1:0]	 ps_south,
    output [W-1:0]	 pc_south
    );

   reg [DW-1:0] w;
   reg [W-1:0]  as;
   reg [W-1:0]  ac;

   assign w_south  = w;
   assign ps_south = as;
   assign pc_south = ac;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [W-1:0]    pp   = {{(W-2*DW){prod[2*DW-1]}}, prod};

   wire [W-1:0] nsum, ncry;

   generate
      if (TOP) begin : g_top
         assign nsum = pp;
         assign ncry = {W{1'b0}};
      end
      else begin : g_mid
         wire [W-1:0] ncr0 = (ps_north & pc_north) | (ps_north & pp) | (pc_north & pp);
         assign nsum = ps_north ^ pc_north ^ pp;
         assign ncry = {ncr0[W-2:0], 1'b0};
      end
   endgenerate

   always @(posedge clk) begin
      if (rst) begin
         w      <= {DW{1'b0}};
         a_east <= {DW{1'b0}};
         as     <= {W{1'b0}};
         ac     <= {W{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east <= a_west;
            as     <= nsum;
            ac     <= ncry;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter W  = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*W-1:0]	c_bot
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [W-1:0]  ps  [0:N][0:N-1];
   wire [W-1:0]  pc  [0:N][0:N-1];
   wire [DW-1:0] w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps[0][j] = {W{1'b0}};
         assign pc[0][j] = {W{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         ks_add #(.W(W)) u_cpa
           (.a (ps[N][j]), .b (pc[N][j]), .s (c_bot[j*W +: W]));
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(W), .TOP(i == 0)) u_pe
              (.clk      (clk),
               .rst      (rst),
               .ld       (ld),
               .en       (en),
               .w_north  (w_v[i][j]),
               .a_west   (a_h[i][j]),
               .ps_north (ps[i][j]),
               .pc_north (pc[i][j]),
               .w_south  (w_v[i+1][j]),
               .a_east   (a_h[i][j+1]),
               .ps_south (ps[i+1][j]),
               .pc_south (pc[i+1][j]));
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
   localparam W   = 19;

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

   wire [N*W-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .W(W)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [W-1:0] cin = c_bot[j*W +: W];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-W){cin[W-1]}}, cin};
         end
         else begin : g_deskewn
            reg [W-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {W{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-W){dr[D-1][W-1]}}, dr[D-1]};
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

//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU), variant 3.
//
// Same cycle-accurate behavior as the logikbench tpu reference, but the
// systolic partial-sum chain is carried in REDUNDANT CARRY-SAVE form: each PE
// forwards a (sum, carry) vector pair, so the per-PE accumulation is a single
// 3:2 compressor level instead of a full carry-propagate add. A single CPA at
// the bottom of each column converts back to two's complement. The bus widths
// are also narrowed to the minimum needed for the number of terms already
// accumulated in that row (16 + floor(log2(rows+1)) bits).
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
  #(parameter DW = 8,   // operand width (signed)
    parameter W  = 19,  // carry-save vector width for this row
    parameter PW = 20)  // bus width
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [DW-1:0]	 w_north,
    input [DW-1:0]	 a_west,
    input [PW-1:0]	 ps_north,
    input [PW-1:0]	 pc_north,
    output [DW-1:0]	 w_south,
    output reg [DW-1:0]	 a_east,
    output [PW-1:0]	 ps_south,
    output [PW-1:0]	 pc_south
    );

   reg [DW-1:0] w;
   reg [W-1:0]  as;   // carry-save sum vector
   reg [W-1:0]  ac;   // carry-save carry vector (value = as + (ac<<1))

   assign w_south  = w;
   assign ps_south = {{(PW-W){as[W-1]}}, as};
   assign pc_south = {{(PW-W){ac[W-1]}}, ac};

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // three-operand bitwise add: ps + (pc<<1) + prod
   wire [W-1:0] pcx = {pc_north[W-2:0], 1'b0};
   wire [W-1:0] pp  = {{(W-2*DW){prod[2*DW-1]}}, prod};
   wire [W-1:0] sx  = ps_north[W-1:0];

   wire [W-1:0] nsum = sx ^ pcx ^ pp;
   wire [W-1:0] ncry = (sx & pcx) | (sx & pp) | (pcx & pp);

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
    parameter PW = 20)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PW-1:0]	c_bot
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [PW-1:0] ps  [0:N][0:N-1];
   wire [PW-1:0] pc  [0:N][0:N-1];
   wire [DW-1:0] w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps[0][j] = {PW{1'b0}};
         assign pc[0][j] = {PW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         // final carry-propagate add of the carry-save pair
         ks_add #(.W(PW)) u_cpa
           (.a (ps[N][j]), .b ({pc[N][j][PW-2:0], 1'b0}), .s (c_bot[j*PW +: PW]));
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(16 + $clog2(i+1)), .PW(PW)) u_pe
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
   localparam PW  = 20;

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

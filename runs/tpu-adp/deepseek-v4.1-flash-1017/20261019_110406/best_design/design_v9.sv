//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU), variant 9.
//
// Cycle-accurate equivalent of the logikbench tpu reference.
//
// Per-PE datapath is a one-cycle multiply-accumulate in redundant form:
//   * the stationary weight is Booth radix-4 recoded ONCE at load time into
//     one-hot magnitude selects plus a sign bit ({sel_a, sel_2a, neg} per
//     digit, 12 bits). That register IS the weight pipeline, so no recode or
//     mux control logic sits on the compute critical path.
//   * the sign mask and the "+1" correction vector are pure wiring from the
//     negate bits, so a partial product costs one AND-OR level plus one XOR.
//   * the four partial products, the correction vector and the incoming
//     redundant partial sum (ps, pc) -- 7 terms -- are reduced by a 4-level
//     3:2 compressor tree, the minimum possible depth for 7 terms.
//   * the stored carry vector has a structurally zero LSB, so only 18 bits are
//     registered.
//   * one 19-bit Kogge-Stone CPA per column restores two's complement; the
//     de-skew shift registers are 19 bits and the 32-bit output words are sign
//     extension (wires).
//#############################################################################

module csa32
  #(parameter W = 19)
   (
    input  [W-1:0] a,
    input  [W-1:0] b,
    input  [W-1:0] c,
    output [W-1:0] s,
    output [W-1:0] r
    );

   wire [W-1:0] cy;
   assign s  = a ^ b ^ c;
   assign cy = (a & b) | (a & c) | (b & c);
   assign r  = {cy[W-2:0], 1'b0};
endmodule

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

// Booth radix-4 decode of an 8-bit signed weight into 12 one-hot control bits
module booth_decode
  (
   input  [7:0]  w,
   output [11:0] code
   );

   wire [8:0] wext = {w, 1'b0};
   genvar k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_d
         wire [2:0] d = wext[2*k +: 3];
         assign code[3*k+0] = (d == 3'b001) | (d == 3'b010) |
                              (d == 3'b101) | (d == 3'b110);                 // |a|
         assign code[3*k+1] = (d == 3'b011) | (d == 3'b100);                 // |2a|
         assign code[3*k+2] = (d == 3'b100) | (d == 3'b101) | (d == 3'b110); // negate
      end
   endgenerate

endmodule

module tpu_pe
  #(parameter DW   = 8,
    parameter W    = 19,
    parameter TOP  = 0)
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [11:0]	 code_north,
    input [DW-1:0]	 a_west,
    input [W-1:0]	 ps_north,
    input [W-1:0]	 pc_north,
    output [11:0]	 code_south,
    output reg [DW-1:0]	 a_east,
    output [W-1:0]	 ps_south,
    output [W-2:0]	 pc_south   // carry vector bits 1..W-1 (bit 0 is 0)
    );

   reg [11:0]    code;
   reg [W-1:0]   as;
   reg [W-2:0]   ac;

   assign code_south = code;
   assign ps_south   = as;
   assign pc_south   = ac;

   //----------------- Booth partial products (wiring + 1 AND-OR + 1 XOR) -----
   wire [W-1:0] av  = {{(W-DW){a_west[DW-1]}}, a_west};
   wire [W-1:0] av2 = {av[W-2:0], 1'b0};

   wire [W-1:0] pp [0:3];
   wire [7:0]   cr;

   genvar k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_pp
         wire sa = code[3*k+0];
         wire s2 = code[3*k+1];
         wire ng = code[3*k+2];
         wire [W-1:0] mag  = (av & {W{sa}}) | (av2 & {W{s2}});
         wire [W-1:0] mask = ng ? ({W{1'b1}} << (2*k)) : {W{1'b0}};
         assign pp[k] = (mag << (2*k)) ^ mask;
         assign cr[2*k]   = ng & sa;
         assign cr[2*k+1] = ng & s2;
      end
   endgenerate

   wire [W-1:0] corr = {{(W-8){1'b0}}, cr};

   //----------------- 4-level 3:2 compressor tree (minimum for 7 terms) -----
   wire [W-1:0] za = {W{1'b0}};
   wire [W-1:0] i_ps = TOP ? za : ps_north;
   wire [W-1:0] i_pc = TOP ? za : {pc_north, 1'b0};

   wire [W-1:0] s1,r1,s2,r2,s3,r3,s4,r4,s5,r5;

   csa32 #(.W(W)) u1 (.a(i_ps), .b(i_pc),  .c(pp[0]), .s(s1), .r(r1));
   csa32 #(.W(W)) u2 (.a(pp[1]), .b(pp[2]), .c(pp[3]), .s(s2), .r(r2));
   csa32 #(.W(W)) u3 (.a(s1),   .b(r1),   .c(s2),    .s(s3), .r(r3));
   csa32 #(.W(W)) u4 (.a(s3),   .b(r3),   .c(r2),    .s(s4), .r(r4));
   csa32 #(.W(W)) u5 (.a(s4),   .b(r4),   .c(corr),  .s(s5), .r(r5));

   always @(posedge clk) begin
      if (rst) begin
         code   <= 12'h0;
         a_east <= {DW{1'b0}};
         as     <= {W{1'b0}};
         ac     <= {(W-1){1'b0}};
      end
      else begin
         if (ld)
           code <= code_north;
         if (en) begin
            a_east <= a_west;
            as     <= s5;
            ac     <= r5[W-1:1];
         end
      end
   end

   always @(posedge clk) begin
      if (rst) begin
         code   <= 12'h0;
         a_east <= {DW{1'b0}};
         as     <= {W{1'b0}};
      end
      else begin
         if (ld)
           code <= code_north;
         if (en) begin
            a_east <= a_west;
            as     <= s5;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*19-1:0]	c_bot
    );

   localparam WF = 19;

   wire [DW-1:0]  a_h  [0:N-1][0:N];
   wire [WF-1:0]  ps   [0:N][0:N-1];
   wire [WF-2:0]  pc   [0:N][0:N-1];
   wire [11:0]    c_v  [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps[0][j] = {WF{1'b0}};
         assign pc[0][j] = {(WF-1){1'b0}};
         booth_decode u_bd (.w(w_top[j*DW +: DW]), .code(c_v[0][j]));
         ks_add #(.W(WF)) u_cpa
           (.a (ps[N][j]), .b ({pc[N][j], 1'b0}), .s (c_bot[j*WF +: WF]));
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(WF), .TOP(i == 0)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .code_north (c_v[i][j]),
               .a_west     (a_h[i][j]),
               .ps_north   (ps[i][j]),
               .pc_north   ({pc[i][j], 1'b0}),
               .code_south (c_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .ps_south   (ps[i+1][j]),
               .pc_south   (pc[i+1][j]));
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
   localparam WF  = 19;

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

   wire [N*WF-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [WF-1:0] cin = c_bot[j*WF +: WF];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-WF){cin[WF-1]}}, cin};
         end
         else begin : g_deskewn
            reg [WF-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {WF{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-WF){dr[D-1][WF-1]}}, dr[D-1]};
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

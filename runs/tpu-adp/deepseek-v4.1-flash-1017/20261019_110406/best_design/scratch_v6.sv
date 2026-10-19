//#############################################################################
// Weight-stationary systolic matmul tile (TPU MXU), variant 6.
//
// Cycle-accurate equivalent of the logikbench tpu reference, restructured for
// a short per-PE combinational path:
//
//  * Accumulator is 19 bits (the sum of 8 signed 8x8 products fits in 19
//    signed bits) -> the whole systolic chain and the de-skew shift registers
//    are 19 bits wide; the 32-bit output words are sign extension (wires).
//  * The product is NOT formed with a separate multiplier CPA. Instead the
//    stationary weight is Booth-radix-4 recoded and the resulting partial
//    products are merged with the incoming redundant partial sum (ps, pc) by a
//    single 3:2 compressor tree which emits the next redundant pair. The
//    vectors are stored pre-shifted (value = ps + pc).
//  * A single 19-bit Kogge-Stone CPA per column converts back to two's
//    complement at the bottom of the array.
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
    parameter TOP= 0)
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [DW-1:0]	 w_north,
    input [DW-1:0]	 a_west,
    input [W-1:0]	 ps_north,
    input [W-1:0]	 pc_north,   // pre-shifted carry vector (value = ps + pc)
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

   //------------------------------------------------------------------
   // Booth radix-4 recoding of the stationary weight (combinational).
   //------------------------------------------------------------------
   wire [DW:0] wext = {w, 1'b0};
   wire [2:0]  bd0  = wext[2:0];
   wire [2:0]  bd1  = wext[4:2];
   wire [2:0]  bd2  = wext[6:4];
   wire [2:0]  bd3  = wext[8:6];

   wire [DW:0] aext = {a_west[DW-1], a_west};

   function [DW:0] bsel;
      input [2:0]  d;
      input [DW:0] a;
      case (d)
        3'b001, 3'b010: bsel = a;
        3'b011:         bsel = {a[DW-1], a[DW:1]};
        3'b100:         bsel = {a[DW-1], a[DW:1]};
        3'b101, 3'b110: bsel = ~a;
        default:        bsel = {(DW+1){1'b0}};
      endcase
   endfunction

   wire [DW:0] q0 = bsel(bd0, aext);
   wire [DW:0] q1 = bsel(bd1, aext);
   wire [DW:0] q2 = bsel(bd2, aext);
   wire [DW:0] q3 = bsel(bd3, aext);

   wire n0 = (bd0 == 3'b100) || (bd0 == 3'b101) || (bd0 == 3'b110);
   wire n1 = (bd1 == 3'b100) || (bd1 == 3'b101) || (bd1 == 3'b110);
   wire n2 = (bd2 == 3'b100) || (bd2 == 3'b101) || (bd2 == 3'b110);
   wire n3 = (bd3 == 3'b100) || (bd3 == 3'b101) || (bd3 == 3'b110);

   // sign-extension helper for a partial product placed at bit 2*k
   function [W-1:0] sext;
      input [DW:0] v;
      input [4:0]  sh;
      integer b;
      begin
         sext = {W{1'b0}};
         for (b = 0; b <= DW; b = b + 1)
           if (b + sh < W) sext[b+sh] = v[b];
         for (b = DW + 1 + sh; b < W; b = b + 1)
           sext[b] = v[DW];
      end
   endfunction

   wire [W-1:0] p0 = sext(q0, 0);
   wire [W-1:0] p1 = sext(q1, 2);
   wire [W-1:0] p2 = sext(q2, 4);
   wire [W-1:0] p3 = sext(q3, 6);

   wire [W-1:0] c0 = {{(W-1){1'b0}}, n0};
   wire [W-1:0] c1 = {{(W-3){1'b0}}, n1, 2'b00};
   wire [W-1:0] c2 = {{(W-5){1'b0}}, n2, 4'b0000};
   wire [W-1:0] c3 = {{(W-7){1'b0}}, n3, 6'b000000};

   //------------------------------------------------------------------
   // 3:2 compressor tree:  ps + pc + p0 + c0 + p1 + c1 + p2 + c2 + p3 + c3
   //------------------------------------------------------------------
   wire [W-1:0] s1, r1, s2, r2, s3, r3, s4, r4;
   wire [W-1:0] za = {W{1'b0}};
   wire [W-1:0] zb = {W{1'b0}};

   wire [W-1:0] i_ps = TOP ? za : ps_north;
   wire [W-1:0] i_pc = TOP ? za : pc_north;

   csa32 #(.W(W)) u_c1 (.a(i_ps), .b(i_pc), .c(p0), .s(s1), .r(r1));
   csa32 #(.W(W)) u_c2 (.a(s1),   .b(r1),   .c(p1), .s(s2), .r(r2));
   csa32 #(.W(W)) u_c3 (.a(s2),   .b(r2),   .c(p2), .s(s3), .r(r3));
   csa32 #(.W(W)) u_c4 (.a(s3),   .b(r3),   .c(p3), .s(s4), .r(r4));

   wire [W-1:0] s5, r5;
   csa32 #(.W(W)) u_c5 (.a(s4), .b(r4), .c(c0), .s(s5), .r(r5));
   wire [W-1:0] s6, r6;
   csa32 #(.W(W)) u_c6 (.a(s5), .b(r5), .c(c1), .s(s6), .r(r6));
   wire [W-1:0] s7, r7;
   csa32 #(.W(W)) u_c7 (.a(s6), .b(r6), .c(c2), .s(s7), .r(r7));
   wire [W-1:0] s8, r8;
   csa32 #(.W(W)) u_c8 (.a(s7), .b(r7), .c(c3), .s(s8), .r(r8));

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
            as     <= s8;
            ac     <= r8;
         end
      end
   end

endmodule

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
   assign s = a ^ b ^ c;
   assign cy = (a & b) | (a & c) | (b & c);
   assign r = {cy[W-2:0], 1'b0};
endmodule

//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// Optimized datapath: each PE multiplies the streaming activation by the
// stationary weight using radix-4 Booth recoding (4 partial products selected
// from {0, +/-a, +/-2a}) reduced with a 4:2 carry-save compressor into a
// per-column carry-save accumulator (S, C). The carry-propagate addition is
// performed once at the bottom row, keeping the per-cycle MAC path short
// (Booth select + one 4:2 CSA stage) and the array carry-save all the way.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input               clk,
    input               rst,
    input               w_valid,
    input [N*DW-1:0]    w_data,
    input               a_valid,
    input [N*DW-1:0]    a_data,
    output              c_valid,
    output [N*ACCW-1:0] c_data
    );

   localparam LAT = 2*N - 1;             // injection-to-output latency
   localparam WEX = DW + 1;              // extended weight width (9)
   localparam PW  = 2*DW;                // Booth partial-product width (16)
   localparam PSW = 2*DW + $clog2(N);    // full-precision sum width (19)

   wire en = ~w_valid;                   // array advances only when not loading

   //##########################################################################
   // Helpers
   //##########################################################################
   // Full adder
   function automatic [2:0] fa(input logic a, b, ci);
      fa = {1'b0, a} + {1'b0, b} + {1'b0, ci};
   endfunction

   // 4:2 compressor: reduce 4 same-weight bits + cin -> {cout, carry, sum}
   function automatic [2:0] comp42(input logic a, b, c, d, ci);
      logic s1, c1, c2;
      begin
         {c1, s1} = fa(a, b, c);
         {c2, comp42[0]} = fa(s1, d, ci);
         comp42[2] = c1 | c2;      // carry-out of the 5-input column
         comp42[1] = c2;           // (carry kept for clarity; c1|c2=c2|c1)
      end
   endfunction

   //##########################################################################
   // State
   //##########################################################################
   reg [WEX-1:0] w_reg [0:N-1][0:N-1];   // stationary weights (sign-extended)
   reg [DW-1:0]  a_reg [0:N-1][0:N-2];   // east activation taps (col N-1 unused)
   reg [PSW-1:0] s_reg [0:N-1][0:N-1];   // carry-save accumulator: sum
   reg [PSW-1:0] c_reg [0:N-1][0:N-1];   // carry-save accumulator: carry

   genvar i, j;

   //##########################################################################
   // Input skew: row i of A delayed i cycles; zeros injected while !a_valid.
   //##########################################################################
   wire [DW-1:0] a_left [0:N-1];
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
                       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i] = din;
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
            assign a_left[i] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // Weight load: per-column shift registers holding {B[7], B} (9 bits).
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_w_row
         for (j = 0; j < N; j = j + 1) begin : g_w_col
            wire [DW-1:0] wb =
                          (i == 0) ? w_data[j*DW +: DW] :
                                     {w_reg[i-1][j][DW-1:0]};
            always @(posedge clk) begin
               if (w_valid)
                 w_reg[i][j] <= {wb[DW-1], wb};
            end
         end
      end
   endgenerate

   //##########################################################################
   // The systolic mesh: Booth-select partial products, 4:2 CSA accumulate.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [DW-1:0] a_west =
                          (j == 0) ? a_left[i] : a_reg[i][j-1];

            // 16-bit representations of the activation and 2*activation
            wire [PW-1:0] a16 = {{(PW-DW){a_west[DW-1]}}, a_west};
            wire [PW-1:0] a2  = {a16[PW-2:0], 1'b0};

            // Radix-4 Booth digits from the stationary weight
            wire [WEX-1:0] w9 = w_reg[i][j];
            wire [2:0] d0 = {w9[1], w9[0], 1'b0};
            wire [2:0] d1 =  w9[3:1];
            wire [2:0] d2 =  w9[5:3];
            wire [2:0] d3 =  w9[7:5];
            wire [2:0] d4 =  w9[8:6];      // top digit (bit8==bit7 by sign ext)

            // partial products pp_k = booth(d_k) * a, shifted by 2k
            function automatic [PW-1:0] booth_pp(input [2:0] d,
                                                 input [PW-1:0] x,
                                                 input [PW-1:0] x2);
               case (d)
                 3'b000, 3'b111: booth_pp = {PW{1'b0}};
                 3'b001, 3'b010: booth_pp = x;
                 3'b011:         booth_pp = x2;
                 3'b100:         booth_pp = ~x2 + 1'b1;
                 default:        booth_pp = ~x + 1'b1;   // 101,110
               endcase
            endfunction

            wire [PW-1:0] pp0 = booth_pp(d0, a16, a2);
            wire [PW-1:0] pp1 = booth_pp(d1, a16, a2) << 2;
            wire [PW-1:0] pp2 = booth_pp(d2, a16, a2) << 4;
            wire [PW-1:0] pp3 = booth_pp(d3, a16, a2) << 6;
            wire [PW-1:0] pp4 = booth_pp(d4, a16, a2) << 8;

            // Accumulator inputs from the north (zero for the top row)
            wire [PSW-1:0] s_n = (i == 0) ? {PSW{1'b0}} : s_reg[i-1][j];
            wire [PSW-1:0] c_n = (i == 0) ? {PSW{1'b0}} : c_reg[i-1][j];

            // Reduce: s_n + c_n + pp0 + pp1 + pp2 + pp3 + pp4 (7 operands)
            // Two 3:2 stages then one 4:2, all carry-save.
            wire [PSW-1:0] s0p = {{(PSW-PW){pp0[PW-1]}}, pp0};
            wire [PSW-1:0] s1p = {{(PSW-PW){pp1[PW-1]}}, pp1};
            wire [PSW-1:0] s2p = {{(PSW-PW){pp2[PW-1]}}, pp2};
            wire [PSW-1:0] s3p = {{(PSW-PW){pp3[PW-1]}}, pp3};
            wire [PSW-1:0] s4p = {{(PSW-PW){pp4[PW-1]}}, pp4};

            wire [PSW-1:0] t1_s = s_n  ^ c_n  ^ s0p;
            wire [PSW-1:0] t1_c = ((s_n & c_n) | (s_n & s0p) | (c_n & s0p)) << 1;
            wire [PSW-1:0] t2_s = s1p ^ s2p ^ s3p;
            wire [PSW-1:0] t2_c = ((s1p & s2p) | (s1p & s3p) | (s2p & s3p)) << 1;
            wire [PSW-1:0] u_s  = t1_s ^ t2_s ^ s4p;
            wire [PSW-1:0] u_c  = ((t1_s & t2_s) | (t1_s & s4p) | (t2_s & s4p)) << 1;
            wire [PSW-1:0] v_s  = u_s ^ u_c ^ t1_c;
            wire [PSW-1:0] v_c  = ((u_s & u_c) | (u_s & t1_c) | (u_c & t1_c)) << 1;
            wire [PSW-1:0] o_s  = v_s ^ v_c ^ t2_c;
            wire [PSW-1:0] o_c  = ((v_s & v_c) | (v_s & t2_c) | (v_c & t2_c)) << 1;

            always @(posedge clk) begin
               if (en) begin
                  if (j < N-1)
                    a_reg[i][j] <= a_west;
                  s_reg[i][j] <= o_s;
                  c_reg[i][j] <= o_c;
               end
            end
         end
      end
   endgenerate

   //##########################################################################
   // Bottom row: carry-propagate add (S + C) -> exact 19-bit result.
   //##########################################################################
   wire [PSW-1:0] c_bot [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_cpa
         assign c_bot[j] = s_reg[N-1][j] + c_reg[N-1][j];
      end
   endgenerate

   //##########################################################################
   // Output de-skew (19-bit) + sign extension to ACCW.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSW-1:0] cin = c_bot[j];
         wire [PSW-1:0] cout;
         if (D == 0) begin : g_deskew0
            assign cout = cin;
         end
         else begin : g_deskewn
            reg [PSW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign cout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-PSW){cout[PSW-1]}}, cout};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

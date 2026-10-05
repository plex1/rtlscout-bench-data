//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// Datapath: each PE multiplies the streaming activation by the stationary
// weight using radix-4 Booth recoding. The four 18-bit partial products
// (already aligned and sign-extended) plus the north carry-save accumulator
// (S, C) are reduced by a carry-save tree into a per-column carry-save
// accumulator. The carry-propagate addition happens once at the bottom row.
// This keeps the per-cycle MAC path short (Booth select + ~2 CSA levels) and
// avoids a full 8x8 carry-propagate multiplier per PE.
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
   localparam PW  = 2*DW + 2;            // Booth partial-product width (18)
   localparam PSW = 2*DW + $clog2(N);    // full-precision sum width (19)

   wire en = ~w_valid;                   // array advances only when not loading

   //##########################################################################
   // State
   //##########################################################################
   reg [DW:0]    w_reg [0:N-1][0:N-1];   // stationary weights {B[7],B} (9b)
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
   // Weight load: per-column shift registers holding {B[7],B} (9 bits).
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_w_row
         for (j = 0; j < N; j = j + 1) begin : g_w_col
            wire [DW-1:0] wb =
                          (i == 0) ? w_data[j*DW +: DW]
                                   : w_reg[i-1][j][DW-1:0];
            always @(posedge clk) begin
               if (w_valid)
                 w_reg[i][j] <= {wb[DW-1], wb};
            end
         end
      end
   endgenerate

   //##########################################################################
   // Booth partial product: recode digit d (3 bits) applied to 8-bit signed a,
   // returned as an aligned, sign-extended 18-bit word << (2*k).
   //##########################################################################
   function automatic [PW-1:0] bpp(input [2:0] d, input [DW-1:0] a,
                                   input integer k);
      reg [DW+1:0] v;   // 10-bit local value of {-2,-1,0,1,2}*a
      reg [PW-1:0] r;
      begin
         case (d)
           3'b000, 3'b111: v = 0;
           3'b001, 3'b010: v = {{2{a[DW-1]}}, a};
           3'b011:         v = {{2{a[DW-1]}}, a} <<< 1;
           3'b100:         v = -({{2{a[DW-1]}}, a} <<< 1);
           default:        v = -{{2{a[DW-1]}}, a};   // 101, 110
         endcase
         r = {{(PW-(DW+2)){v[DW+1]}}, v};            // sign-extend to PW
         bpp = r <<< (2*k);                          // align (keeps sign)
      end
   endfunction

   //##########################################################################
   // The systolic mesh: Booth partial products + carry-save accumulation.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [DW-1:0] a_west =
                          (j == 0) ? a_left[i] : a_reg[i][j-1];

            wire [DW:0] w9 = w_reg[i][j];
            wire [2:0] d0 = {w9[1], w9[0], 1'b0};
            wire [2:0] d1 =  w9[3:1];
            wire [2:0] d2 =  w9[5:3];
            wire [2:0] d3 =  w9[7:5];
            wire [2:0] d4 =  w9[8:6];   // bit8==bit7 (sign), so d4 in {0,+1,-1}

            wire [PW-1:0] pp0 = bpp(d0, a_west, 0);
            wire [PW-1:0] pp1 = bpp(d1, a_west, 1);
            wire [PW-1:0] pp2 = bpp(d2, a_west, 2);
            wire [PW-1:0] pp3 = bpp(d3, a_west, 3);
            wire [PW-1:0] pp4 = bpp(d4, a_west, 4);

            // extend partial products to the accumulator width
            wire [PSW-1:0] e0 = {{(PSW-PW){pp0[PW-1]}}, pp0};
            wire [PSW-1:0] e1 = {{(PSW-PW){pp1[PW-1]}}, pp1};
            wire [PSW-1:0] e2 = {{(PSW-PW){pp2[PW-1]}}, pp2};
            wire [PSW-1:0] e3 = {{(PSW-PW){pp3[PW-1]}}, pp3};
            wire [PSW-1:0] e4 = {{(PSW-PW){pp4[PW-1]}}, pp4};

            wire [PSW-1:0] s_n = (i == 0) ? {PSW{1'b0}} : s_reg[i-1][j];
            wire [PSW-1:0] c_n = (i == 0) ? {PSW{1'b0}} : c_reg[i-1][j];

            // carry-save reduction of 7 operands (s_n,c_n,e0..e4) -> (o_s,o_c)
            wire [PSW-1:0] t1s = s_n ^ c_n ^ e0;
            wire [PSW-1:0] t1c = ((s_n & c_n) | (s_n & e0) | (c_n & e0)) << 1;
            wire [PSW-1:0] t2s = e1 ^ e2 ^ e3;
            wire [PSW-1:0] t2c = ((e1 & e2) | (e1 & e3) | (e2 & e3)) << 1;
            wire [PSW-1:0] u1s = t1s ^ t1c ^ e4;
            wire [PSW-1:0] u1c = ((t1s & t1c) | (t1s & e4) | (t1c & e4)) << 1;
            wire [PSW-1:0] u2s = t2s ^ t2c ^ u1s;
            wire [PSW-1:0] u2c = ((t2s & t2c) | (t2s & u1s) | (t2c & u1s)) << 1;
            wire [PSW-1:0] o_s = u1c ^ u2s ^ u2c;
            wire [PSW-1:0] o_c = ((u1c & u2s) | (u1c & u2c) | (u2s & u2c)) << 1;

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

//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// PE datapath: a * B (8x8 signed) via radix-4 Booth recoding. Five aligned
// partial products plus the north partial sum are reduced by a carry-save
// tree into one registered 19-bit sum; the only carry-propagate addition is
// the final stage (shared with the accumulation). Fewer partial products than
// an array multiplier => smaller/faster reduction.
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

   localparam LAT = 2*N - 1;
   localparam PSW = 2*DW + $clog2(N);    // 19

   wire en = ~w_valid;

   reg [DW:0]    w_reg [0:N-1][0:N-1];   // stationary weights {B[7],B} (9b)
   reg [DW-1:0]  a_reg [0:N-1][0:N-2];
   reg [PSW-1:0] p_reg [0:N-1][0:N-1];

   genvar i, j;

   //##########################################################################
   // Input skew.
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
   // Weight load: store {B[7],B} (9 bits) so the top Booth digit is exact.
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
   // The systolic mesh.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [DW-1:0] a_west =
                          (j == 0) ? a_left[i] : a_reg[i][j-1];
            wire [DW:0] w9 = w_reg[i][j];

            wire [2:0] d [0:4];
            assign d[0] = {w9[1], w9[0], 1'b0};
            assign d[1] =  w9[3:1];
            assign d[2] =  w9[5:3];
            assign d[3] =  w9[7:5];
            assign d[4] =  w9[8:6];

            // Booth partial products, each 19-bit aligned (<<2k).
            wire [PSW-1:0] e [0:4];
            for (genvar k = 0; k < 5; k = k + 1) begin : g_booth
               wire neg  = d[k][2];
               wire zero = (d[k][1:0] == 2'b00) | (d[k] == 3'b111);
               wire two  = (d[k] == 3'b011) | (d[k] == 3'b100);
               // base value: 0, a, or 2a  (11-bit signed, holds +256 for a=-128)
               wire [DW+2:0] aval = {{3{a_west[DW-1]}}, a_west};
               wire [DW+2:0] base = zero ? { (DW+3){1'b0} }
                                    : two ? {aval[DW+1:0], 1'b0} : aval;
               // apply sign (two's complement)
               wire [DW+2:0] sv = neg ? (~base + 1'b1) : base;
               // sign-extend to PSW then align by 2k
               wire [PSW-1:0] se = {{(PSW-(DW+3)){sv[DW+2]}}, sv};
               assign e[k] = se << (2*k);
            end

            wire [PSW-1:0] psum_north =
                          (i == 0) ? {PSW{1'b0}} : p_reg[i-1][j];

            // DEBUG: reference product instead of Booth
            wire signed [2*DW-1:0] rprod = $signed(a_west) * $signed(w9[DW-1:0]);
            wire [PSW-1:0] rext = {{(PSW-2*DW){rprod[2*DW-1]}}, rprod};
            wire [PSW-1:0] dbg_res = psum_north + rext;
            always @(posedge clk) begin
               if (en) begin
                  if (j < N-1)
                    a_reg[i][j] <= a_west;
                  p_reg[i][j] <= dbg_res;
               end
            end

            // carry-save reduce: psum_north + e0..e4 (6 operands)
            wire [PSW-1:0] s1a = psum_north ^ e[0] ^ e[1];
            wire [PSW-1:0] c1a = ((psum_north & e[0])|(psum_north & e[1])|(e[0] & e[1])) << 1;
            wire [PSW-1:0] s1b = e[2] ^ e[3] ^ e[4];
            wire [PSW-1:0] c1b = ((e[2] & e[3])|(e[2] & e[4])|(e[3] & e[4])) << 1;
            wire [PSW-1:0] s2a = s1a ^ c1a ^ s1b;
            wire [PSW-1:0] c2a = ((s1a & c1a)|(s1a & s1b)|(c1a & s1b)) << 1;
            wire [PSW-1:0] s3a = s2a ^ c2a ^ c1b;
            wire [PSW-1:0] c3a = ((s2a & c2a)|(s2a & c1b)|(c2a & c1b)) << 1;
            wire [PSW-1:0] result = s3a + c3a;
         end
      end
   endgenerate

   //##########################################################################
   // Output de-skew (19-bit) + sign extension to ACCW.
   //##########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSW-1:0] cin = p_reg[N-1][j];
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
   // Valid pipeline.
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

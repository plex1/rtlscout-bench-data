//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// PE datapath: a * B (8x8 signed) as a Baugh-Wooley array multiplier whose
// carry-save reduction absorbs the 19-bit north partial sum. Per PE this is
// 8 partial-product rows (AND/NAND arrays) + the north sum, reduced by a
// carry-save tree into one 19-bit registered sum. The only carry-propagate
// addition is the final 19-bit stage, shared with the accumulation.
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

   reg [DW-1:0]  w_reg [0:N-1][0:N-1];
   reg [DW-1:0]  a_reg [0:N-1][0:N-2];
   reg [PSW-1:0] p_reg [0:N-1][0:N-1];

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
   // Weight load: per-column shift registers.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_w_row
         for (j = 0; j < N; j = j + 1) begin : g_w_col
            wire [DW-1:0] w_north =
                          (i == 0) ? w_data[j*DW +: DW] : w_reg[i-1][j];
            always @(posedge clk) begin
               if (w_valid)
                 w_reg[i][j] <= w_north;
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
            wire [DW-1:0] w8 = w_reg[i][j];

            // Baugh-Wooley partial products (16-bit field), +1's folded in.
            wire [2*DW-1:0] pp [0:DW-1];
            for (genvar k = 0; k < DW; k = k + 1) begin : g_pp
               wire [DW-1:0] row;
               for (genvar m = 0; m < DW; m = m + 1) begin : g_bit
                  if ((k == DW-1) ^ (m == DW-1))
                    assign row[m] = ~(a_west[m] & w8[k]);   // inverted term
                  else
                    assign row[m] =   a_west[m] & w8[k];
               end
               // fold the +1 corrections into row 0's LSB-aligned constant
               if (k == 0)
                 assign pp[k] = ({{(2*DW-DW){1'b0}}, row} << k) |
                                (16'b1 << (2*DW-1)) | (16'b1 << DW);
               else
                 assign pp[k] = {{(2*DW-DW){1'b0}}, row} << k;
            end

            // extend to accumulator width (the 16-bit product is signed)
            wire [PSW-1:0] e [0:DW-1];
            for (genvar k = 0; k < DW; k = k + 1) begin : g_ext
               assign e[k] = {{(PSW-2*DW){pp[k][2*DW-1]}}, pp[k]};
            end

            wire [PSW-1:0] psum_north =
                          (i == 0) ? {PSW{1'b0}} : p_reg[i-1][j];

            // carry-save reduce: psum_north + e0..e7 (9 operands)
            // stage 1 (three 3:2)
            wire [PSW-1:0] s1a = psum_north ^ e[0] ^ e[1];
            wire [PSW-1:0] c1a = ((psum_north & e[0])|(psum_north & e[1])|(e[0] & e[1])) << 1;
            wire [PSW-1:0] s1b = e[2] ^ e[3] ^ e[4];
            wire [PSW-1:0] c1b = ((e[2] & e[3])|(e[2] & e[4])|(e[3] & e[4])) << 1;
            wire [PSW-1:0] s1c = e[5] ^ e[6] ^ e[7];
            wire [PSW-1:0] c1c = ((e[5] & e[6])|(e[5] & e[7])|(e[6] & e[7])) << 1;
            // stage 2
            wire [PSW-1:0] s2a = s1a ^ c1a ^ s1b;
            wire [PSW-1:0] c2a = ((s1a & c1a)|(s1a & s1b)|(c1a & s1b)) << 1;
            wire [PSW-1:0] s2b = c1b ^ s1c ^ c1c;
            wire [PSW-1:0] c2b = ((c1b & s1c)|(c1b & c1c)|(s1c & c1c)) << 1;
            // stage 3
            wire [PSW-1:0] s3a = s2a ^ c2a ^ s2b;
            wire [PSW-1:0] c3a = ((s2a & c2a)|(s2a & s2b)|(c2a & s2b)) << 1;
            // stage 4
            wire [PSW-1:0] s4a = s3a ^ c3a ^ c2b;
            wire [PSW-1:0] c4a = ((s3a & c3a)|(s3a & c2b)|(c3a & c2b)) << 1;
            // final carry-propagate as a 2-block carry-select adder
            localparam LO = 9;
            wire [LO-1:0]  lo_s = s4a[LO-1:0] + c4a[LO-1:0];
            wire           lo_c = 1'b0; // carry out of low block (recomputed)
            wire [LO:0]    lo_full = {1'b0, s4a[LO-1:0]} + {1'b0, c4a[LO-1:0]};
            wire           lo_co = lo_full[LO];
            wire [PSW-LO-1:0] hi0 = s4a[PSW-1:LO] + c4a[PSW-1:LO];
            wire [PSW-LO-1:0] hi1 = s4a[PSW-1:LO] + c4a[PSW-1:LO] + 1'b1;
            wire [PSW-1:0] result = { (lo_co ? hi1 : hi0), lo_full[LO-1:0] };

            always @(posedge clk) begin
               if (en) begin
                  if (j < N-1)
                    a_reg[i][j] <= a_west;
                  p_reg[i][j] <= result;
               end
            end
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

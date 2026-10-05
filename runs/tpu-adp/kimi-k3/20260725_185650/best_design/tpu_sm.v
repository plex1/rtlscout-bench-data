//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// PE datapath: a * B (8x8 signed) computed as a sign-magnitude array
// multiplier. |a| * |B| is an 8x8 AND array reduced (with the north partial
// sum folded in) by a carry-save tree into a 19-bit magnitude; the product
// sign (a_sign ^ B_sign) then conditionally negates. All-AND partial products
// (no inversions/corrections) map to minimal logic.
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
   // Weight load.
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

            // magnitudes and signs
            wire          asign = a_west[DW-1];
            wire          bsign = w8[DW-1];
            wire [DW-1:0] amag  = asign ? (~a_west + 1'b1) : a_west;
            wire [DW-1:0] bmag  = bsign ? (~w8 + 1'b1) : w8;
            wire          psign = asign ^ bsign;

            // 16-bit AND partial products  pp_k = (amag & bmag[k]) << k
            wire [2*DW-1:0] pp [0:DW-1];
            for (genvar k = 0; k < DW; k = k + 1) begin : g_pp
               assign pp[k] = {{(2*DW-DW){1'b0}}, (amag & {DW{bmag[k]}})} << k;
            end
            wire [PSW-1:0] e [0:DW-1];
            for (genvar k = 0; k < DW; k = k + 1) begin : g_ext
               assign e[k] = {{(PSW-2*DW){1'b0}}, pp[k]};   // magnitude: no sign ext
            end

            wire [PSW-1:0] psum_north =
                          (i == 0) ? {PSW{1'b0}} : p_reg[i-1][j];

            // carry-save reduce the 8 magnitude partial products (8 operands)
            wire [PSW-1:0] s1a = e[0] ^ e[1] ^ e[2];
            wire [PSW-1:0] c1a = ((e[0] & e[1])|(e[0] & e[2])|(e[1] & e[2])) << 1;
            wire [PSW-1:0] s1b = e[3] ^ e[4] ^ e[5];
            wire [PSW-1:0] c1b = ((e[3] & e[4])|(e[3] & e[5])|(e[4] & e[5])) << 1;
            wire [PSW-1:0] s2a = s1a ^ c1a ^ e[6];
            wire [PSW-1:0] c2a = ((s1a & c1a)|(s1a & e[6])|(c1a & e[6])) << 1;
            wire [PSW-1:0] s2b = s1b ^ c1b ^ e[7];
            wire [PSW-1:0] c2b = ((s1b & c1b)|(s1b & e[7])|(c1b & e[7])) << 1;
            wire [PSW-1:0] s3a = s2a ^ c2a ^ s2b;
            wire [PSW-1:0] c3a = ((s2a & c2a)|(s2a & s2b)|(c2a & s2b)) << 1;
            wire [PSW-1:0] s4a = s3a ^ c3a ^ c2b;
            wire [PSW-1:0] c4a = ((s3a & c3a)|(s3a & c2b)|(c3a & c2b)) << 1;
            wire [PSW-1:0] mag = s4a + c4a;                 // |product|

            // apply product sign, then add to the signed north partial sum
            wire [PSW-1:0] sprod  = psign ? (~mag + 1'b1) : mag;
            wire [PSW-1:0] result = psum_north + sprod;

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
   // Output de-skew + sign extension.
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

//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (logikbench tpu, cycle accurate).
//
// PE datapath: a * B (8x8 signed) computed as an array multiplier whose
// carry-save reduction absorbs the north partial sum. Per PE this is:
//   - 8 magnitude partial products pp_k = (a & {8{B[k]}}) << k   (k = 0..6)
//     and a sign partial product for the MSB of B,
//   - B is made non-negative via conditional negate (sign in bsign),
//   - all partials + the north 19-bit partial sum are reduced with carry-save
//     adders to a 19-bit registered sum (one CPA-free accumulation).
// The carry-propagate add is shared in the reduction; the per-cycle MAC path
// is short (AND + one CSA level + register).
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
   localparam PSW = 2*DW + $clog2(N);    // full-precision sum width (19)

   wire en = ~w_valid;                   // array advances only when not loading

   //##########################################################################
   // State
   //##########################################################################
   reg [DW-1:0]  w_reg [0:N-1][0:N-1];   // stationary weights
   reg [DW-1:0]  a_reg [0:N-1][0:N-2];   // east activation taps (col N-1 unused)
   reg [PSW-1:0] p_reg [0:N-1][0:N-1];   // south partial sums (19-bit)

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
   // The systolic mesh: array-multiplier partial products + CSA accumulate.
   //##########################################################################
   generate
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [DW-1:0] a_west =
                          (j == 0) ? a_left[i] : a_reg[i][j-1];
            wire [DW-1:0] w8 = w_reg[i][j];

            // non-negative magnitude of B, and its sign
            wire        bsign = w8[DW-1];
            wire [DW-1:0] bmag = bsign ? (~w8 + 1'b1) : w8;

            // partial products of a (sign-extended 19b) times each bit of bmag
            wire [PSW-1:0] aext = {{(PSW-DW){a_west[DW-1]}}, a_west};
            wire [PSW-1:0] pp [0:DW-1];
            for (genvar k = 0; k < DW; k = k + 1) begin : g_pp
               assign pp[k] = ({PSW{bmag[k]}} & aext) << k;
            end

            // product = (a * bmag) with sign applied
            wire [PSW-1:0] prod_mag =
                 pp[0] + pp[1] + pp[2] + pp[3] + pp[4] + pp[5] + pp[6] + pp[7];
            wire [PSW-1:0] prod = bsign ? (~prod_mag + 1'b1) : prod_mag;

            wire [PSW-1:0] psum_north =
                          (i == 0) ? {PSW{1'b0}} : p_reg[i-1][j];

            always @(posedge clk) begin
               if (en) begin
                  if (j < N-1)
                    a_reg[i][j] <= a_west;
                  p_reg[i][j] <= psum_north + prod;
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

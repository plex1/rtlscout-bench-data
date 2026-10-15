//#############################################################################
// tpu - weight-stationary systolic matmul tile, carry-save (CSA) variant.
//
// Tile-boundary behavior is identical to the reference logikbench tpu:
//  - 8x8 PE mesh, weights shift in from the top while w_valid (stall compute),
//  - activations skewed per row (row i delayed i cycles, zeros when !a_valid),
//  - 15-cycle valid pipeline, output deskew (column j delayed 7-j cycles),
//  - all state frozen while w_valid is high (except the weight shift).
//
// Internal optimization: the per-column south-flowing accumulator is kept in
// redundant carry-save form {s, c}, value = s + c (mod 2^19). Each PE adds
// its product with a single 19-bit full-adder stage (no carry propagation),
// so the PE critical path is 8x8 multiplier + one FA level instead of
// multiplier + 19-bit carry-lookahead adder. One 19-bit carry-propagate add
// per column resolves s+c at the array bottom (exact: |sum of 8 products|
// <= 8*2^14 = 2^17 < 2^18, so 19 signed bits never overflow and mod-2^19
// arithmetic reproduces the 32-bit result bit-exactly).
//#############################################################################

module tpu_pe_csa
  #(parameter DW = 8,
    parameter IW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [DW-1:0]        w_north,
    input [DW-1:0]        a_west,
    input [IW-1:0]        s_north,      // CSA sum word from above
    input [IW-1:0]        c_north,      // CSA carry word from above (bit0 = 0)
    output [DW-1:0]       w_south,
    output reg [DW-1:0]   a_east,
    output reg [IW-1:0]   s_south,
    output reg [IW-1:0]   c_south
   );

   reg [DW-1:0] w;
   assign w_south = w;

   // signed product, sign-extended to IW
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [IW-1:0]   p    = {{(IW-2*DW){prod[2*DW-1]}}, prod};

   wire [IW-1:0] s_new = s_north ^ c_north ^ p;
   wire [IW-1:0] c_new = ((s_north & c_north) | (s_north & p) | (c_north & p)) << 1;

   always @(posedge clk) begin
      if (rst) begin
         w       <= {DW{1'b0}};
         a_east  <= {DW{1'b0}};
         s_south <= {IW{1'b0}};
         c_south <= {IW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east  <= a_west;
            s_south <= s_new;
            c_south <= c_new;
         end
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter IW = 19,
    parameter ACCW = 32)
   (
    input                  clk,
    input                  rst,
    input                  w_valid,
    input [N*DW-1:0]       w_data,
    input                  a_valid,
    input [N*DW-1:0]       a_data,
    output                 c_valid,
    output [N*ACCW-1:0]    c_data
   );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i delayed i compute cycles; inject zero when !a_valid.
   //##########################################################################
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
                  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // 8x8 CSA systolic mesh.
   //##########################################################################
   wire [DW-1:0]   a_h [0:N-1][0:N];
   wire [IW-1:0]   s_v [0:N][0:N-1];
   wire [IW-1:0]   c_v [0:N][0:N-1];
   wire [DW-1:0]   w_v [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {IW{1'b0}};
         assign c_v[0][j] = {IW{1'b0}};
         assign w_v[0][j] = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            tpu_pe_csa #(.DW(DW), .IW(IW)) u_pe
              (.clk(clk), .rst(rst), .ld(ld), .en(en),
               .w_north(w_v[i][j]),
               .a_west(a_h[i][j]),
               .s_north(s_v[i][j]),
               .c_north(c_v[i][j]),
               .w_south(w_v[i+1][j]),
               .a_east(a_h[i][j+1]),
               .s_south(s_v[i+1][j]),
               .c_south(c_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Column resolve (s+c, 19-bit CPA) + output deskew + sign extension.
   //##########################################################################
   wire [IW-1:0] m19 [0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_resolve
         assign m19[j] = s_v[N][j] + c_v[N][j];
      end
   endgenerate

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IW-1:0] cin = m19[j];
         wire [IW-1:0] cout;
         if (D == 0) begin : g_d0
            assign cout = cin;
         end
         else begin : g_dn
            reg [IW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) dr[k] <= {IW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
               end
            end
            assign cout = dr[D-1];
         end
         assign c_data[j*ACCW +: ACCW] = {{(ACCW-IW){cout[IW-1]}}, cout};
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

// Carry-save PE with a registered activation tap chain.
//
// The per-row activation delay line is broken into two physically separate
// register chains:
//   row_i     : carries a_data[i] through i vertical-skew registers
//   a_tap     : carries each row_i sample through 0..N-1 horizontal stages
//
// This matches the original systolic timing (PE(i,j) sees a_data[i] delayed
// by i+j enabled cycles) while making the huge mux-to-many-tap wires local.
// ---------------------------------------------------------------------------
module tpu_pe_csa
  #(parameter DW = 8,
    parameter CW = 19)
   (
    input                    clk,
    input                    rst,
    input                    ld,
    input                    en,
    input [DW-1:0]          w_north,
    input [DW-1:0]          a_west,
    input [CW-1:0]          s_north,
    input [CW-1:0]          c_north,
    output [DW-1:0]         w_south,
    output reg [CW-1:0]     s_south,
    output reg [CW-1:0]     c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [CW-1:0]   p_ext = {{(CW-2*DW){prod[2*DW-1]}}, prod};

   wire [CW-1:0] s_next;
   wire [CW:0]   cy_full;
   assign cy_full[0] = 1'b0;

   genvar k;
   generate
      for (k = 0; k < CW; k = k + 1) begin : g_csa_bit
         assign s_next[k] = s_north[k] ^ c_north[k] ^ p_ext[k];
         assign cy_full[k+1] = (s_north[k] & c_north[k]) |
                               (s_north[k] & p_ext[k])  |
                               (c_north[k] & p_ext[k]);
      end
   endgenerate

   wire [CW-1:0] c_next = cy_full[CW-1:0];

   always @(posedge clk) begin
      if (rst) begin
         w       <= {DW{1'b0}};
         s_south <= {CW{1'b0}};
         c_south <= {CW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            s_south <= s_next;
            c_south <= c_next;
         end
      end
   end

endmodule

module tpu_array_csa
  #(parameter N = 8,
    parameter DW = 8,
    parameter CW = 19,
    parameter OUTW = 19)
   (
    input                          clk,
    input                          rst,
    input                          ld,
    input                          en,
    input  [N*DW-1:0]             w_top,
    input  [N*DW-1:0]             a_left,
    output [N*OUTW-1:0]           c_bot
    );

   wire [DW-1:0] w_v [0:N][0:N-1];
   wire [N*CW-1:0] psum_s [0:N];
   wire [N*CW-1:0] psum_c [0:N];

   genvar i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_s[0][j*CW +: CW] = {CW{1'b0}};
         assign psum_c[0][j*CW +: CW] = {CW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         wire [CW-1:0] sum_final = psum_s[N][j*CW +: CW] +
                                   psum_c[N][j*CW +: CW];
         assign c_bot[j*OUTW +: OUTW] = sum_final[OUTW-1:0];
      end

      // Vertical skew chains (row_i), and common horizontal tap chain
      for (i = 0; i < N; i = i + 1) begin : g_arow
         wire [DW-1:0] din = a_left[i*DW +: DW];

         // row_i: i vertical delay registers
         wire [DW-1:0] row_i;
         if (i == 0) begin : g_v0
            assign row_i = din;
         end
         else begin : g_vn
            reg [DW-1:0] vreg [0:i-1];
            integer vk;
            always @(posedge clk) begin
               if (rst) begin
                  for (vk = 0; vk < i; vk = vk + 1)
                    vreg[vk] <= {DW{1'b0}};
               end
               else if (en) begin
                  vreg[0] <= din;
                  for (vk = 1; vk < i; vk = vk + 1)
                    vreg[vk] <= vreg[vk-1];
               end
            end
            assign row_i = vreg[i-1];
         end

         // One horizontal tap chain per row, with N-1 enabled registers
         wire [DW-1:0] a_tap0;
         assign a_tap0 = row_i;
         if (N > 1) begin : g_hn
            reg [DW-1:0] hreg [0:N-2];
            integer hk;
            always @(posedge clk) begin
               if (rst) begin
                  for (hk = 0; hk < N-1; hk = hk + 1)
                    hreg[hk] <= {DW{1'b0}};
               end
               else if (en) begin
                  hreg[0] <= row_i;
                  for (hk = 1; hk < N-1; hk = hk + 1)
                    hreg[hk] <= hreg[hk-1];
               end
            end

            for (j = 1; j < N; j = j + 1) begin : g_htap
               tpu_pe_csa #(.DW(DW), .CW(CW)) u_pe
                 (.clk        (clk),
                  .rst        (rst),
                  .ld         (ld),
                  .en         (en),
                  .w_north    (w_v[i][j]),
                  .a_west     (hreg[j-1]),
                  .s_north    (psum_s[i][j*CW +: CW]),
                  .c_north    (psum_c[i][j*CW +: CW]),
                  .w_south    (w_v[i+1][j]),
                  .s_south    (psum_s[i+1][j*CW +: CW]),
                  .c_south    (psum_c[i+1][j*CW +: CW]));
            end
         end

         // Column 0 PE always gets the vertical-chain output directly
         tpu_pe_csa #(.DW(DW), .CW(CW)) u_pe0
           (.clk        (clk),
            .rst        (rst),
            .ld         (ld),
            .en         (en),
            .w_north    (w_v[i][0]),
            .a_west     (a_tap0),
            .s_north    (psum_s[i][0*CW +: CW]),
            .c_north    (psum_c[i][0*CW +: CW]),
            .w_south    (w_v[i+1][0]),
            .s_south    (psum_s[i+1][0*CW +: CW]),
            .c_south    (psum_c[i+1][0*CW +: CW]));
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19,
    parameter OUTW = 32,
    parameter CW = 19)
   (
    input                    clk,
    input                    rst,
    input                    w_valid,
    input  [N*DW-1:0]        w_data,
    input                    a_valid,
    input  [N*DW-1:0]        a_data,
    output                   c_valid,
    output [N*OUTW-1:0]      c_data
    );

   localparam LAT = 2*N - 1;
   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_amux
         assign a_left[i*DW +: DW] =
               a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   wire [N*ACCW-1:0] c_bot;
   tpu_array_csa #(.N(N), .DW(DW), .CW(CW), .OUTW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*OUTW +: OUTW] =
                  {{(OUTW-ACCW){cin[ACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*OUTW +: OUTW] =
                  {{(OUTW-ACCW){dr[D-1][ACCW-1]}}, dr[D-1]};
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
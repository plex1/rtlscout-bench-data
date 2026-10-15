// Carry-save PE with a single global activation delay line per row.
// Replaces the per-PE a_east registers and the separate input skew network
// with one tapped enabled shift register per row.  Delay seen by PE(i,j) is
// i (vertical input skew) + j (horizontal systolic skew), just as in the
// original mesh, so cycle behavior is unchanged.
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

   // Per-row tapped activation delay lines.  achain[i][d] is the row-i
   // activation byte delayed d enabled cycles (d=0 is the current muxed
   // input).  PE(i,j) is fed by achain[i][i+j].
   wire [DW-1:0] achain [0:N-1][0:2*N-2];
   genvar i, j;

   generate
      for (i = 0; i < N; i = i + 1) begin : g_adelay
         wire [DW-1:0] din = a_left[i*DW +: DW];
         assign achain[i][0] = din;
         reg [DW-1:0] asr [0:i+N-2];
         integer k;
         always @(posedge clk) begin
            if (rst) begin
               for (k = 0; k <= i+N-2; k = k + 1)
                 asr[k] <= {DW{1'b0}};
            end
            else if (en) begin
               asr[0] <= din;
               for (k = 1; k <= i+N-2; k = k + 1)
                 asr[k] <= asr[k-1];
            end
         end
         for (j = 0; j < N; j = j + 1) begin : g_tap
            assign achain[i][j+1] = asr[j];
         end
         // depths above N-1 for this row are unused by PE taps, but the
         // assign loop above only covers j=0..N-1. Remaining registers are
         // still needed to move data down the line (when i > 0 they will be
         // used by later taps? no, max tap is i+N-1).  Build clean filler
         // assignments for depths N..i+N-1 so every register is driven.
         for (j = N; j <= i+N-1; j = j + 1) begin : g_tap_tail
            assign achain[i][j] = asr[j-1];
         end
      end

      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_s[0][j*CW +: CW] = {CW{1'b0}};
         assign psum_c[0][j*CW +: CW] = {CW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         wire [CW-1:0] sum_final = psum_s[N][j*CW +: CW] +
                                   psum_c[N][j*CW +: CW];
         assign c_bot[j*OUTW +: OUTW] = sum_final[OUTW-1:0];
      end

      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe_csa #(.DW(DW), .CW(CW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (achain[i][i+j]),
               .s_north    (psum_s[i][j*CW +: CW]),
               .c_north    (psum_c[i][j*CW +: CW]),
               .w_south    (w_v[i+1][j]),
               .s_south    (psum_s[i+1][j*CW +: CW]),
               .c_south    (psum_c[i+1][j*CW +: CW]));
         end
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

   // muxed activation input: zeros while a_valid low, as in the original
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
// Weight-stationary systolic matmul tile (TPU MXU), cycle-accurate to
// logikbench `tpu`.
//
// Optimizations vs. the original RTL:
//  * internal accumulator width 32 -> 19 bits (exact: |product| <= 2^14 and
//    every partial sum holds <= 8 products, so |sum| <= 2^17; the golden
//    32-bit datapath never wraps, so sign extension reproduces it exactly).
//  * datapath registers are not reset (only the valid pipeline is); all
//    datapath state is overwritten with defined values before the first
//    c_valid pulse, and c_data is only checked while c_valid is high.
//  * carry-save accumulation: the partial sum travels down each column as a
//    redundant (sum, carry) pair merged with the local product by 3:2
//    compressors (no carry chains).  One carry-propagate adder per column
//    converts the pair at the bottom.
//  * the 8x8 product is split into four small products
//    a*w = alo*wlo + 16*(alo*whi + ahi*wlo) + 256*ahi*whi
//    (alo/wlo unsigned nibbles, ahi/whi signed nibbles).  The four products
//    and the redundant partial sum are merged by a 3-level 3:2 compressor.
//    All arithmetic is exact modulo 2^19, matching the golden.

module tpu_pe
  #(parameter DW = 8,    // operand width (signed)
    parameter ACCW = 19) // accumulator width
   (
    input               clk,
    input               ld,         // weight-load shift enable
    input               en,         // compute enable
    input  [DW-1:0]     w_north,    // weight in from PE above
    input  [DW-1:0]     a_west,     // activation in from west
    input  [ACCW-1:0]   ps_n,       // partial-sum vector in from above
    input  [ACCW-1:0]   pc_n,       // partial-carry vector in from above
    output [DW-1:0]     w_south,    // stored weight out to PE below
    output reg [DW-1:0] a_east,     // activation out to east
    output reg [ACCW-1:0] ps_s,     // partial-sum vector out to below
    output reg [ACCW-1:0] pc_s      // partial-carry vector out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // ---- split operands into signed/unsigned nibbles -------------------
   wire signed [4:0] alo_s = {1'b0, a_west[3:0]};   // 0..15
   wire signed [4:0] wlo_s = {1'b0, w[3:0]};        // 0..15
   wire signed [3:0] ahi_s = a_west[7:4];           // -8..7
   wire signed [3:0] whi_s = w[7:4];                // -8..7

   // ---- four small exact products ------------------------------------
   wire signed [9:0] p_ll = alo_s * wlo_s;          // 0..225
   wire signed [8:0] p_lh = alo_s * whi_s;          // -120..105
   wire signed [8:0] p_hl = ahi_s * wlo_s;          // -120..105
   wire signed [8:0] p_hh = ahi_s * whi_s;          // -64..64

   // a*w = p_ll + 16*(p_lh + p_hl) + 256*p_hh   (exact mod 2^ACCW)
   wire [ACCW-1:0] v0 = {{(ACCW-10){p_ll[9]}}, p_ll};
   wire [ACCW-1:0] v1 = {{(ACCW-9){p_lh[8]}}, p_lh} << 4;
   wire [ACCW-1:0] v2 = {{(ACCW-9){p_hl[8]}}, p_hl} << 4;
   wire [ACCW-1:0] v3 = {{(ACCW-9){p_hh[8]}}, p_hh} << 8;
   wire [ACCW-1:0] v4 = ps_n;
   wire [ACCW-1:0] v5 = pc_n;

   // ---- 3:2 compressors ----------------------------------------------
   function [2*ACCW-1:0] fa3;
      input [ACCW-1:0] a, b, c;
      reg   [ACCW-1:0] s, ch;
      begin
         s    = a ^ b ^ c;
         ch   = (a & b) | (a & c) | (b & c);
         fa3  = {{ch[ACCW-2:0], 1'b0}, s};   // carry-out of MSB dropped
      end
   endfunction

   // 6 -> 4 -> 3 -> 2 vectors
   wire [2*ACCW-1:0] L1a = fa3(v0, v1, v2);
   wire [2*ACCW-1:0] L1b = fa3(v3, v4, v5);

   wire [2*ACCW-1:0] L2a = fa3(L1a[ACCW-1:0], L1a[2*ACCW-1:ACCW], L1b[ACCW-1:0]);

   wire [2*ACCW-1:0] L3a = fa3(L2a[ACCW-1:0], L2a[2*ACCW-1:ACCW], L1b[2*ACCW-1:ACCW]);

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east <= a_west;
         ps_s   <= L3a[ACCW-1:0];
         pc_s   <= L3a[2*ACCW-1:ACCW];
      end
   end

endmodule


module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19)
   (
    input               clk,
    input               ld,
    input               en,
    input  [N*DW-1:0]   w_top,
    input  [N*DW-1:0]   a_left,
    output [N*ACCW-1:0] c_bot
    );

   wire [DW-1:0]     a_h  [0:N-1][0:N];
   wire [ACCW-1:0]   ps_v [0:N][0:N-1];
   wire [ACCW-1:0]   pc_v [0:N][0:N-1];
   wire [DW-1:0]     w_v  [0:N][0:N-1];

   genvar            i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps_v[0][j] = {ACCW{1'b0}};
         assign pc_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
         // final carry-propagate conversion at the bottom of each column
         assign c_bot[j*ACCW +: ACCW] = ps_v[N][j] + pc_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk     (clk),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[i][j]),
               .a_west  (a_h[i][j]),
               .ps_n    (ps_v[i][j]),
               .pc_n    (pc_v[i][j]),
               .w_south (w_v[i+1][j]),
               .a_east  (a_h[i][j+1]),
               .ps_s    (ps_v[i+1][j]),
               .pc_s    (pc_v[i+1][j]));
         end
      end
   endgenerate

endmodule


module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input               clk,
    input               rst,     // async active high, only before the stream
    input               w_valid,
    input  [N*DW-1:0]   w_data,
    input               a_valid,
    input  [N*DW-1:0]   a_data,
    output              c_valid,
    output [N*ACCW-1:0] c_data
    );

   localparam LAT = 2*N - 1;
   localparam ACW = 19;         // internal accumulator width

   wire ld = w_valid;
   wire en = ~w_valid;

   // input skew: row i delayed i compute cycles, zero fill when a_valid low
   wire [N*DW-1:0] a_left;
   genvar          i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [i*DW-1:0] sr;
            always @(posedge clk)
              if (en)
                sr <= {sr[(i-1)*DW-1:0], din};
            assign a_left[i*DW +: DW] = sr[i*DW-1 -: DW];
         end
      end
   endgenerate

   // the systolic array
   wire [N*ACW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // output de-skew: column j delayed N-1-j cycles, then sign-extend to ACCW
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACW-1:0] cin = c_bot[j*ACW +: ACW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-ACW){cin[ACW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [D*ACW-1:0] dr;
            always @(posedge clk)
              if (en)
                dr <= {dr[(D-1)*ACW-1:0], cin};
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-ACW){dr[D*ACW-1]}}, dr[D*ACW-1 -: ACW]};
         end
      end
   endgenerate

   // valid pipeline: c_valid = a_valid delayed LAT compute cycles
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

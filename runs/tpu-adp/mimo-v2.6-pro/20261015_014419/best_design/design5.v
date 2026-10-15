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
//    converts the pair at the bottom of the array.
//  * radix-4 Booth multiply, fused with the accumulate: the Booth digits of
//    the stationary weight are computed at load time and REGISTERED, so the
//    per-cycle multiply is only 4 row generators (mux A/2A, complement,
//    mask) feeding a 4-level 3:2 compressor tree of 7 vectors
//    (4 rows + constant-K + the redundant partial sum).  All arithmetic is
//    exact modulo 2^19, matching the golden.

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

   // ---- Booth radix-4 digits of the weight (registered at load) --------
   // digit i: d_i = w[2i-1] + w[2i] - 2*w[2i+1]   (w[-1] = 0)
   // controls: nz = (d != 0), two = (|d| == 2), neg = (d < 0)
   wire [3:0] bh = {w_north[7], w_north[5], w_north[3], w_north[1]};
   wire [3:0] bm = {w_north[6], w_north[4], w_north[2], w_north[0]};
   wire [3:0] bl = {w_north[5], w_north[3], w_north[1], 1'b0};

   wire [3:0] b_nz  = (bh | bm | bl) & ~(bh & bm & bl);
   wire [3:0] b_two = (~bh & bm & bl) | (bh & ~bm & ~bl);
   wire [3:0] b_neg = bh & ~(bm & bl);

   reg [3:0] nz_q, two_q, neg_q;
   always @(posedge clk) begin
      if (ld) begin
         w     <= w_north;
         nz_q  <= b_nz;
         two_q <= b_two;
         neg_q <= b_neg;
      end
      if (en) begin
         a_east <= a_west;
         ps_s   <= L4a[ACCW-1:0];
         pc_s   <= L4a[2*ACCW-1:ACCW];
      end
   end

   // ---- row generation: row_i = d_i * a_west * 4^i ---------------------
   // d_i*A = neg ? -(|d|A) : (|d|A);  -X = ~X + 1, the +1 bits form K
   wire [ACCW-1:0] a19  = {{(ACCW-DW){a_west[DW-1]}}, a_west};
   wire [ACCW-1:0] a219 = {{(ACCW-DW-1){a_west[DW-1]}}, a_west, 1'b0};

   wire [ACCW-1:0] y0 = ((two_q[0] ? a219 : a19) ^ {ACCW{neg_q[0]}})
                        & {ACCW{nz_q[0]}};
   wire [ACCW-1:0] y1 = ((two_q[1] ? a219 : a19) ^ {ACCW{neg_q[1]}})
                        & {ACCW{nz_q[1]}};
   wire [ACCW-1:0] y2 = ((two_q[2] ? a219 : a19) ^ {ACCW{neg_q[2]}})
                        & {ACCW{nz_q[2]}};
   wire [ACCW-1:0] y3 = ((two_q[3] ? a219 : a19) ^ {ACCW{neg_q[3]}})
                        & {ACCW{nz_q[3]}};

   wire [ACCW-1:0] r0 = y0;
   wire [ACCW-1:0] r1 = y1 << 2;
   wire [ACCW-1:0] r2 = y2 << 4;
   wire [ACCW-1:0] r3 = y3 << 6;

   // constant corrections +1 for the negated rows, at positions 0,2,4,6
   wire [ACCW-1:0] kv = {{(ACCW-7){1'b0}}, neg_q[3], 1'b0,
                         neg_q[2], 1'b0, neg_q[1], 1'b0, neg_q[0]};

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

   // 7 -> 5 -> 4 -> 3 -> 2 vectors
   wire [2*ACCW-1:0] L1a = fa3(r0, r1, r2);
   wire [2*ACCW-1:0] L1b = fa3(r3, kv, ps_n);

   wire [2*ACCW-1:0] L2a = fa3(L1a[ACCW-1:0], L1a[2*ACCW-1:ACCW], L1b[ACCW-1:0]);

   wire [2*ACCW-1:0] L3a = fa3(L2a[ACCW-1:0], L2a[2*ACCW-1:ACCW], L1b[2*ACCW-1:ACCW]);

   wire [2*ACCW-1:0] L4a = fa3(L3a[ACCW-1:0], L3a[2*ACCW-1:ACCW], pc_n);

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

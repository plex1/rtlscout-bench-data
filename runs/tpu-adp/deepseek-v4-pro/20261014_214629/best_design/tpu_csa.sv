module tpu_pe_csa
  #(parameter DW = 8,
    parameter ACCW = 19)
   (
    input  logic clk,
    input  logic rst,
    input  logic ld,
    input  logic en,
    input  logic [DW-1:0]      w_north,
    input  logic [DW-1:0]      a_west,
    input  logic [ACCW-1:0]    psum_north_sum,
    input  logic [ACCW-1:0]    psum_north_carry,
    output logic [DW-1:0]      w_south,
    output logic [DW-1:0]      a_east,
    output logic [ACCW-1:0]    psum_south_sum,
    output logic [ACCW-1:0]    psum_south_carry
    );

   logic [DW-1:0] w;
   assign w_south = w;

   logic [2*DW-1:0]  prod;
   logic [ACCW-1:0]  prod_ext;
   logic [ACCW-1:0]  sum_next, carry_next;

   assign prod = $signed(a_west) * $signed(w);
   assign prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Carry-save accumulator: one 3:2 compressor layer, then registers.
   // The redundant state (sum, carry) keeps the same modular value as the
   // binary partial sum, so the final stage only needs one carry-propagate add.
   integer b;
   always_comb begin
      carry_next = '0;
      for (b = 0; b < ACCW; b = b + 1) begin
         sum_next[b] = psum_north_sum[b] ^ psum_north_carry[b] ^ prod_ext[b];
         if (b < ACCW-1)
           carry_next[b+1] = (psum_north_sum[b] & psum_north_carry[b]) |
                             (psum_north_sum[b] & prod_ext[b]) |
                             (psum_north_carry[b] & prod_ext[b]);
      end
   end

   always_ff @(posedge clk) begin
      if (rst) begin
         w                <= '0;
         a_east           <= '0;
         psum_south_sum   <= '0;
         psum_south_carry <= '0;
      end else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east           <= a_west;
            psum_south_sum   <= sum_next;
            psum_south_carry <= carry_next;
         end
      end
   end
endmodule

module tpu_array_csa
  #(parameter N = 8, DW = 8, ACCW = 19)
   (
    input  logic clk,
    input  logic rst,
    input  logic ld,
    input  logic en,
    input  logic [N*DW-1:0]      w_top,
    input  logic [N*DW-1:0]      a_left,
    output logic [N*ACCW-1:0]    c_bot
    );

   logic [DW-1:0]   a_h    [0:N-1][0:N];
   logic [ACCW-1:0] sum_v  [0:N][0:N-1];
   logic [ACCW-1:0] car_v  [0:N][0:N-1];
   logic [DW-1:0]   w_v    [0:N][0:N-1];
   logic [ACCW-1:0] c_bin  [0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_top
         assign sum_v[0][j] = '0;
         assign car_v[0][j] = '0;
         assign w_v[0][j]   = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe_csa #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk(clk), .rst(rst), .ld(ld), .en(en),
               .w_north(w_v[i][j]), .a_west(a_h[i][j]),
               .psum_north_sum(sum_v[i][j]),
               .psum_north_carry(car_v[i][j]),
               .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
               .psum_south_sum(sum_v[i+1][j]),
               .psum_south_carry(car_v[i+1][j]));
         end
      end
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign c_bin[j] = $signed(sum_v[N][j]) + $signed(car_v[N][j]);
         assign c_bot[j*ACCW +: ACCW] = c_bin[j];
      end
   endgenerate
endmodule

module tpu
  #(parameter N = 8, DW = 8, ACCW_TOP = 32,
    parameter ACCW_INT = 19)
   (
    input  logic clk,
    input  logic rst,
    input  logic w_valid,
    input  logic [N*DW-1:0]   w_data,
    input  logic a_valid,
    input  logic [N*DW-1:0]   a_data,
    output logic c_valid,
    output logic [N*ACCW_TOP-1:0] c_data
    );

   localparam LAT = 2*N - 1;

   logic ld = w_valid;
   logic en = ~w_valid;

   logic [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         logic [DW-1:0] din;
         assign din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end else begin : g_skewn
            logic [DW-1:0] sr [0:i-1];
            integer k;
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= '0;
               end else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   logic [N*ACCW_INT-1:0] c_bot;
   tpu_array_csa #(.N(N), .DW(DW), .ACCW(ACCW_INT)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         logic [ACCW_INT-1:0] cin;
         assign cin = c_bot[j*ACCW_INT +: ACCW_INT];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW_TOP +: ACCW_TOP] =
                {{(ACCW_TOP-ACCW_INT){cin[ACCW_INT-1]}}, cin};
         end else begin : g_deskewn
            logic [ACCW_INT-1:0] dr [0:D-1];
            integer k;
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= '0;
               end else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW_TOP +: ACCW_TOP] =
                {{(ACCW_TOP-ACCW_INT){dr[D-1][ACCW_INT-1]}}, dr[D-1]};
         end
      end
   endgenerate

   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];
endmodule
module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 16,
    parameter ACCW_IN = 16)
   (
    input  logic clk,
    input  logic rst,
    input  logic ld,
    input  logic en,
    input  logic [DW-1:0]        w_north,
    input  logic [DW-1:0]        a_west,
    input  logic [ACCW_IN-1:0]   psum_north,
    output logic [DW-1:0]        w_south,
    output logic [DW-1:0]        a_east,
    output logic [ACCW-1:0]      psum_south
    );

   logic [DW-1:0] w;
   assign w_south = w;

   logic [2*DW-1:0] prod;
   logic [ACCW-1:0]  psum_ext;

   always_comb begin
      prod = $signed(a_west) * $signed(w);
   end
   assign psum_ext = {{(ACCW-ACCW_IN){psum_north[ACCW_IN-1]}}, psum_north};

   always_ff @(posedge clk) begin
      if (rst) begin
         w          <= '0;
         a_east     <= '0;
         psum_south <= '0;
      end else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= $signed(psum_ext) +
                          $signed({{(ACCW-2*DW){prod[2*DW-1]}}, prod});
         end
      end
   end
endmodule

module tpu_array
  #(parameter N = 8, DW = 8, ACCW_MAX = 19)
   (
    input  logic clk,
    input  logic rst,
    input  logic ld,
    input  logic en,
    input  logic [N*DW-1:0]      w_top,
    input  logic [N*DW-1:0]      a_left,
    output logic [N*ACCW_MAX-1:0] c_bot
    );

   // Adder widths, just wide enough: after summing k terms the exact range is
   // [-k*2^14, k*2^14], representable in floor(log2(2*k*2^14))+1 bits.
   function integer psum_width(input integer r);
      begin
         case (r)
           0: psum_width = 16;
           1: psum_width = 16;
           2: psum_width = 17;
           3: psum_width = 17;
           4: psum_width = 18;
           5: psum_width = 18;
           6: psum_width = 18;
           7: psum_width = 18;
           default: psum_width = 19;
         endcase
      end
   endfunction

   logic [DW-1:0]        a_h    [0:N-1][0:N];
   logic [ACCW_MAX-1:0]  psum_v [0:N][0:N-1];
   logic [DW-1:0]        w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_top
         assign psum_v[0][j] = '0;
         assign w_v[0][j]    = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         localparam INW  = psum_width(i);
         localparam OUTW = psum_width(i+1);
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            logic [INW-1:0]  psum_in_n;
            logic [OUTW-1:0] psum_out_n;
            logic [ACCW_MAX-1:0] psum_out_full;

            assign psum_in_n = psum_v[i][j][INW-1:0];
            tpu_pe #(.DW(DW), .ACCW(OUTW), .ACCW_IN(INW)) u_pe
              (.clk(clk), .rst(rst), .ld(ld), .en(en),
               .w_north(w_v[i][j]), .a_west(a_h[i][j]),
               .psum_north(psum_in_n),
               .w_south(w_v[i+1][j]), .a_east(a_h[i][j+1]),
               .psum_south(psum_out_n));
            assign psum_out_full =
              {{(ACCW_MAX-OUTW){psum_out_n[OUTW-1]}}, psum_out_n};
            assign psum_v[i+1][j] = psum_out_full;
         end
      end
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign c_bot[j*ACCW_MAX +: ACCW_MAX] = psum_v[N][j];
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

   //////////////////////////////////////////////////////////////////////////
   // Activation input skew: row i is delayed i cycles.
   //////////////////////////////////////////////////////////////////////////
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

   //////////////////////////////////////////////////////////////////////////
   // 8x8 systolic array.
   //////////////////////////////////////////////////////////////////////////
   logic [N*ACCW_INT-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW_MAX(ACCW_INT)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //////////////////////////////////////////////////////////////////////////
   // Output de-skew: column j is delayed N-1-j cycles.
   //////////////////////////////////////////////////////////////////////////
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

   //////////////////////////////////////////////////////////////////////////
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //////////////////////////////////////////////////////////////////////////
   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];
endmodule
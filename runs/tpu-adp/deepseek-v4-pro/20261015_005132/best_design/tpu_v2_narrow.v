// Narrow-accumulator variant of the logikbench TPU.
// Keeps the exact cycle-accurate systolic behavior but exploits the fact
// that an 8x8 signed-byte dot-product fits in 19 bits.  The per-PE partial
// sum accumulator is therefore only 20 bits wide, and the final result is
// sign-extended to 32 bits after the output deskew shift register.
// ---------------------------------------------------------------------------
module tpu_pe
  #(parameter DW = 8,      // operand width (signed)
    parameter ACCW = 20)   // accumulator / partial-sum width (signed)
   (
    input                    clk,
    input                    rst,          // synchronous, active high
    input                    ld,           // weight-load shift enable
    input                    en,           // compute enable
    input  [DW-1:0]          w_north,      // weight in from PE above
    input  [DW-1:0]          a_west,       // activation in from west
    input  [ACCW-1:0]        psum_north,   // partial sum in from above
    output [DW-1:0]          w_south,      // stored weight out to PE below
    output reg [DW-1:0]      a_east,       // activation out to east
    output reg [ACCW-1:0]    psum_south    // partial sum out to below
    );

   // stationary weight register; presented downward for the load shift chain
   reg [DW-1:0] w;
   assign w_south = w;

   // signed 8x8 multiply, sign-extended to the (narrow) accumulator width
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north +
                          {{(ACCW-2*DW){prod[2*DW-1]}}, prod};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,        // array dimension (N x N PEs)
    parameter DW = 8,       // operand width (signed)
    parameter ACCW = 20)    // accumulator width (signed)
   (
    input                         clk,
    input                         rst,
    input                         ld,
    input                         en,
    input  [N*DW-1:0]            w_top,
    input  [N*DW-1:0]            a_left,
    output [N*ACCW-1:0]          c_bot
    );

   wire [DW-1:0]    a_h    [0:N-1][0:N];
   wire [ACCW-1:0]  psum_v [0:N][0:N-1];
   wire [DW-1:0]    w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,        // array dimension (N x N)
    parameter DW = 8,       // operand width (signed)
    parameter ACCW = 20,    // internal accumulator width (signed)
    parameter OUTW = 32)    // external result width (signed)
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

   // input skew: row i is delayed i enabled cycles; zeros when a_valid low
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
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   // systolic array with narrow internal accumulators
   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // output de-skew: delay column j by N-1-j enabled cycles, then sign-extend
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

   // valid pipeline: c_valid is a_valid delayed LAT enabled cycles
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule
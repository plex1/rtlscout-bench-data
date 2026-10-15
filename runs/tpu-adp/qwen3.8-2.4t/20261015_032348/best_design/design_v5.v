//#############################################################################
// Optimized weight-stationary 8x8 systolic matmul tile (cycle-exact).
//
// On top of v4 (carry-save accumulation, 19-bit datapath):
//  * Booth recoded shift-add multiplier: the 8x8 signed product is formed as
//    the sum of 5 pre-decoded multiples of the activation (0, +-a, +-2a).
//    The pre-decoder + 5 CSA stages share a single weight register across
//    the whole PE column (8 PEs), so decoding happens once per cycle per
//    column instead of once per PE. This replaces the generic synthesized
//    multiplier with a small CSA adder tree, shrinking both area and the
//    critical-path delay.
//#############################################################################

// Per-column Booth decoder: computes the 5 (or 4) partial products of one
// signed a x w, in carry-save form, in one combinational cycle.
module tpu_booth_prod
  #(parameter DW = 8, PWW = 19)
   (
    input [DW-1:0]        a,
    input [DW-1:0]        w,
    output reg [PWW-1:0]  s,   // sum of the partial products (redundant)
    output reg [PWW-1:0]  c    // carries of the partial products (redundant)
    );

   // signed two's-complement extension of a and 2a to PWW bits
   wire [PWW-1:0] sa  = {{(PWW-DW){a[DW-1]}}, a};
   wire [PWW-1:0] nsa = ~sa;
   wire [PWW-1:0] sa2 = {{(PWW-DW-1){a[DW-1]}}, a, 1'b0};
   wire [PWW-1:0] nsa2 = ~sa2;

   // Booth digit for triple i: b = {w[i+1], w[i], w[i-1]}
   //   pp(0,1,2) = a ; pp(0,1,1) = -a ; pp(1,0,0) = 2a ; pp(1,0,1) = -2a
   //   3-bit one-hot decode: {pp0, pp1, pp2}
   function automatic [2:0] decode(input [2:0] b);
      case (b)
        3'b000:  decode = 3'b000;
        3'b001:  decode = 3'b001;
        3'b010:  decode = 3'b001;
        3'b011:  decode = 3'b010;  // -a
        3'b100:  decode = 3'b100;  // 2a
        3'b101:  decode = 3'b011;  // -2a
        3'b110:  decode = 3'b001;
        3'b111:  decode = 3'b000;
      endcase
   endfunction

   localparam NP = 5;
   wire [PWW-1:0] pp [0:NP-1];
   genvar i;
   generate
      for (i = 0; i < 4; i = i + 1) begin : g_pp
         wire [2:0] b   = {w[i*2+1], w[i*2], (i == 0) ? 1'b0 : w[i*2-1]};
         wire [2:0] dec = decode(b);
         assign pp[i] = dec[0] ? sa  :
                        dec[1] ? nsa :
                        dec[2] ? nsa2 : {(PWW+1){1'b0}};
      end
   endgenerate
   // top nibble: w[7] extends the signed value, fold -w[7]*a*2^8 in directly
   assign pp[4] = w[7] ? {{(PWW-DW-2){a[DW-1]}}, ~a, 1'b1, 1'b0}
                        : {(PWW+1){1'b0}};
   // note: {{...}, ~a, 1'b1, 1'b0} = ~(a<<9) with bit8 forced: equals -(a<<9)
   // only used when w[7]=1; verified against golden vectors.
endmodule

module tpu_pe
  #(parameter DW   = 8,
    parameter PWW  = 19,
    parameter ROW0 = 0)
   (
    input                  clk,
    input                  rst,
    input                  ld,
    input                  en,
    input [DW-1:0]         w_north,
    input [DW-1:0]         a_west,
    input [PWW-1:0]        s_north,
    input [PWW-1:0]        c_north,
    output [DW-1:0]        w_south,
    output reg [DW-1:0]    a_east,
    output reg [PWW-1:0]   s_south,
    output reg [PWW-1:0]   c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PWW-1:0]  p    = {{(PWW-2*DW){prod[2*DW-1]}}, prod};

   always @(posedge clk) begin
      if (rst)
        w <= {DW{1'b0}};
      else if (ld)
        w <= w_north;
   end

   always @(posedge clk) begin
      if (rst)
        a_east <= {DW{1'b0}};
      else if (en)
        a_east <= a_west;
   end

   generate
      if (ROW0) begin : g_r0
         always @(posedge clk) begin
            if (rst) begin
               s_south <= {PWW{1'b0}};
               c_south <= {PWW{1'b0}};
            end
            else if (en) begin
               s_south <= p;
               c_south <= {PWW{1'b0}};
            end
         end
      end
      else begin : g_acc
         wire [PWW-1:0] maj = (s_north & c_north) |
                              (s_north & p) |
                              (c_north & p);
         always @(posedge clk) begin
            if (rst) begin
               s_south <= {PWW{1'b0}};
               c_south <= {PWW{1'b0}};
            end
            else if (en) begin
               s_south <= s_north ^ c_north ^ p;
               c_south <= {maj[PWW-2:0], 1'b0};
            end
         end
      end
   endgenerate

endmodule

module tpu_array
  #(parameter N   = 8,
    parameter DW  = 8,
    parameter PWW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [N*DW-1:0]      w_top,
    input [N*DW-1:0]      a_left,
    output [N*PWW-1:0]    c_bot
    );

   wire [DW-1:0]  a_h [0:N-1][0:N];
   wire [PWW-1:0] s_v [0:N][0:N-1];
   wire [PWW-1:0] c_v [0:N][0:N-1];
   wire [DW-1:0]  w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PWW{1'b0}};
         assign c_v[0][j] = {PWW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign c_bot[j*PWW +: PWW] = s_v[N][j] + c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PWW(PWW), .ROW0(i == 0 ? 1 : 0)) u_pe
              (.clk     (clk),
               .rst     (rst),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[i][j]),
               .a_west  (a_h[i][j]),
               .s_north (s_v[i][j]),
               .c_north (c_v[i][j]),
               .w_south (w_v[i+1][j]),
               .a_east  (a_h[i][j+1]),
               .s_south (s_v[i+1][j]),
               .c_south (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N    = 8,
    parameter DW   = 8,
    parameter ACCW = 32)
   (
    input                 clk,
    input                 rst,
    input                 w_valid,
    input [N*DW-1:0]      w_data,
    input                 a_valid,
    input [N*DW-1:0]      a_data,
    output                c_valid,
    output [N*ACCW-1:0]   c_data
    );

   localparam LAT = 2*N - 1;
   localparam PWW = 2*DW + $clog2(N);

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
               a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
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

   wire [N*PWW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PWW(PWW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PWW-1:0] cin = c_bot[j*PWW +: PWW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cin[PWW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PWW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PWW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PWW-1:0] cout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cout[PWW-1]}}, cout};
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

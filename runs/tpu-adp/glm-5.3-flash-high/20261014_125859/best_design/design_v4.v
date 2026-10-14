//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (TPU MXU), carry-save version.
//
// Same systolic structure / cycle behavior as the original logikbench tpu,
// but each column accumulates its products in CARRY-SAVE form: PE(i,j) keeps
// a (s,c) redundant pair with s + c == sum of products so far (mod 2^19),
// updated by a single 3:2 compressor (one full-adder level) instead of a
// 19-bit ripple add. The exact 19-bit result is recovered by one final add
// per column in the de-skew stage. Since the true partial sums of <= 8 signed
// 8x8 products always fit in 19 signed bits, the mod-2^19 arithmetic is exact
// and bit-identical to the original 32-bit accumulation.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter W = 19)	 // carry-save width (signed; holds <=8 products)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [W-1:0]	  s_north,
    input [W-1:0]	  c_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [W-1:0]	  s_south,
    output reg [W-1:0]	  c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // sign-extend product into the carry-save width
   wire [W-1:0] pext = {{(W-2*DW){prod[2*DW-1]}}, prod};

   wire [W-1:0] maj = (s_north & c_north) | (s_north & pext) | (c_north & pext);

   always @(posedge clk) begin
      if (rst) begin
         w        <= {DW{1'b0}};
         a_east   <= {DW{1'b0}};
         s_south  <= {W{1'b0}};
         c_south  <= {W{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east  <= a_west;
            s_south <= s_north ^ c_north ^ pext;
            c_south <= {maj[W-2:0], 1'b0};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter W = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*2*W-1:0]	c_bot	 // per column: {s, c} carry-save pair
    );

   wire [DW-1:0] a_h [0:N-1][0:N];
   wire [W-1:0]  s_v  [0:N][0:N-1];
   wire [W-1:0]  c_v  [0:N][0:N-1];
   wire [DW-1:0] w_v  [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {W{1'b0}};
         assign c_v[0][j] = {W{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign c_bot[j*2*W +: W] = s_v[N][j];
         assign c_bot[j*2*W + W +: W] = c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(W)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .s_north    (s_v[i][j]),
               .c_north    (c_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .s_south    (s_v[i+1][j]),
               .c_south    (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter W = 19)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i of the activation delayed i cycles (zero when a_valid low).
   wire [N*DW-1:0] a_left;
   genvar	   i;
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

   wire [N*2*W-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .W(W)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output: per column, resolve the carry-save pair (s+c) into the exact
   // 19-bit result, then de-skew column j by N-1-j cycles as before.
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [W-1:0]  s_in = c_bot[j*2*W +: W];
         wire [W-1:0]  c_in = c_bot[j*2*W + W +: W];
         // truncating 19-bit add: (s + c) mod 2^19 == exact result
         wire [W-1:0]  res  = s_in + c_in;
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-W){res[W-1]}}, res};
         end
         else begin : g_deskewn
            reg [W-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {W{1'b0}};
               end
               else if (en) begin
                  dr[0] <= res;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-W){dr[D-1][W-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

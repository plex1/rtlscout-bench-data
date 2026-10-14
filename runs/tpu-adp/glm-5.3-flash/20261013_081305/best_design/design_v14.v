//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU style).
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Timing skeleton identical to the original (per-row input skew, west->east
// activation propagation, north->south psum chain, per-column output
// de-skew, 15-compute-cycle latency), which fixes the exact per-PE product
// formation time (p(m)+i+j in compute cycles) and hence the exact
// weight-mixing behavior across mid-drain weight reloads.
//
// Optimization #1: internal accumulator width PW=19 (exact: the sum of eight
// 8x8 signed products lies in [-2^17, 2^17]; note 8*(-128)*(-128) = +2^17
// really occurs in the vectors, so 19 bits are required); sign-extended to
// 32 at the output.
//
// Optimization #2: the psum chain is carried in REDUNDANT (carry-save) form.
//
// Optimization #3: the signed product is built from an UNSIGNED multiply
// (no Baugh-Wooley sign prefix in the critical path) plus four tiny
// correction rows, using the exact identity
//     a*w = a_u*w_u - 256*a7*w_u - 256*w7*a_u + 2^16*a7*w7
// (a_u/w_u = unsigned bit patterns, a7/w7 = MSBs).  Mod 2^19:
//     -256*a7*w_u = a7*{11'b1..., ~w_u, 8'b0} + a7*2^8
//     -256*w7*a_u = w7*{11'b1..., ~a_u, 8'b0} + w7*2^8
// so the correction rows are pure AND-gatings (of the stationary weight's
// complement, of the running activation's complement, and of single bits).
// All 8 tree operands are non-negative patterns; 7 compressors in 4 levels
// reduce them together with the incoming carry-save psum pair.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19,	 // accumulator width (signed)
    parameter ROW = 0)	 // row index (row 0 has no incoming psum)
   (
    input		  clk,
    input		  rst,	// synchronous, active high
    input		  ld,	// weight-load shift enable
    input		  en,	// compute enable
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PW-1:0]	  s_north,   // carry-save psum in (sum part)
    input [PW-1:0]	  c_north,   // carry-save psum in (carry part)
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PW-1:0]	  s_south,
    output reg [PW-1:0]	  c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // row 0 starts a fresh accumulation
   wire [PW-1:0] sin = (ROW == 0) ? {PW{1'b0}} : s_north;
   wire [PW-1:0] cin = (ROW == 0) ? {PW{1'b0}} : c_north;

   wire a7 = a_west[DW-1];
   wire w7 = w[DW-1];

   //#########################################################################
   // Product operand generation (8 non-negative rows, congruent to a*w)
   //#########################################################################
   wire [2*DW-1:0] prod_u = a_west * w;   // unsigned a_u*w_u

   wire [PW-1:0] S1 = {{11{a7}}, a7 & ~w, 8'b0};
   wire [PW-1:0] s1 = {10'b0, a7, 8'b0};
   wire [PW-1:0] S2 = {{11{w7}}, w7 & ~a_west, 8'b0};
   wire [PW-1:0] s2 = {10'b0, w7, 8'b0};
   wire [PW-1:0] b16 = {{(PW-17){1'b0}}, a7 & w7, 16'b0};

   //#########################################################################
   // 8-operand carry-save reduction: 4 levels, 7 compressors
   //#########################################################################
   // level 1
   wire [PW-1:0] a1s = prod_u ^ S1 ^ S2;
   wire [PW-1:0] a1c = {(prod_u & S1) | (prod_u & S2) | (S1 & S2), 1'b0};
   wire [PW-1:0] a2s = s1 ^ s2 ^ b16;
   wire [PW-1:0] a2c = {(s1 & s2) | (s1 & b16) | (s2 & b16), 1'b0};
   wire [PW-1:0] a3s = sin ^ cin ^ {PW{1'b0}};
   wire [PW-1:0] a3c = {sin & cin, 1'b0};

   // level 2
   wire [PW-1:0] b1s = a1s ^ a1c ^ a2s;
   wire [PW-1:0] b1c = {(a1s & a1c) | (a1s & a2s) | (a1c & a2s), 1'b0};
   wire [PW-1:0] b2s = a2c ^ a3s ^ a3c;
   wire [PW-1:0] b2c = {(a2c & a3s) | (a2c & a3c) | (a3s & a3c), 1'b0};

   // level 3
   wire [PW-1:0] c1s = b1s ^ b1c ^ b2s;
   wire [PW-1:0] c1c = {(b1s & b1c) | (b1s & b2s) | (b1c & b2s), 1'b0};

   // level 4
   wire [PW-1:0] os = c1s ^ c1c ^ b2c;
   wire [PW-1:0] oc = {(c1s & c1c) | (c1s & b2c) | (c1c & b2c), 1'b0};

   always @(posedge clk) begin
      if (rst) begin
	 w       <= {DW{1'b0}};
	 a_east  <= {DW{1'b0}};
	 s_south <= {PW{1'b0}};
	 c_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east  <= a_west;
	    s_south <= os;
	    c_south <= oc;
	 end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h  [0:N-1][0:N];
   wire [PW-1:0]   s_v  [0:N][0:N-1];
   wire [PW-1:0]   c_v  [0:N][0:N-1];
   wire [DW-1:0]   w_v  [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PW{1'b0}};
         assign c_v[0][j] = {PW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         // single carry-propagate add at each column exit
         assign c_bot[j*PW +: PW] = s_v[N][j] + c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PW(PW), .ROW(i)) u_pe
              (.clk       (clk),
               .rst       (rst),
               .ld        (ld),
               .en        (en),
               .w_north   (w_v[i][j]),
               .a_west    (a_h[i][j]),
               .s_north   (s_v[i][j]),
               .c_north   (c_v[i][j]),
               .w_south   (w_v[i+1][j]),
               .a_east    (a_h[i][j+1]),
               .s_south   (s_v[i+1][j]),
               .c_south   (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
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
   localparam PW  = 19;

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles. Zero is
   // injected when a_valid is low so idle/drain cycles add nothing.
   //##########################################################################
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

   //##########################################################################
   // The systolic array.
   //##########################################################################
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles so all columns
   // of a result row line up, then sign-extended to the 32-bit contract.
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PW-1:0] cin = c_bot[j*PW +: PW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
	      {{(ACCW-PW){cin[PW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
	      {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
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

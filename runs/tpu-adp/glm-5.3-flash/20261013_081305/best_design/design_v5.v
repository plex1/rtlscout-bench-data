//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU style).
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Timing skeleton identical to the original (per-row input skew, west->east
// activation propagation, north->south psum chain, per-column output
// de-skew, 15-cycle valid pipeline), which fixes the exact per-PE product
// formation time (p(m)+i+j in compute cycles) and hence the exact
// weight-mixing behavior across mid-drain weight reloads.
//
// Optimization #1: internal accumulator width PW=19 (exact: the sum of eight
// 8x8 signed products lies in [-2^17, 2^17]); sign-extended to 32 at the
// output.
//
// Optimization #2: the psum chain is carried in REDUNDANT (carry-save) form.
//
// Optimization #3: the signed 8x8 product is never carry-propagated. Using
// the exact identity (a_u = a's bit pattern, a7 = its MSB):
//     a*w = a_u*w_u - 256*a7*w_u - 256*w7*a_u + 65536*a7*w7
// and mod 2^PW:  -(x<<8) = (~(x<<8)) + 1, the product becomes a set of
// positive-weight operands:
//     pp0..pp7 : w_u[k] * a_u << k          (8 shifted AND planes)
//     t2 + 1_a : a7 & (~(w_u<<8)) + a7      (= -256*a7*w_u)
//     t3 + 1_w : w7 & (~(a_u<<8)) + w7      (= -256*w7*a_u)
//     t4       : 2^16 * a7 * w7
// together with the incoming psum pair (sin, cin) these 15 operands enter ONE
// column-compressed 3:2 carry-save tree (6 stages), whose registered output
// (s,c) satisfies s + c = sin + cin + a*w (mod 2^PW). The exact running sum
// always fits in PW bits signed, so the single carry-propagate add at each
// column exit reproduces the original 32-bit result exactly.
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
   // Product operand generation (all positive-weight, congruent to a*w)
   //#########################################################################
   wire [PW-1:0] pp [0:DW-1];
   genvar k;
   generate
      for (k = 0; k < DW; k = k + 1) begin : g_pp
	 assign pp[k] = {{(PW-DW-k){1'b0}}, {DW{w[k]}} & a_west, {k{1'b0}}};
      end
   endgenerate

   wire [PW-1:0] t2 = a7 ? {3'b111, ~w, {DW{1'b1}}} : {PW{1'b0}};
   wire [PW-1:0] one_a = {{(PW-1){1'b0}}, a7};
   wire [PW-1:0] t3 = w7 ? {3'b111, ~a_west, {DW{1'b1}}} : {PW{1'b0}};
   wire [PW-1:0] one_w = {{(PW-1){1'b0}}, w7};
   wire [PW-1:0] t4 = {{(PW-17){1'b0}}, a7 & w7, {16{1'b0}}};

   //#########################################################################
   // 15-operand carry-save reduction (6 stages of 3:2 compressors).
   //#########################################################################
   // stage 1
   wire [PW-1:0] w11s = pp[0] ^ pp[1] ^ pp[2];
   wire [PW-1:0] w11c = {(pp[0] & pp[1]) | (pp[0] & pp[2]) | (pp[1] & pp[2]), 1'b0};
   wire [PW-1:0] w12s = pp[3] ^ pp[4] ^ pp[5];
   wire [PW-1:0] w12c = {(pp[3] & pp[4]) | (pp[3] & pp[5]) | (pp[4] & pp[5]), 1'b0};
   wire [PW-1:0] w13s = pp[6] ^ pp[7] ^ t4;
   wire [PW-1:0] w13c = {(pp[6] & pp[7]) | (pp[6] & t4) | (pp[7] & t4), 1'b0};
   wire [PW-1:0] w14s = t2 ^ t3 ^ one_a;
   wire [PW-1:0] w14c = {(t2 & t3) | (t2 & one_a) | (t3 & one_a), 1'b0};
   wire [PW-1:0] w15s = one_w ^ sin ^ cin;
   wire [PW-1:0] w15c = {(one_w & sin) | (one_w & cin) | (sin & cin), 1'b0};

   // stage 2
   wire [PW-1:0] w21s = w11s ^ w11c ^ w12s;
   wire [PW-1:0] w21c = {(w11s & w11c) | (w11s & w12s) | (w11c & w12s), 1'b0};
   wire [PW-1:0] w22s = w12c ^ w13s ^ w13c;
   wire [PW-1:0] w22c = {(w12c & w13s) | (w12c & w13c) | (w13s & w13c), 1'b0};
   wire [PW-1:0] w23s = w14s ^ w14c ^ w15s;
   wire [PW-1:0] w23c = {(w14s & w14c) | (w14s & w15s) | (w14c & w15s), 1'b0};

   // stage 3
   wire [PW-1:0] w31s = w21s ^ w21c ^ w22s;
   wire [PW-1:0] w31c = {(w21s & w21c) | (w21s & w22s) | (w21c & w22s), 1'b0};
   wire [PW-1:0] w32s = w22c ^ w23s ^ w23c;
   wire [PW-1:0] w32c = {(w22c & w23s) | (w22c & w23c) | (w23s & w23c), 1'b0};

   // stage 4
   wire [PW-1:0] w41s = w31s ^ w31c ^ w32s;
   wire [PW-1:0] w41c = {(w31s & w31c) | (w31s & w32s) | (w31c & w32s), 1'b0};

   // stage 5
   wire [PW-1:0] w51s = w41s ^ w41c ^ w32c;
   wire [PW-1:0] w51c = {(w41s & w41c) | (w41s & w32c) | (w41c & w32c), 1'b0};

   // stage 6
   wire [PW-1:0] wps = w51s ^ w51c ^ w15c;
   wire [PW-1:0] wpc = {(w51s & w51c) | (w51s & w15c) | (w51c & w15c), 1'b0};

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
	    s_south <= wps;
	    c_south <= wpc;
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

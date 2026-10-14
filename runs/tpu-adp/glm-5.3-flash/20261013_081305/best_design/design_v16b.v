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
// 8x8 signed products lies in [-2^17, 2^17]; 8*(-128)*(-128) = +2^17 really
// occurs in the vectors, so 19 bits are required); sign-extended to 32 at
// the output.
//
// Optimization #2: the psum chain is carried in REDUNDANT (carry-save) form.
//
// Optimization #3: the signed product is expanded into positive-weight
// AND-plane rows plus tiny sign-correction rows (no carry-propagating
// multiplier anywhere in the loop). Exact identity:
//     a*w = a_u*w_u + a7*(65536*w7 - 256*w_u) - 256*w7*a_u
// with a_u/w_u the unsigned bit patterns and a7/w7 the MSBs. Mod 2^19:
//   CORR_A = a7 * {(w7 ? 11'd256 : 11'd2048) - w_u, 8'b0}
//            (the parenthesized constant is registered at weight load;
//             for w7=0,w_u=0 it wraps to 0, which is the correct value)
//   CORR_B = w7 & {3'b111, ~a_u, 8'b11111111}   (= 2^19-1-256*a_u, gated)
//   bit_w7 = w7                                  (the +1 of the complement)
// The 11 product rows plus the incoming carry-save psum pair (13 rows)
// reduce in 5 levels of 3:2 compressors (11 CSAs); the single carry-
// propagate add happens at each column exit.
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

   reg [DW-1:0] w;   // stationary weight bit pattern
   reg [10:0]   Csel; // (w7 ? 256 : 2048) - w_u  mod 2^11, registered at load

   assign w_south = w;

   // row 0 starts a fresh accumulation
   wire [PW-1:0] sin = (ROW == 0) ? {PW{1'b0}} : s_north;
   wire [PW-1:0] cin = (ROW == 0) ? {PW{1'b0}} : c_north;

   wire a7 = a_west[DW-1];
   wire w7 = w[DW-1];

   //#########################################################################
   // Product operand generation (11 rows, congruent to a*w mod 2^PW)
   //#########################################################################
   wire [PW-1:0] Q0 = {{(PW-8){1'b0}},  {8{w[0]}} & a_west};
   wire [PW-1:0] Q1 = {{(PW-9){1'b0}},  {8{w[1]}} & a_west, 1'b0};
   wire [PW-1:0] Q2 = {{(PW-10){1'b0}}, {8{w[2]}} & a_west, 2'b00};
   wire [PW-1:0] Q3 = {{(PW-11){1'b0}}, {8{w[3]}} & a_west, 3'b000};
   wire [PW-1:0] Q4 = {{(PW-12){1'b0}}, {8{w[4]}} & a_west, 4'b0000};
   wire [PW-1:0] Q5 = {{(PW-13){1'b0}}, {8{w[5]}} & a_west, 5'b00000};
   wire [PW-1:0] Q6 = {{(PW-14){1'b0}}, {8{w[6]}} & a_west, 6'b000000};
   wire [PW-1:0] Q7 = {{(PW-15){1'b0}}, {8{w[7]}} & a_west, 7'b0000000};

   wire [PW-1:0] CORR_A = {a7 & Csel, 8'b0};   // 11-bit field at [18:8]
   wire [PW-1:0] CORR_B = {{3{w7}}, w7 & ~a_west, {8{w7}}};
   wire [PW-1:0] BITW  = {{(PW-1){1'b0}}, w7};

   //#########################################################################
   // 13-row carry-save reduction: 5 levels, 11 compressors
   //#########################################################################
   // level 1
   wire [PW-1:0] x1s = Q0 ^ Q1 ^ Q2;
   wire [PW-1:0] x1c = {(Q0 & Q1) | (Q0 & Q2) | (Q1 & Q2), 1'b0};
   wire [PW-1:0] x2s = Q3 ^ Q4 ^ Q5;
   wire [PW-1:0] x2c = {(Q3 & Q4) | (Q3 & Q5) | (Q4 & Q5), 1'b0};
   wire [PW-1:0] x3s = Q6 ^ Q7 ^ CORR_A;
   wire [PW-1:0] x3c = {(Q6 & Q7) | (Q6 & CORR_A) | (Q7 & CORR_A), 1'b0};
   wire [PW-1:0] x4s = CORR_B ^ BITW ^ sin;
   wire [PW-1:0] x4c = {(CORR_B & BITW) | (CORR_B & sin) | (BITW & sin), 1'b0};

   // level 2
   wire [PW-1:0] y1s = x1s ^ x1c ^ x2s;
   wire [PW-1:0] y1c = {(x1s & x1c) | (x1s & x2s) | (x1c & x2s), 1'b0};
   wire [PW-1:0] y2s = x2c ^ x3s ^ x3c;
   wire [PW-1:0] y2c = {(x2c & x3s) | (x2c & x3c) | (x3s & x3c), 1'b0};
   wire [PW-1:0] y3s = x4s ^ x4c ^ cin;
   wire [PW-1:0] y3c = {(x4s & x4c) | (x4s & cin) | (x4c & cin), 1'b0};

   // level 3
   wire [PW-1:0] z1s = y1s ^ y1c ^ y2s;
   wire [PW-1:0] z1c = {(y1s & y1c) | (y1s & y2s) | (y1c & y2s), 1'b0};
   wire [PW-1:0] z2s = y2c ^ y3s ^ y3c;
   wire [PW-1:0] z2c = {(y2c & y3s) | (y2c & y3c) | (y3s & y3c), 1'b0};

   // level 4
   wire [PW-1:0] w1s = z1s ^ z1c ^ z2s;
   wire [PW-1:0] w1c = {(z1s & z1c) | (z1s & z2s) | (z1c & z2s), 1'b0};

   // level 5
   wire [PW-1:0] os = w1s ^ w1c ^ z2c;
   wire [PW-1:0] oc = {(w1s & w1c) | (w1s & z2c) | (w1c & z2c), 1'b0};

   //#########################################################################
   // Registers
   //#########################################################################
   always @(posedge clk) begin
      if (rst) begin
	 w       <= {DW{1'b0}};
	 Csel    <= 11'b0;
	 a_east  <= {DW{1'b0}};
	 s_south <= {PW{1'b0}};
	 c_south <= {PW{1'b0}};
      end
      else begin
	 if (ld) begin
	    w    <= w_north;
	    Csel <= (w_north[DW-1] ? 11'd256 : 11'd2048)
		    - {3'b000, w_north};
	 end
	 if (en) begin
	    a_east  <= a_west;
	    s_south <= os;
	    c_south <= oc;
	 end
      end
   end

   // Csel update shares the ld cycle (small 11-bit subtract, off the
   // compute critical path)

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

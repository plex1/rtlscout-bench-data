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
// Optimization #3: the signed 8x8 product is built from 9 operands and
// reduced together with the incoming psum pair (sin,cin) in 5 levels of
// 3:2 compressors (9 compressor cells), with no carry propagation anywhere
// in the loop. Exact identities used (a_u/w_u = unsigned bit pattern,
// a7/w7 = MSB):
//     a*w = (a_u - 128*a7)*(w_u - 128*w7)
//         = a_u*w_u - 128*a7*w_u - 128*w7*a_u + 16384*a7*w7
//     -128*a7*w_u = -a7*(w_u<<7) = a7*( (~(w_u<<7)) + 1 )  (mod 2^PW)
// The negated row is a 10-bit shifted-inverted pattern; its +1 correction
// is a single gated bit; the a7w7 term reduces to one gated bit at 2^14
// because 16384 = 2^14 (not 2^16: (a_u-128*a7) already shifts by 7).
// Free mod-2^19 sign extension of the positive-weight rows a_u<<j makes
// them usable as single 19-bit operands without any extra logic.
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
   // Product operand generation (9 operands, all congruent to a*w mod 2^PW)
   //#########################################################################
   wire [PW-1:0] q0 = {{8{1'b0}}, w[0] & a_west};                    // bit 8
   wire [PW-1:0] q1 = {{7{1'b0}}, w[1] & a_west, 1'b0};              // bit 9
   wire [PW-1:0] q2 = {{6{1'b0}}, w[2] & a_west, 2'b00};             // bit 10
   wire [PW-1:0] q3 = {{5{1'b0}}, w[3] & a_west, 3'b000};            // bit 11
   wire [PW-1:0] q4 = {{4{1'b0}}, w[4] & a_west, 4'b0000};           // bit 12
   wire [PW-1:0] q5 = {{3{1'b0}}, w[5] & a_west, 5'b00000};          // bit 13
   wire [PW-1:0] q6 = {{2{1'b0}}, w[6] & a_west, 6'b000000};         // bit 14
   // negated row: -a7*(w_u<<7) = a7*(~w_u<<7)+a7
   wire [PW-1:0] q7 = a7 ? {6'b111111, ~w, 1'b1} : {PW{1'b0}};       // bits 7..17
   wire [PW-1:0] q8 = {{(PW-15){1'b0}}, w7 & a7, 14'b00000000000000};// bit 14
   wire [PW-1:0] q9 = {{(PW-1){1'b0}}, w7};			 // bit 0

   //#########################################################################
   // 11-operand carry-save reduction: 5 levels, 9 compressors.
   //#########################################################################
   // level 1 (5 compressors)
   wire [PW-1:0] u1s = q0 ^ q1 ^ q2;
   wire [PW-1:0] u1c = {(q0 & q1) | (q0 & q2) | (q1 & q2), 1'b0};
   wire [PW-1:0] u2s = q3 ^ q4 ^ q5;
   wire [PW-1:0] u2c = {(q3 & q4) | (q3 & q5) | (q4 & q5), 1'b0};
   wire [PW-1:0] u3s = q6 ^ q7 ^ q8;
   wire [PW-1:0] u3c = {(q6 & q7) | (q6 & q8) | (q7 & q8), 1'b0};
   wire [PW-1:0] u4s = q9 ^ sin ^ cin;
   wire [PW-1:0] u4c = {(q9 & sin) | (q9 & cin) | (sin & cin), 1'b0};

   // level 2 (2 compressors)
   wire [PW-1:0] v1s = u1s ^ u1c ^ u2s;
   wire [PW-1:0] v1c = {(u1s & u1c) | (u1s & u2s) | (u1c & u2s), 1'b0};
   wire [PW-1:0] v2s = u2c ^ u3s ^ u3c;
   wire [PW-1:0] v2c = {(u2c & u3s) | (u2c & u3c) | (u3s & u3c), 1'b0};

   // level 3 (1 compressor)
   wire [PW-1:0] z1s = v1s ^ v1c ^ v2s;
   wire [PW-1:0] z1c = {(v1s & v1c) | (v1s & v2s) | (v1c & v2s), 1'b0};

   // level 4 (1 compressor)
   wire [PW-1:0] z2s = z1s ^ z1c ^ v2c;
   wire [PW-1:0] z2c = {(z1s & z1c) | (z1s & v2c) | (z1c & v2c), 1'b0};

   // level 5 (1 compressor)
   wire [PW-1:0] z3s = z2s ^ z2c ^ u4s;
   wire [PW-1:0] z3c = {(z2s & z2c) | (z2s & u4s) | (z2c & u4s), 1'b0};

   // level 6 (1 compressor) - outputs
   wire [PW-1:0] os = z3s ^ z3c ^ u4c;
   wire [PW-1:0] oc = {(z3s & z3c) | (z3s & u4c) | (z3c & u4c), 1'b0};

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

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
// Optimization #1: internal accumulator width PW=18 (exact: the sum of eight
// 8x8 signed products lies in [-130048, 129032] which fits the 18-bit signed
// range [-131072, 131071]); sign-extended to 32 at the output.
//
// Optimization #2: the psum chain is carried in REDUNDANT (carry-save) form:
// every CSA level preserves s + c = true sum (mod 2^PW), and the single
// carry-propagate add happens at each column exit.
//
// Optimization #3: the signed 8x8 product is expanded into positive-weight
// AND-plane rows plus tiny sign-correction rows, so no carry-propagating
// multiplier exists anywhere in the compute loop. Exact identity:
//   a*w = a_u*w_u - 256*a7*w_u - 256*w7*a_u + 65536*a7*w7
// with a_u/w_u the unsigned bit patterns and a7/w7 the MSBs. Mod 2^PW:
//   S1 = a7*(2^18 - 256*w_u) = a7 & {K1, 8'b0},  K1 = (2^10 - w_u) mod 2^10
//        (pure AND-gating of a load-time-registered constant)
//   S2 = w7*(2^18 - 256*a_u) = {2'b11, ~a_u, 8'b0} & {10{w7}} (w7 stationary)
//   b8 = w7 (the +2^8 inside 2^18 - 256*a_u), b16 = a7&w7 at bit 16
// All rows are non-negative patterns (free zero extension) or AND-gated
// constants; the 14 tree operands reduce in 5 carry-save levels (12 CSAs).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 18,	 // accumulator width (signed)
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
   reg [9:0]    K1;  // (2^10 - w_u) mod 2^10, registered at load

   assign w_south = w;

   // row 0 starts a fresh accumulation
   wire [PW-1:0] sin = (ROW == 0) ? {PW{1'b0}} : s_north;
   wire [PW-1:0] cin = (ROW == 0) ? {PW{1'b0}} : c_north;

   wire a7 = a_west[DW-1];
   wire w7 = w[DW-1];

   //#########################################################################
   // Product operand generation (congruent to a*w mod 2^PW)
   //#########################################################################
   // 8 AND planes: a_u * w_u
   wire [PW-1:0] P0 = {{(PW-8){1'b0}}, {8{w[0]}} & a_west};
   wire [PW-1:0] P1 = {{(PW-9){1'b0}}, {8{w[1]}} & a_west, 1'b0};
   wire [PW-1:0] P2 = {{(PW-10){1'b0}}, {8{w[2]}} & a_west, 2'b00};
   wire [PW-1:0] P3 = {{(PW-11){1'b0}}, {8{w[3]}} & a_west, 3'b000};
   wire [PW-1:0] P4 = {{(PW-12){1'b0}}, {8{w[4]}} & a_west, 4'b0000};
   wire [PW-1:0] P5 = {{(PW-13){1'b0}}, {8{w[5]}} & a_west, 5'b00000};
   wire [PW-1:0] P6 = {{(PW-14){1'b0}}, {8{w[6]}} & a_west, 6'b000000};
   wire [PW-1:0] P7 = {{(PW-15){1'b0}}, {8{w[7]}} & a_west, 7'b0000000};

   // S1 = -256*a7*w_u = a7 & {K1, 8'b0}
   wire [PW-1:0] S1 = {10{a7}} & {K1, 8'b0};
   // S2 = -256*w7*a_u = w7 & {2'b11, ~a_u, 8'b0}
   wire [PW-1:0] S2 = {w7, w7, {8{w7}} & ~a_west, 8'b0};
   // single-bit corrections
   wire [PW-1:0] B8  = {{(PW-9){1'b0}}, w7, 8'b0};
   wire [PW-1:0] B16 = {{(PW-17){1'b0}}, a7 & w7, 16'b0};

   //#########################################################################
   // 14-operand carry-save reduction: 5 levels, 12 compressors
   //#########################################################################
   // level 1
   wire [PW-1:0] f1s = P0 ^ P1 ^ P2;
   wire [PW-1:0] f1c = {(P0 & P1) | (P0 & P2) | (P1 & P2), 1'b0};
   wire [PW-1:0] f2s = P3 ^ P4 ^ P5;
   wire [PW-1:0] f2c = {(P3 & P4) | (P3 & P5) | (P4 & P5), 1'b0};
   wire [PW-1:0] f3s = P6 ^ P7 ^ S1;
   wire [PW-1:0] f3c = {(P6 & P7) | (P6 & S1) | (P7 & S1), 1'b0};
   wire [PW-1:0] f4s = S2 ^ B8 ^ B16;
   wire [PW-1:0] f4c = {(S2 & B8) | (S2 & B16) | (B8 & B16), 1'b0};

   // level 2
   wire [PW-1:0] g1s = f1s ^ f1c ^ f2s;
   wire [PW-1:0] g1c = {(f1s & f1c) | (f1s & f2s) | (f1c & f2s), 1'b0};
   wire [PW-1:0] g2s = f2c ^ f3s ^ f3c;
   wire [PW-1:0] g2c = {(f2c & f3s) | (f2c & f3c) | (f3s & f3c), 1'b0};
   wire [PW-1:0] g3s = f4s ^ f4c ^ sin;
   wire [PW-1:0] g3c = {(f4s & f4c) | (f4s & sin) | (f4c & sin), 1'b0};

   // level 3
   wire [PW-1:0] h1s = g1s ^ g1c ^ g2s;
   wire [PW-1:0] h1c = {(g1s & g1c) | (g1s & g2s) | (g1c & g2s), 1'b0};
   wire [PW-1:0] h2s = g2c ^ g3s ^ g3c;
   wire [PW-1:0] h2c = {(g2c & g3s) | (g2c & g3c) | (g3s & g3c), 1'b0};

   // level 4
   wire [PW-1:0] i1s = h1s ^ h1c ^ h2s;
   wire [PW-1:0] i1c = {(h1s & h1c) | (h1s & h2s) | (h1c & h2s), 1'b0};
   wire [PW-1:0] i2s = h2c ^ cin ^ {PW{1'b0}};
   wire [PW-1:0] i2c = {h2c & cin, 1'b0};

   // level 5
   wire [PW-1:0] os = i1s ^ i1c ^ i2s;
   wire [PW-1:0] oc = {(i1s & i1c) | (i1s & i2s) | (i1c & i2s), 1'b0};

   //#########################################################################
   // Registers
   //#########################################################################
   wire [DW-1:0] l_inv = ~w_north;
   always @(posedge clk) begin
      if (rst) begin
	 w       <= {DW{1'b0}};
	 K1      <= 10'b0;
	 a_east  <= {DW{1'b0}};
	 s_south <= {PW{1'b0}};
	 c_south <= {PW{1'b0}};
      end
      else begin
	 if (ld) begin
	    w  <= w_north;
	    K1 <= {2'b00, l_inv[7:0]} + 10'd1;   // (2^10 - w_u) mod 2^10
	 end
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
    parameter PW = 18)
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
   localparam PW  = 18;

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

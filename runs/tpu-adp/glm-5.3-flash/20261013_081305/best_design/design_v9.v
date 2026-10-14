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
// Optimization #2: the psum chain is carried in REDUNDANT (carry-save) form:
// per PE a single level of 3:2 compression, with one carry-propagate add at
// each column exit.
//
// Optimization #3: the signed 8x8 product is formed by radix-4 (Booth)
// recoding of the STATIONARY weight, so no carry-propagating multiplier
// exists anywhere in the loop:
//   w = sum_k d_k * 2^(2k),  d_k in {-2,-1,0,1,2}  (exact for signed w,
//   with the implicit sign bit w[8] = w[7]).
//   a*w = a_u*w - 256*a7*w   (a = a_u - 256*a7, a_u = unsigned bit pattern)
//       = sum_k d_k*a_u*2^(2k)  +  a7*(2^19 - 256*w_u)   (mod 2^19)
//   The a7 term is a7 * {negw, 8'b0} with negw = (2^11 - w_u) mod 2^11
//   (registered at weight load), plus the single 2^16*a7*w7 correction bit.
//   Each Booth row is a 10-bit field at position 2k: X = |d_k|*a_u (shifted
//   by d_k's MSB), complemented when d_k < 0; the per-row complement
//   constant ((2^2k - 2^(2k+10))*neg_k) is pre-summed at weight load into
//   the registered constant Creg. So the tree has 9 operands and 4 levels.
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
   reg [10:0]   negw;      // (2^11 - w_u) mod 2^11, registered at load
   reg [3:0]    negr;      // Booth sign digits, registered at load
   reg [3:0]    mag1r;     // |d_k| != 0
   reg [3:0]    mag2r;     // |d_k| == 2
   reg [PW-1:0] Creg;      // summed complement constants, registered at load

   assign w_south = w;

   // row 0 starts a fresh accumulation
   wire [PW-1:0] sin = (ROW == 0) ? {PW{1'b0}} : s_north;
   wire [PW-1:0] cin = (ROW == 0) ? {PW{1'b0}} : c_north;

   wire a7 = a_west[DW-1];
   wire [DW-1:0] a_u = a_west;

   //#########################################################################
   // Booth digits of the incoming weight (load path), k = 0..3
   // d_k = -2*w[2k+1] + w[2k] + w[2k-1]  (w[8] := w[7], w[-1] := 0)
   // mag1 = y^z, mag2 = ~mag1 & (x^y), neg = x & ~(y & z)  with x,y,z =
   // w[2k+1], w[2k], w[2k-1]
   //#########################################################################
   wire [3:0] l_m1, l_m2, l_ng;
   assign l_m1[0] = w_north[0];
   assign l_m1[1] = w_north[2] ^ w_north[1];
   assign l_m1[2] = w_north[4] ^ w_north[3];
   assign l_m1[3] = w_north[6] ^ w_north[5];
   assign l_m2[0] = ~l_m1[0] & (w_north[1] ^ w_north[0]);
   assign l_m2[1] = ~l_m1[1] & (w_north[3] ^ w_north[2]);
   assign l_m2[2] = ~l_m1[2] & (w_north[5] ^ w_north[4]);
   assign l_m2[3] = ~l_m1[3] & (w_north[7] ^ w_north[6]);
   assign l_ng[0] = w_north[1] & ~w_north[0];
   assign l_ng[1] = w_north[3] & ~(w_north[2] & w_north[1]);
   assign l_ng[2] = w_north[5] & ~(w_north[4] & w_north[3]);
   assign l_ng[3] = w_north[7] & ~(w_north[6] & w_north[5]);

   always @(posedge clk) begin
      if (rst) begin
	 w      <= {DW{1'b0}};
	 negw   <= 11'b0;
	 negr   <= 4'b0;
	 mag1r  <= 4'b0;
	 mag2r  <= 4'b0;
	 Creg   <= {PW{1'b0}};
      end
      else if (ld) begin
	 w     <= w_north;
	 negw  <= 11'd0 - w_north;   // (2^11 - w_u) mod 2^11
	 negr  <= l_ng;
	 mag1r <= l_m1;
	 mag2r <= l_m2;
	 Creg  <= (l_ng[0] ? -19'd1023  : 19'd0)
		+ (l_ng[1] ? -19'd4092  : 19'd0)
		+ (l_ng[2] ? -19'd16368 : 19'd0)
		+ (l_ng[3] ? -19'd65472 : 19'd0);
      end
   end

   //#########################################################################
   // Booth partial rows: 10-bit field at position 2k
   //#########################################################################
   wire [9:0] X0 = mag1r[0] ? (mag2r[0] ? {1'b0, a_u, 1'b0}
					: {2'b00, a_u}) : 10'b0;
   wire [9:0] X1 = mag1r[1] ? (mag2r[1] ? {1'b0, a_u, 1'b0}
					: {2'b00, a_u}) : 10'b0;
   wire [9:0] X2 = mag1r[2] ? (mag2r[2] ? {1'b0, a_u, 1'b0}
					: {2'b00, a_u}) : 10'b0;
   wire [9:0] X3 = mag1r[3] ? (mag2r[3] ? {1'b0, a_u, 1'b0}
					: {2'b00, a_u}) : 10'b0;

   wire [PW-1:0] R0 = {{(PW-10){1'b0}}, negr[0] ? ~X0 : X0};
   wire [PW-1:0] R1 = {{(PW-12){1'b0}}, negr[1] ? ~X1 : X1, 2'b00};
   wire [PW-1:0] R2 = {{(PW-14){1'b0}}, negr[2] ? ~X2 : X2, 4'b0000};
   wire [PW-1:0] R3 = {{(PW-16){1'b0}}, negr[3] ? ~X3 : X3, 6'b000000};

   // a7 correction row: a7 * (2^19 - 256*w_u) = a7 * {negw, 8'b0}
   wire [PW-1:0] N1 = a7 ? {negw, 8'b0} : {PW{1'b0}};
   // 2^16 * a7 * w7
   wire [PW-1:0] T4 = {{(PW-17){1'b0}}, a7 & w[DW-1], 16'b0};

   //#########################################################################
   // 9-operand carry-save reduction: 4 levels, 7 compressors
   //#########################################################################
   wire [PW-1:0] a1s = R0 ^ R1 ^ Creg;
   wire [PW-1:0] a1c = {(R0 & R1) | (R0 & Creg) | (R1 & Creg), 1'b0};
   wire [PW-1:0] a2s = R2 ^ R3 ^ N1;
   wire [PW-1:0] a2c = {(R2 & R3) | (R2 & N1) | (R3 & N1), 1'b0};
   wire [PW-1:0] a3s = T4 ^ sin ^ cin;
   wire [PW-1:0] a3c = {(T4 & sin) | (T4 & cin) | (sin & cin), 1'b0};

   wire [PW-1:0] b1s = a1s ^ a1c ^ a2s;
   wire [PW-1:0] b1c = {(a1s & a1c) | (a1s & a2s) | (a1c & a2s), 1'b0};
   wire [PW-1:0] b2s = a2c ^ a3s ^ a3c;
   wire [PW-1:0] b2c = {(a2c & a3s) | (a2c & a3c) | (a3s & a3c), 1'b0};

   wire [PW-1:0] c1s = b1s ^ b1c ^ b2s;
   wire [PW-1:0] c1c = {(b1s & b1c) | (b1s & b2s) | (b1c & b2s), 1'b0};

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

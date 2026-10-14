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
// Optimization #3: the signed 8x8 product is formed by radix-4 Booth
// recoding of the STATIONARY weight, so no carry-propagating multiplier
// exists anywhere in the compute loop. All Booth control signals (digit
// magnitude/sign decodes, the 2^11-w constant, the negative-row count) are
// pure COMBINATIONAL functions of the single registered weight word - the
// weight is stationary during compute, so no extra state is needed:
//   w = sum_k d_k*2^(2k), d_k in {-2,-1,0,1,2}  (w[8] := w[7], w[-1] := 0)
//   a*w = (a_u - 256*a7)*w = sum_k d_k*a_u*2^(2k) - 256*a7*w
// Each Booth row is a 10-bit field |d_k|*a_u placed at bit 2k, complemented
// (full 19-bit word) when d_k < 0; a full-word complement contributes
// exactly -1 (mod 2^PW), so the per-row corrections collapse into the 3-bit
// sum corr3 = sum_k ng_k. The -256*a7*w term is a7 * {negw, 8'b0} with
// negw = (2^11 - w) mod 2^11 (sign-extended w), which also absorbs the
// 2^16*a7*w7 piece since -256*a7*w = a7*(2^19 - 256*w) mod 2^19 for the
// SIGNED w. Total: 8 tree operands -> 4 CSA levels, 6 compressors.
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

   //#####################################################################
   // Booth digit decode of the STATIONARY weight (combinational)
   // d_k = -2*x + y + z, x=w[2k+1], y=w[2k], z=w[2k-1] (z_0 = 0)
   // mag1 = y^z, mag2 = ~mag1 & (x^y), ng = x & ~(y & z)
   //#####################################################################
   wire [3:0] b_m1, b_m2, b_ng;
   assign b_m1[0] = w[0];
   assign b_m2[0] = ~w[0] & (w[1] ^ w[0]);
   assign b_ng[0] = w[1];
   assign b_m1[1] = w[2] ^ w[1];
   assign b_m2[1] = ~b_m1[1] & (w[3] ^ w[2]);
   assign b_ng[1] = w[3] & ~(w[2] & w[1]);
   assign b_m1[2] = w[4] ^ w[3];
   assign b_m2[2] = ~b_m1[2] & (w[5] ^ w[4]);
   assign b_ng[2] = w[5] & ~(w[4] & w[3]);
   assign b_m1[3] = w[6] ^ w[5];
   assign b_m2[3] = ~b_m1[3] & (w[7] ^ w[6]);
   assign b_ng[3] = w[7] & ~(w[6] & w[5]);

   wire [10:0] negw  = 11'd0 - {{3{w[DW-1]}}, w};
   wire [2:0]  corr3 = {2'b00, b_ng[0]} + {2'b00, b_ng[1]}
		     + {2'b00, b_ng[2]} + {2'b00, b_ng[3]};

   //#####################################################################
   // Booth partial rows (runtime)
   //#####################################################################
   wire a7 = a_west[DW-1];

   wire [9:0] X0 = b_m2[0] ? {1'b0, a_west, 1'b0}
			   : (b_m1[0] ? {2'b00, a_west} : 10'b0);
   wire [9:0] X1 = b_m2[1] ? {1'b0, a_west, 1'b0}
			   : (b_m1[1] ? {2'b00, a_west} : 10'b0);
   wire [9:0] X2 = b_m2[2] ? {1'b0, a_west, 1'b0}
			   : (b_m1[2] ? {2'b00, a_west} : 10'b0);
   wire [9:0] X3 = b_m2[3] ? {1'b0, a_west, 1'b0}
			   : (b_m1[3] ? {2'b00, a_west} : 10'b0);

   wire [PW-1:0] Z0 = {{(PW-10){1'b0}}, X0};
   wire [PW-1:0] Z1 = {{(PW-12){1'b0}}, X1, 2'b00};
   wire [PW-1:0] Z2 = {{(PW-14){1'b0}}, X2, 4'b0000};
   wire [PW-1:0] Z3 = {{(PW-16){1'b0}}, X3, 6'b000000};

   wire [PW-1:0] R0 = b_ng[0] ? ~Z0 : Z0;
   wire [PW-1:0] R1 = b_ng[1] ? ~Z1 : Z1;
   wire [PW-1:0] R2 = b_ng[2] ? ~Z2 : Z2;
   wire [PW-1:0] R3 = b_ng[3] ? ~Z3 : Z3;

   wire [PW-1:0] N1 = a7 ? {negw, 8'b0} : {PW{1'b0}};
   wire [PW-1:0] CR = {{(PW-3){1'b0}}, corr3};

   // row 0 starts a fresh accumulation
   wire [PW-1:0] sin = (ROW == 0) ? {PW{1'b0}} : s_north;
   wire [PW-1:0] cin = (ROW == 0) ? {PW{1'b0}} : c_north;

   //#####################################################################
   // 8-operand carry-save reduction: 4 levels, 6 compressors
   //#####################################################################
   wire [PW-1:0] p1s = R0 ^ R1 ^ R2;
   wire [PW-1:0] p1c = {(R0 & R1) | (R0 & R2) | (R1 & R2), 1'b0};
   wire [PW-1:0] p2s = R3 ^ CR ^ N1;
   wire [PW-1:0] p2c = {(R3 & CR) | (R3 & N1) | (CR & N1), 1'b0};

   wire [PW-1:0] q1s = p1s ^ p1c ^ p2s;
   wire [PW-1:0] q1c = {(p1s & p1c) | (p1s & p2s) | (p1c & p2s), 1'b0};
   wire [PW-1:0] q2s = p2c ^ sin ^ cin;
   wire [PW-1:0] q2c = {(p2c & sin) | (p2c & cin) | (sin & cin), 1'b0};

   wire [PW-1:0] r1s = q1s ^ q1c ^ q2s;
   wire [PW-1:0] r1c = {(q1s & q1c) | (q1s & q2s) | (q1c & q2s), 1'b0};

   wire [PW-1:0] os = r1s ^ r1c ^ q2c;
   wire [PW-1:0] oc = {(r1s & r1c) | (r1s & q2c) | (r1c & q2c), 1'b0};

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

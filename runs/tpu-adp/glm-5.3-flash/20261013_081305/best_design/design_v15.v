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
// Optimization #2: the per-PE accumulate add is an explicit Sklansky
// parallel-prefix (carry-lookahead) adder: G = a&b, P = a|b, 5 prefix levels
// for 19 bits, sum = (a^b) ^ carry. This keeps the flop count minimal (one
// 19-bit accumulator per PE, no redundant carry register) while cutting the
// carry chain from ~19 ripple stages to 5 prefix levels.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // accumulator width (signed)
   (
    input		  clk,
    input		  rst,	// synchronous, active high
    input		  ld,	// weight-load shift enable
    input		  en,	// compute enable
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PW-1:0]	  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PW-1:0]   p    = {{(PW-2*DW){prod[2*DW-1]}}, prod};

   //#####################################################################
   // 19-bit Sklansky parallel-prefix adder: psum_north + p
   //#####################################################################
   wire [PW-1:0] GI = psum_north & p;
   wire [PW-1:0] PI = psum_north | p;
   wire [PW-1:0] PX = psum_north ^ p;

   // prefix state per level: GT[l*PW + i], PT[l*PW + i], l = 0..5
   wire [(5+1)*PW-1:0] GT, PT;
   generate
      for (genvar i = 0; i < PW; i = i + 1) begin : g_l0
	 assign GT[i]  = GI[i];
	 assign PT[i]  = PI[i];
      end
      for (genvar l = 1; l <= 5; l = l + 1) begin : g_lvl
	 for (genvar k = 0; k < (1 << (5 - l)); k = k + 1) begin : g_grp
	    for (genvar i = 0; i < (1 << (l - 1)); i = i + 1) begin : g_bit
	       if (k * (1 << l) + i < PW) begin : g_white
		  assign GT[l*PW + k*(1<<l) + i] =
		    GT[(l-1)*PW + k*(1<<l) + i];
		  assign PT[l*PW + k*(1<<l) + i] =
		    PT[(l-1)*PW + k*(1<<l) + i];
	       end
	       if (k * (1 << l) + (1 << (l - 1)) + i < PW) begin : g_black
		  assign GT[l*PW + k*(1<<l) + (1<<(l-1)) + i] =
		    GT[(l-1)*PW + k*(1<<l) + (1<<(l-1)) + i] |
		    (PT[(l-1)*PW + k*(1<<l) + (1<<(l-1)) + i] &
		     GT[(l-1)*PW + k*(1<<l) + (1<<(l-1)) - 1]);
		  assign PT[l*PW + k*(1<<l) + (1<<(l-1)) + i] =
		    PT[(l-1)*PW + k*(1<<l) + (1<<(l-1)) + i] &
		    PT[(l-1)*PW + k*(1<<l) + (1<<(l-1)) - 1];
	       end
	    end
	 end
      end
   endgenerate

   wire [PW-1:0] GO = GT[5*PW +: PW];
   wire [PW-1:0] sum = PX ^ {GO[PW-2:0], 1'b0};

   always @(posedge clk) begin
      if (rst) begin
	 w          <= {DW{1'b0}};
	 a_east     <= {DW{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east     <= a_west;
	    psum_south <= sum;
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

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [PW-1:0]   psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*PW +: PW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PW(PW)) u_pe
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

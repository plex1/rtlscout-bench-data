//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Schedule (identical to the original, required for cycle accuracy):
//   product A[m][i]*B[i][j] is formed at compute edge E(m)+i+j, where E(m) is
//   the compute edge that accepted activation row m.  Row i's activation is
//   delayed i cycles (input skew) and then flows east one hop per cycle; the
//   partial sums flow south; column j's finished result is delayed 7-j cycles
//   (output de-skew) so the whole result row leaves together, 2N-1 = 15
//   compute cycles after its A row was accepted.
//
// Optimization vs. the original: the accumulators are PW = 2*DW + clog2(N) =
// 19 bits instead of 32.  Every partial sum is a sum of at most N signed
// DWxDW products, |.| <= N*2^(2DW-2) = 131072 = 2^17, which fits in 19 signed
// bits, so the 19-bit datapath is exact; results are sign-extended to 32 bits
// at the outputs only.  This shrinks 64 accumulator registers/adders and the
// 28 de-skew registers with zero behavior change.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // accumulator / partial-sum width (signed)
   (
    input		      clk,
    input		      rst,	// synchronous, active high
    input		      ld,	// weight-load shift enable
    input		      en,	// compute enable
    input [DW-1:0]	      w_north,	// weight in from PE above
    input [DW-1:0]	      a_west,	// activation in from west
    input [PW-1:0]	      psum_north,// partial sum in from above
    output [DW-1:0]	      w_south,	// stored weight out to PE below
    output reg [DW-1:0]	      a_east,	// activation out to east
    output reg [PW-1:0]	      psum_south // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
	 w	   <= {DW{1'b0}};
	 a_east	   <= {DW{1'b0}};
	 psum_south <= {PW{1'b0}};
      end
      else begin
	 if (ld)
	   w <= w_north;
	 if (en) begin
	    a_east	   <= a_west;
	    psum_south <= psum_north +
			  {{(PW-2*DW){prod[2*DW-1]}}, prod};
	 end
      end
   end

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		      clk,
    input		      rst,	 // synchronous, active high
    input		      w_valid,	 // weight-row valid (load phase)
    input [N*DW-1:0]	      w_data,	 // weight row, packed { col N-1 ... col 0 }
    input		      a_valid,	 // activation-row valid (compute phase)
    input [N*DW-1:0]	      a_data,	 // activation row, packed { col N-1 ... col 0 }
    output		      c_valid,	 // result-row valid
    output [N*ACCW-1:0]	      c_data	 // result row, packed { col N-1 ... col 0 }
    );

   localparam PW  = 2*DW + $clog2(N);
   localparam LAT = 2*N - 1;

   wire en = ~w_valid;

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles.  Zero is
   // injected while a_valid is low so idle/drain cycles add nothing.
   //##########################################################################
   wire [N*DW-1:0] a_left;

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer	    k;
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

   genvar	   gi, gj;
   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [PW-1:0]	  psum_v [0:N][0:N-1];
   wire [DW-1:0]	  w_v    [0:N][0:N-1];

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_west
	 assign a_h[gi][0] = a_left[gi*DW +: DW];
      end
      for (gj = 0; gj < N; gj = gj + 1) begin : g_edge
	 assign psum_v[0][gj] = {PW{1'b0}};
	 assign w_v[0][gj]    = w_data[gj*DW +: DW];
	 assign c_bot[gj*PW +: PW] = psum_v[N][gj];
      end
      for (gi = 0; gi < N; gi = gi + 1) begin : g_pe_row
	 for (gj = 0; gj < N; gj = gj + 1) begin : g_pe_col
	    tpu_pe #(.DW(DW), .PW(PW)) u_pe
	      (.clk(clk), .rst(rst), .ld(w_valid), .en(en),
	       .w_north(w_v[gi][gj]), .a_west(a_h[gi][gj]),
	       .psum_north(psum_v[gi][gj]),
	       .w_south(w_v[gi+1][gj]), .a_east(a_h[gi][gj+1]),
	       .psum_south(psum_v[gi+1][gj]));
	 end
      end
   endgenerate

   //##########################################################################
   // Output de-skew: column j is delayed N-1-j compute cycles so all columns
   // of a result row line up.  Values are sign-extended to ACCW here.
   //##########################################################################
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_deskew
	 localparam D = N - 1 - gj;
	 wire [PW-1:0] cin = c_bot[gj*PW +: PW];
	 if (D == 0) begin : g_deskew0
	    assign c_data[gj*ACCW +: ACCW] =
	      {{(ACCW-PW){cin[PW-1]}}, cin};
	 end
	 else begin : g_deskewn
	    reg [PW-1:0] dr [0:D-1];
	    integer	     k;
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
	    assign c_data[gj*ACCW +: ACCW] =
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

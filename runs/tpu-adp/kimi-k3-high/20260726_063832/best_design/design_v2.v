//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) - optimized.
//
// Optimizations vs reference:
//  - Internal partial-sum width reduced to 19 bits: 8 signed 8x8 products
//    (each in [-16256, 16384]) sum to [-130048, 131072], which fits in
//    signed 19 bits with margin. Results are sign-extended to 32 bits at
//    the output (free). Shrinks adders, psum regs, and deskew regs and
//    shortens the critical carry chain.
//  - Datapath resets removed (only the c_valid pipeline keeps one). Reset
//    is asserted only before streaming; any stale register state drains out
//    of the array before the first cycle c_valid can possibly be high, so
//    the resets are functionally unobservable.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PSW = 19)	 // internal partial-sum width (signed)
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PSW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PSW-1:0]  psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // signed 8x8 product, sign-extended to partial-sum width
   wire signed [PSW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east     <= a_west;
         psum_south <= psum_north + prod;
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PSW = 19)
   (
    input		clk,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*PSW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [PSW-1:0]  psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PSW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*PSW +: PSW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PSW(PSW)) u_pe
              (.clk        (clk),
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
    parameter ACCW = 32,
    parameter PSW = 19)	// internal accumulator width (sign-extended to ACCW)
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

   //##########################################################################
   // Input skew: row i delayed i cycles. Zeros injected when a_valid is low.
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
               if (en) begin
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
   // The systolic array (19-bit internal partial sums).
   //##########################################################################
   wire [N*PSW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PSW(PSW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew on 19-bit values, sign-extended to ACCW afterwards.
   //##########################################################################
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PSW-1:0] cin = c_bot[j*PSW +: PSW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-PSW){cin[PSW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PSW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-PSW){dr[D-1][PSW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   // (The only reset in the design: c_valid must be deterministic from
   // cycle zero since it is checked every cycle.)
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

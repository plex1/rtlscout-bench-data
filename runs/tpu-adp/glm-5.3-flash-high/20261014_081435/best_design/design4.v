// Optimized TPU MXU tile (logikbench tpu, cycle-accurate) - FMA PE variant.
//
// Same structure/timing as the golden logikbench RTL (PE mesh, per-row input
// skew, eastward activation shift, southward psum accumulation, per-column
// output de-skew), with:
//   * Internal accumulators narrowed to 19 bits (max |sum| = 8*16384 = 131072
//     fits exactly in 19 signed bits); sign-extend to 32 at the pins only.
//   * FMA-style PE: the stationary weight multiply is expanded into shifted
//     partial products of the activation (sign-extended to 19 bits), which are
//     compressed together with psum_in in a single adder tree instead of a
//     multiplier CPA cascaded into an accumulator CPA.
//   * Reset only on the c_valid pipeline; datapath registers are naturally
//     flushed by the trace (zero-weight/zero-A cycles after reset).

module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19)
   (
    input		  clk,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // activation sign-extended to accumulator width
   wire signed [ACCW-1:0] ax = {{(ACCW-DW){a_west[DW-1]}}, a_west};

   // psum + a*w as one compression tree (mod 2^ACCW, true result fits ACCW)
   wire [ACCW-1:0] pp0 = w[0] ? ax : {ACCW{1'b0}};		// +ax
   wire [ACCW-1:0] pp1 = w[1] ? (ax <<< 1) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp2 = w[2] ? (ax <<< 2) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp3 = w[3] ? (ax <<< 3) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp4 = w[4] ? (ax <<< 4) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp5 = w[5] ? (ax <<< 5) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp6 = w[6] ? (ax <<< 6) : {ACCW{1'b0}};
   wire [ACCW-1:0] pp7 = w[7] ? (~((ax <<< 7)) + 1'b1) : {ACCW{1'b0}}; // -ax<<7

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east     <= a_west;
         psum_south <= psum_north + pp0 + pp1 + pp2 + pp3 +
                       pp4 + pp5 + pp6 + pp7;
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,   // external result width
    parameter PACCW = 19)  // internal accumulator width
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

   // Input skew: row i delayed i cycles; zero injected when a_valid low.
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

   // PE mesh
   wire [DW-1:0]    a_h    [0:N-1][0:N];
   wire [PACCW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]    w_v    [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PACCW{1'b0}};
         assign w_v[0][j]    = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(PACCW)) u_pe
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

   // Output de-skew: column j delayed N-1-j cycles; sign-extend to 32 bits.
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = psum_v[N][j];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PACCW){cin[PACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PACCW){dr[D-1][PACCW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

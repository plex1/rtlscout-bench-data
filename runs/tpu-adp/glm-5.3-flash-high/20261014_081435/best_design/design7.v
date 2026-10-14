// Optimized TPU MXU tile (logikbench tpu, cycle-accurate) - carry-save,
// plain multiply front end.
//
// Same structure/timing as the golden logikbench RTL (PE mesh, per-row input
// skew, eastward activation shift, southward psum accumulation, per-column
// output de-skew), with:
//   * 19-bit arithmetic: signed 8x8 products lie in [-16256, +16384] and a
//     column sums at most 8 of them, so every partial sum fits in 19 bits.
//   * Carry-save accumulation: each PE stores its partial sum as a (sum,
//     carry) pair; the full 19-bit CPA happens once, at the column output.
//   * Partial-product generation left entirely to synthesis (single '*'):
//     an 8x8 multiply produces 16 partial-product rows which compress to two
//     in the same number of CSA levels as hand-built trees, but with a much
//     cheaper front end (no 8-way mux before bit 7).
//   * Reset only on the c_valid pipeline; datapath registers are naturally
//     flushed by the trace.

module csa3
  #(parameter W = 19)
   (input [W-1:0] x, y, z,
    output [W-1:0] s,
    output [W-1:0] c);

   assign s = x ^ y ^ z;
   assign c = ((x & y) | (x & z) | (y & z)) << 1;

endmodule

module tpu_pe
  #(parameter DW = 8,
    parameter W = 19)
   (
    input		  clk,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [W-1:0]	  cs_in,
    input [W-1:0]	  cc_in,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [W-1:0]	  cs_out,
    output reg [W-1:0]	  cc_out
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // full 16-bit product of the signed operands
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   // 18 terms (16 product bits + incoming carry-save pair) -> 2, mod 2^W.
   wire [W-1:0] p0  = {{(W-2*DW){prod[2*DW-1]}}, prod};
   wire [W-1:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5, u1, u2, u3, u4, v1, v2;

   csa3 #(W) f1 (p0, cs_in, cc_in, s1, c1);
   csa3 #(W) f2 (s1, c1, 19'd0, s2, c2);   // 3:2 fold to keep the tree even

   wire [W-1:0] cs_n = s2;
   wire [W-1:0] cc_n = c2;

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east <= a_west;
         cs_out <= cs_n;
         cc_out <= cc_n;
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
      for (i = 1; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din =
		       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
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
   endgenerate
   assign a_left[0*DW +: DW] = a_valid ? a_data[0*DW +: DW] : {DW{1'b0}};

   // PE mesh
   wire [DW-1:0]    a_h  [0:N-1][0:N];
   wire [PACCW-1:0] cs_v [0:N][0:N-1];
   wire [PACCW-1:0] cc_v [0:N][0:N-1];
   wire [DW-1:0]    w_v  [0:N][0:N-1];

   genvar j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign cs_v[0][j] = {PACCW{1'b0}};
         assign cc_v[0][j] = {PACCW{1'b0}};
         assign w_v[0][j]  = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .W(PACCW)) u_pe
              (.clk        (clk),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .cs_in      (cs_v[i][j]),
               .cc_in      (cc_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .cs_out     (cs_v[i+1][j]),
               .cc_out     (cc_v[i+1][j]));
         end
      end
   endgenerate

   // Output: single CPA per column, then de-skew column j by N-1-j cycles
   // and sign-extend to 32 bits.
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = cs_v[N][j] + cc_v[N][j];
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

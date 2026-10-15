// Optimized TPU MXU, v5:
//  - Per-row tapered partial-sum widths. Row i accumulates at most i+1
//    products (the top edge injects 0 every cycle), and each product has
//    magnitude <= 2^14, so |psum at row i| <= (i+1)*2^14. A signed PW-bit
//    register covers magnitudes up to 2^(PW-1); equality is fine because
//    adding a further product of either sign stays representable:
//      row 0      -> 16 bits (1 product,  max 16384)
//      row 1      -> 17 bits (2 products, max 32768)
//      rows 2-3   -> 18 bits (<=4 products, max 65536)
//      rows 4-7   -> 19 bits (<=8 products, max 131072)
//    Widths are non-decreasing going south, so boundaries connect bit-for-bit.
//  - Generic Kogge-Stone prefix adders for the accumulation.
//  - No datapath resets (sim starts zero; values always loaded/flushed before
//    observation). Only vpipe keeps a reset.

module ksadd #(parameter W = 19) (
    input  [W-1:0] A,
    input  [W-1:0] B,
    output [W-1:0] S
);
   // flattened levels: G[l*W + k], P[l*W + k], levels 0..5
   wire [6*W-1:0] G, P;
   assign G[W-1:0] = A & B;
   assign P[W-1:0] = A ^ B;
   genvar l, k;
   generate
      for (l = 1; l <= 5; l = l + 1) begin : lvl
         for (k = 0; k < W; k = k + 1) begin : b
            if (k >= (1 << l)) begin : merge
               assign G[l*W+k] = G[(l-1)*W+k] |
                                 (P[(l-1)*W+k] & G[(l-1)*W+k-(1<<l)]);
               assign P[l*W+k] = P[(l-1)*W+k] & P[(l-1)*W+k-(1<<l)];
            end
            else begin : pass
               assign G[l*W+k] = G[(l-1)*W+k];
               assign P[l*W+k] = P[(l-1)*W+k];
            end
         end
      end
   endgenerate
   assign S = P[W-1:0] ^ {G[5*W + W-2 : 5*W], 1'b0};
endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PW = 19)	 // this PE's partial-sum width (signed)
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [PW-1:0]	  psum_north, // partial sum in from above
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output reg [PW-1:0]	  psum_south  // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PW-1:0] prod_x = {{(PW-2*DW){prod[2*DW-1]}}, prod};

   wire [PW-1:0] nxt;
   ksadd #(.W(PW)) u_add (.A(psum_north), .B(prod_x), .S(nxt));

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east     <= a_west;
         psum_south <= nxt;
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8)	 // operand width (signed)
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [151:0]	c_bot	// column results, one per column (19b each)
    );

   // tapered widths per row
   function integer pwf;
      input integer i;
      begin
         case (i)
            0: pwf = 16;
            1: pwf = 17;
            2,3: pwf = 18;
            default: pwf = 19;
         endcase
      end
   endfunction

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [18:0]     psum_v [0:N][0:N-1]; // extra top bits are unused (0)
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = 19'd0;
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*19 +: 19] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PW(pwf(i))) u_pe
              (.clk        (clk),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j][pwf(i)-1:0]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (psum_v[i+1][j][pwf(i)-1:0]));
         end
      end
      // sign-extend the unused upper bits of each psum net
      for (i = 1; i <= N; i = i + 1) begin : g_tie
         for (j = 0; j < N; j = j + 1) begin : g_tiec
            if (pwf(i-1) < 19) begin : tiehi
               assign psum_v[i][j][18:pwf(i-1)] =
                      {19-pwf(i-1){psum_v[i][j][pwf(i-1)-1]}};
            end
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32, // result width (signed), at the port
    parameter PW = 19)	 // output partial-sum width
   (
    input		clk,
    input		rst,	 // synchronous, active high
    input		w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]	w_data,	 // weight row, packed { col N-1 ... col 0 }
    input		a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]	a_data,	 // activation row, packed { col N-1 ... col 0 }
    output		c_valid, // result-row valid
    output [N*ACCW-1:0]	c_data	 // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i delayed i compute cycles; zeros injected when !a_valid.
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

   // Systolic array with tapered partial sums, 19-bit output.
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW)) u_array
     (.clk(clk), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output de-skew at PW bits, then sign-extend to ACCW (pure wiring).
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
               if (en) begin
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

   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

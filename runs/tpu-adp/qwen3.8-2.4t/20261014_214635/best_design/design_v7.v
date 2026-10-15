// Optimized TPU MXU, v7:
//  Sign-magnitude PE front-end + tapered internal widths + KS adder w/ cin.
//
//  PE structure:
//   - weight stored as magnitude wm (8b) + sign ws, precomputed at load time
//     (removes the per-cycle weight abs/negate from the critical path).
//   - activation registered in two's complement; magnitude computed
//     combinationally (small 8-bit incrementer).
//   - unsigned multiplier: mag = a_mag * wm; product sign s = as ^ ws.
//   - accumulate: new = ps + prod with prod = s ? ~mag+1 : mag. Implemented
//     as one PW-bit Kogge-Stone adder: B = mag ^ s, cin = s (the +1 fold-in).
//     Bits of B above 15 are just s (since mag[15]=0), pure wiring.
//   - tapered psum register widths per row (proven bounds, see v6).

module ksadd #(parameter W = 16) (
    input  [W-1:0] A,
    input  [W-1:0] B,
    input	       cin,
    output [W-1:0] S
);
   // Kogge-Stone prefix with offsets 1,2,4,... ; 5 levels cover W <= 19.
   localparam L = (W > 16) ? 5 : 4;
   wire [(L+1)*W-1:0] G, P;
   assign G[W-1:0] = A & B;
   assign P[W-1:0] = A ^ B;
   genvar l, k;
   generate
      for (l = 1; l <= L; l = l + 1) begin : lvl
         // bit 0: fold cin in at level 1, pass through afterwards
         if (l == 1) begin : b0f
            assign G[W+0] = G[0] | (P[0] & cin);
            assign P[W+0] = P[0];
         end
         else begin : b0p
            assign G[l*W+0] = G[(l-1)*W+0];
            assign P[l*W+0] = P[(l-1)*W+0];
         end
         for (k = 1; k < W; k = k + 1) begin : b
            if (k >= (1 << (l-1))) begin : merge
               assign G[l*W+k] = G[(l-1)*W+k] |
                                 (P[(l-1)*W+k] & G[(l-1)*W+k-(1<<(l-1))]);
               assign P[l*W+k] = P[(l-1)*W+k] & P[(l-1)*W+k-(1<<(l-1))];
            end
            else begin : pass
               assign G[l*W+k] = G[(l-1)*W+k];
               assign P[l*W+k] = P[(l-1)*W+k];
            end
         end
      end
   endgenerate
   assign S[0] = P[0] ^ cin;
   // S[i] = P0[i] ^ carry_in[i] = P0[i] ^ G_final[i-1]
   assign S[W-1:1] = P[W-1:1] ^ G[L*W + W-2 : L*W];
endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter AW = 19,	 // inter-row (port) partial-sum width
    parameter PW = 19)	 // internal (register/adder) partial-sum width
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above (2's comp)
    input [DW-1:0]	  a_west,     // activation in from west
    input [AW-1:0]	  psum_north, // partial sum in from above (AW bits)
    output [DW-1:0]	  w_south,    // stored weight out to PE below (2's comp)
    output reg [DW-1:0]	  a_east,     // activation out to east
    output [AW-1:0]	  psum_south  // partial sum out to below (AW bits)
    );

   // stationary weight, stored sign-magnitude (precomputed at load)
   reg [DW-1:0] wm;	// magnitude
   reg	       ws;	// sign
   wire [DW-1:0] w_mag_in = w_north[DW-1] ? (~w_north + 1'b1) : w_north;
   assign w_south = ws ? (~wm + 1'b1) : wm; // back to two's complement

   // activation magnitude, straight off a_west (matches original timing:
   // product is computed from the PRE-registered activation)
   wire [DW-1:0] a_mag = a_west[DW-1] ? (~a_west + 1'b1) : a_west;
   wire	       as = a_west[DW-1];

   // unsigned product magnitude + sign
   wire [2*DW-1:0] mag = a_mag * wm;
   wire	       s = as ^ ws;

   // accumulation: ps + (s ? ~mag+1 : mag), via KS adder with cin = s
   wire [PW-1:0] ps_in = psum_north[PW-1:0];
   wire [PW-1:0] b_op  = {{(PW-16){s}}, mag ^ {16{s}}};
   wire [PW-1:0] nxt;
   ksadd #(.W(PW)) u_add (.A(ps_in), .B(b_op), .cin(s), .S(nxt));

   reg [PW-1:0] ps;
   always @(posedge clk) begin
      if (ld) begin
         wm <= w_mag_in;
         ws <= w_north[DW-1];
      end
      if (en) begin
         a_east <= a_west;
         ps     <= nxt;
      end
   end

   assign psum_south = {{(AW-PW){ps[PW-1]}}, ps};

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter AW = 19)	 // inter-row partial-sum width
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*AW-1:0]	c_bot	// column results, one per column (AW bits)
    );

   // per-row internal width (non-decreasing going south)
   function integer pwf;
      input integer i;
      begin
         case (i)
            0:      pwf = 16;
            1:      pwf = 17;
            2:      pwf = 17;
            3,4,5,6: pwf = 18;
            default: pwf = 19;
         endcase
      end
   endfunction

   wire [DW-1:0]  a_h    [0:N-1][0:N];
   wire [AW-1:0]  psum_v [0:N][0:N-1];
   wire [DW-1:0]  w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*AW +: AW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW), .PW(pwf(i))) u_pe
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

   // Systolic array with tapered partial sums, AW-bit output.
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .AW(PW)) u_array
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

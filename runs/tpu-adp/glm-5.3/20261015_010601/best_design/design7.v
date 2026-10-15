//#############################################################################
// Cycle-exact rewrite of the logikbench TPU (weight-stationary MXU tile).
//
// Same product schedule as the original array: product (i,j) of the row
// accepted at compute step t is computed at step t+i+j with the weight in
// force at that step (this is what makes mid-drain weight reloads work).
//
// Delay optimization: each cell computes  S' + 2C' = S + 2C + a*w  with a
// radix-4 Booth recoder feeding a 4-level carry-save tree, so a cycle never
// contains a carry-propagate adder.  The products' partial products are
// summed straight into the accumulator's redundant form; the only CPAs in
// the design are the eight 19-bit final column adders, which sit in the
// de-skew slack (and combinationally on the output for column 7).
//
// Booth details: groups (a[2i+1],a[2i],a[2i-1]) select PP in
// {0,+w,+2w,~w,~2w}; a negative selection contributes ~X, so a correction
// constant K = sum_i neg_i<<2i is added as a seventh tree operand.  All
// arithmetic is modulo 2^19, which is exact because |column sum| <= 2^17.
//#############################################################################

module tpu_cell (
   input		 clk,
   input		 en,
   input [7:0]		 a,
   input [7:0]		 w,
   input [18:0]		 Sin,
   input [18:0]		 Cin,
   output reg [18:0]	 Sout,
   output reg [18:0]	 Cout
   );

   wire [18:0] ws = {{11{w[7]}}, w};
   wire [18:0] w2 = {ws[17:0], 1'b0};

   wire [2:0] b [0:3];
   assign b[0] = {a[1], a[0], 1'b0};
   assign b[1] = {a[3], a[2], a[1]};
   assign b[2] = {a[5], a[4], a[3]};
   assign b[3] = {a[7], a[6], a[5]};

   wire [18:0] T [0:3];
   wire [3:0]	 neg;

   genvar g;
   generate
      for (g = 0; g < 4; g = g + 1) begin : g_pp
	 wire b2 = b[g][2], b1 = b[g][1], b0 = b[g][0];
	 // digit = -2*b2 + b1 + b0 :  +1 for 001/010, +2 for 011,
	 // -1 for 101/110, -2 for 100, 0 for 000/111
	 wire p1 = ~b2 & (b1 ^ b0);	  // +w
	 wire p2 = ~b2 & b1 & b0;	  // +2w
	 wire m1 = b2 & (b1 ^ b0);	  // -w
	 wire m2 = b2 & ~b1 & ~b0;	  // -2w
	 assign neg[g] = b2 & ~(b1 & b0);
	 wire [18:0] pp = p2 ? w2 : m2 ? ~w2 : p1 ? ws : m1 ? ~ws : 19'd0;
	 assign T[g] = pp <<< (2*g);
      end
   endgenerate

   wire [18:0] K = {11'd0, neg[3], 1'b0, neg[2], 1'b0, neg[1], 1'b0, neg[0]};

   wire [18:0] y0 = {Cin[17:0], 1'b0};

   // level 1
   wire [18:0] a0 = Sin ^ y0 ^ T[0];
   wire [18:0] a1 = (Sin & y0) | (Sin & T[0]) | (y0 & T[0]);
   wire [18:0] b0t = T[1] ^ T[2] ^ T[3];
   wire [18:0] b1t = (T[1] & T[2]) | (T[1] & T[3]) | (T[2] & T[3]);
   wire [18:0] ya1 = {a1[17:0], 1'b0};
   wire [18:0] yb1 = {b1t[17:0], 1'b0};

   // level 2
   wire [18:0] c0 = a0 ^ ya1 ^ b0t;
   wire [18:0] c1 = (a0 & ya1) | (a0 & b0t) | (ya1 & b0t);
   wire [18:0] yc1 = {c1[17:0], 1'b0};

   // level 3
   wire [18:0] d0 = c0 ^ yc1 ^ yb1;
   wire [18:0] d1 = (c0 & yc1) | (c0 & yb1) | (yc1 & yb1);
   wire [18:0] yd1 = {d1[17:0], 1'b0};

   // level 4
   wire [18:0] e0 = d0 ^ yd1 ^ K;
   wire [18:0] e1 = (d0 & yd1) | (d0 & K) | (yd1 & K);

   always @(posedge clk) begin
      if (en) begin
	 Sout <= e0;
	 Cout <= e1;
      end
   end

endmodule

//#############################################################################

module tpu (
   input		 clk,
   input		 rst,
   input		 w_valid,
   input [63:0]		 w_data,
   input		 a_valid,
   input [63:0]		 a_data,
   output		 c_valid,
   output [255:0]	 c_data
   );

   wire en = ~w_valid;

   // weight storage: row i at W[i*64 +: 64], byte j at + j*8
   reg [511:0] W;
   always @(posedge clk) begin
      if (w_valid) W <= {W[447:0], w_data};
   end

   // activation skew chains: tap for (i,j) = element i delayed (i+j)
   wire [511:0] tap;

   genvar i, j;
   generate
      for (i = 0; i < 8; i = i + 1) begin : g_arow
	 localparam L = i + 7;
	 reg [8*L-1:0] sr;
	 wire [7:0]    din = a_valid ? a_data[i*8 +: 8] : 8'd0;
	 always @(posedge clk) begin
	    if (en)	 sr <= {sr[8*L-9:0], din};
	 end
	 for (j = 0; j < 8; j = j + 1) begin : g_atap
	    if (i + j == 0) begin : g_z
	       assign tap[7:0] = din;
	    end
	    else begin : g_n
	       assign tap[(i*8+j)*8 +: 8] = sr[8*(i+j-1) +: 8];
	    end
	 end
      end
   endgenerate

   // per-column carry-save accumulation chains
   generate
      for (j = 0; j < 8; j = j + 1) begin : g_col
	 wire [18:0] S [0:7];
	 wire [18:0] C [0:7];

	 for (i = 0; i < 8; i = i + 1) begin : g_pe
	    if (i == 0) begin : g_p0
	       tpu_cell u_cell
		 (.clk(clk), .en(en),
		  .a(tap[j*8 +: 8]), .w(W[j*8 +: 8]),
		  .Sin(19'd0), .Cin(19'd0),
		  .Sout(S[0]), .Cout(C[0]));
	    end
	    else begin : g_pn
	       tpu_cell u_cell
		 (.clk(clk), .en(en),
		  .a(tap[(i*8+j)*8 +: 8]), .w(W[i*64 + j*8 +: 8]),
		  .Sin(S[i-1]), .Cin(C[i-1]),
		  .Sout(S[i]), .Cout(C[i]));
	    end
	 end

	 // single carry-propagate adder per column (carry vector weighted x2)
	 wire signed [18:0] colsum = S[7] + {C[7][17:0], 1'b0};

	 localparam D = 7 - j;
	 if (D == 0) begin : g_d0
	    assign c_data[j*32 +: 32] = {{13{colsum[18]}}, colsum};
	 end
	 else begin : g_dn
	    reg signed [18:0] d [0:D-1];
	    integer t;
	    always @(posedge clk) begin
	       if (en) begin
		  d[0] <= colsum;
		  for (t = 1; t < D; t = t + 1)
		    d[t] <= d[t-1];
	       end
	    end
	    assign c_data[j*32 +: 32] = {{13{d[D-1][18]}}, d[D-1]};
	 end
      end
   endgenerate

   reg [14:0] vp;
   always @(posedge clk) begin
      if (rst)	  vp <= 15'd0;
      else if (en) vp <= {vp[13:0], a_valid};
   end
   assign c_valid = vp[14];

endmodule

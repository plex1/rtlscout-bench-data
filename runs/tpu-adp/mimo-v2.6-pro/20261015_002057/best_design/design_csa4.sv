// TPU MXU tile (weight-stationary systolic array), cycle-accurate to the
// logikbench 'tpu' benchmark. Optimized datapath:
//   * accumulators 19 bits (|C| <= 8*16384 = 131072)
//   * PE product as an explicit Baugh-Wooley partial-product matrix (64 AND
//     terms, sign corrections folded into row 7 and a constant)
//   * column partial sums carried south in 4 redundant words: each PE absorbs
//     its 8 product rows into the 4 words with only 3 carry-save levels (8
//     rows + 4 words -> 4 words); one 19-bit CPA per column (plus two CSAs)
//     converts to binary at the array bottom.
module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter AW = 19)	 // accumulator width (signed)
   (
    input		    clk,
    input		    rst,	// synchronous, active high
    input		    w_valid,
    input [N*DW-1:0]	    w_data,
    input		    a_valid,
    input [N*DW-1:0]	    a_data,
    output		    c_valid,
    output [N*32-1:0]	    c_data
    );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   // 3:2 carry-save compressor (mod 2^AW; true values always fit in AW bits)
   function [AW-1:0] csa_s;
      input [AW-1:0] x, y, z;
      begin
	 csa_s = x ^ y ^ z;
      end
   endfunction
   function [AW-1:0] csa_c;
      input [AW-1:0] x, y, z;
      begin
	 csa_c = ((x & y) | (y & z) | (x & z)) << 1;
      end
   endfunction

   wire [N*DW-1:0] a_left;
   wire [DW-1:0] a_h [0:N-1][0:N];   // west -> east activations
   wire [DW-1:0] w_v [0:N][0:N-1];   // north -> south weight load chain
   wire [AW-1:0] p0 [0:N][0:N-1];    // north -> south partial sums (4 words)
   wire [AW-1:0] p1 [0:N][0:N-1];
   wire [AW-1:0] p2 [0:N][0:N-1];
   wire [AW-1:0] p3 [0:N][0:N-1];

   genvar i, j, gp, gb;
   generate
      // ---------------- input skew: row i delayed i compute cycles --------
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer k;
	    always @(posedge clk) begin
	       if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < i; k = k + 1) sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[i*DW +: DW] = sr[i-1];
	 end
      end

      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign p0[0][j] = {AW{1'b0}};
	 assign p1[0][j] = {AW{1'b0}};
	 assign p2[0][j] = {AW{1'b0}};
	 assign p3[0][j] = {AW{1'b0}};
	 assign w_v[0][j] = w_data[j*DW +: DW];
      end

      // ---------------- PE mesh ------------------------------------------
      for (i = 0; i < N; i = i + 1) begin : g_row
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    reg [DW-1:0] wreg;
	    reg [DW-1:0] areg;
	    reg [AW-1:0] q0, q1, q2, q3;

	    wire [DW-1:0] a_w = a_h[i][j];
	    wire [AW-1:0] t0 = p0[i][j];
	    wire [AW-1:0] t1 = p1[i][j];
	    wire [AW-1:0] t2 = p2[i][j];
	    wire [AW-1:0] t3 = p3[i][j];

	    // Baugh-Wooley partial products: 8x8 AND matrix, sign corrections
	    // folded into row 7 and the constant K.
	    //   prod == sum_p (row_p << p) + K   (mod 2^AW)
	    wire [7:0] prow [0:7];
	    for (gp = 0; gp < 8; gp = gp + 1) begin : g_pp
	       for (gb = 0; gb < 8; gb = gb + 1) begin : g_pb
		  if ((gp < 7) && (gb < 7))
		    assign prow[gp][gb] = a_w[gp] & wreg[gb];
		  else if (gp < 7)
		    assign prow[gp][gb] = ~(a_w[gp] & wreg[7]);
		  else if (gb < 7)
		    assign prow[gp][gb] = ~(a_w[7] & wreg[gb]);
		  else
		    assign prow[gp][gb] = a_w[7] & wreg[7];
	       end
	    end
	    wire [AW-1:0] r0 = {{(AW-8){1'b0}}, prow[0]};
	    wire [AW-1:0] r1 = {{(AW-8){1'b0}}, prow[1]} << 1;
	    wire [AW-1:0] r2 = {{(AW-8){1'b0}}, prow[2]} << 2;
	    wire [AW-1:0] r3 = {{(AW-8){1'b0}}, prow[3]} << 3;
	    wire [AW-1:0] r4 = {{(AW-8){1'b0}}, prow[4]} << 4;
	    wire [AW-1:0] r5 = {{(AW-8){1'b0}}, prow[5]} << 5;
	    wire [AW-1:0] r6 = {{(AW-8){1'b0}}, prow[6]} << 6;
	    wire [AW-1:0] r7 = {{(AW-8){1'b0}}, prow[7]} << 7;
	    wire [AW-1:0] cr = 19'h78100;	 // constant correction

	    // carry-save tree: 8 rows + 4 incoming words -> 4 words (3 levels)
	    wire [AW-1:0] s1a = csa_s(r0, r1, r2);
	    wire [AW-1:0] c1a = csa_c(r0, r1, r2);
	    wire [AW-1:0] s1b = csa_s(r3, r4, r5);
	    wire [AW-1:0] c1b = csa_c(r3, r4, r5);
	    wire [AW-1:0] s1c = csa_s(r6, r7, cr);
	    wire [AW-1:0] c1c = csa_c(r6, r7, cr);
	    wire [AW-1:0] s1d = csa_s(t0, t1, t2);
	    wire [AW-1:0] c1d = csa_c(t0, t1, t2);
	    wire [AW-1:0] s2a = csa_s(s1a, c1a, s1b);
	    wire [AW-1:0] c2a = csa_c(s1a, c1a, s1b);
	    wire [AW-1:0] s2b = csa_s(c1b, s1c, c1c);
	    wire [AW-1:0] c2b = csa_c(c1b, s1c, c1c);
	    wire [AW-1:0] s2c = csa_s(s1d, c1d, t3);
	    wire [AW-1:0] c2c = csa_c(s1d, c1d, t3);
	    wire [AW-1:0] s3a = csa_s(s2a, c2a, s2b);
	    wire [AW-1:0] c3a = csa_c(s2a, c2a, s2b);
	    wire [AW-1:0] s3b = csa_s(c2b, s2c, c2c);
	    wire [AW-1:0] c3b = csa_c(c2b, s2c, c2c);

	    always @(posedge clk) begin
	       if (ld) wreg <= w_v[i][j];
	       if (en) begin
		  areg <= a_w;
		  q0 <= s3a;
		  q1 <= c3a;
		  q2 <= s3b;
		  q3 <= c3b;
	       end
	    end

	    assign w_v[i+1][j] = wreg;
	    assign a_h[i][j+1] = areg;
	    assign p0[i+1][j] = q0;
	    assign p1[i+1][j] = q1;
	    assign p2[i+1][j] = q2;
	    assign p3[i+1][j] = q3;
	 end
      end

      // ---------------- column bottom: reduce 4 words -> binary -----------
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [AW-1:0] bs = csa_s(p0[N][j], p1[N][j], p2[N][j]);
	 wire [AW-1:0] bc = csa_c(p0[N][j], p1[N][j], p2[N][j]);
	 wire [AW-1:0] cs = csa_s(bs, bc, p3[N][j]);
	 wire [AW-1:0] cc = csa_c(bs, bc, p3[N][j]);
	 wire [AW-1:0] cin = cs + cc;
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*32 +: 32] = {{(32-AW){cin[AW-1]}}, cin};
	 end
	 else begin : g_deskewn
	    reg [AW-1:0] dr [0:D-1];
	    integer k;
	    always @(posedge clk) begin
	       if (en) begin
		  dr[0] <= cin;
		  for (k = 1; k < D; k = k + 1) dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*32 +: 32] = {{(32-AW){dr[D-1][AW-1]}}, dr[D-1]};
	 end
      end
   endgenerate

   // ---------------- valid pipeline --------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

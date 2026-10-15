// TPU MXU tile (weight-stationary systolic array), cycle-accurate to the
// logikbench 'tpu' benchmark. Optimized datapath:
//   * accumulators narrowed to 19 bits (|C| <= 8*16384 = 131072)
//   * PE product built as explicit shift-add partial rows
//   * column partial sums carried in carry-save form (two vectors) down the
//     column: the PE critical path is row muxes + CSA levels only; a single
//     19-bit CPA per column converts to binary at the array bottom.
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

   // 3:2 carry-save compressor (result mod 2^AW; true values always fit)
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
   wire [AW-1:0] p_s [0:N][0:N-1];   // north -> south partial sums (save)
   wire [AW-1:0] p_c [0:N][0:N-1];   // north -> south partial sums (carry)

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
	 assign p_s[0][j] = {AW{1'b0}};
	 assign p_c[0][j] = {AW{1'b0}};
	 assign w_v[0][j] = w_data[j*DW +: DW];
      end

      // ---------------- PE mesh ------------------------------------------
      for (i = 0; i < N; i = i + 1) begin : g_row
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    reg [DW-1:0] wreg;
	    reg [DW-1:0] areg;
	    reg [AW-1:0] sreg;
	    reg [AW-1:0] creg;

	    wire [DW-1:0] a_raw = a_h[i][j];
	    // zero the product rows while the array is stalled (w_valid)
	    wire [DW-1:0] a_w  = a_raw & {DW{~ld}};
	    wire [AW-1:0] s_in = p_s[i][j];
	    wire [AW-1:0] c_in = p_c[i][j];

	    // Baugh-Wooley partial products: 8 x 8 AND matrix, sign corrections
	    // folded into row 7 and the constant K.  prod == sum_p (row_p << p) + K
	    // (mod 2^AW), rows are plain 8-bit vectors (no sign extension).
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
	    wire [AW-1:0] cr = 491776;		 // constant correction

	    // carry-save tree: eleven operands reduce to two words
	    wire [AW-1:0] s1a = csa_s(r0, r1, r2);
	    wire [AW-1:0] c1a = csa_c(r0, r1, r2);
	    wire [AW-1:0] s1b = csa_s(r3, r4, r5);
	    wire [AW-1:0] c1b = csa_c(r3, r4, r5);
	    wire [AW-1:0] s1c = csa_s(r6, r7, s_in);
	    wire [AW-1:0] c1c = csa_c(r6, r7, s_in);
	    wire [AW-1:0] s2a = csa_s(s1a, c1a, s1b);
	    wire [AW-1:0] c2a = csa_c(s1a, c1a, s1b);
	    wire [AW-1:0] s2b = csa_s(c1b, s1c, c1c);
	    wire [AW-1:0] c2b = csa_c(c1b, s1c, c1c);
	    wire [AW-1:0] s3a = csa_s(s2a, c2a, s2b);
	    wire [AW-1:0] c3a = csa_c(s2a, c2a, s2b);
	    wire [AW-1:0] s3b = csa_s(c2b, c_in, cr);
	    wire [AW-1:0] c3b = csa_c(c2b, c_in, cr);
	    wire [AW-1:0] s4a = csa_s(s3a, c3a, s3b);
	    wire [AW-1:0] c4a = csa_c(s3a, c3a, s3b);
	    wire [AW-1:0] s5a = csa_s(s4a, c4a, c3b);
	    wire [AW-1:0] c5a = csa_c(s4a, c4a, c3b);

	    always @(posedge clk) begin
	       if (ld) wreg <= w_v[i][j];
	       if (en) begin
		  areg <= a_raw;
		  sreg <= s5a;
		  creg <= c5a;
	       end
	    end

	    assign w_v[i+1][j] = wreg;
	    assign a_h[i][j+1] = areg;
	    assign p_s[i+1][j] = sreg;
	    assign p_c[i+1][j] = creg;
	 end
      end

      // ---------------- column bottom: CPA + output de-skew ---------------
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [AW-1:0] cin = p_s[N][j] + p_c[N][j];
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

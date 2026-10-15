// TPU MXU tile (weight-stationary systolic array), cycle-accurate to the
// logikbench 'tpu' benchmark. Optimized datapath:
//   * accumulators narrowed to 19 bits (|C| <= 8*16384 = 131072)
//   * PE product built as explicit shift-add partial rows and fused with the
//     accumulation through a carry-save tree: the PE critical path is
//     row muxes + 5 CSA levels + one 19-bit CPA (no multiplier CPA + adder).
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
   wire [AW-1:0] p_v [0:N][0:N-1];   // north -> south partial sums

   genvar i, j;
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
	       if (rst) begin
		  for (k = 0; k < i; k = k + 1) sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
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
	 assign p_v[0][j] = {AW{1'b0}};
	 assign w_v[0][j] = w_data[j*DW +: DW];
      end

      // ---------------- PE mesh ------------------------------------------
      for (i = 0; i < N; i = i + 1) begin : g_row
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    reg [DW-1:0] wreg;
	    reg [DW-1:0] areg;
	    reg [AW-1:0] preg;

	    wire [DW-1:0] a_w  = a_h[i][j];
	    wire [AW-1:0] p_in = p_v[i][j];

	    // shift-add partial rows of a_w * wreg (two's complement)
	    wire [AW-1:0] aa = {{(AW-DW){a_w[DW-1]}}, a_w};
	    wire [AW-1:0] r0 = wreg[0] ? aa       : {AW{1'b0}};
	    wire [AW-1:0] r1 = wreg[1] ? (aa << 1) : {AW{1'b0}};
	    wire [AW-1:0] r2 = wreg[2] ? (aa << 2) : {AW{1'b0}};
	    wire [AW-1:0] r3 = wreg[3] ? (aa << 3) : {AW{1'b0}};
	    wire [AW-1:0] r4 = wreg[4] ? (aa << 4) : {AW{1'b0}};
	    wire [AW-1:0] r5 = wreg[5] ? (aa << 5) : {AW{1'b0}};
	    wire [AW-1:0] r6 = wreg[6] ? (aa << 6) : {AW{1'b0}};
	    wire [AW-1:0] r7 = wreg[7] ? ~(aa << 7) : {AW{1'b0}};
	    wire [AW-1:0] cr = {{(AW-1){1'b0}}, wreg[7]};

	    // carry-save tree: 10 operands -> 2
	    wire [AW-1:0] s1a = csa_s(r0, r1, r2);
	    wire [AW-1:0] c1a = csa_c(r0, r1, r2);
	    wire [AW-1:0] s1b = csa_s(r3, r4, r5);
	    wire [AW-1:0] c1b = csa_c(r3, r4, r5);
	    wire [AW-1:0] s1c = csa_s(r6, r7, cr);
	    wire [AW-1:0] c1c = csa_c(r6, r7, cr);
	    wire [AW-1:0] s2a = csa_s(s1a, c1a, s1b);
	    wire [AW-1:0] c2a = csa_c(s1a, c1a, s1b);
	    wire [AW-1:0] s2b = csa_s(c1b, s1c, c1c);
	    wire [AW-1:0] c2b = csa_c(c1b, s1c, c1c);
	    wire [AW-1:0] s3a = csa_s(s2a, c2a, s2b);
	    wire [AW-1:0] c3a = csa_c(s2a, c2a, s2b);
	    wire [AW-1:0] s4a = csa_s(s3a, c3a, c2b);
	    wire [AW-1:0] c4a = csa_c(s3a, c3a, c2b);
	    wire [AW-1:0] s5a = csa_s(s4a, c4a, p_in);
	    wire [AW-1:0] c5a = csa_c(s4a, c4a, p_in);

	    wire [AW-1:0] psum_next = s5a + c5a;

	    always @(posedge clk) begin
	       if (rst) begin
		  wreg  <= {DW{1'b0}};
		  areg  <= {DW{1'b0}};
		  preg  <= {AW{1'b0}};
	       end
	       else begin
		  if (ld) wreg <= w_v[i][j];
		  if (en) begin
		     areg <= a_w;
		     preg <= psum_next;
		  end
	       end
	    end

	    assign w_v[i+1][j] = wreg;
	    assign a_h[i][j+1] = areg;
	    assign p_v[i+1][j] = preg;
	 end
      end

      // ---------------- output de-skew ------------------------------------
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [AW-1:0] cin = p_v[N][j];
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*32 +: 32] = {{(32-AW){cin[AW-1]}}, cin};
	 end
	 else begin : g_deskewn
	    reg [AW-1:0] dr [0:D-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < D; k = k + 1) dr[k] <= {AW{1'b0}};
	       end
	       else if (en) begin
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

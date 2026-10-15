// Optimized TPU MXU, v9: carry-save (redundant) accumulation.
//
//  Each PE accumulates into a redundant pair (cs, ss) with cs + ss = partial
//  sum, using a 3:2 compressor per bit:
//      cs' = cs ^ ss ^ prod
//      ss' = maj(cs, ss, prod) << 1
//  No carry propagation inside the array: the accumulator is ~2 gate levels
//  after the multiplier. The partial-sum bounds from v6 apply to cs+ss, and
//  loose induction bounds |cs|,|ss| <= 8*2^15 = 2^18, so CW=21-bit vectors are
//  safe everywhere. A single carry-propagate adder per column converts the
//  bottom-row pair to the final 21-bit sum; deskew + sign extension follow.
//  All timing/latency/stall behavior is identical (same register stages).

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter CW = 21)	 // redundant-vector width (signed)
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [CW-1:0]	  cs_north,   // carry-save vector in
    input [CW-1:0]	  ss_north,   // save-sum vector in
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output reg [CW-1:0]	  cs_south,   // carry-save vector out
    output reg [CW-1:0]	  ss_south    // save-sum vector out
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [CW-1:0] px = {{(CW-2*DW){prod[2*DW-1]}}, prod};

   wire [CW-1:0] nxt_cs = cs_north ^ ss_north ^ px;
   wire [CW-1:0] nxt_ss = ((cs_north & ss_north) |
                           (cs_north & px) |
                           (ss_north & px)) << 1;

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east   <= a_west;
         cs_south <= nxt_cs;
         ss_south <= nxt_ss;
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter CW = 21,	 // redundant-vector width
    parameter AW = 21)	 // output width per column (cpa result)
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*AW-1:0]	c_bot	// column results, one per column
    );

   wire [DW-1:0]  a_h   [0:N-1][0:N];
   wire [CW-1:0]  cs_v  [0:N][0:N-1];
   wire [CW-1:0]  ss_v  [0:N][0:N-1];
   wire [DW-1:0]  w_v   [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign cs_v[0][j] = {CW{1'b0}};
         assign ss_v[0][j] = {CW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
         // final carry-propagate add per column (pure combinational here)
         assign c_bot[j*AW +: AW] = cs_v[N][j] + ss_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .CW(CW)) u_pe
              (.clk        (clk),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .cs_north   (cs_v[i][j]),
               .ss_north   (ss_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .cs_south   (cs_v[i+1][j]),
               .ss_south   (ss_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32, // result width (signed), at the port
    parameter PW = 21)	 // bottom-row sum width
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

   // Systolic array with carry-save accumulation.
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .CW(PW), .AW(PW)) u_array
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

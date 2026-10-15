// Optimized TPU MXU, v11: tapered-width carry-save accumulation.
//
//  Same redundant (cs,ss) accumulation as v9, but with per-row tapered vector
//  widths matching the partial-sum magnitude bound at each row:
//    row i partial sum magnitude <= (i+1)*2^14, so the redundant pair only
//    needs to hold the sum mod 2^CW with CW bits of signed range:
//      row 0 -> 16, row 1 -> 17, rows 2-3 -> 18, rows 4-7 -> 19.
//  Each PE keeps narrow (CW-bit) cs/ss registers and outputs full 19-bit
//  sign-extended vectors (pure wiring), keeping inter-row nets uniform.
//  Final per-column CPA at the bottom: cs + ss (19-bit), then deskew.
//  Modular arithmetic keeps the invariant exact; tapering is safe because
//  widths are non-decreasing going south.

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter AW = 19,	 // inter-row (port) vector width
    parameter CW = 19)	 // internal (register) vector width
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [AW-1:0]	  cs_north,   // carry-save vector in
    input [AW-1:0]	  ss_north,   // save-sum vector in
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output [AW-1:0]	  cs_south,   // carry-save vector out
    output [AW-1:0]	  ss_south    // save-sum vector out
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [CW-1:0] px = {{(CW-2*DW){prod[2*DW-1]}}, prod};

   wire [CW-1:0] cn = cs_north[CW-1:0];
   wire [CW-1:0] sn = ss_north[CW-1:0];
   wire [CW-1:0] nxt_cs = cn ^ sn ^ px;
   wire [CW-1:0] nxt_ss = ((cn & sn) | (cn & px) | (sn & px)) << 1;

   reg [CW-1:0] cs, ss;
   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east <= a_west;
         cs     <= nxt_cs;
         ss     <= nxt_ss;
      end
   end

   assign cs_south = {{(AW-CW){cs[CW-1]}}, cs};
   assign ss_south = {{(AW-CW){ss[CW-1]}}, ss};

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter AW = 19)	 // inter-row vector width
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*AW-1:0]	c_bot	// column results, one per column
    );

   function integer cwf;
      input integer i;
      begin
         case (i)
            0:      cwf = 16;
            1:      cwf = 17;
            2,3:    cwf = 18;
            default: cwf = 19;
         endcase
      end
   endfunction

   wire [DW-1:0]  a_h   [0:N-1][0:N];
   wire [AW-1:0]  cs_v  [0:N][0:N-1];
   wire [AW-1:0]  ss_v  [0:N][0:N-1];
   wire [DW-1:0]  w_v   [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign cs_v[0][j] = {AW{1'b0}};
         assign ss_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
         // final carry-propagate add per column; 19-bit result is exact
         assign c_bot[j*AW +: AW] = cs_v[N][j] + ss_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW), .CW(cwf(i))) u_pe
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
    parameter PW = 19)	 // bottom-row sum width
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

   // Systolic array with tapered carry-save accumulation.
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

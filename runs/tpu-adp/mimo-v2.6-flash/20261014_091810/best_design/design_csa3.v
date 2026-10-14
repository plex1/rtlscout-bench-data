//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- variant:
// radix-4 Booth weights + THREE-LEVEL carry-save accumulation (full fold).
//
// All six addends (redundant accumulator c, s and the four shifted Booth
// product terms t0, t1<<2, t2<<4, t3<<6) are reduced by pure full-adder
// levels -- no carry-propagate addition remains in the multiply-accumulate
// cone:
//     L1 (parallel): CSA(c, s, t0) -> (x1, y1)
//                    CSA(t1', t2', t3') -> (x2, y2)
//     L2:            CSA(x1, y1, x2) -> (x3, y3)
//     L3:            CSA(y2, x3, y3) -> (s4, c4)
// Value identity: s4 + c4 = c + s + t0 + t1' + t2' + t3' (exact).
// Component widths grow by 3 bits per row (exact, no wrap).
// De-skew stores the exact 19-bit column result; reset only on the valid
// pipeline; term timing u = i+j (exact under mid-drain reloads).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter EW = 12,	 // booth-encoded weight width
    parameter CWI = 16,	 // carry-in width
    parameter SWI = 16,	 // sum-in width
    parameter CWO = 16,	 // carry-out width
    parameter SWO = 16)	 // sum-out width
   (
    input		      clk,
    input		      rst,	      // (used only by the valid pipe)
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [EW-1:0]	      w_north,    // encoded weight in from PE above
    input [DW-1:0]	      a_west,     // activation in from west
    input [CWI-1:0]	      c_north,    // carry partial sum in
    input [SWI-1:0]	      s_north,    // sum partial sum in
    output [EW-1:0]	      w_south,    // encoded weight out to PE below
    output reg [DW-1:0]	      a_east,     // activation out to east
    output reg [CWO-1:0]      c_south,    // carry partial sum out
    output reg [SWO-1:0]      s_south     // sum partial sum out
    );

   // operand widths through the three CSA levels
   localparam M0 = (CWI > SWI) ? ((CWI > 2*DW) ? CWI : 2*DW)
			       : ((SWI > 2*DW) ? SWI : 2*DW);
   localparam M1 = M0 + 1;
   localparam M2 = M0 + 2;

   reg [EW-1:0] w;
   assign w_south = w;

   // signed value of one Booth digit d in {0, +1, +2, -1, -2} times a
   function signed [15:0] dmul;
      input [2:0]	       d;
      input [DW-1:0]	       a;
      reg signed [15:0]       as;
      begin
	 as = {{(16-DW){a[DW-1]}}, a};
	 case (d)
	   3'b000:  dmul = 16'sd0;
	   3'b001:  dmul = as;
	   3'b010:  dmul = {as[14:0], 1'b0};
	   3'b111:  dmul = -as;
	   3'b110:  dmul = -{as[14:0], 1'b0};
	   default: dmul = 16'sd0;
	 endcase
      end
   endfunction

   wire signed [15:0] t0 = dmul(w[2:0],   a_west);
   wire signed [15:0] t1 = dmul(w[5:3],   a_west);
   wire signed [15:0] t2 = dmul(w[8:6],   a_west);
   wire signed [15:0] t3 = dmul(w[11:9],  a_west);

   // shifted terms (each still exact in 16 signed bits: |t|<=256)
   wire [15:0]	      t1s = {t1[13:0], 2'b00};
   wire [15:0]	      t2s = {t2[11:0], 4'b0000};
   wire [15:0]	      t3s = {t3[9:0],  6'b000000};

   // ---- level 1 (two CSAs in parallel) -----------------------------------
   wire signed [M0-1:0] ce = $signed(c_north);
   wire signed [M0-1:0] se = $signed(s_north);
   wire signed [M0-1:0] p0 = $signed(t0);
   wire [M0-1:0]	       x1 = ce ^ se ^ p0;
   wire [M0:0]	       y1 = {((ce & se) | (se & p0) | (ce & p0)), 1'b0};

   wire signed [M0-1:0] q1 = $signed(t1s);
   wire signed [M0-1:0] q2 = $signed(t2s);
   wire signed [M0-1:0] q3 = $signed(t3s);
   wire [M0-1:0]	       x2 = q1 ^ q2 ^ q3;
   wire [M0:0]	       y2 = {((q1 & q2) | (q2 & q3) | (q1 & q3)), 1'b0};

   // ---- level 2 ----------------------------------------------------------
   wire signed [M1-1:0] x1e = $signed(x1);
   wire signed [M1-1:0] y1e = $signed(y1);
   wire signed [M1-1:0] x2e = $signed(x2);
   wire [M1-1:0]	       x3 = x1e ^ y1e ^ x2e;
   wire [M1:0]	       y3 = {((x1e & y1e) | (y1e & x2e) | (x1e & x2e)), 1'b0};

   // ---- level 3 ----------------------------------------------------------
   wire signed [M2-1:0] y2e = $signed(y2);
   wire signed [M2-1:0] x3e = $signed(x3);
   wire signed [M2-1:0] y3e = $signed(y3);
   wire [M2-1:0]	       s4 = y2e ^ x3e ^ y3e;
   wire [M2:0]	       c4 = {((y2e & x3e) | (x3e & y3e) | (y2e & y3e)), 1'b0};

   always @(posedge clk) begin
      if (ld)
	w <= w_north;
      if (en) begin
	 a_east  <= a_west;
	 s_south <= s4[SWO-1:0];
	 c_south <= c4[CWO-1:0];
      end
   end

endmodule

//#############################################################################

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input		clk,
    input		rst,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*ACCW-1:0]	c_bot	// column results, one per column
    );

   localparam NDIG = DW / 2;
   localparam EW   = 3 * NDIG;

   // width of the carry register stored at row i (iterative, no recursion)
   function integer wcw;
      input integer i;
      integer	       k, cw, sw, m;
      begin
	 cw = 2*DW + 3;	// row 0: M0 = 2*DW -> CWO = M0 + 3
	 sw = 2*DW + 2;		// row 0: SWO = M0 + 2
	 for (k = 1; k <= i; k = k + 1) begin
	    m = (cw > sw) ? cw : sw;
	    if (m < 2*DW) m = 2*DW;
	    sw = m + 2;
	    cw = m + 3;
	 end
	 wcw = cw;
      end
   endfunction

   function integer wsw;
      input integer i;
      integer	       k, cw, sw, m;
      begin
	 cw = 2*DW + 3;
	 sw = 2*DW + 2;
	 for (k = 1; k <= i; k = k + 1) begin
	    m = (cw > sw) ? cw : sw;
	    if (m < 2*DW) m = 2*DW;
	    sw = m + 2;
	    cw = m + 3;
	 end
	 wsw = sw;
      end
   endfunction

   // slot width must hold the widest component (row N-1's carry, up to 40
   // bits for N=8) -- the slots are wires, so extra bits cost nothing
   localparam SLOTW = wcw(N - 1);

   // radix-4 Booth encoding of one signed weight byte
   function [EW-1:0] booth_enc;
      input [7:0]  w;
      integer      k, idx;
      reg	       b2, b1, b0;
      reg [2:0]      d;
      begin
	 for (k = 0; k < NDIG; k = k + 1) begin
	    b2 = w[2*k+1];
	    b1 = w[2*k];
	    idx = (k == 0) ? 0 : (2*k - 1);
	    b0 = (k == 0) ? 1'b0 : w[idx];
	    case ({b2, b1, b0})
	      3'b000:  d = 3'b000; //  0
	      3'b001:  d = 3'b001; // +1
	      3'b010:  d = 3'b001; // +1
	      3'b011:  d = 3'b010; // +2
	      3'b100:  d = 3'b110; // -2
	      3'b101:  d = 3'b111; // -1
	      3'b110:  d = 3'b111; // -1
	      3'b111:  d = 3'b000; //  0
	      default: d = 3'b000;
	    endcase
	    booth_enc[3*k +: 3] = d;
	 end
      end
   endfunction

   wire [((N+1)*N*SLOTW)-1:0] c_v;
   wire [((N+1)*N*SLOTW)-1:0] s_v;
   wire [DW-1:0]	       a_h    [0:N-1][0:N];
   wire [EW-1:0]	       w_v    [0:N][0:N-1];

   assign c_v[0 +: N*SLOTW] = {(N*SLOTW){1'b0}};
   assign s_v[0 +: N*SLOTW] = {(N*SLOTW){1'b0}};

   genvar		       jb;
   generate
      for (jb = 0; jb < N; jb = jb + 1) begin : g_bot
	 wire [SLOTW-1:0] cv = c_v[(N*N+jb)*SLOTW +: SLOTW];
	 wire [SLOTW-1:0] sv = s_v[(N*N+jb)*SLOTW +: SLOTW];
	 wire [SLOTW-1:0] sum = cv + sv;   // exact: |V| < 2^17 << 2^SLOTW
	 assign c_bot[jb*ACCW +: ACCW] = sum[ACCW-1:0];
      end
   endgenerate

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign w_v[0][j] = booth_enc(w_top[j*DW +: DW]);
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 localparam CWI = (i == 0) ? 2*DW : wcw(i - 1);
	 localparam SWI = (i == 0) ? 2*DW : wsw(i - 1);
	 localparam CWO = wcw(i);
	 localparam SWO = wsw(i);
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    wire [CWI-1:0]  c_north = c_v[(i*N+j)*SLOTW +: CWI];
	    wire [SWI-1:0]  s_north = s_v[(i*N+j)*SLOTW +: SWI];
	    wire [CWO-1:0]  c_south;
	    wire [SWO-1:0]  s_south;
	    tpu_pe #(.DW(DW), .EW(EW),
		     .CWI(CWI), .SWI(SWI), .CWO(CWO), .SWO(SWO)) u_pe
	      (.clk        (clk),
	       .rst        (rst),
	       .ld         (ld),
	       .en         (en),
	       .w_north    (w_v[i][j]),
	       .a_west     (a_h[i][j]),
	       .c_north    (c_north),
	       .s_north    (s_north),
	       .w_south    (w_v[i+1][j]),
	       .a_east     (a_h[i][j+1]),
	       .c_south    (c_south),
	       .s_south    (s_south));
	    assign c_v[((i+1)*N+j)*SLOTW +: SLOTW] = {{(SLOTW-CWO){1'b0}}, c_south};
	    assign s_v[((i+1)*N+j)*SLOTW +: SLOTW] = {{(SLOTW-SWO){1'b0}}, s_south};
	 end
      end
   endgenerate

endmodule

//#############################################################################

module tpu
  #(parameter N = 8,	 // array dimension (N x N)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 32) // result / accumulator width (signed)
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

   function integer pwidth;
      input integer k;
      integer	       w, pmax, pmin;
      begin
	 pmax = (1 << (2*DW - 2));
	 pmin = -((1 << (DW-1)) * ((1 << (DW-1)) - 1));
	 w = 2*DW;
	 while ((k*pmax > ((1 << (w-1)) - 1)) || (k*pmin < -((1 << (w-1)))))
	   w = w + 1;
	 pwidth = w;
      end
   endfunction

   localparam RESW = pwidth(N);

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar	  i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
	 wire [DW-1:0] din =
		       a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
	 if (i == 0) begin : g_skew0
	    assign a_left[i*DW +: DW] = din;
	 end
	 else begin : g_skewn
	    reg [DW-1:0] sr [0:i-1];
	    integer	       k;
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

   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
	 localparam D = N - 1 - j;
	 wire [RESW-1:0] cin = c_bot[j*ACCW +: RESW];
	 if (D == 0) begin : g_deskew0
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-RESW){cin[RESW-1]}}, cin};
	 end
	 else begin : g_deskewn
	    reg [RESW-1:0] dr [0:D-1];
	    integer	    k;
	    always @(posedge clk) begin
	       if (en) begin
		  dr[0] <= cin;
		  for (k = 1; k < D; k = k + 1)
		    dr[k] <= dr[k-1];
	       end
	    end
	    assign c_data[j*ACCW +: ACCW] =
		    {{(ACCW-RESW){dr[D-1][RESW-1]}}, dr[D-1]};
	 end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

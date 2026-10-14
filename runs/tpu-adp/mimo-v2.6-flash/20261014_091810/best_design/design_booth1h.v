//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) -- variant:
// ONE-HOT radix-4 Booth digit storage (5 bits per digit, digit in
// {0,+1,+2,-1,-2}), encoded once at the load input.  The multiply cone then
// selects term k with a plain AND-OR (no binary decode).  Everything else
// matches design_booth.v: exact-width partial sums / de-skew, reset only on
// the valid pipeline, anti-diagonal term timing u = i+j.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter EW = 20,	 // one-hot booth digit width (4 digits x 5)
    parameter PIN = 16,	 // psum_north width (signed)
    parameter POUT = 16) // psum_south width (signed)
   (
    input		      clk,
    input		      rst,	      // (used only by the valid pipe)
    input		      ld,	      // weight-load shift enable
    input		      en,	      // compute enable
    input [EW-1:0]	      w_north,    // encoded weight in from PE above
    input [DW-1:0]	      a_west,     // activation in from west
    input [PIN-1:0]	      psum_north, // partial sum in from above
    output [EW-1:0]	      w_south,    // encoded weight out to PE below
    output reg [DW-1:0]	      a_east,     // activation out to east
    output reg [POUT-1:0]     psum_south  // partial sum out to below
    );

   reg [EW-1:0] w;
   assign w_south = w;

   // term k = digit_k * a, selected by the one-hot digit line
   // one-hot order per digit: {p2, p1, z, m1, m2}  (values 0,+1,+2,-1,-2)
   function signed [15:0] dmul;
      input [4:0]	       oh;
      input [DW-1:0]	       a;
      reg signed [15:0]       as;
      begin
	 as = {{(16-DW){a[DW-1]}}, a};
	 dmul = (oh[4] ? {as[14:0], 1'b0} : 16'sd0) |
		(oh[3] ? as              : 16'sd0) |
		(oh[1] ? -as             : 16'sd0) |
		(oh[0] ? -{as[14:0], 1'b0} : 16'sd0);
      end
   endfunction

   wire signed [15:0] t0 = dmul(w[4:0],    a_west);
   wire signed [15:0] t1 = dmul(w[9:5],    a_west);
   wire signed [15:0] t2 = dmul(w[14:10],  a_west);
   wire signed [15:0] t3 = dmul(w[19:15],  a_west);

   wire [15:0]	      lo = t0 + {t1[13:0], 2'b00};
   wire [15:0]	      hi = {t2[11:0], 4'b0000} + {t3[9:0], 6'b000000};
   wire [2*DW-1:0]    prod = lo + hi;

   always @(posedge clk) begin
      if (ld)
	w <= w_north;
      if (en) begin
	 a_east     <= a_west;
	 psum_south <= $signed(psum_north) + $signed(prod);
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

   localparam NDIG = DW / 2;	   // booth digits per byte
   localparam EW   = 5 * NDIG;	   // one-hot digit width (20 for DW=8)

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

   // one-hot radix-4 Booth encoding: per digit order {p2,p1,z,m1,m2}
   function [EW-1:0] booth_enc;
      input [7:0]  w;
      integer      k, idx;
      reg	       b2, b1, b0;
      reg [4:0]      oh;
      begin
	 for (k = 0; k < NDIG; k = k + 1) begin
	    b2 = w[2*k+1];
	    b1 = w[2*k];
	    idx = (k == 0) ? 0 : (2*k - 1);
	    b0 = (k == 0) ? 1'b0 : w[idx];
	    case ({b2, b1, b0})
	      3'b000:  oh = 5'b00100; //  0
	      3'b001:  oh = 5'b01000; // +1
	      3'b010:  oh = 5'b01000; // +1
	      3'b011:  oh = 5'b10000; // +2
	      3'b100:  oh = 5'b00001; // -2
	      3'b101:  oh = 5'b00010; // -1
	      3'b110:  oh = 5'b00010; // -1
	      3'b111:  oh = 5'b00100; //  0
	      default: oh = 5'b00100;
	    endcase
	    booth_enc[5*k +: 5] = oh;
	 end
      end
   endfunction

   wire [((N+1)*N*ACCW)-1:0] psum_v;
   wire [DW-1:0]	       a_h    [0:N-1][0:N];
   wire [EW-1:0]	       w_v    [0:N][0:N-1];

   assign psum_v[0 +: N*ACCW] = {(N*ACCW){1'b0}};
   assign c_bot = psum_v[(N*N*ACCW) +: N*ACCW];

   genvar		       i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
	 assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
	 assign w_v[0][j] = booth_enc(w_top[j*DW +: DW]);
      end
      for (i = 0; i < N; i = i + 1) begin : g_row
	 localparam PWI = pwidth(i);
	 localparam PWO = pwidth(i + 1);
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    wire [PWI-1:0]  psum_north = psum_v[(i*N+j)*ACCW +: PWI];
	    wire [PWO-1:0]  psum_south;
	    tpu_pe #(.DW(DW), .EW(EW), .PIN(PWI), .POUT(PWO)) u_pe
	      (.clk        (clk),
	       .rst        (rst),
	       .ld         (ld),
	       .en         (en),
	       .w_north    (w_v[i][j]),
	       .a_west     (a_h[i][j]),
	       .psum_north (psum_north),
	       .w_south    (w_v[i+1][j]),
	       .a_east     (a_h[i][j+1]),
	       .psum_south (psum_south));
	    assign psum_v[((i+1)*N+j)*ACCW +: ACCW] = {{(ACCW-PWO){1'b0}}, psum_south};
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

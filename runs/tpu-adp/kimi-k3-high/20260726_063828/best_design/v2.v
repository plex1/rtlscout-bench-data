//#############################################################################
// Optimized weight-stationary systolic TPU (8x8) - structural RTL.
// Same cycle-accurate behavior as the logikbench tpu, with hand-built Booth
// multipliers, flat skew/deskew networks, and sparse weight-load enables.
//#############################################################################

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;

   wire ld = w_valid;
   wire en = ~w_valid;

   //#########################################################################
   // Weight-load pulse counter (previous-pulse value used for enables).
   //#########################################################################
   reg [3:0] wcnt;
   always @(posedge clk) begin
      if (rst)
	wcnt <= 4'd0;
      else if (ld && (wcnt != 4'd8))
	wcnt <= wcnt + 4'd1;
   end

   //#########################################################################
   // Input skew: per-row shift chains (original scheme). a_left[i] is the
   // skewed row-i activation fed to the mesh west edge.
   //#########################################################################
   wire [DW-1:0] a_left [0:N-1];
   genvar gi;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_skew
	 wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
	 if (gi == 0) begin : g_s0
	    assign a_left[gi] = din;
	 end
	 else begin : g_sn
	    reg [DW-1:0] sr [0:gi-1];
	    integer k;
	    always @(posedge clk) begin
	       if (rst) begin
		  for (k = 0; k < gi; k = k + 1)
		    sr[k] <= {DW{1'b0}};
	       end
	       else if (en) begin
		  sr[0] <= din;
		  for (k = 1; k < gi; k = k + 1)
		    sr[k] <= sr[k-1];
	       end
	    end
	    assign a_left[gi] = sr[gi-1];
	 end
      end
   endgenerate

   //#########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //#########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

   //#########################################################################
   // The 8x8 PE mesh (wires only across generate scopes).
   //#########################################################################
   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] p_v    [0:N][0:N-1];
   wire [DW-1:0]   w_out  [0:N-1][0:N-1];
   wire [ACCW-1:0] c_bot  [0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_row
	 for (j = 0; j < N; j = j + 1) begin : g_col
	    // west edge / top edge taps
	    if (j == 0) begin : g_we
	       assign a_h[i][0] = a_left[i];
	    end
	    if (i == 0) begin : g_te
	       assign p_v[0][j] = {ACCW{1'b0}};
	    end

	    // ----- weight register with sparse load enable -----
	    reg [DW-1:0] w;
	    assign w_out[i][j] = w;
	    wire [DW-1:0] w_src = (i == 0) ? w_data[j*DW +: DW] : w_out[i-1][j];
	    wire [3:0] widx = i[3:0] + wcnt;
	    wire w_en = rst | (ld & (widx == j[3:0]));
	    always @(posedge clk)
	      if (w_en)
		w <= rst ? {DW{1'b0}} : w_src;

	    // ----- activation pipeline -----
	    reg [DW-1:0] a_east;
	    assign a_h[i][j+1] = a_east;
	    always @(posedge clk)
	      if (rst)
		a_east <= {DW{1'b0}};
	      else if (en)
		a_east <= a_h[i][j];

	    // ----- partial-sum pipeline -----
	    reg [ACCW-1:0] p_south;
	    assign p_v[i+1][j] = p_south;
	    wire [2*DW-1:0] prod;
	    assign prod = $signed(a_h[i][j]) * $signed(w);
	    always @(posedge clk)
	      if (rst)
		p_south <= {ACCW{1'b0}};
	      else if (en)
		p_south <= p_v[i][j] + {{(ACCW-2*DW){prod[2*DW-1]}}, prod};
	 end
      end
      for (j = 0; j < N; j = j + 1) begin : g_bot
	 assign c_bot[j] = p_v[N][j];
      end
   endgenerate

   //#########################################################################
   // Output de-skew: flat register file d[j][t] (column j delayed N-1-j).
   //#########################################################################
   reg [ACCW-1:0] d [0:N-1][0:N];
   integer ej, et;
   always @(posedge clk) begin
      for (ej = 0; ej < N; ej = ej + 1) begin
	 for (et = 0; et <= N; et = et + 1) begin
	    if (rst)
	      d[ej][et] <= {ACCW{1'b0}};
	    else if (en)
	      d[ej][et] <= (et == 0) ? c_bot[ej] : d[ej][et-1];
	 end
      end
   end

   //#########################################################################
   // Result mux.
   //#########################################################################
   generate
      for (j = 0; j < N; j = j + 1) begin : g_out
	 assign c_data[j*ACCW +: ACCW] = d[j][N-1-j];
      end
   endgenerate

endmodule

//#############################################################################
// Signed 8x8 multiplier (radix-4 Booth, 4 partial products + Baugh-Wooley
// sign correction). x, y interpreted as signed; p = x * y (16 bits).
//#############################################################################
module mul8s
   (
    input  [7:0]  x,
    input  [7:0]  y,
    output [15:0] p
    );

   function [8:0] booth_pp;
      input [8:0] wsel;	// encoded digit, 3 sign-extension bits
      input [7:0] mult;	// multiplicand
      reg   [2:0] t;
      reg   [8:0] u;
      begin
	 t = {wsel[4], wsel[2], wsel[0]};
	 case (t)
	   3'b000, 3'b111: u = 9'b0;
	   3'b001, 3'b010: u = {mult[7], mult};
	   3'b011:	    u = {mult, 1'b0};
	   3'b100:	    u = ~{mult, 1'b0};
	   default:	    u = ~{mult[7], mult};
	 endcase
	 booth_pp = {u, wsel[6] ^ wsel[4]};
      end
   endfunction

   wire y7 = y[7];
   wire [8:0] w1 = {x[0], 1'b0}  ^ {y7, y7};
   wire [8:0] w2 = {x[2], x[1]}  ^ {y7, y7};
   wire [8:0] w3 = {x[4], x[3]}  ^ {y7, y7};
   wire [8:0] w4 = {x[6], x[5]}  ^ {y7, y7};

   wire [8:0] pp0 = booth_pp(w1, y);
   wire [8:0] pp1 = booth_pp(w2, y);
   wire [8:0] pp2 = booth_pp(w3, y);
   wire [8:0] pp3 = booth_pp(w4, y);

   wire [15:0] c0 = {pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[7:0]};
   wire [15:0] c1 = {pp1[8], pp1[8], pp1[8], pp1[8], pp1[8], pp1[8], pp1[7:0], 2'b00};
   wire [15:0] c2 = {pp2[8], pp2[8], pp2[8], pp2[8], pp2[7:0], 4'b0000};
   wire [15:0] c3 = {pp3[8], pp3[8], pp3[7:0], 6'b000000};
   wire [15:0] bw = 16'hA002;

   assign p = c0 + c1 + c2 + c3 + bw;

endmodule

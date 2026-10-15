// Optimized TPU MXU, v4:
//  v2 base (PW=19 psums, no datapath resets) + the PE accumulation adder
//  written as an explicit Kogge-Stone parallel-prefix adder (5 levels for
//  19 bits) instead of relying on the tool's default ripple structure.

module ksadd19 (
    input  [18:0] A,
    input  [18:0] B,
    output [18:0] S
);
   // Kogge-Stone prefix: S = A + B, 19 bits, offsets 1,2,4,8,16
   wire [18:0] G0 = A & B;
   wire [18:0] P0 = A ^ B;
   genvar k;
   generate
      // level 1: offset 1
      wire [18:0] G1, P1;
      assign G1[0] = G0[0];  assign P1[0] = P0[0];
      for (k = 1; k < 19; k = k + 1) begin : gl1
         assign G1[k] = G0[k] | (P0[k] & G0[k-1]);
         assign P1[k] = P0[k] & P0[k-1];
      end
      // level 2: offset 2
      wire [18:0] G2, P2;
      assign G2[1:0] = G1[1:0];  assign P2[1:0] = P1[1:0];
      for (k = 2; k < 19; k = k + 1) begin : gl2
         assign G2[k] = G1[k] | (P1[k] & G1[k-2]);
         assign P2[k] = P1[k] & P1[k-2];
      end
      // level 3: offset 4
      wire [18:0] G3, P3;
      assign G3[3:0] = G2[3:0];  assign P3[3:0] = P2[3:0];
      for (k = 4; k < 19; k = k + 1) begin : gl3
         assign G3[k] = G2[k] | (P2[k] & G2[k-4]);
         assign P3[k] = P2[k] & P2[k-4];
      end
      // level 4: offset 8
      wire [18:0] G4, P4;
      assign G4[7:0] = G3[7:0];  assign P4[7:0] = P3[7:0];
      for (k = 8; k < 19; k = k + 1) begin : gl4
         assign G4[k] = G3[k] | (P3[k] & G3[k-8]);
         assign P4[k] = P3[k] & P3[k-8];
      end
      // level 5: offset 16
      wire [18:0] G5, P5;
      assign G5[15:0] = G4[15:0];  assign P5[15:0] = P4[15:0];
      for (k = 16; k < 19; k = k + 1) begin : gl5
         assign G5[k] = G4[k] | (P4[k] & G4[k-16]);
         assign P5[k] = P4[k] & P4[k-16];
      end
   endgenerate
   assign S = P0 ^ {G5[17:0], 1'b0};
endmodule

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 19) // internal partial-sum width (signed)
   (
    input		  clk,
    input		  ld,	      // weight-load shift enable
    input		  en,	      // compute enable
    input [DW-1:0]	  w_north,    // weight in from PE above
    input [DW-1:0]	  a_west,     // activation in from west
    input [ACCW-1:0]	  psum_north, // partial sum in from above
    output [DW-1:0]	  w_south,    // stored weight out to PE below
    output reg [DW-1:0]	  a_east,     // activation out to east
    output reg [ACCW-1:0] psum_south  // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_x = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] nxt;
   ksadd19 u_add (.A(psum_north), .B(prod_x), .S(nxt));

   always @(posedge clk) begin
      if (ld)
        w <= w_north;
      if (en) begin
         a_east     <= a_west;
         psum_south <= nxt;
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,	 // array dimension (N x N PEs)
    parameter DW = 8,	 // operand width (signed)
    parameter ACCW = 19) // internal partial-sum width (signed)
   (
    input		clk,
    input		ld,	// weight-load shift enable
    input		en,	// compute enable
    input [N*DW-1:0]	w_top,	// weight inputs, one per column (load)
    input [N*DW-1:0]	a_left,	// activation inputs, one per row
    output [N*ACCW-1:0]	c_bot	// column results, one per column
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
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
    parameter PW = 19)	 // internal partial-sum width
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

   // Systolic array with narrow (PW-bit) partial sums.
   wire [N*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(PW)) u_array
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

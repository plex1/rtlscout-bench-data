// Carry-save systolic TPU. Each PE replaces its multiply-accumulate CPA with a
// 3:2 compressor: the partial sum travels south as a (S, C) carry-save pair.
// S_out = S ^ C ^ P ; C_out = MAJ(S, C, P)  =>  S + 2C is preserved mod 2^19.
// The exact 8-product sum fits in 19 signed bits, so mod-2^19 carry-save
// arithmetic is lossless. One final CPA per column (after the de-skew delay)
// reconstructs the 19-bit result, sign-extended to 32 at the pins.

module tpu_pe
  #(parameter DW = 8,
    parameter CW = 19)   // carry-save pair width
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [CW-1:0]	  s_north,
    input [CW-1:0]	  c_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [CW-1:0]	  s_south,
    output reg [CW-1:0]	  c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [CW-1:0]   pext = {{(CW-2*DW){prod[2*DW-1]}}, prod};
   wire [CW-1:0]   maj  = (s_north & c_north) |
			  (s_north & pext) | (c_north & pext);

   always @(posedge clk) begin
      if (rst) begin
	 w	  <= {DW{1'b0}};
	 a_east	  <= {DW{1'b0}};
	 s_south  <= {CW{1'b0}};
	 c_south  <= {CW{1'b0}};
      end
      else begin
	 if (ld) w <= w_north;
	 if (en) begin
	    a_east  <= a_west;
	    // 3:2 compressor: S' + C' = S + C + P (mod 2^CW)
	    s_south <= s_north ^ c_north ^ pext;
	    c_south <= {maj[CW-2:0], 1'b0};
	 end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter CW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*CW-1:0]	s_bot,
    output [N*CW-1:0]	c_bot
    );

   wire [DW-1:0] a_h    [0:N-1][0:N];
   wire [CW-1:0] s_v    [0:N][0:N-1];
   wire [CW-1:0] c_v    [0:N][0:N-1];
   wire [DW-1:0] w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {CW{1'b0}};   // fresh accumulation at the top
         assign c_v[0][j] = {CW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign s_bot[j*CW +: CW] = s_v[N][j];
         assign c_bot[j*CW +: CW] = c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .CW(CW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .s_north    (s_v[i][j]),
               .c_north    (c_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .s_south    (s_v[i+1][j]),
               .c_south    (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,
    parameter IACCW = 19)
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

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar	   j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_skew
         wire [DW-1:0] din =
		       a_valid ? a_data[j*DW +: DW] : {DW{1'b0}};
         if (j == 0) begin : g_skew0
            assign a_left[j*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:j-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < j; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < j; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[j*DW +: DW] = sr[j-1];
         end
      end
   endgenerate

   wire [N*IACCW-1:0] s_bot, c_bot;
   tpu_array #(.N(N), .DW(DW), .CW(IACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .s_bot(s_bot), .c_bot(c_bot));

   // Output de-skew: column j delayed N-1-j cycles; then one final CPA
   // reconstructs the value from the carry-save pair: V = S + 2*C (mod 2^19).
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [IACCW-1:0] sin = s_bot[j*IACCW +: IACCW];
         wire [IACCW-1:0] cin = c_bot[j*IACCW +: IACCW];
         if (D == 0) begin : g_deskew0
            wire [IACCW-1:0] v = sin + cin;
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-IACCW){v[IACCW-1]}}, v};
         end
         else begin : g_deskewn
            // D delay stages on the carry-save pair, then the final CPA
            reg [IACCW-1:0] dsr [0:D-1];
            reg [IACCW-1:0] dcr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1) begin
                     dsr[k] <= {IACCW{1'b0}};
                     dcr[k] <= {IACCW{1'b0}};
                  end
               end
               else if (en) begin
                  dsr[0] <= sin;
                  dcr[0] <= cin;
                  for (k = 1; k < D; k = k + 1) begin
                     dsr[k] <= dsr[k-1];
                     dcr[k] <= dcr[k-1];
                  end
               end
            end
            wire [IACCW-1:0] v = dsr[D-1] + dcr[D-1];
            assign c_data[j*ACCW +: ACCW] =
		   {{(ACCW-IACCW){v[IACCW-1]}}, v};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   integer kk;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

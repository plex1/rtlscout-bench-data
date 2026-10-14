//#############################################################################
// Weight-stationary systolic 8x8 matmul tile (TPU MXU), width-optimized v3.
//
// Same systolic structure as the original logikbench tpu (cycle-accurate,
// including mid-drain weight reloads), with narrowed arithmetic:
//
//  - Each signed 8x8 product lies in [-2^14, 2^14]. Row i of the array
//    accumulates exactly (i+1) products, so PE row i's partial sums lie in
//    [-(i+1)*2^14, (i+1)*2^14] and need only
//        PSW(i) = 16 + ceil(log2(i+1))  bits signed
//    (row0:16, row1:17, row2:17, row3:18, row4:18, row5:18, row6:18,
//     row7:19) instead of 32. The product itself is 16 bits, so upper rows
//    sign-extend it into the row's partial-sum width.
//  - De-skew registers use the narrow bottom-row width (19 bits).
//#############################################################################

module tpu_pe
  #(parameter DW = 8,	 // operand width (signed)
    parameter PSW = 19)  // partial-sum width (signed)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PSW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PSW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {PSW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= psum_north +
                          {{(PSW-2*DW){prod[2*DW-1]}}, prod};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*19-1:0]	c_bot	 // bottom-row partial sums (19 bits each)
    );

   // partial-sum width per row (row i accumulates i+1 products)
   function integer psw;
      input integer row;
      begin
	 psw = 16;
	 if (row >= 1) psw = 17;
	 if (row >= 3) psw = 18;
	 if (row >= 7) psw = 19;
      end
   endfunction

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [18:0]     psum_v [0:N][0:N-1]; // max width 19, slices use low bits
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = 19'd0;
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*19 +: 19] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            wire [psw(i)-1:0] sout;
            tpu_pe #(.DW(DW), .PSW(psw(i))) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_north (psum_v[i][j][psw(i)-1:0]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south (sout));
            // sign-extend this row's partial sum into the 19-bit net
            assign psum_v[i+1][j] = {{(19-psw(i)){sout[psw(i)-1]}}, sout};
         end
      end
   endgenerate

endmodule

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

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew: row i of the activation delayed i cycles (zero when a_valid low).
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
               if (rst) begin
                  for (k = 0; k < i; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   wire [N*19-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output de-skew: column j delayed N-1-j cycles (19-bit registers).
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [18:0] cin = c_bot[j*19 +: 19];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-19){cin[18]}}, cin};
         end
         else begin : g_deskewn
            reg [18:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= 19'd0;
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-19){dr[D-1][18]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// Optimized version 12: carry-save accumulation in PE
// Instead of psum_south = psum_north + prod (one big adder),
// use carry-save: store sum and carry separately, combine at the bottom.
// This reduces the critical path per PE from (multiply + add) to (multiply + CSA)
module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_n,    // sum in
    input [ACCW-1:0]	  pcarry_n,  // carry in
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_s,    // sum out
    output reg [ACCW-1:0] pcarry_s   // carry out
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};
   
   // Carry-save: psum_s = psum_n ^ prod_ext ^ pcarry_n
   //             pcarry_s = (psum_n & prod_ext) | (psum_n & pcarry_n) | (prod_ext & pcarry_n)
   // This is a 3:2 compressor (full adder)
   wire [ACCW-1:0] a = psum_n;
   wire [ACCW-1:0] b = prod_ext;
   wire [ACCW-1:0] c = pcarry_n;
   
   wire [ACCW-1:0] cs = a ^ b ^ c;                           // sum
   wire [ACCW-1:0] cc = (a & b) | (a & c) | (b & c);          // carry (unshifted, carry to next stage)
   // The carry should be shifted left by 1 in a proper carry-save, but since
   // we accumulate the carry into the next PE, we need to shift it.
   // Actually in a proper CSA accumulator chain, the carry from stage k 
   // feeds into stage k+1 shifted by 1. But the psum chain is spatial (PE to PE),
   // not temporal. So we can't easily do carry-save across PEs.
   // 
   // Let me instead use a simpler approach: just use the original add but
   // try a different accumulator width.

   always @(posedge clk or posedge rst) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_s     <= {ACCW{1'b0}};
         pcarry_s   <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_s     <= cs;
            pcarry_s   <= {cc[ACCW-2:0], 1'b0};  // shift carry left
         end
      end
   end

endmodule

// This approach has a fundamental problem: the carry from one PE needs to 
// be added in the next PE, but the carry-save representation needs the final
// adder at the end. Also the carry shift means the carry from PE(i,j) goes to
// PE(i+1,j) but shifted by 1 bit, which means PE(i+1,j) needs to handle an
// extra carry input. This doesn't work cleanly in a systolic array.
// 
// Let me abandon this and go back to the simpler approach.

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*ACCW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [ACCW-1:0] pcarry_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {ACCW{1'b0}};
         assign pcarry_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j] + pcarry_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_h[i][j]),
               .psum_n     (psum_v[i][j]),
               .pcarry_n   (pcarry_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_s     (psum_v[i+1][j]),
               .pcarry_s   (pcarry_v[i+1][j]));
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
   localparam PACCW = 19;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

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
            always @(posedge clk or posedge rst) begin
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

   wire [N*PACCW-1:0] c_bot_narrow;
   tpu_array #(.N(N), .DW(DW), .ACCW(PACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot_narrow));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] cin = c_bot_narrow[j*PACCW +: PACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = 
               {{(ACCW-PACCW){cin[PACCW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk or posedge rst) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PACCW-1:0] dout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] = 
               {{(ACCW-PACCW){dout[PACCW-1]}}, dout};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk or posedge rst) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

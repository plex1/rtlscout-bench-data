// PE with fast prefix adder for the accumulate
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
    input [ACCW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Prefix adder: compute sum = psum_north + prod_ext
   wire [ACCW-1:0] GI, PI, PT;
   wire [ACCW-1:0] GO, PO;

   assign GI = psum_north & prod_ext;
   assign PI = psum_north | prod_ext;
   assign PT = psum_north ^ prod_ext;

   // Sklansky prefix tree
   wire [ACCW-1:0] g0, p0, g1, p1, g2, p2, g3, p3, g4, p4;
   assign g0 = GI;
   assign p0 = PI;

   // Level 1: span 2
   genvar b;
   generate
      for (b = 0; b < ACCW; b = b + 1) begin : lvl1
         if (b >= 1) begin
            assign g1[b] = GI[b] | (PI[b] & GI[b-1]);
            assign p1[b] = PI[b] & PI[b-1];
         end else begin
            assign g1[b] = GI[b];
            assign p1[b] = PI[b];
         end
      end
      // Level 2: span 4
      for (b = 0; b < ACCW; b = b + 1) begin : lvl2
         if (b >= 3) begin
            assign g2[b] = g1[b] | (p1[b] & g1[b-2]);
            assign p2[b] = p1[b] & p1[b-2];
         end else if (b >= 1) begin
            assign g2[b] = g1[b] | (p1[b] & GI[b-1]);
            assign p2[b] = p1[b] & PI[b-1];
         end else begin
            assign g2[b] = GI[b];
            assign p2[b] = PI[b];
         end
      end
      // Level 3: span 8
      for (b = 0; b < ACCW; b = b + 1) begin : lvl3
         if (b >= 7) begin
            assign g3[b] = g2[b] | (p2[b] & g2[b-7]);
            assign p3[b] = p2[b] & p2[b-7];
         end else if (b >= 3) begin
            assign g3[b] = g2[b] | (p2[b] & g1[b-3]);
            assign p3[b] = p2[b] & p1[b-3];
         end else if (b >= 1) begin
            assign g3[b] = g2[b] | (p2[b] & GI[b-1]);
            assign p3[b] = p2[b] & PI[b-1];
         end else begin
            assign g3[b] = g2[b];
            assign p3[b] = p2[b];
         end
      end
      // Level 4: span 16 (for 19 bits)
      for (b = 0; b < ACCW; b = b + 1) begin : lvl4
         if (b >= 15) begin
            assign g4[b] = g3[b] | (p3[b] & g3[b-15]);
            assign p4[b] = p3[b] & p3[b-15];
         end else begin
            assign g4[b] = g3[b];
            assign p4[b] = p3[b];
         end
      end
   endgenerate

   wire [ACCW-1:0] sum = PT ^ {g4[ACCW-2:0], 1'b0};

   always @(posedge clk or posedge rst) begin
      if (rst) begin
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= sum;
         end
      end
   end

endmodule

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
               .rst        (rst),
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
  #(parameter N = 8,
    parameter DW = 8)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*32-1:0]	c_data
    );

   localparam LAT = 2*N - 1;
   localparam ACCW = 19;

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
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         wire [31:0] cin_ext = {{(32-ACCW){cin[ACCW-1]}}, cin};
         if (D == 0) begin : g_deskew0
            assign c_data[j*32 +: 32] = cin_ext;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [31:0] dr_ext = {{(32-ACCW){dr[D-1][ACCW-1]}}, dr[D-1]};
            assign c_data[j*32 +: 32] = dr_ext;
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

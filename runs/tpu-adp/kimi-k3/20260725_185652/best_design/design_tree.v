//#############################################################################
// Weight-stationary 8x8 matmul tile, cycle-compatible with logikbench tpu.
// Optimization: with the weight held stationary per column, the column dot
// product is computed as 8 parallel multiplies feeding a balanced adder tree
// (depth 3) instead of an 8-stage serial MAC chain. The column result is then
// re-registered per row to reproduce the systolic pipeline timing.
//#############################################################################

module tpu_tree
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		     clk,
    input		     rst,
    input		     ld,	      // weight-load shift enable
    input		     en,	      // compute enable
    input [N*DW-1:0]     w_top,	      // weight inputs, one per column
    input [N*DW-1:0]     a_left,	      // activation inputs, one per row
    output [N*ACCW-1:0]  c_bot	      // column results (registered per row)
    );

   genvar i, j;

   // ------------------------------------------------------------------
   // Weight registers: per-column vertical shift chain (load from top).
   // w[i][j] holds B[i][j] after N load cycles.
   // ------------------------------------------------------------------
   wire [DW-1:0] w [0:N][0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_wedge
         assign w[0][j] = w_top[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_wrow
         for (j = 0; j < N; j = j + 1) begin : g_wcol
            reg [DW-1:0] wr;
            always @(posedge clk) begin
               if (rst)
                 wr <= {DW{1'b0}};
               else if (ld)
                 wr <= w[i][j];
            end
            assign w[i+1][j] = wr;
         end
      end
   endgenerate

   // ------------------------------------------------------------------
   // Column datapath: 8 parallel signed multiplies on the west-edge
   // (skewed) activation + balanced adder tree. The column-j dot product is
   // registered once (colr), then re-registered (N-1)+j times to reproduce
   // the systolic bottom-row timing so the output de-skew lines up.
   // ------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         // partial products (combinational) on the west-edge activation
         wire [2*DW-1:0] pp [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_pp
            assign pp[i] = $signed(a_left[i*DW +: DW]) * $signed(w[i+1][j]);
         end
         // balanced adder tree over 8 products, in ACCW precision
         wire [ACCW-1:0] s [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_sext
            assign s[i] = {{(ACCW-2*DW){pp[i][2*DW-1]}}, pp[i]};
         end
         // level 1: 4 adds
         wire [ACCW-1:0] l1 [0:3];
         for (i = 0; i < 4; i = i + 1) begin : g_l1
            assign l1[i] = s[2*i] + s[2*i+1];
         end
         // level 2: 2 adds
         wire [ACCW-1:0] l2 [0:1];
         assign l2[0] = l1[0] + l1[1];
         assign l2[1] = l1[2] + l1[3];
         // level 3: final
         wire [ACCW-1:0] dot = l2[0] + l2[1];

         // register the dot product (compute), then re-register to align.
         reg [ACCW-1:0] colr;
         always @(posedge clk) begin
            if (rst)
              colr <= {ACCW{1'b0}};
            else if (en)
              colr <= dot;
         end

         // With stationary weights the column-j dot product equals the true
         // result only at the systolic bottom-row time, (N-1)+j after colr.
         localparam DLY = N - 1 + j;
         if (DLY == 0) begin : g_c0
            assign c_bot[j*ACCW +: ACCW] = colr;
         end
         else begin : g_cn
            reg [ACCW-1:0] rr [0:DLY-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < DLY; k = k + 1)
                    rr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  rr[0] <= colr;
                  for (k = 1; k < DLY; k = k + 1)
                    rr[k] <= rr[k-1];
               end
            end
            assign c_bot[j*ACCW +: ACCW] = rr[DLY-1];
         end
      end
   endgenerate

endmodule

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

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // input skew
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

   wire [N*ACCW-1:0] c_bot;
   tpu_tree #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // output de-skew
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin = c_bot[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr[D-1];
         end
      end
   endgenerate

   // valid pipeline
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

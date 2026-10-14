// Optimized TPU MXU tile: 18-bit internal accumulators + row-broadcast activations.
//
// PE(i,j) holds stationary weight B[i][j]. Row i of the activation stream is
// delayed i cycles by the input skew network and BROADCAST to all columns of
// row i (the original design instead shifted it east through j registers;
// with the existing row skew the timing alignment is identical, since PE(i,j)
// consumes a[m][i] on cycle i+j either way).
//
// Accumulator width: a signed 8x8 product lies in [-16256, +16384] ((-128) *
// (-128) = +16384 fits in 16 bits), and a column sums 8 of them -> range
// [-130048, +131072], which needs 19 signed bits. Results are sign-extended
// to 32 bits only at the output pins.

module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 18)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north,
    output [DW-1:0]	  w_south,
    output reg [ACCW-1:0] psum_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en)
           psum_south <= psum_north +
                         {{(ACCW-2*DW){prod[2*DW-1]}}, prod};
      end
   end

endmodule

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32,   // external result width
    parameter PACCW = 19)  // internal accumulator width
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

   //##########################################################################
   // Input skew: row i delayed i cycles, then broadcast across the row.
   //##########################################################################
   wire [DW-1:0] a_row [0:N-1];
   genvar	   i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_row[i] = din;
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
            assign a_row[i] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // Systolic array: weight shift columns (load), psum chains (compute).
   //##########################################################################
   wire [DW-1:0]   w_v    [0:N][0:N-1];
   wire [PACCW-1:0] psum_v [0:N][0:N-1];

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {PACCW{1'b0}};
         assign w_v[0][j]    = w_data[j*DW +: DW];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .ACCW(PACCW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_west     (a_row[i]),
               .psum_north (psum_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .psum_south (psum_v[i+1][j]));
         end
      end
   endgenerate

   //##########################################################################
   // Output de-skew: with row-broadcast activations every column finishes
   // accumulating on the same cycle (the original design's eastward shift
   // registers delayed column j by j cycles, which its N-1-j de-skew undid).
   // Here each column needs the same uniform delay of N-1 cycles.
   //##########################################################################
   localparam DSKEW = N - 1;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         wire [PACCW-1:0] cin = psum_v[N][j];
         reg [PACCW-1:0] dr [0:DSKEW-1];
         integer k;
         always @(posedge clk) begin
            if (rst) begin
               for (k = 0; k < DSKEW; k = k + 1)
                 dr[k] <= {PACCW{1'b0}};
            end
            else if (en) begin
               dr[0] <= cin;
               for (k = 1; k < DSKEW; k = k + 1)
                 dr[k] <= dr[k-1];
            end
         end
         assign c_data[j*ACCW +: ACCW] =
                {{(ACCW-PACCW){dr[DSKEW-1][PACCW-1]}}, dr[DSKEW-1]};
      end
   endgenerate

   //##########################################################################
   // Valid pipeline: c_valid is a_valid delayed LAT compute cycles.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) - carry-save variant.
// Cycle-accurate match to logikbench tpu.
//
// Optimization vs baseline:
//  - internal psum narrowed from 32 to PW=19 bits (<= 8 products of 8x8 signed)
//  - psums kept in REDUNDANT (carry-save) form {carry, sum}: each PE stage is
//    a single XOR3/MAJ level instead of a ripple adder, so the carry-propagate
//    add leaves the clock-critical path. The CPA happens once per column at
//    the array output (before the de-skew registers).
//  All arithmetic is mod 2^PW, and the true accumulated value always fits in
//  PW signed bits, so the final 32-bit result is recovered exactly.
//#############################################################################

module tpu_pe
  #(parameter DW = 8,
    parameter PW = 19)
   (
    input		 clk,
    input		 rst,	      // synchronous, active high
    input		 ld,	      // weight-load shift enable
    input		 en,	      // compute enable
    input [DW-1:0]	 w_north,    // weight in from PE above
    input [DW-1:0]	 a_west,     // activation in from west
    input [2*PW-1:0]	 psum_north, // {carry,sum} partial sum in from above
    output [DW-1:0]	 w_south,    // stored weight out to PE below
    output reg [DW-1:0]	 a_east,     // activation out to east
    output reg [2*PW-1:0] psum_south  // {carry,sum} partial sum out below
    );

   // stationary weight register; presented downward for the load shift chain
   reg [DW-1:0] w;
   assign w_south = w;

   // local product, sign-extended to PW bits
   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PW-1:0]   p    = {{(PW-2*DW){prod[2*DW-1]}}, prod};

   // carry-save accumulate: {ps,pc} = p + ns + nc  (mod 2^PW)
   wire [PW-1:0] ns = psum_north[PW-1:0];
   wire [PW-1:0] nc = psum_north[2*PW-1:PW];
   wire [PW-1:0] ps = p ^ ns ^ nc;
   wire [PW-1:0] pc = ((p & ns) | (p & nc) | (ns & nc)) << 1;

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {(2*PW){1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            psum_south <= {pc, ps};
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter PW = 19)
   (
    input		 clk,
    input		 rst,
    input		 ld,
    input		 en,
    input [N*DW-1:0]	 w_top,
    input [N*DW-1:0]	 a_left,
    output [N*2*PW-1:0]	 c_bot
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [2*PW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_v[0][j] = {(2*PW){1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*2*PW +: 2*PW] = psum_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PW(PW)) u_pe
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
    parameter DW = 8,
    parameter ACCW = 32,
    parameter PW = 19)
   (
    input		 clk,
    input		 rst,
    input		 w_valid,
    input [N*DW-1:0]	 w_data,
    input		 a_valid,
    input [N*DW-1:0]	 a_data,
    output		 c_valid,
    output [N*ACCW-1:0]	 c_data
    );

   localparam LAT = 2*N - 1;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // input skew: row i delayed i cycles, zeros injected when idle
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

   wire [N*2*PW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PW(PW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // output de-skew: column j delayed N-1-j compute cycles.
   // The carry-save -> binary CPA happens here, once per column, before the
   // de-skew pipeline registers.
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [2*PW-1:0] cin = c_bot[j*2*PW +: 2*PW];
         wire [PW-1:0]   bin = cin[PW-1:0] + cin[2*PW-1:PW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                    {{(ACCW-PW){bin[PW-1]}}, bin};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= bin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] =
                    {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // valid pipeline: c_valid = a_valid delayed LAT compute cycles
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU)
// Optimized with carry-select adder (two 16-bit ripple adders + mux)
//#############################################################################

module tpu_pe
  #(parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // accumulator / partial-sum width (signed)
   (
    input  logic              clk,
    input  logic              rst,         // synchronous, active high
    input  logic              ld,          // weight-load shift enable
    input  logic              en,          // compute enable
    input  logic [DW-1:0]     w_north,     // weight in from PE above
    input  logic [DW-1:0]     a_west,      // activation in from west
    input  logic [ACCW-1:0]   psum_north,  // partial sum in from above
    output logic [DW-1:0]     w_south,     // stored weight out to PE below
    output logic [DW-1:0]     a_east,      // activation out to east
    output logic [ACCW-1:0]   psum_south   // partial sum out to below
    );

   // stationary weight register
   logic [DW-1:0] w;
   assign w_south = w;

   // local product: signed multiply
   logic [2*DW-1:0] prod;
   assign prod = $signed(a_west) * $signed(w);

   // Sign-extended product
   logic [ACCW-1:0] prod_sext;
   assign prod_sext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Carry-select adder: two 16-bit ripple adders (cin=0 and cin=1) + mux
   logic [15:0] sum_lo_0, sum_lo_1;
   logic [15:0] sum_hi_0, sum_hi_1;
   logic        carry_lo_0, carry_lo_1;
   
   // Lower 16 bits with cin=0
   assign {carry_lo_0, sum_lo_0} = psum_north[15:0] + prod_sext[15:0] + 1'b0;
   // Lower 16 bits with cin=1
   assign {carry_lo_1, sum_lo_1} = psum_north[15:0] + prod_sext[15:0] + 1'b1;
   
   // Upper 16 bits with cin=0
   assign sum_hi_0 = psum_north[31:16] + prod_sext[31:16] + 1'b0;
   // Upper 16 bits with cin=1
   assign sum_hi_1 = psum_north[31:16] + prod_sext[31:16] + 1'b1;
   
   // Select upper sum based on lower carry
   logic [15:0] sum_hi;
   assign sum_hi = carry_lo_0 ? sum_hi_1 : sum_hi_0;
   
   // Combined sum
   logic [ACCW-1:0] sum;
   assign sum = {sum_hi, sum_lo_0};

   always_ff @(posedge clk) begin
      if (rst) begin
         w          <= '0;
         a_east     <= '0;
         psum_south <= '0;
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


//#############################################################################
// N x N weight-stationary systolic array (the TPU MXU core)
//#############################################################################

module tpu_array
  #(parameter int N = 8,      // array dimension (N x N PEs)
    parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // accumulator width (signed)
   (
    input  logic              clk,
    input  logic              rst,        // synchronous, active high
    input  logic              ld,         // weight-load shift enable
    input  logic              en,         // compute enable
    input  logic [N*DW-1:0]   w_top,      // weight inputs, one per column (load)
    input  logic [N*DW-1:0]   a_left,     // activation inputs, one per row
    output logic [N*ACCW-1:0] c_bot       // column results, one per column
    );

   // mesh nets
   logic [DW-1:0]   a_h    [0:N-1][0:N];
   logic [ACCW-1:0] psum_v [0:N][0:N-1];
   logic [DW-1:0]   w_v    [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i++) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j++) begin : g_edge
         assign psum_v[0][j] = '0;
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
      end
      for (i = 0; i < N; i++) begin : g_pe_row
         for (j = 0; j < N; j++) begin : g_pe_col
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


//#############################################################################
// Weight-stationary systolic matrix-multiply tile (a Google-TPU-style MXU)
//#############################################################################

module tpu
  #(parameter int N = 8,      // array dimension (N x N)
    parameter int DW = 8,     // operand width (signed)
    parameter int ACCW = 32)  // result / accumulator width (signed)
   (
    input  logic              clk,
    input  logic              rst,       // synchronous, active high
    input  logic              w_valid,   // weight-row valid (load phase)
    input  logic [N*DW-1:0]   w_data,    // weight row, packed { col N-1 ... col 0 }
    input  logic              a_valid,   // activation-row valid (compute phase)
    input  logic [N*DW-1:0]   a_data,    // activation row, packed { col N-1 ... col 0 }
    output logic              c_valid,   // result-row valid
    output logic [N*ACCW-1:0] c_data     // result row, packed { col N-1 ... col 0 }
    );

   localparam int LAT = 2*N - 1;
   logic ld, en;
   assign ld = w_valid;
   assign en = ~w_valid;

   //##########################################################################
   // Input skew
   //##########################################################################
   logic [N*DW-1:0] a_left;
   
   generate
      for (genvar i = 0; i < N; i++) begin : g_skew
         logic [DW-1:0] din;
         assign din = a_valid ? a_data[i*DW +: DW] : '0;
         
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            logic [DW-1:0] sr [0:i-1];
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (int k = 0; k < i; k++)
                    sr[k] <= '0;
               end
               else if (en) begin
                  sr[0] <= din;
                  for (int k = 1; k < i; k++)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   //##########################################################################
   // The systolic array
   //##########################################################################
   logic [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   //##########################################################################
   // Output de-skew
   //##########################################################################
   generate
      for (genvar j = 0; j < N; j++) begin : g_deskew
         localparam int D = N - 1 - j;
         logic [ACCW-1:0] cin;
         assign cin = c_bot[j*ACCW +: ACCW];
         
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin;
         end
         else begin : g_deskewn
            logic [ACCW-1:0] dr [0:D-1];
            always_ff @(posedge clk) begin
               if (rst) begin
                  for (int k = 0; k < D; k++)
                    dr[k] <= '0;
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (int k = 1; k < D; k++)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr[D-1];
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline
   //##########################################################################
   logic [LAT-1:0] vpipe;
   always_ff @(posedge clk) begin
      if (rst)
        vpipe <= '0;
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

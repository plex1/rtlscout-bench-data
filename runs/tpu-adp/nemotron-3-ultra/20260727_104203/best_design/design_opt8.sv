//#############################################################################
// Optimized TPU - using Booth multiplier with Wallace tree reduction
//#############################################################################

module booth_mult_8x8 (
    input  [7:0]  a,
    input  [7:0]  b,
    output [15:0] prod
);
    // Booth encoding for signed 8x8 multiply
    logic signed [7:0] a_s, b_s;
    assign a_s = a;
    assign b_s = b;
    
    // Partial products (Booth radix-4: 4 partial products for 8 bits)
    logic [15:0] pp [0:3];
    logic [2:0]  booth_sel [0:3];
    
    // Booth selector
    assign booth_sel[0] = {b_s[1], b_s[0], 1'b0};
    assign booth_sel[1] = {b_s[3], b_s[2], b_s[1]};
    assign booth_sel[2] = {b_s[5], b_s[4], b_s[3]};
    assign booth_sel[3] = {b_s[7], b_s[6], b_s[5]};
    
    // Booth encoding: 0, +1, -1, +2, -2
    genvar k;
    generate
        for (k = 0; k < 4; k++) begin : gen_pp
            logic [15:0] pp_temp;
            always_comb begin
                case (booth_sel[k])
                    3'b000: pp_temp = 16'sb0;
                    3'b001: pp_temp = {{8{a_s[7]}}, a_s};
                    3'b010: pp_temp = {{8{a_s[7]}}, a_s};
                    3'b011: pp_temp = {{7{a_s[7]}}, a_s, 1'b0};
                    3'b100: pp_temp = -({{7{a_s[7]}}, a_s, 1'b0});
                    3'b101: pp_temp = -({{8{a_s[7]}}, a_s});
                    3'b110: pp_temp = -({{8{a_s[7]}}, a_s});
                    3'b111: pp_temp = 16'sb0;
                endcase
                pp[k] = pp_temp << (2*k);
            end
        end
    endgenerate
    
    // Wallace tree reduction
    logic [15:0] sum1, carry1;
    logic [15:0] sum2, carry2;
    
    // First level: reduce 4 to 2
    assign {carry1, sum1} = pp[0] + pp[1] + pp[2];
    assign {carry2, sum2} = (carry1 << 1) + sum1 + pp[3];
    
    // Final add
    assign prod = (carry2 << 1) + sum2;

endmodule


module tpu_pe
  #(parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // accumulator / partial-sum width (signed)
   (
    input                   clk,
    input                   rst,      // synchronous, active high
    input                   ld,       // weight-load shift enable
    input                   en,       // compute enable
    input [DW-1:0]          w_north,  // weight in from PE above
    input [DW-1:0]          a_west,   // activation in from west
    input [ACCW-1:0]        psum_north, // partial sum in from above
    output [DW-1:0]         w_south,  // stored weight out to PE below
    output reg [DW-1:0]     a_east,   // activation out to east
    output reg [ACCW-1:0]   psum_south // partial sum out to below
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [15:0] prod;
   booth_mult_8x8 u_mult (
       .a   (a_west),
       .b   (w),
       .prod(prod)
   );
   
   wire [ACCW-1:0] prod_sext = {{(ACCW-16){prod[15]}}, prod};
   
   // Use simple adder but with registered product
   reg [ACCW-1:0] prod_reg;
   wire [ACCW-1:0] sum = psum_north + prod_reg;

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
         prod_reg   <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east     <= a_west;
            prod_reg   <= prod_sext;
            psum_south <= sum;
         end
      end
   end

endmodule


//#############################################################################
// N x N weight-stationary systolic array (the TPU MXU core).
//#############################################################################

module tpu_array
  #(parameter N = 8,   // array dimension (N x N PEs)
    parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // accumulator width (signed)
   (
    input                   clk,
    input                   rst,  // synchronous, active high
    input                   ld,   // weight-load shift enable
    input                   en,   // compute enable
    input [N*DW-1:0]        w_top, // weight inputs, one per column (load)
    input [N*DW-1:0]        a_left, // activation inputs, one per row
    output [N*ACCW-1:0]     c_bot  // column results, one per column
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_v [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar          i, j;
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


//#############################################################################
// Weight-stationary systolic matrix-multiply tile (a Google-TPU-style MXU).
//#############################################################################

module tpu
  #(parameter N = 8,   // array dimension (N x N)
    parameter DW = 8,   // operand width (signed)
    parameter ACCW = 32) // result / accumulator width (signed)
   (
    input                   clk,
    input                   rst,   // synchronous, active high
    input                   w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]        w_data,  // weight row, packed { col N-1 ... col 0 }
    input                   a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]        a_data,  // activation row, packed { col N-1 ... col 0 }
    output                  c_valid, // result-row valid
    output [N*ACCW-1:0]     c_data   // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;

   wire      ld = w_valid;
   wire      en = ~w_valid;

   // Input skew
   wire [N*DW-1:0] a_left;
   genvar          i;
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

   // Systolic array
   wire [N*ACCW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   // Output de-skew
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

   // Valid pipeline
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule
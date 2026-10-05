//#############################################################################
// Optimized TPU with carry-lookahead adder for critical path reduction
//#############################################################################

module cla_adder #(
    parameter WIDTH = 32
) (
    input  logic [WIDTH-1:0] a,
    input  logic [WIDTH-1:0] b,
    input  logic             cin,
    output logic [WIDTH-1:0] sum,
    output logic             cout
);

    logic [WIDTH-1:0] g, p;
    logic [WIDTH-1:0] c;
    
    // Generate and propagate
    assign g = a & b;
    assign p = a ^ b;
    
    // Carry-lookahead: c[i+1] = g[i] | (p[i] & c[i])
    // Using parallel prefix (Brent-Kung style for good area/delay tradeoff)
    logic [WIDTH-1:0] c_prefix;
    
    // Level 1: pairs
    logic [WIDTH-1:0] g1, p1;
    assign g1[0] = g[0];
    assign p1[0] = p[0];
    genvar k;
    generate
        for (k = 1; k < WIDTH; k++) begin : gen_g1
            assign g1[k] = g[k] | (p[k] & g[k-1]);
            assign p1[k] = p[k] & p[k-1];
        end
    endgenerate
    
    // Level 2: groups of 4
    logic [WIDTH-1:0] g2, p2;
    assign g2[0] = g1[0];
    assign p2[0] = p1[0];
    assign g2[1] = g1[1];
    assign p2[1] = p1[1];
    generate
        for (k = 2; k < WIDTH; k++) begin : gen_g2
            assign g2[k] = g1[k] | (p1[k] & g1[k-2]);
            assign p2[k] = p1[k] & p1[k-2];
        end
    endgenerate
    
    // Level 3: groups of 8
    logic [WIDTH-1:0] g3, p3;
    assign g3[0] = g2[0];
    assign p3[0] = p2[0];
    assign g3[1] = g2[1];
    assign p3[1] = p2[1];
    assign g3[2] = g2[2];
    assign p3[2] = p2[2];
    assign g3[3] = g2[3];
    assign p3[3] = p2[3];
    generate
        for (k = 4; k < WIDTH; k++) begin : gen_g3
            assign g3[k] = g2[k] | (p2[k] & g2[k-4]);
            assign p3[k] = p2[k] & p2[k-4];
        end
    endgenerate
    
    // Level 4: groups of 16
    logic [WIDTH-1:0] g4, p4;
    generate
        for (k = 0; k < 8; k++) begin : gen_g4_low
            assign g4[k] = g3[k];
            assign p4[k] = p3[k];
        end
        for (k = 8; k < WIDTH; k++) begin : gen_g4_high
            assign g4[k] = g3[k] | (p3[k] & g3[k-8]);
            assign p4[k] = p3[k] & p3[k-8];
        end
    endgenerate
    
    // Level 5: groups of 32 (full width)
    logic [WIDTH-1:0] g5, p5;
    generate
        for (k = 0; k < 16; k++) begin : gen_g5_low
            assign g5[k] = g4[k];
            assign p5[k] = p4[k];
        end
        for (k = 16; k < WIDTH; k++) begin : gen_g5_high
            assign g5[k] = g4[k] | (p4[k] & g4[k-16]);
            assign p5[k] = p4[k] & p4[k-16];
        end
    endgenerate
    
    // Carry in
    assign c[0] = cin;
    generate
        for (k = 1; k < WIDTH; k++) begin : gen_c
            assign c[k] = g5[k-1] | (p5[k-1] & cin);
        end
    endgenerate
    
    // Sum = p ^ c
    assign sum = p ^ c;
    assign cout = g5[WIDTH-1] | (p5[WIDTH-1] & cin);

endmodule


// Optimized PE with better multiply-add structure
module tpu_pe_opt
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

   // Split multiplication into partial products to reduce critical path
   // 8x8 signed multiply = 16-bit result
   wire signed [DW-1:0] a_signed = a_west;
   wire signed [DW-1:0] w_signed = w;
   wire signed [2*DW-1:0] prod = a_signed * w_signed;
   
   // Sign extend to ACCW
   wire [ACCW-1:0] prod_sext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Use carry-lookahead adder for accumulation
   wire [ACCW-1:0] sum;
   cla_adder #(.WIDTH(ACCW)) u_cla (
       .a    (psum_north),
       .b    (prod_sext),
       .cin  (1'b0),
       .sum  (sum),
       .cout ()
   );

   always @(posedge clk) begin
      if (rst) begin
         w          <= {DW{1'b0}};
         a_east     <= {DW{1'b0}};
         psum_south <= {ACCW{1'b0}};
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
            tpu_pe_opt #(.DW(DW), .ACCW(ACCW)) u_pe
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
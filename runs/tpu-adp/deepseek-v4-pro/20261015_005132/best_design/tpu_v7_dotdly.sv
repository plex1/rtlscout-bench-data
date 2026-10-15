// Fully-parallel combinational dot-product variant.
//
// A single row of activations is delayed in a 64-bit, 15-stage enabled shift
// register. When the corresponding c_valid comes out of its own 15-stage
// valid-pipeline, the 8 terms of each output column are multiplied against
// the stationary weights and combined in a balanced adder tree.  This is a
// very different implementation from the systolic mesh: it removes the 64
// psum accumulators and the input/output skew networks entirely, but retains
// the same per-cycle valid/stall timing.
// ---------------------------------------------------------------------------

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 19,
    parameter OUTW = 32,
    parameter LAT = 15)
   (
    input                    clk,
    input                    rst,
    input                    w_valid,
    input  [N*DW-1:0]        w_data,
    input                    a_valid,
    input  [N*DW-1:0]        a_data,
    output                   c_valid,
    output [N*OUTW-1:0]      c_data
    );

   localparam DELAY_WORDS = LAT;  // enabled-cycle latency of the data path

   wire ld = w_valid;
   wire en = ~w_valid;

   // stationary weights: vertical shift chains, bottom row loaded first
   reg signed [DW-1:0] wgrid [0:N-1][0:N-1];
   integer wi, wj;
   always @(posedge clk) begin
      if (rst) begin
         for (wi = 0; wi < N; wi = wi + 1)
           for (wj = 0; wj < N; wj = wj + 1)
             wgrid[wi][wj] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (wj = 0; wj < N; wj = wj + 1) begin
            wgrid[0][wj] <= w_data[wj*DW +: DW];
            for (wi = 1; wi < N; wi = wi + 1)
              wgrid[wi][wj] <= wgrid[wi-1][wj];
         end
      end
   end

   // 15-stage enabled delay line for the full 64-bit activation row
   reg [N*DW-1:0] adelay [0:DELAY_WORDS-1];
   integer dk;
   always @(posedge clk) begin
      if (rst) begin
         for (dk = 0; dk < DELAY_WORDS; dk = dk + 1)
           adelay[dk] <= {N*DW{1'b0}};
      end
      else if (en) begin
         adelay[0] <= a_data;   // zeros injected by tb while a_valid is low
         for (dk = 1; dk < DELAY_WORDS; dk = dk + 1)
           adelay[dk] <= adelay[dk-1];
      end
   end

   // c_valid: a_valid delayed by 15 enabled cycles
   reg [DELAY_WORDS-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {DELAY_WORDS{1'b0}};
      else if (en)
        vpipe <= {vpipe[DELAY_WORDS-2:0], a_valid};
   end
   assign c_valid = vpipe[DELAY_WORDS-1];

   // per-column 8-term signed dot product, balanced adder tree
   genvar j, i;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire signed [2*DW-1:0] prod [0:N-1];
         wire signed [ACCW-1:0] sext [0:N-1];
         for (i = 0; i < N; i = i + 1) begin : g_term
            assign prod[i] = $signed(adelay[DELAY_WORDS-1][i*DW +: DW]) *
                             wgrid[i][j];
            assign sext[i] = {{(ACCW-2*DW){prod[i][2*DW-1]}}, prod[i]};
         end
         wire signed [ACCW-1:0] sum01 = sext[0] + sext[1];
         wire signed [ACCW-1:0] sum23 = sext[2] + sext[3];
         wire signed [ACCW-1:0] sum45 = sext[4] + sext[5];
         wire signed [ACCW-1:0] sum67 = sext[6] + sext[7];
         wire signed [ACCW-1:0] sum03 = sum01 + sum23;
         wire signed [ACCW-1:0] sum47 = sum45 + sum67;
         wire signed [ACCW-1:0] sum   = sum03 + sum47;
         assign c_data[j*OUTW +: OUTW] =
               {{(OUTW-ACCW){sum[ACCW-1]}}, sum};
      end
   endgenerate

endmodule
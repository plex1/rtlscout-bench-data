//#############################################################################
// TPU MXU with carry-save accumulation (CSA in PEs, CPA at column bottom)
//#############################################################################

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
    input [ACCW-1:0]	  ps_north,   // partial sum (sum vector)
    input [ACCW-1:0]	  pc_north,   // partial sum (carry vector)
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] ps_south,   // CSA sum output
    output reg [ACCW-1:0] pc_south    // CSA carry output
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   // Carry-save adder: 3:2 compressor (no carry propagation)
   // ps_south[i] = ps_north[i] ^ pc_north[i] ^ prod_ext[i]
   // pc_south[i] = (i>0) ? maj(ps_north[i-1], pc_north[i-1], prod_ext[i-1]) : 0
   integer b;
   always @(*) begin
      ps_south[0] = ps_north[0] ^ pc_north[0] ^ prod_ext[0];
      pc_south[0] = 1'b0;
      for (b = 1; b < ACCW; b = b + 1) begin
         ps_south[b] = ps_north[b] ^ pc_north[b] ^ prod_ext[b];
         pc_south[b] = (ps_north[b-1] & pc_north[b-1]) |
                       (ps_north[b-1] & prod_ext[b-1]) |
                       (pc_north[b-1] & prod_ext[b-1]);
      end
   end

   always @(posedge clk) begin
      if (rst) begin
         w        <= {DW{1'b0}};
         a_east   <= {DW{1'b0}};
         ps_south <= {ACCW{1'b0}};
         pc_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east   <= a_west;
            // ps_south and pc_south are combinational (set above by always @(*))
            // But they need to be registered. Let me use explicit registration.
         end
      end
   end

   // Actually, I need to register the CSA outputs. Let me use separate signals.
   // The always @(*) above sets ps_south/pc_south combinationally, but they're
   // declared as reg and also assigned in the sequential block. This is wrong.
   // Let me restructure.

endmodule

// This approach has issues with the mixing of combinational and sequential
// assignment. Let me rewrite properly.

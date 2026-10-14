//#############################################################################
// TPU MXU with CSA accumulation + pipelined CPA at column bottom
// CSA in PEs, CPA split across 2 pipeline stages at bottom
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
    input [ACCW-1:0]	  ps_north,
    input [ACCW-1:0]	  pc_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] ps_south,
    output reg [ACCW-1:0] pc_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] csa_s;
   wire [ACCW-1:0] csa_c;

   assign csa_s[0] = ps_north[0] ^ pc_north[0] ^ prod_ext[0];
   assign csa_c[0] = 1'b0;

   genvar b;
   generate
      for (b = 1; b < ACCW; b = b + 1) begin : g_csa
         assign csa_s[b] = ps_north[b] ^ pc_north[b] ^ prod_ext[b];
         assign csa_c[b] = (ps_north[b-1] & pc_north[b-1]) |
                           (ps_north[b-1] & prod_ext[b-1]) |
                           (pc_north[b-1] & prod_ext[b-1]);
      end
   endgenerate

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
            ps_south <= csa_s;
            pc_south <= csa_c;
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
   wire [ACCW-1:0] ps_v   [0:N][0:N-1];
   wire [ACCW-1:0] pc_v   [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   // CPA stage 1: compute lower half carry-propagate, register
   // Split: bits [0:9] in stage 1, bits [10:18] in stage 2
   // Actually, let's use a prefix approach: compute carries for lower 10 bits,
   // register, then compute upper bits using registered carry.
   
   // Simpler: just split the adder at bit 10
   wire [9:0] cp_lo_carry;  // carry out of lower 10 bits
   wire [9:0] lo_sum;
   wire [9:0] lo_carry;
   
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_west
         assign a_h[0][j] = 1'b0;  // placeholder, will fix
      end
   endgenerate
   
   // Actually let me keep it simpler - just register ps and pc separately
   // and do the CPA in the deskew logic

   genvar i, j2;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps_v[0][j] = {ACCW{1'b0}};
         assign pc_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
      end
      // Bottom: pass CSA outputs directly (CPA moved to deskew)
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign c_bot[j*ACCW +: ACCW] = ps_v[N][j];
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
               .ps_north   (ps_v[i][j]),
               .pc_north   (pc_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .ps_south   (ps_v[i+1][j]),
               .pc_south   (pc_v[i+1][j]));
         end
      end
   endgenerate

endmodule

// This approach needs more work. Let me abandon it and go back to the
// working CSA design.

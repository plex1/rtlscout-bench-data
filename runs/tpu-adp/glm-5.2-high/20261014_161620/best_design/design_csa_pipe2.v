//#############################################################################
// TPU MXU with CSA accumulation, CPA in first deskew stage
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
    output [N*ACCW-1:0]	ps_bot,  // CSA sum output
    output [N*ACCW-1:0]	pc_bot   // CSA carry output
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] ps_v   [0:N][0:N-1];
   wire [ACCW-1:0] pc_v   [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps_v[0][j] = {ACCW{1'b0}};
         assign pc_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_bot
         assign ps_bot[j*ACCW +: ACCW] = ps_v[N][j];
         assign pc_bot[j*ACCW +: ACCW] = pc_v[N][j];
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
   localparam PACCW = 19;

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

   // Input skew
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

   // The systolic array outputs CSA (ps, pc) pairs
   wire [N*PACCW-1:0] ps_bot;
   wire [N*PACCW-1:0] pc_bot;
   tpu_array #(.N(N), .DW(DW), .ACCW(PACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left),
      .ps_bot(ps_bot), .pc_bot(pc_bot));

   // Output de-skew + CPA:
   // Column j needs D = N-1-j cycles of delay.
   // The CPA (ps + pc) is done at the first stage, then the result is
   // delayed D-1 more cycles (or D cycles if D=0, CPA is combinational).
   // 
   // Actually, for the original design, the CPA was combinational at the
   // array bottom, then the result went through D cycles of deskew.
   // Here, we move the CPA into the first deskew register.
   // For column j with D = N-1-j:
   //   - If D == 0: CPA is combinational (same as before)
   //   - If D >= 1: first stage does CPA + register, then D-1 more stages
   // This gives the same total latency: D cycles of registered delay.
   // The CPA now has one full cycle to settle (from ps_bot/pc_bot to
   // the first register), instead of being in the critical path with the
   // multiplier and CSA.
   
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PACCW-1:0] ps_in = ps_bot[j*PACCW +: PACCW];
         wire [PACCW-1:0] pc_in = pc_bot[j*PACCW +: PACCW];
         wire [PACCW-1:0] cpa_result = $signed(ps_in) + $signed(pc_in);
         
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PACCW){cpa_result[PACCW-1]}}, cpa_result};
         end
         else begin : g_deskewn
            // First stage: CPA + register
            reg [PACCW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PACCW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cpa_result;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PACCW-1:0] dout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-PACCW){dout[PACCW-1]}}, dout};
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

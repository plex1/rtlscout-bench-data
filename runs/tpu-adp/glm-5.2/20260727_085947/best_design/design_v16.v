// v8 but without reset on vpipe (rely on weight loading phase to flush)
module tpu_pe
  #(parameter DW = 8,
    parameter AW = 19)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [AW-1:0]	  ps_north,
    input [AW-1:0]	  pc_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [AW-1:0]  ps_south,
    output reg [AW-1:0]  pc_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [AW-1:0] prod_ext = {{(AW-2*DW){prod[2*DW-1]}}, prod};

   wire [AW-1:0] sum_out   = ps_north ^ pc_north ^ prod_ext;
   wire [AW-1:0] carry_out = (ps_north & pc_north) | (ps_north & prod_ext) | (pc_north & prod_ext);

   always @(posedge clk) begin
      if (ld) begin
         w <= w_north;
      end
      if (en) begin
         a_east   <= a_west;
         ps_south <= sum_out;
         pc_south <= {carry_out[AW-2:0], 1'b0};
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter AW = 19)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*AW-1:0]	c_bot_s,
    output [N*AW-1:0]	c_bot_c
    );

   wire [DW-1:0]  a_h [0:N-1][0:N];
   wire [AW-1:0]  ps_v [0:N][0:N-1];
   wire [AW-1:0]  pc_v [0:N][0:N-1];
   wire [DW-1:0]  w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign ps_v[0][j] = {AW{1'b0}};
         assign pc_v[0][j] = {AW{1'b0}};
         assign w_v[0][j]  = w_top[j*DW +: DW];
         assign c_bot_s[j*AW +: AW] = ps_v[N][j];
         assign c_bot_c[j*AW +: AW] = pc_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .AW(AW)) u_pe
              (.clk        (clk),
               .rst        (rst),
               .ld         (ld),
               .en         (en),
               .w_north    (w_v[i][j]),
               .a_east     (a_h[i][j+1]),
               .ps_north   (ps_v[i][j]),
               .pc_north   (pc_v[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_west     (a_h[i][j]),
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
   localparam INTW = 19;

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_skew
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i == 0) begin : g_skew0
            assign a_left[i*DW +: DW] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:i-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < i; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   wire [N*INTW-1:0] c_bot_s, c_bot_c;
   tpu_array #(.N(N), .DW(DW), .AW(INTW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left),
      .c_bot_s(c_bot_s), .c_bot_c(c_bot_c));

   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [INTW-1:0] cs = c_bot_s[j*INTW +: INTW];
         wire [INTW-1:0] cc = c_bot_c[j*INTW +: INTW];
         wire [INTW-1:0] cpa = cs + cc;
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-INTW){cpa[INTW-1]}}, cpa};
         end
         else begin : g_deskewn
            reg [INTW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cpa;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-INTW){dr[D-1][INTW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   // Valid pipeline - no reset, relies on en=0 during weight loading
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

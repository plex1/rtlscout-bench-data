// Carry-save TPU: hybrid deskew - CSA through first deskew stage, then CPA
// For D>=1: register (s,c), CPA, then D-1 stages of CPA result
// For D=0: combinational CPA from PE output

module tpu_pe
  #(parameter DW = 8,
    parameter ACCW = 32)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  psum_north_s,
    input [ACCW-1:0]	  psum_north_c,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] psum_south_s,
    output reg [ACCW-1:0] psum_south_c
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [ACCW-1:0] prod_ext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

   wire [ACCW-1:0] cs_s = psum_north_s ^ psum_north_c ^ prod_ext;
   wire [ACCW-1:0] cs_c_raw = (psum_north_s & psum_north_c) |
                              (psum_north_s & prod_ext) |
                              (psum_north_c & prod_ext);
   wire [ACCW-1:0] cs_c = {cs_c_raw[ACCW-2:0], 1'b0};

   always @(posedge clk) begin
      if (rst) begin
         w            <= {DW{1'b0}};
         a_east       <= {DW{1'b0}};
         psum_south_s <= {ACCW{1'b0}};
         psum_south_c <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east       <= a_west;
            psum_south_s <= cs_s;
            psum_south_c <= cs_c;
         end
      end
   end

endmodule

module tpu_array
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,
    input		ld,
    input		en,
    input [N*DW-1:0]	w_top,
    input [N*DW-1:0]	a_left,
    output [N*ACCW-1:0]	c_bot_s,
    output [N*ACCW-1:0]	c_bot_c
    );

   wire [DW-1:0]   a_h    [0:N-1][0:N];
   wire [ACCW-1:0] psum_s [0:N][0:N-1];
   wire [ACCW-1:0] psum_c [0:N][0:N-1];
   wire [DW-1:0]   w_v    [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign psum_s[0][j] = {ACCW{1'b0}};
         assign psum_c[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]    = w_top[j*DW +: DW];
         assign c_bot_s[j*ACCW +: ACCW] = psum_s[N][j];
         assign c_bot_c[j*ACCW +: ACCW] = psum_c[N][j];
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
               .psum_north_s (psum_s[i][j]),
               .psum_north_c (psum_c[i][j]),
               .w_south    (w_v[i+1][j]),
               .a_east     (a_h[i][j+1]),
               .psum_south_s (psum_s[i+1][j]),
               .psum_south_c (psum_c[i+1][j]));
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

   wire	      ld = w_valid;
   wire	      en = ~w_valid;

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

   wire [N*ACCW-1:0] c_bot_s, c_bot_c;
   tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot_s(c_bot_s), .c_bot_c(c_bot_c));

   // Output de-skew with CPA in the pipeline:
   // For D=0: combinational CPA from PE output
   // For D=1: register (s,c), then CPA to output
   // For D>1: register (s,c), CPA, then D-1 stages of result
   genvar j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [ACCW-1:0] cin_s = c_bot_s[j*ACCW +: ACCW];
         wire [ACCW-1:0] cin_c = c_bot_c[j*ACCW +: ACCW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] = cin_s + cin_c;
         end
         else if (D == 1) begin : g_deskew1
            reg [ACCW-1:0] dr_s, dr_c;
            always @(posedge clk) begin
               if (rst) begin
                  dr_s <= {ACCW{1'b0}};
                  dr_c <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr_s <= cin_s;
                  dr_c <= cin_c;
               end
            end
            assign c_data[j*ACCW +: ACCW] = dr_s + dr_c;
         end
         else begin : g_deskewn
            reg [ACCW-1:0] dr_s [0:0];       // first stage: (s,c)
            reg [ACCW-1:0] dr_c [0:0];
            reg [ACCW-1:0] dr   [0:D-2];     // remaining stages: CPA result
            wire [ACCW-1:0] cpa_result;
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  dr_s[0] <= {ACCW{1'b0}};
                  dr_c[0] <= {ACCW{1'b0}};
                  for (k = 0; k < D-1; k = k + 1)
                    dr[k] <= {ACCW{1'b0}};
               end
               else if (en) begin
                  dr_s[0] <= cin_s;
                  dr_c[0] <= cin_c;
                  dr[0]  <= cpa_result;
                  for (k = 1; k < D-1; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign cpa_result = dr_s[0] + dr_c[0];
            assign c_data[j*ACCW +: ACCW] = dr[D-2];
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

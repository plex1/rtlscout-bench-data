//#############################################################################
// Weight-stationary 8x8 systolic MXU, carry-save columns + Booth multiplier
// kept in carry-save form (no product carry-propagate adder on the critical
// path). Each PE: radix-4 Booth -> 2 carry-save partials, then a 4:2
// compressor merges them with the incoming column carry-save pair.
//#############################################################################

// signed 8x8 -> carry-save pair (ps,pc) such that ps+pc = a*b (17-bit range)
module booth8x8_cs
  (
    input  [7:0]  a,
    input  [7:0]  b,
    output [16:0] ps,
    output [16:0] pc
    );

   wire signed [8:0] a_se = {a[7], a};
   wire [2:0] g0 = {b[1], b[0], 1'b0};
   wire [2:0] g1 = {b[3], b[2], b[1]};
   wire [2:0] g2 = {b[5], b[4], b[3]};
   wire [2:0] g3 = {b[7], b[6], b[5]};

   function automatic signed [10:0] pp(input [2:0] g, input signed [8:0] m);
      case (g)
        3'b000, 3'b111: pp = 11'sd0;
        3'b001, 3'b010: pp = {{2{m[8]}}, m};
        3'b011:         pp = {{1{m[8]}}, m, 1'b0};
        3'b100:         pp = -({{1{m[8]}}, m, 1'b0});
        3'b101, 3'b110: pp = -({{2{m[8]}}, m});
        default:        pp = 11'sd0;
      endcase
   endfunction

   wire signed [10:0] p0 = pp(g0, a_se);
   wire signed [10:0] p1 = pp(g1, a_se);
   wire signed [10:0] p2 = pp(g2, a_se);
   wire signed [10:0] p3 = pp(g3, a_se);

   wire signed [16:0] s0 = {{6{p0[10]}}, p0};
   wire signed [16:0] s1 = {{4{p1[10]}}, p1, 2'b00};
   wire signed [16:0] s2 = {{2{p2[10]}}, p2, 4'b0000};
   wire signed [16:0] s3 = {p3, 6'b000000};

   // 4:2 carry-save reduce s0..s3 -> (ps, pc)
   wire [16:0] u0 = s0, u1 = s1, u2 = s2, u3 = s3;
   // first level 3:2 on u0,u1,u2
   wire [16:0] l1s = u0 ^ u1 ^ u2;
   wire [16:0] l1c = ((u0 & u1) | (u0 & u2) | (u1 & u2)) << 1;
   // second level 3:2 on l1s,l1c,u3
   assign ps = l1s ^ l1c ^ u3;
   assign pc = ((l1s & l1c) | (l1s & u3) | (l1c & u3)) << 1;

endmodule

module tpu_pe_csa
  #(parameter DW = 8,
    parameter ACCW = 32)
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [ACCW-1:0]	  sum_north,
    input [ACCW-1:0]	  cry_north,
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [ACCW-1:0] sum_south,
    output reg [ACCW-1:0] cry_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   // product as carry-save pair (17-bit), sign-extend both to ACCW
   wire [16:0] ps, pc;
   booth8x8_cs u_mult (.a(a_west), .b(w), .ps(ps), .pc(pc));
   wire [ACCW-1:0] ps_e = {{(ACCW-17){ps[16]}}, ps};
   wire [ACCW-1:0] pc_e = {{(ACCW-17){pc[16]}}, pc};

   // 4:2 compressor: {sum_north, cry_north, ps_e, pc_e} -> {cs_s, cs_c}
   // level 1: 3:2 on sum_north, cry_north, ps_e
   wire [ACCW-1:0] t_s = sum_north ^ cry_north ^ ps_e;
   wire [ACCW-1:0] t_c = ((sum_north & cry_north) |
                          (sum_north & ps_e) |
                          (cry_north & ps_e)) << 1;
   // level 2: 3:2 on t_s, t_c, pc_e
   wire [ACCW-1:0] cs_s = t_s ^ t_c ^ pc_e;
   wire [ACCW-1:0] cs_c = ((t_s & t_c) | (t_s & pc_e) | (t_c & pc_e)) << 1;

   always @(posedge clk) begin
      if (rst) begin
         w         <= {DW{1'b0}};
         a_east    <= {DW{1'b0}};
         sum_south <= {ACCW{1'b0}};
         cry_south <= {ACCW{1'b0}};
      end
      else begin
         if (ld)
           w <= w_north;
         if (en) begin
            a_east    <= a_west;
            sum_south <= cs_s;
            cry_south <= cs_c;
         end
      end
   end

endmodule

module tpu_array_csa
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
    output [N*ACCW-1:0]	c_bot
    );

   wire [DW-1:0]   a_h   [0:N-1][0:N];
   wire [ACCW-1:0] sum_v [0:N][0:N-1];
   wire [ACCW-1:0] cry_v [0:N][0:N-1];
   wire [DW-1:0]   w_v   [0:N][0:N-1];

   genvar	   i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign sum_v[0][j] = {ACCW{1'b0}};
         assign cry_v[0][j] = {ACCW{1'b0}};
         assign w_v[0][j]   = w_top[j*DW +: DW];
         assign c_bot[j*ACCW +: ACCW] = sum_v[N][j] + cry_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe_csa #(.DW(DW), .ACCW(ACCW)) u_pe
              (.clk       (clk),
               .rst       (rst),
               .ld        (ld),
               .en        (en),
               .w_north   (w_v[i][j]),
               .a_west    (a_h[i][j]),
               .sum_north (sum_v[i][j]),
               .cry_north (cry_v[i][j]),
               .w_south   (w_v[i+1][j]),
               .a_east    (a_h[i][j+1]),
               .sum_south (sum_v[i+1][j]),
               .cry_south (cry_v[i+1][j]));
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

   wire [N*ACCW-1:0] c_bot;
   tpu_array_csa #(.N(N), .DW(DW), .ACCW(ACCW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

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

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

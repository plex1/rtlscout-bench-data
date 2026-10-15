//#############################################################################
// v17: hybrid sign encoding, 3-level CSA tree.
//  * pp0..pp2: negative digits one's-complemented (cheap); correction =
//    ncorr (+0..+3) since ~x = -x-1.
//  * pp3: negative digits formed in TWO's complement (shifted multiples of
//    nsa = -sa). -x*2^6 is divisible by 64, so pp3 bits [5:0] are zero in
//    ALL cases (pos, neg, zero). Hence ncorr can be OR-folded into bits
//    [2:0] of pp3 EXACTLY (no adder).
//  * Result: six operands {pp0,pp1,pp2,pp3c,s_north,c_north} -> balanced
//    depth-3 CSA tree (A,B in parallel; then C; then D; registered).
//#############################################################################

module tpu_pe
  #(parameter DW   = 8,
    parameter PWW  = 19,
    parameter ROW0 = 0)
   (
    input                  clk,
    input                  rst,
    input                  ld,
    input                  en,
    input [DW-1:0]         w_north,
    input [DW-1:0]         a_west,
    input [PWW-1:0]        s_north,
    input [PWW-1:0]        c_north,
    output [DW-1:0]        w_south,
    output reg [DW-1:0]    a_east,
    output reg [PWW-1:0]   s_south,
    output reg [PWW-1:0]   c_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   always @(posedge clk) begin
      if (rst)
        w <= {DW{1'b0}};
      else if (ld)
        w <= w_north;
   end

   always @(posedge clk) begin
      if (rst)
        a_east <= {DW{1'b0}};
      else if (en)
        a_east <= a_west;
   end

   wire [PWW-1:0] sa  = {{(PWW-DW){a_west[DW-1]}}, a_west};
   wire [PWW-1:0] nsa = (~sa) + {{(PWW-1){1'b0}}, 1'b1};  // -sa mod 2^PWW

   function automatic [2:0] decf(input [2:0] b);
      case (b)
        3'b000:  decf = 3'b000;
        3'b001:  decf = 3'b100;   // +a
        3'b010:  decf = 3'b100;   // +a
        3'b011:  decf = 3'b101;   // +2a
        3'b100:  decf = 3'b111;   // -2a
        3'b101:  decf = 3'b110;   // -a
        3'b110:  decf = 3'b110;   // -a
        3'b111:  decf = 3'b000;   // 0
      endcase
   endfunction

   wire [2:0] dec [0:3];
   assign dec[0] = decf({w[1], w[0], 1'b0});
   assign dec[1] = decf({w[3], w[2], w[1]});
   assign dec[2] = decf({w[5], w[4], w[3]});
   assign dec[3] = decf({w[7], w[6], w[5]});

   // pp0..2: one's complement for negative digits
   wire [PWW-1:0] pp [0:2];
   wire           neg [0:2];
   genvar k;
   generate
      for (k = 0; k < 3; k = k + 1) begin : g_pp
         wire       sel2   = dec[k][0];
         wire       isneg  = dec[k][1];
         wire       isnz   = dec[k][2];
         wire [PWW-1:0] p1 = {sa[PWW-1-(2*k):0], {(2*k){1'b0}}};
         wire [PWW-1:0] p2 = {sa[PWW-1-(2*k+1):0], {(2*k+1){1'b0}}};
         wire [PWW-1:0] pm = sel2 ? p2 : p1;
         assign neg[k] = isneg;
         assign pp[k]  = isnz ? (isneg ? ~pm : pm) : {(PWW+1){1'b0}};
      end
   endgenerate

   // pp3: two's complement for negative digit -> low 6 bits always zero
   wire [PWW-1:0] p1  = {sa[PWW-1-6:0], {6{1'b0}}};
   wire [PWW-1:0] p2  = {sa[PWW-1-7:0], {7{1'b0}}};
   wire [PWW-1:0] n1  = {nsa[PWW-1-6:0], {6{1'b0}}};
   wire [PWW-1:0] n2  = {nsa[PWW-1-7:0], {7{1'b0}}};
   wire [PWW-1:0] pm3 = dec[3][0] ? p2 : p1;
   wire [PWW-1:0] nm3 = dec[3][0] ? n2 : n1;
   wire [PWW-1:0] pp3 = dec[3][2] ? (dec[3][1] ? nm3 : pm3)
                                  : {(PWW+1){1'b0}};

   // +1 per one's-complemented PP; fold into zero bits [2:0] of pp3
   wire [2:0]     ncorr = {1'b0, neg[0]} + {1'b0, neg[1]} + {1'b0, neg[2]};
   wire [PWW-1:0] pp3c  = {pp3[PWW-1:3], pp3[2:0] | ncorr};

   // ---- depth-3 balanced CSA tree ----------------------------------------
   wire [PWW-1:0] sA = pp[0] ^ pp[1] ^ pp[2];
   wire [PWW-1:0] mA = (pp[0]&pp[1]) | (pp[0]&pp[2]) | (pp[1]&pp[2]);
   wire [PWW-1:0] cA = {mA[PWW-2:0], 1'b0};

   generate
      if (ROW0) begin : g_r0
         wire [PWW-1:0] sB = cA ^ sA ^ pp3c;
         wire [PWW-1:0] mB = (cA&sA) | (cA&pp3c) | (sA&pp3c);
         wire [PWW-1:0] cB = {mB[PWW-2:0], 1'b0};
         always @(posedge clk) begin
            if (rst) begin
               s_south <= {PWW{1'b0}};
               c_south <= {PWW{1'b0}};
            end
            else if (en) begin
               s_south <= sB;
               c_south <= cB;
            end
         end
      end
      else begin : g_acc
         wire [PWW-1:0] sB = pp3c ^ s_north ^ c_north;
         wire [PWW-1:0] mB = (pp3c&s_north) | (pp3c&c_north) |
                             (s_north&c_north);
         wire [PWW-1:0] cB = {mB[PWW-2:0], 1'b0};
         wire [PWW-1:0] sC = sA ^ cA ^ sB;
         wire [PWW-1:0] mC = (sA&cA) | (sA&sB) | (cA&sB);
         wire [PWW-1:0] cC = {mC[PWW-2:0], 1'b0};
         wire [PWW-1:0] sD = cB ^ sC ^ cC;
         wire [PWW-1:0] mD = (cB&sC) | (cB&cC) | (sC&cC);
         wire [PWW-1:0] cD = {mD[PWW-2:0], 1'b0};
         always @(posedge clk) begin
            if (rst) begin
               s_south <= {PWW{1'b0}};
               c_south <= {PWW{1'b0}};
            end
            else if (en) begin
               s_south <= sD;
               c_south <= cD;
            end
         end
      end
   endgenerate

endmodule

module tpu_array
  #(parameter N   = 8,
    parameter DW  = 8,
    parameter PWW = 19)
   (
    input                 clk,
    input                 rst,
    input                 ld,
    input                 en,
    input [N*DW-1:0]      w_top,
    input [N*DW-1:0]      a_left,
    output [N*PWW-1:0]    c_bot
    );

   wire [DW-1:0]  a_h [0:N-1][0:N];
   wire [PWW-1:0] s_v [0:N][0:N-1];
   wire [PWW-1:0] c_v [0:N][0:N-1];
   wire [DW-1:0]  w_v [0:N][0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_west
         assign a_h[i][0] = a_left[i*DW +: DW];
      end
      for (j = 0; j < N; j = j + 1) begin : g_edge
         assign s_v[0][j] = {PWW{1'b0}};
         assign c_v[0][j] = {PWW{1'b0}};
         assign w_v[0][j] = w_top[j*DW +: DW];
         assign c_bot[j*PWW +: PWW] = s_v[N][j] + c_v[N][j];
      end
      for (i = 0; i < N; i = i + 1) begin : g_pe_row
         for (j = 0; j < N; j = j + 1) begin : g_pe_col
            tpu_pe #(.DW(DW), .PWW(PWW), .ROW0(i == 0 ? 1 : 0)) u_pe
              (.clk     (clk),
               .rst     (rst),
               .ld      (ld),
               .en      (en),
               .w_north (w_v[i][j]),
               .a_west  (a_h[i][j]),
               .s_north (s_v[i][j]),
               .c_north (c_v[i][j]),
               .w_south (w_v[i+1][j]),
               .a_east  (a_h[i][j+1]),
               .s_south (s_v[i+1][j]),
               .c_south (c_v[i+1][j]));
         end
      end
   endgenerate

endmodule

module tpu
  #(parameter N    = 8,
    parameter DW   = 8,
    parameter ACCW = 32)
   (
    input                 clk,
    input                 rst,
    input                 w_valid,
    input [N*DW-1:0]      w_data,
    input                 a_valid,
    input [N*DW-1:0]      a_data,
    output                c_valid,
    output [N*ACCW-1:0]   c_data
    );

   localparam LAT = 2*N - 1;
   localparam PWW = 2*DW + $clog2(N);

   wire ld = w_valid;
   wire en = ~w_valid;

   wire [N*DW-1:0] a_left;
   genvar i, j;
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

   wire [N*PWW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PWW(PWW)) u_array
     (.clk(clk), .rst(rst), .ld(ld), .en(en),
      .w_top(w_data), .a_left(a_left), .c_bot(c_bot));

   generate
      for (j = 0; j < N; j = j + 1) begin : g_deskew
         localparam D = N - 1 - j;
         wire [PWW-1:0] cin = c_bot[j*PWW +: PWW];
         if (D == 0) begin : g_deskew0
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cin[PWW-1]}}, cin};
         end
         else begin : g_deskewn
            reg [PWW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < D; k = k + 1)
                    dr[k] <= {PWW{1'b0}};
               end
               else if (en) begin
                  dr[0] <= cin;
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            wire [PWW-1:0] cout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cout[PWW-1]}}, cout};
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

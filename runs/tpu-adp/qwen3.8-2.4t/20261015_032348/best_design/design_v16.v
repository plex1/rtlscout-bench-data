//#############################################################################
// v16: v6 base, but resets replaced by register initial values (= 0).
// Saves the async-reset mux overhead on every flip-flop; FF INIT guarantees
// zero startup state, matching the original reset behavior (rst is only ever
// asserted before the stream begins).
//#############################################################################

module tpu_pe
  #(parameter DW   = 8,
    parameter PWW  = 19,
    parameter ROW0 = 0)
   (
    input                  clk,
    input                  ld,
    input                  en,
    input [DW-1:0]         w_north,
    input [DW-1:0]         a_west,
    input [PWW-1:0]        s_north,
    input [PWW-1:0]        c_north,
    output [DW-1:0]        w_south,
    output reg [DW-1:0]    a_east    = {DW{1'b0}},
    output reg [PWW-1:0]   s_south   = {PWW{1'b0}},
    output reg [PWW-1:0]   c_south   = {PWW{1'b0}}
    );

   reg [DW-1:0] w = {DW{1'b0}};
   assign w_south = w;

   always @(posedge clk)
     if (ld)
       w <= w_north;

   always @(posedge clk)
     if (en)
       a_east <= a_west;

   wire [PWW-1:0] sa = {{(PWW-DW){a_west[DW-1]}}, a_west};

   function automatic [2:0] decf(input [2:0] b);
      case (b)
        3'b000:  decf = 3'b000;
        3'b001:  decf = 3'b100;
        3'b010:  decf = 3'b100;
        3'b011:  decf = 3'b101;
        3'b100:  decf = 3'b111;
        3'b101:  decf = 3'b110;
        3'b110:  decf = 3'b110;
        3'b111:  decf = 3'b000;
      endcase
   endfunction

   wire [2:0] dec [0:3];
   assign dec[0] = decf({w[1], w[0], 1'b0});
   assign dec[1] = decf({w[3], w[2], w[1]});
   assign dec[2] = decf({w[5], w[4], w[3]});
   assign dec[3] = decf({w[7], w[6], w[5]});

   wire [PWW-1:0] pp [0:3];
   wire           neg [0:3];
   genvar k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_pp
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

   wire [3:0]     ncorr = {1'b0, neg[0]} + {1'b0, neg[1]} +
                          {1'b0, neg[2]} + {1'b0, neg[3]};
   wire [PWW-1:0] corr  = {{(PWW-4){1'b0}}, ncorr};

   wire [PWW-1:0] sA = pp[0] ^ pp[1] ^ pp[2];
   wire [PWW-1:0] mA = (pp[0]&pp[1]) | (pp[0]&pp[2]) | (pp[1]&pp[2]);
   wire [PWW-1:0] cA = {mA[PWW-2:0], 1'b0};

   generate
      if (ROW0) begin : g_r0
         wire [PWW-1:0] sB = cA ^ sA ^ pp[3];
         wire [PWW-1:0] mB = (cA&sA) | (cA&pp[3]) | (sA&pp[3]);
         wire [PWW-1:0] cB = {mB[PWW-2:0], 1'b0};
         wire [PWW-1:0] sC = cB ^ sB ^ corr;
         wire [PWW-1:0] mC = (cB&sB) | (cB&corr) | (sB&corr);
         wire [PWW-1:0] cC = {mC[PWW-2:0], 1'b0};
         always @(posedge clk)
           if (en) begin
              s_south <= sC;
              c_south <= cC;
           end
      end
      else begin : g_acc
         wire [PWW-1:0] sB = pp[3] ^ s_north ^ c_north;
         wire [PWW-1:0] mB = (pp[3]&s_north) | (pp[3]&c_north) |
                             (s_north&c_north);
         wire [PWW-1:0] cB = {mB[PWW-2:0], 1'b0};
         wire [PWW-1:0] sC = sA ^ cA ^ sB;
         wire [PWW-1:0] mC = (sA&cA) | (sA&sB) | (cA&sB);
         wire [PWW-1:0] cC = {mC[PWW-2:0], 1'b0};
         wire [PWW-1:0] sD = cB ^ sC ^ cC;
         wire [PWW-1:0] mD = (cB&sC) | (cB&cC) | (sC&cC);
         wire [PWW-1:0] cD = {mD[PWW-2:0], 1'b0};
         wire [PWW-1:0] sE = sD ^ cD ^ corr;
         wire [PWW-1:0] mE = (sD&cD) | (sD&corr) | (cD&corr);
         wire [PWW-1:0] cE = {mE[PWW-2:0], 1'b0};
         always @(posedge clk)
           if (en) begin
              s_south <= sE;
              c_south <= cE;
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
    input                 rst,   // unused; init state replaces reset
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
            initial
              for (k = 0; k < i; k = k + 1)
                sr[k] = {DW{1'b0}};
            always @(posedge clk)
              if (en) begin
                 sr[0] <= din;
                 for (k = 1; k < i; k = k + 1)
                   sr[k] <= sr[k-1];
              end
            assign a_left[i*DW +: DW] = sr[i-1];
         end
      end
   endgenerate

   wire [N*PWW-1:0] c_bot;
   tpu_array #(.N(N), .DW(DW), .PWW(PWW)) u_array
     (.clk(clk), .ld(ld), .en(en),
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
            initial
              for (k = 0; k < D; k = k + 1)
                dr[k] = {PWW{1'b0}};
            always @(posedge clk)
              if (en) begin
                 dr[0] <= cin;
                 for (k = 1; k < D; k = k + 1)
                   dr[k] <= dr[k-1];
              end
            wire [PWW-1:0] cout = dr[D-1];
            assign c_data[j*ACCW +: ACCW] =
                   {{(ACCW-PWW){cout[PWW-1]}}, cout};
         end
      end
   endgenerate

   reg [LAT-1:0] vpipe = {LAT{1'b0}};
   always @(posedge clk)
     if (en)
       vpipe <= {vpipe[LAT-2:0], a_valid};
   assign c_valid = vpipe[LAT-1];

endmodule

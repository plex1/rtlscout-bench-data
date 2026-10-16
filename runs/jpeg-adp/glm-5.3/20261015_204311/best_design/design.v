//-----------------------------------------------------------------
// jpeg_idct_2d — optimized 8x8 2-D IDCT (row pass / corner-turn / column pass)
//-----------------------------------------------------------------
module jpeg_idct_2d
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 15:0]  inport_data0_i
    ,input  [ 15:0]  inport_data1_i
    ,input  [ 15:0]  inport_data2_i
    ,input  [ 15:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);
localparam XW = 20;   // x-pass output width (|x_out| < 2^19 provably)
localparam YW = 18;   // y-pass output width (|y_out| < 2^17)

wire              x_valid_w;
wire [XW-1:0]     x_data_w;
wire [ 5:0]       x_idx_w;
wire              x_accept_w;

jpeg_idct_core #(.OUT_SHIFT(11), .IN_WIDTH(16), .OUT_WIDTH(XW), .COL_MODE(1'b0))
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(inport_data0_i)
    ,.inport_data1_i(inport_data1_i)
    ,.inport_data2_i(inport_data2_i)
    ,.inport_data3_i(inport_data3_i)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

wire              t_valid_w;
wire [XW-1:0]     t_data0_w;
wire [XW-1:0]     t_data1_w;
wire [XW-1:0]     t_data2_w;
wire [XW-1:0]     t_data3_w;
wire [ 2:0]       t_idx_w;

jpeg_idct_transpose_wrap #(.W(XW))
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o(x_accept_w)
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
    ,.outport_ready_i(1'b1)
);

wire              y_valid_w;
wire [YW-1:0]     y_data_w;

jpeg_idct_core #(.OUT_SHIFT(15), .IN_WIDTH(XW), .OUT_WIDTH(YW), .COL_MODE(1'b1))
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(t_valid_w)
    ,.inport_data0_i(t_data0_w)
    ,.inport_data1_i(t_data1_w)
    ,.inport_data2_i(t_data2_w)
    ,.inport_data3_i(t_data3_w)
    ,.inport_idx_i(t_idx_w)
    ,.outport_valid_o(y_valid_w)
    ,.outport_data_o(y_data_w)
    ,.outport_idx_o(outport_idx_o)
);

assign outport_valid_o = y_valid_w;
assign outport_data_o  = {{(32-YW){y_data_w[YW-1]}}, y_data_w};

endmodule
//-----------------------------------------------------------------
// Two-cycle constant-multiplier slot:  a*b split into 4 partial products
//   a = ah*1024 + al ,  b = bh*64 + bl
//   p = ah*bh<<16 + ah*bl<<10 + al*bh<<6 + al*bl   (exact, then mod 2^32)
//-----------------------------------------------------------------
module idct_mul_slot
#(parameter AW = 20)
(
     input             clk
    ,input             rst
    ,input  signed [AW-1:0] a_i
    ,input  signed [12:0]   b_i
    ,output reg  signed [31:0] p_o
);
    localparam HIB = AW - 10;
    localparam PW  = AW + 13;
    localparam MAXW = (PW > 32) ? PW : 32;

    wire signed [HIB-1:0] ah = a_i[AW-1:10];
    wire signed [10:0]    al = {1'b0, a_i[9:0]};
    wire signed [6:0]     bh = b_i[12:6];
    wire signed [6:0]     bl = {1'b0, b_i[5:0]};

    wire signed [HIB+6:0] m_hh = ah * bh;
    wire signed [HIB+6:0] m_hl = ah * bl;
    wire signed [17:0]    m_lh = al * bh;
    wire signed [17:0]    m_ll = al * bl;

    reg  signed [HIB+6:0] m_hh_q, m_hl_q;
    reg  signed [17:0]    m_lh_q, m_ll_q;

    always @(posedge clk)
    if (rst)
    begin
        m_hh_q <= 0; m_hl_q <= 0; m_lh_q <= 0; m_ll_q <= 0;
    end
    else
    begin
        m_hh_q <= m_hh; m_hl_q <= m_hl; m_lh_q <= m_lh; m_ll_q <= m_ll;
    end

    wire signed [PW-1:0] s_hh = m_hh_q;
    wire signed [PW-1:0] s_hl = m_hl_q;
    wire signed [PW-1:0] s_lh = m_lh_q;
    wire signed [PW-1:0] s_ll = m_ll_q;
    wire signed [MAXW-1:0] p_sgn = (s_hh <<< 16) + (s_hl <<< 10) + (s_lh <<< 6) + s_ll;
    wire signed [31:0] p_val = p_sgn[31:0];

    always @(posedge clk)
    if (rst)
        p_o <= 0;
    else
        p_o <= p_val;
endmodule
//-----------------------------------------------------------------
// x * 181 mod 2^32  (carry-save tree over the 5 set bits of 181)
//-----------------------------------------------------------------
function [31:0] mul181c;
    input [31:0] x;
    reg [31:0] a, b, c, d, e;
    reg [31:0] s1, c1, s2, c2, s3, c3;
    begin
        a = x;       b = x << 7; c = x << 5; d = x << 4; e = x << 2;
        s1 = a ^ b ^ c;            c1 = ((a & b) | ((a ^ b) & c)) << 1;
        s2 = s1 ^ c1 ^ d;          c2 = ((s1 & c1) | ((s1 ^ c1) & d)) << 1;
        s3 = s2 ^ c2 ^ e;          c3 = ((s2 & c2) | ((s2 ^ c2) & e)) << 1;
        mul181c = s3 + c3;
    end
endfunction
//-----------------------------------------------------------------
// Time-shared 1-D IDCT pass
//-----------------------------------------------------------------
module jpeg_idct_core
#(
     parameter OUT_SHIFT  = 11
    ,parameter IN_WIDTH   = 16
    ,parameter OUT_WIDTH  = 20
    ,parameter COL_MODE   = 0
)
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [IN_WIDTH-1:0] inport_data0_i
    ,input  [IN_WIDTH-1:0] inport_data1_i
    ,input  [IN_WIDTH-1:0] inport_data2_i
    ,input  [IN_WIDTH-1:0] inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [OUT_WIDTH-1:0] outport_data_o
    ,output [  5:0]  outport_idx_o
);
localparam AW = IN_WIDTH + 3;   // room for i0 = a0 +/- a4 without overflow

localparam signed [12:0] C1_16 = 13'sd4017;
localparam signed [12:0] C2_16 = 13'sd3784;
localparam signed [12:0] C3_16 = 13'sd3406;
localparam signed [12:0] C4_16 = 13'sd2896;
localparam signed [12:0] C5_16 = 13'sd2276;
localparam signed [12:0] C6_16 = 13'sd1567;
localparam signed [12:0] C7_16 = 13'sd799;

wire signed [AW-1:0] block_in_0_1 = $signed(inport_data0_i);
wire signed [AW-1:0] block_in_2_3 = $signed(inport_data1_i);
wire signed [AW-1:0] block_in_4_5 = $signed(inport_data2_i);
wire signed [AW-1:0] block_in_6_7 = $signed(inport_data3_i);

reg signed [AW-1:0] i0;
reg signed [AW-1:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
reg signed [12:0]   mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

always @(posedge clk)
if (rst)
begin
    i0 <= 0;
    mul0_a <= 0; mul1_a <= 0; mul2_a <= 0; mul3_a <= 0; mul4_a <= 0;
    mul0_b <= 0; mul1_b <= 0; mul2_b <= 0; mul3_b <= 0; mul4_b <= 0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= C2_16;
        mul1_a <= block_in_6_7; mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1; mul0_b <= C1_16;
        mul1_a <= block_in_6_7; mul1_b <= C7_16;
        mul2_a <= block_in_4_5; mul2_b <= C5_16;
        mul3_a <= block_in_2_3; mul3_b <= C3_16;
        mul4_a <= i0;           mul4_b <= C4_16;
    end
    3'd2:
    begin
        i0 <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3; mul0_b <= C6_16;
        mul1_a <= block_in_6_7; mul1_b <= C2_16;
        mul4_a <= i0;           mul4_b <= C4_16;
    end
    default: ;
    endcase
end

wire signed [31:0] mul0, mul1, mul2, mul3, mul4;

idct_mul_slot #(.AW(AW)) u_ms0 (.clk(clk), .rst(rst), .a_i(mul0_a), .b_i(mul0_b), .p_o(mul0));
idct_mul_slot #(.AW(AW)) u_ms1 (.clk(clk), .rst(rst), .a_i(mul1_a), .b_i(mul1_b), .p_o(mul1));
idct_mul_slot #(.AW(AW)) u_ms2 (.clk(clk), .rst(rst), .a_i(mul2_a), .b_i(mul2_b), .p_o(mul2));
idct_mul_slot #(.AW(AW)) u_ms3 (.clk(clk), .rst(rst), .a_i(mul3_a), .b_i(mul3_b), .p_o(mul3));
idct_mul_slot #(.AW(AW)) u_ms4 (.clk(clk), .rst(rst), .a_i(mul4_a), .b_i(mul4_b), .p_o(mul4));

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
always @(posedge clk) if (rst) begin out_stg0_valid_q <= 0; out_stg0_idx_q <= 0; end
else begin out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i; end

reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
always @(posedge clk) if (rst) begin out_stg1_valid_q <= 0; out_stg1_idx_q <= 0; end
else begin out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q; end

reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
always @(posedge clk) if (rst) begin out_stg2_valid_q <= 0; out_stg2_idx_q <= 0; end
else begin out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q; end

reg signed [31:0] o_s5, o_s6, o_s7, o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] m181u, m181v;   // o_t6*181, o_t5*181  (mod 2^32)
reg signed [31:0] o_q5, o_q6;     // (o_t6 -/+ o_t5)*181  (mod 2^32)

always @(posedge clk)
if (rst)
begin
    o_s5<=0; o_s6<=0; o_s7<=0; o_t0<=0; o_t1<=0; o_t2<=0; o_t3<=0; o_t4<=0;
    o_t5<=0; o_t6<=0; m181u<=0; m181v<=0; o_q5<=0; o_q6<=0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0: begin o_t3 <= mul0 + mul1; end
    3'd1: begin o_s7 <= mul0 + mul1; o_s6 <= mul2 + mul3; o_t0 <= mul4; end
    3'd2: begin o_t0 <= o_t0 + o_t3; o_t3 <= o_t0 - o_t3; o_t7 <= o_s6 + o_s7; end
    3'd3: begin o_t4 <= (mul0 - mul1) + (mul2 - mul3); end
    3'd4: begin o_t0 <= mul0 - mul1; o_s5 <= mul2 - mul3; end
    3'd5: begin o_t3 <= mul0 - mul1; o_t4 <= mul4;
                 o_t5 <= o_t0 - o_s5; o_t6 <= o_s7 - o_s6; end
    3'd6: begin o_t1 <= o_t4 + o_t3; o_t2 <= o_t4 - o_t3;
                 m181u <= mul181c(o_t6); m181v <= mul181c(o_t5); end
    default: begin
        o_q5 <= m181u - m181v;   // (o_t6 - o_t5) * 181
        o_q6 <= m181u + m181v;   // (o_t5 + o_t6) * 181
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;
always @(posedge clk) if (rst) begin out_stg3_valid_q <= 0; out_stg3_idx_q <= 0; end
else begin out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q; end

// q>>8 with truncation toward zero:  (q>>>8) + (q<0 && q[7:0]!=0)
wire signed [31:0] q5h = o_q5 >>> 8;
wire               q5i = o_q5[31] & (|o_q5[7:0]);
wire signed [31:0] q6h = o_q6 >>> 8;
wire               q6i = o_q6[31] & (|o_q6[7:0]);

reg signed [OUT_WIDTH-1:0] block_out[0:7];
reg signed [OUT_WIDTH-1:0] block_out_tmp;
always @(posedge clk)
if (rst)
begin
    block_out[0]<=0; block_out[1]<=0; block_out[2]<=0; block_out[3]<=0; block_out[4]<=0;
    block_out[5]<=0; block_out[6]<=0; block_out[7]<=0; block_out_tmp<=0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end
    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;
    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + q5h + q5i) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - q5h - q5i) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + q6h + q6i) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - q6h - q6i) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;
always @(posedge clk) if (rst) valid_q <= 8'b0;
else if (img_start_i) valid_q <= 8'b0;
else valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;
always @(posedge clk) if (rst) ptr_q <= 6'd0;
else if (img_start_i) ptr_q <= 6'd0;
else if (outport_valid_o) ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
generate
if (COL_MODE) begin : gen_col
    assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
end else begin : gen_row
    assign outport_idx_o = ptr_q;
end
endgenerate
endmodule
//-----------------------------------------------------------------
// Corner-turn double buffer (4 banks, one per output lane)
//-----------------------------------------------------------------
module jpeg_idct_transpose_wrap
#(parameter W = 20)
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [W-1:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [W-1:0]  outport_data0_o
    ,output [W-1:0]  outport_data1_o
    ,output [W-1:0]  outport_data2_o
    ,output [W-1:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);
reg        blk_wr_q, blk_rd_q;
reg [1:0]  block_ready_q;
reg [5:0]  rd_idx_q;
reg [W-1:0] rd0_q, rd1_q, rd2_q, rd3_q;
reg [3:0]  addr_nq;      // RAM address (within bank) for the next read
localparam STATE_IDLE=2'd0, STATE_SETUP=2'd1, STATE_ACTIVE=2'd2;
reg [1:0]  state_q;

wire        wr = inport_valid_i & inport_accept_o;
wire [4:0]  wr_addr = {blk_wr_q, inport_idx_i[3:0]};
wire [5:0]  rd_next = (state_q == STATE_ACTIVE) ? (rd_idx_q + 6'd1) : 6'd0;
wire [4:0]  rd_addr = {blk_rd_q, addr_nq};

reg [W-1:0] ram0 [0:31];
reg [W-1:0] ram1 [0:31];
reg [W-1:0] ram2 [0:31];
reg [W-1:0] ram3 [0:31];

always @(posedge clk)
begin
    if (wr && inport_idx_i[5:4] == 2'd0) ram0[wr_addr] <= inport_data_i;
    if (wr && inport_idx_i[5:4] == 2'd1) ram1[wr_addr] <= inport_data_i;
    if (wr && inport_idx_i[5:4] == 2'd2) ram2[wr_addr] <= inport_data_i;
    if (wr && inport_idx_i[5:4] == 2'd3) ram3[wr_addr] <= inport_data_i;
    rd0_q <= ram0[rd_addr];
    rd1_q <= ram1[rd_addr];
    rd2_q <= ram2[rd_addr];
    rd3_q <= ram3[rd_addr];
    addr_nq <= {rd_next[0] ^ rd_next[2], rd_next[5:3]};
end

assign inport_accept_o = ~block_ready_q[blk_wr_q];
assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];
assign outport_data0_o = rd0_q;
assign outport_data1_o = rd1_q;
assign outport_data2_o = rd2_q;
assign outport_data3_o = rd3_q;

always @(posedge clk)
if (rst)
begin
    blk_wr_q <= 0; blk_rd_q <= 0; block_ready_q <= 0; rd_idx_q <= 0; state_q <= STATE_IDLE;
end
else if (img_start_i)
begin
    block_ready_q <= 0; blk_wr_q <= 0; blk_rd_q <= 0; rd_idx_q <= 0; state_q <= STATE_IDLE;
end
else
begin
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o)
    begin
        block_ready_q[blk_wr_q] <= 1'b1;
        blk_wr_q <= ~blk_wr_q;
    end
    if (outport_valid_o && rd_idx_q == 6'd63)
    begin
        block_ready_q[blk_rd_q] <= 1'b0;
        blk_rd_q <= ~blk_rd_q;
    end
    case (state_q)
    STATE_IDLE: if (block_ready_q[blk_rd_q] && outport_ready_i) state_q <= STATE_SETUP;
    STATE_SETUP: state_q <= STATE_ACTIVE;
    STATE_ACTIVE: if (rd_idx_q == 6'd63) state_q <= STATE_IDLE;
    default: ;
    endcase
    if (state_q == STATE_ACTIVE) rd_idx_q <= rd_idx_q + 6'd1;
end
endmodule

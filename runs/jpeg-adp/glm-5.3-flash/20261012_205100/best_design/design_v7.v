//-----------------------------------------------------------------
// jpeg_idct_2d — full 8x8 2-D IDCT (core_jpeg): x-pass -> double-buffered
// transpose -> y-pass. v7: two-stage multipliers with per-array b-register
// copies (fanout split) + exact signed truncating (x*181)/256.
// All arithmetic bit-exact vs the original for every possible input.
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

wire          idct_x_valid_w;
wire [ 19:0]  idct_x_data_w;
wire [  5:0]  idct_x_idx_w;

jpeg_idct_x
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
    ,.outport_valid_o(idct_x_valid_w)
    ,.outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [ 19:0]  transpose_data0_w;
wire [ 19:0]  transpose_data1_w;
wire [ 19:0]  transpose_data2_w;
wire [ 19:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(idct_x_valid_w)
    ,.inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w)
    ,.inport_accept_o()
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
    ,.outport_ready_i(1'b1)
);

jpeg_idct_y
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(transpose_valid_w)
    ,.inport_data0_i(transpose_data0_w)
    ,.inport_data1_i(transpose_data1_w)
    ,.inport_data2_i(transpose_data2_w)
    ,.inport_data3_i(transpose_data3_w)
    ,.inport_idx_i(transpose_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//-----------------------------------------------------------------
// 1-D IDCT, row pass (input coefficients 16-bit signed)
//-----------------------------------------------------------------
module jpeg_idct_x
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
    ,output [ 19:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam OUT_SHIFT = 11;

localparam signed [12:0] C1_16 = 13'sd4017;
localparam signed [12:0] C2_16 = 13'sd3784;
localparam signed [12:0] C3_16 = 13'sd3406;
localparam signed [12:0] C4_16 = 13'sd2896;
localparam signed [12:0] C5_16 = 13'sd2276;
localparam signed [12:0] C6_16 = 13'sd1567;
localparam signed [12:0] C7_16 = 13'sd799;

wire signed [15:0] d0 = inport_data0_i;
wire signed [15:0] d2 = inport_data1_i;
wire signed [15:0] d4 = inport_data2_i;
wire signed [15:0] d6 = inport_data3_i;

reg signed [16:0] i0;
reg signed [15:0] mul0_a, mul1_a, mul2_a, mul3_a;
reg signed [16:0] mul4_a;
// two copies of each constant register: one feeds the lo partial-product
// array, the other the hi array (halves the fanout on each net)
reg signed [12:0] mul0_bl, mul0_bh;
reg signed [12:0] mul1_bl, mul1_bh;
reg signed [12:0] mul2_bl, mul2_bh;
reg signed [12:0] mul3_bl, mul3_bh;
reg signed [12:0] mul4_bl, mul4_bh;

always @ (posedge clk)
if (rst)
begin
    i0 <= 17'b0;
    mul0_a <= 16'b0; mul0_bl <= 13'b0; mul0_bh <= 13'b0;
    mul1_a <= 16'b0; mul1_bl <= 13'b0; mul1_bh <= 13'b0;
    mul2_a <= 16'b0; mul2_bl <= 13'b0; mul2_bh <= 13'b0;
    mul3_a <= 16'b0; mul3_bl <= 13'b0; mul3_bh <= 13'b0;
    mul4_a <= 17'b0; mul4_bl <= 13'b0; mul4_bh <= 13'b0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= d0 + d4;
        mul0_a <= d2; mul0_bl <= C2_16; mul0_bh <= C2_16;
        mul1_a <= d6; mul1_bl <= C6_16; mul1_bh <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= d0; mul0_bl <= C1_16; mul0_bh <= C1_16;
        mul1_a <= d6; mul1_bl <= C7_16; mul1_bh <= C7_16;
        mul2_a <= d4; mul2_bl <= C5_16; mul2_bh <= C5_16;
        mul3_a <= d2; mul3_bl <= C3_16; mul3_bh <= C3_16;
        mul4_a <= i0; mul4_bl <= C4_16; mul4_bh <= C4_16;
    end
    3'd2:
    begin
        i0     <= d0 - d4;
    end
    3'd3:
    begin
        mul0_a <= d0; mul0_bl <= C7_16; mul0_bh <= C7_16;
        mul1_a <= d6; mul1_bl <= C1_16; mul1_bh <= C1_16;
        mul2_a <= d4; mul2_bl <= C3_16; mul2_bh <= C3_16;
        mul3_a <= d2; mul3_bl <= C5_16; mul3_bh <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= d0; mul0_bl <= C7_16; mul0_bh <= C7_16;
        mul1_a <= d6; mul1_bl <= C1_16; mul1_bh <= C1_16;
        mul2_a <= d4; mul2_bl <= C3_16; mul2_bh <= C3_16;
        mul3_a <= d2; mul3_bl <= C5_16; mul3_bh <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= d2; mul0_bl <= C6_16; mul0_bh <= C6_16;
        mul1_a <= d6; mul1_bl <= C2_16; mul1_bh <= C2_16;
        mul4_a <= i0; mul4_bl <= C4_16; mul4_bh <= C4_16;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

// 2-stage split of each product a*b (a = 16/17-bit signed, b = 13-bit const):
// stage 1 registers lo = a[7:0]*b and hi = a[15:8]*b (hi signed);
// stage 2 forms lo + (hi << 8) = a*b exactly (no wrap at these widths).
wire signed [20:0] p0lo = $signed({1'b0, mul0_a[7:0]}) * mul0_bl;
wire signed [20:0] p0hi = $signed(mul0_a[15:8]) * mul0_bh;
wire signed [20:0] p1lo = $signed({1'b0, mul1_a[7:0]}) * mul1_bl;
wire signed [20:0] p1hi = $signed(mul1_a[15:8]) * mul1_bh;
wire signed [20:0] p2lo = $signed({1'b0, mul2_a[7:0]}) * mul2_bl;
wire signed [20:0] p2hi = $signed(mul2_a[15:8]) * mul2_bh;
wire signed [20:0] p3lo = $signed({1'b0, mul3_a[7:0]}) * mul3_bl;
wire signed [20:0] p3hi = $signed(mul3_a[15:8]) * mul3_bh;
wire signed [20:0] p4lo = $signed({1'b0, mul4_a[7:0]}) * mul4_bl;
wire signed [20:0] p4hi = $signed(mul4_a[16:8]) * mul4_bh;

reg signed [20:0] lo0_q, lo1_q, lo2_q, lo3_q, lo4_q;
reg signed [20:0] hi0_q, hi1_q, hi2_q, hi3_q, hi4_q;

reg signed [27:0] mul0, mul1, mul2, mul3;
reg signed [28:0] mul4;

always @ (posedge clk)
if (rst)
begin
    lo0_q <= 0; lo1_q <= 0; lo2_q <= 0; lo3_q <= 0; lo4_q <= 0;
    hi0_q <= 0; hi1_q <= 0; hi2_q <= 0; hi3_q <= 0; hi4_q <= 0;
    mul0 <= 0; mul1 <= 0; mul2 <= 0; mul3 <= 0; mul4 <= 0;
end
else
begin
    lo0_q <= p0lo; hi0_q <= p0hi;
    lo1_q <= p1lo; hi1_q <= p1hi;
    lo2_q <= p2lo; hi2_q <= p2hi;
    lo3_q <= p3lo; hi3_q <= p3hi;
    lo4_q <= p4lo; hi4_q <= p4hi;
    mul0 <= lo0_q + (hi0_q <<< 8);
    mul1 <= lo1_q + (hi1_q <<< 8);
    mul2 <= lo2_q + (hi2_q <<< 8);
    mul3 <= lo3_q + (hi3_q <<< 8);
    mul4 <= lo4_q + (hi4_q <<< 8);
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q;
end

// all x-pass internals bounded by 570.7M < 2^30 -> 31-bit signed
reg signed [30:0] o_s5, o_s6, o_s7;
reg signed [30:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [30:0] o_t6_5, o_t5_6;

//--- 1/sqrt(2) scaling: d = trunc((x*181 mod 2^32)/256) (signed, exact) ------
// p = (x*181) mod 2^32 ; d = (p >>> 8) + (p[31] && p[7:0]!=0)
// p[7:0] = (x[7:0]*181) mod 256 (181 odd) -> zero detect computed early.
wire signed [31:0] px5 = o_t6_5 * 181;
wire signed [31:0] px6 = o_t5_6 * 181;
wire        nzx5 = ((o_t6_5[7:0] * 8'd181) != 8'd0);
wire        nzx6 = ((o_t5_6[7:0] * 8'd181) != 8'd0);

reg signed [23:0] s5hi_q, s6hi_q;   // p >>> 8 (fits 24-bit signed)
reg signed [31:0] f5_q, f6_q;       // truncation fixup (0 or +1), signed

always @ (posedge clk)
if (rst)
begin
    s5hi_q <= 24'b0; s6hi_q <= 24'b0; f5_q <= 32'b0; f6_q <= 32'b0;
end
else if (out_stg2_idx_q == 3'd7)
begin
    s5hi_q <= px5[31:8];
    s6hi_q <= px6[31:8];
    f5_q   <= px5[31] & nzx5;
    f6_q   <= px6[31] & nzx6;
end

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0; o_t0 <= 0; o_t1 <= 0; o_t2 <= 0;
    o_t3 <= 0; o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0: begin o_t3 <= mul0 + mul1; end
    3'd1: begin o_s7 <= mul0 + mul1; o_s6 <= mul2 + mul3; o_t0 <= mul4; end
    3'd2: begin o_t0 <= o_t0 + o_t3; o_t3 <= o_t0 - o_t3; o_t7 <= o_s6 + o_s7; end
    3'd3: begin o_t4 <= (mul0 - mul1) + (mul2 - mul3); end
    3'd4: begin o_t0 <= mul0 - mul1; o_s5 <= mul2 - mul3; end
    3'd5: begin o_t3 <= mul0 - mul1; o_t4 <= mul4; o_t5 <= o_t0 - o_s5; o_t6 <= o_s7 - o_s6; end
    3'd6: begin o_t1 <= o_t4 + o_t3; o_t2 <= o_t4 - o_t3; o_t6_5 <= o_t6 - o_t5; o_t5_6 <= o_t5 + o_t6; end
    default: ;
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
end
else
begin
    out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q;
end

// |out| <= 346243 < 2^19 -> 20-bit exact
reg signed [19:0] block_out [0:7];
reg signed [19:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end
else if (out_stg3_valid_q)
begin
    /* verilator lint_off WIDTH */
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
        block_out[2] <= ((o_t2 + s5hi_q + f5_q) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - s5hi_q - f5_q) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + s6hi_q + f6_q) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - s6hi_q - f6_q) >>> OUT_SHIFT);
    end
    /* verilator lint_on WIDTH */
end

reg [7:0] valid_q;
reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = ptr_q;

endmodule

//-----------------------------------------------------------------
// Double-buffered corner turn (4 lanes x 2 halves x 16 entries)
// Read path: one-hot decode + AND-OR tree (low fanout, shallow).
//-----------------------------------------------------------------
module jt_ram
#(
     parameter W = 20
)
(
     input           clk
    ,input           wr_en
    ,input           block_wr
    ,input  [  3:0]  wr_addr
    ,input  [W-1:0]  wr_data
    ,input           block_rd
    ,input  [  3:0]  rd_addr
    ,output reg [W-1:0]  rd_data
);

reg [W-1:0] mem0 [0:15];
reg [W-1:0] mem1 [0:15];

// one-hot decode of the read address
reg [15:0] oh;
integer i;
always @ *
begin
    for (i = 0; i < 16; i = i + 1)
        oh[i] = (rd_addr == i);
end

wire [15:0] oh0 = oh & {16{~block_rd}};
wire [15:0] oh1 = oh & {16{ block_rd}};

// gated words: exactly one of the 32 is (possibly) non-zero
reg [W-1:0] gt [0:31];
reg [W-1:0] t1 [0:15];
reg [W-1:0] t2 [0:7];
reg [W-1:0] t3 [0:3];
reg [W-1:0] t4 [0:1];
reg [W-1:0] acc;
integer k;
always @ *
begin
    for (k = 0; k < 16; k = k + 1)
    begin
        gt[k]    = mem0[k] & {W{oh0[k]}};
        gt[k+16] = mem1[k] & {W{oh1[k]}};
    end
    for (k = 0; k < 16; k = k + 1) t1[k] = gt[2*k] | gt[2*k+1];
    for (k = 0; k < 8;  k = k + 1) t2[k] = t1[2*k] | t1[2*k+1];
    for (k = 0; k < 4;  k = k + 1) t3[k] = t2[2*k] | t2[2*k+1];
    for (k = 0; k < 2;  k = k + 1) t4[k] = t3[2*k] | t3[2*k+1];
    acc = t4[0] | t4[1];
end

always @ (posedge clk)
begin
    if (wr_en)
    begin
        if (block_wr)
            mem1[wr_addr] <= wr_data;
        else
            mem0[wr_addr] <= wr_data;
    end

    rd_data <= acc;
end

endmodule

module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 19:0]  outport_data0_o
    ,output [ 19:0]  outport_data1_o
    ,output [ 19:0]  outport_data2_o
    ,output [ 19:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;
reg [3:0]   rda0_q, rda1_q, rda2_q, rda3_q;

wire [3:0]  wr_addr_w = inport_idx_i[3:0];

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

jt_ram
#(.W(20))
u_ram0
(
     .clk(clk)
    ,.wr_en(wr0_w)
    ,.block_wr(block_wr_q)
    ,.wr_addr(wr_addr_w)
    ,.wr_data(inport_data_i)
    ,.block_rd(block_rd_q)
    ,.rd_addr(rda0_q)
    ,.rd_data(outport_data0_o)
);

jt_ram
#(.W(20))
u_ram1
(
     .clk(clk)
    ,.wr_en(wr1_w)
    ,.block_wr(block_wr_q)
    ,.wr_addr(wr_addr_w)
    ,.wr_data(inport_data_i)
    ,.block_rd(block_rd_q)
    ,.rd_addr(rda1_q)
    ,.rd_data(outport_data1_o)
);

jt_ram
#(.W(20))
u_ram2
(
     .clk(clk)
    ,.wr_en(wr2_w)
    ,.block_wr(block_wr_q)
    ,.wr_addr(wr_addr_w)
    ,.wr_data(inport_data_i)
    ,.block_rd(block_rd_q)
    ,.rd_addr(rda2_q)
    ,.rd_data(outport_data2_o)
);

jt_ram
#(.W(20))
u_ram3
(
     .clk(clk)
    ,.wr_en(wr3_w)
    ,.block_wr(block_wr_q)
    ,.wr_addr(wr_addr_w)
    ,.wr_data(inport_data_i)
    ,.block_rd(block_rd_q)
    ,.rd_addr(rda3_q)
    ,.rd_data(outport_data3_o)
);

reg [1:0] block_ready_q;

always @ (posedge clk)
if (rst)
begin
    block_ready_q <= 2'b0;
    block_wr_q    <= 1'b0;
    block_rd_q    <= 1'b0;
end
else if (img_start_i)
begin
    block_ready_q <= 2'b0;
    block_wr_q    <= 1'b0;
    block_rd_q    <= 1'b0;
end
else
begin
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o)
    begin
        block_ready_q[block_wr_q] <= 1'b1;
        block_wr_q                <= ~block_wr_q;
    end

    if (outport_valid_o && rd_idx_q[5:0] == 6'd63)
    begin
        block_ready_q[block_rd_q] <= 1'b0;
        block_rd_q <= ~block_rd_q;
    end
end

assign inport_accept_o = ~block_ready_q[block_wr_q];

localparam STATE_W           = 2;
localparam STATE_IDLE        = 2'd0;
localparam STATE_SETUP       = 2'd1;
localparam STATE_ACTIVE      = 2'd2;

reg [STATE_W-1:0] state_q;
reg [STATE_W-1:0] next_state_r;

always @ *
begin
    next_state_r = state_q;
    case (state_q)
    STATE_IDLE:
        if (block_ready_q[block_rd_q] && outport_ready_i)
            next_state_r = STATE_SETUP;
    STATE_SETUP:
        next_state_r = STATE_ACTIVE;
    STATE_ACTIVE:
        if (outport_valid_o && rd_idx_q == 6'd63)
            next_state_r = STATE_IDLE;
    default: ;
    endcase
    if (img_start_i)
        next_state_r = STATE_IDLE;
end

always @ (posedge clk)
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @ (posedge clk)
if (rst)
    rd_idx_q <= 6'b0;
else if (img_start_i)
    rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE)
    rd_idx_q <= rd_idx_q + 6'd1;

// read address (identical schedule to original), replicated per RAM to
// cut select-net fanout
reg [3:0] rda_next;

always @ *
begin
    case (state_q)
    STATE_IDLE:
        rda_next = 4'd0;
    STATE_SETUP:
        rda_next = 4'd8;
    STATE_ACTIVE:
    begin
        case (rd_idx_q[2:0])
        3'd0: rda_next = rd_addr_q - 4'd8;
        3'd1: rda_next = rd_addr_q + 4'd8;
        3'd2: rda_next = rd_addr_q;
        3'd3: rda_next = rd_addr_q - 4'd8;
        3'd4: rda_next = rd_addr_q + 4'd8;
        3'd5: rda_next = rd_addr_q - 4'd8;
        3'd6: rda_next = rd_addr_q + 4'd1;
        3'd7: rda_next = rd_addr_q + 4'd8;
        endcase
    end
    default:
        rda_next = rd_addr_q;
    endcase
end

always @ (posedge clk)
if (rst)
begin
    rd_addr_q <= 4'b0;
    rda0_q <= 4'b0;
    rda1_q <= 4'b0;
    rda2_q <= 4'b0;
    rda3_q <= 4'b0;
end
else
begin
    rd_addr_q <= rda_next;
    rda0_q <= rda_next;
    rda1_q <= rda_next;
    rda2_q <= rda_next;
    rda3_q <= rda_next;
end

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

endmodule

//-----------------------------------------------------------------
// 1-D IDCT, column pass (inputs are the 20-bit x-pass results)
//-----------------------------------------------------------------
module jpeg_idct_y
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data0_i
    ,input  [ 19:0]  inport_data1_i
    ,input  [ 19:0]  inport_data2_i
    ,input  [ 19:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam OUT_SHIFT = 15;

localparam signed [12:0] C1_16 = 13'sd4017;
localparam signed [12:0] C2_16 = 13'sd3784;
localparam signed [12:0] C3_16 = 13'sd3406;
localparam signed [12:0] C4_16 = 13'sd2896;
localparam signed [12:0] C5_16 = 13'sd2276;
localparam signed [12:0] C6_16 = 13'sd1567;
localparam signed [12:0] C7_16 = 13'sd799;

wire signed [19:0] block_in_0_1 = inport_data0_i;
wire signed [19:0] block_in_2_3 = inport_data1_i;
wire signed [19:0] block_in_4_5 = inport_data2_i;
wire signed [19:0] block_in_6_7 = inport_data3_i;

reg signed [20:0] i0;
reg signed [19:0] mul0_a, mul1_a, mul2_a, mul3_a;
reg signed [20:0] mul4_a;
reg signed [12:0] mul0_bl, mul0_bh;
reg signed [12:0] mul1_bl, mul1_bh;
reg signed [12:0] mul2_bl, mul2_bh;
reg signed [12:0] mul3_bl, mul3_bh;
reg signed [12:0] mul4_bl, mul4_bh;

always @ (posedge clk)
if (rst)
begin
    i0 <= 21'b0;
    mul0_a <= 20'b0; mul0_bl <= 13'b0; mul0_bh <= 13'b0;
    mul1_a <= 20'b0; mul1_bl <= 13'b0; mul1_bh <= 13'b0;
    mul2_a <= 20'b0; mul2_bl <= 13'b0; mul2_bh <= 13'b0;
    mul3_a <= 20'b0; mul3_bl <= 13'b0; mul3_bh <= 13'b0;
    mul4_a <= 21'b0; mul4_bl <= 13'b0; mul4_bh <= 13'b0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_bl <= C2_16; mul0_bh <= C2_16;
        mul1_a <= block_in_6_7; mul1_bl <= C6_16; mul1_bh <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1; mul0_bl <= C1_16; mul0_bh <= C1_16;
        mul1_a <= block_in_6_7; mul1_bl <= C7_16; mul1_bh <= C7_16;
        mul2_a <= block_in_4_5; mul2_bl <= C5_16; mul2_bh <= C5_16;
        mul3_a <= block_in_2_3; mul3_bl <= C3_16; mul3_bh <= C3_16;
        mul4_a <= i0;           mul4_bl <= C4_16; mul4_bh <= C4_16;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1; mul0_bl <= C7_16; mul0_bh <= C7_16;
        mul1_a <= block_in_6_7; mul1_bl <= C1_16; mul1_bh <= C1_16;
        mul2_a <= block_in_4_5; mul2_bl <= C3_16; mul2_bh <= C3_16;
        mul3_a <= block_in_2_3; mul3_bl <= C5_16; mul3_bh <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1; mul0_bl <= C7_16; mul0_bh <= C7_16;
        mul1_a <= block_in_6_7; mul1_bl <= C1_16; mul1_bh <= C1_16;
        mul2_a <= block_in_4_5; mul2_bl <= C3_16; mul2_bh <= C3_16;
        mul3_a <= block_in_2_3; mul3_bl <= C5_16; mul3_bh <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3; mul0_bl <= C6_16; mul0_bh <= C6_16;
        mul1_a <= block_in_6_7; mul1_bl <= C2_16; mul1_bh <= C2_16;
        mul4_a <= i0;           mul4_bl <= C4_16; mul4_bh <= C4_16;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

// 2-stage split of each product a*b, taken mod 2^32 exactly as the original:
// stage 1 registers lo = a[9:0]*b (zero-extended) and hi = a[19:10]*b (signed)
// (mul4: hi = a[20:10]); stage 2 forms lo + (hi << 10) in 32-bit wrapping
// context — identical low-32 bits to the original 32x32 multiply.
wire signed [22:0] p0lo = $signed({1'b0, mul0_a[9:0]}) * mul0_bl;
wire signed [22:0] p0hi = $signed(mul0_a[19:10]) * mul0_bh;
wire signed [22:0] p1lo = $signed({1'b0, mul1_a[9:0]}) * mul1_bl;
wire signed [22:0] p1hi = $signed(mul1_a[19:10]) * mul1_bh;
wire signed [22:0] p2lo = $signed({1'b0, mul2_a[9:0]}) * mul2_bl;
wire signed [22:0] p2hi = $signed(mul2_a[19:10]) * mul2_bh;
wire signed [22:0] p3lo = $signed({1'b0, mul3_a[9:0]}) * mul3_bl;
wire signed [22:0] p3hi = $signed(mul3_a[19:10]) * mul3_bh;
wire signed [22:0] p4lo = $signed({1'b0, mul4_a[9:0]}) * mul4_bl;
wire signed [22:0] p4hi = $signed(mul4_a[20:10]) * mul4_bh;

reg signed [22:0] lo0_q, lo1_q, lo2_q, lo3_q, lo4_q;
reg signed [22:0] hi0_q, hi1_q, hi2_q, hi3_q;
reg signed [22:0] hi4_q;

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

always @ (posedge clk)
if (rst)
begin
    lo0_q <= 0; lo1_q <= 0; lo2_q <= 0; lo3_q <= 0; lo4_q <= 0;
    hi0_q <= 0; hi1_q <= 0; hi2_q <= 0; hi3_q <= 0; hi4_q <= 0;
    mul0 <= 0; mul1 <= 0; mul2 <= 0; mul3 <= 0; mul4 <= 0;
end
else
begin
    lo0_q <= p0lo; hi0_q <= p0hi;
    lo1_q <= p1lo; hi1_q <= p1hi;
    lo2_q <= p2lo; hi2_q <= p2hi;
    lo3_q <= p3lo; hi3_q <= p3hi;
    lo4_q <= p4lo; hi4_q <= p4hi;
    /* verilator lint_off WIDTH */
    mul0 <= lo0_q + (hi0_q <<< 10);
    mul1 <= lo1_q + (hi1_q <<< 10);
    mul2 <= lo2_q + (hi2_q <<< 10);
    mul3 <= lo3_q + (hi3_q <<< 10);
    mul4 <= lo4_q + (hi4_q <<< 10);
    /* verilator lint_on WIDTH */
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q;
end

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

//--- 1/sqrt(2) scaling: d = trunc((x*181 mod 2^32)/256) (signed, exact) ------
// p = (x*181) mod 2^32 ; d = (p >>> 8) + (p[31] && p[7:0]!=0)
// p[7:0] = (x[7:0]*181) mod 256 (181 odd) -> zero detect computed early.
wire signed [31:0] p5 = o_t6_5 * 181;
wire signed [31:0] p6 = o_t5_6 * 181;
wire        nz5 = ((o_t6_5[7:0] * 8'd181) != 8'd0);
wire        nz6 = ((o_t5_6[7:0] * 8'd181) != 8'd0);

reg signed [23:0] s5hi_q, s6hi_q;   // p >>> 8 (fits 24-bit signed)
reg signed [31:0] f5_q, f6_q;       // truncation fixup (0 or +1), signed

always @ (posedge clk)
if (rst)
begin
    s5hi_q <= 24'b0; s6hi_q <= 24'b0; f5_q <= 32'b0; f6_q <= 32'b0;
end
else if (out_stg2_idx_q == 3'd7)
begin
    s5hi_q <= p5[31:8];
    s6hi_q <= p6[31:8];
    f5_q   <= p5[31] & nz5;
    f6_q   <= p6[31] & nz6;
end

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0; o_t0 <= 0; o_t1 <= 0; o_t2 <= 0;
    o_t3 <= 0; o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0: begin o_t3 <= mul0 + mul1; end
    3'd1: begin o_s7 <= mul0 + mul1; o_s6 <= mul2 + mul3; o_t0 <= mul4; end
    3'd2: begin o_t0 <= o_t0 + o_t3; o_t3 <= o_t0 - o_t3; o_t7 <= o_s6 + o_s7; end
    3'd3: begin o_t4 <= (mul0 - mul1) + (mul2 - mul3); end
    3'd4: begin o_t0 <= mul0 - mul1; o_s5 <= mul2 - mul3; end
    3'd5: begin o_t3 <= mul0 - mul1; o_t4 <= mul4; o_t5 <= o_t0 - o_s5; o_t6 <= o_s7 - o_s6; end
    3'd6: begin o_t1 <= o_t4 + o_t3; o_t2 <= o_t4 - o_t3; o_t6_5 <= o_t6 - o_t5; o_t5_6 <= o_t5 + o_t6; end
    default: ;
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
end
else
begin
    out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q;
end

// (wrapped 32-bit sum) >> 15 always fits 17-bit signed
reg signed [16:0] block_out [0:7];
reg signed [16:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end
else if (out_stg3_valid_q)
begin
    /* verilator lint_off WIDTH */
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
        block_out[2] <= ((o_t2 + s5hi_q + f5_q) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - s5hi_q - f5_q) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + s6hi_q + f6_q) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - s6hi_q - f6_q) >>> OUT_SHIFT);
    end
    /* verilator lint_on WIDTH */
end

reg [7:0] valid_q;
reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

wire signed [16:0] out_sel = block_out[ptr_q[2:0]];

assign outport_valid_o = valid_q[6];
assign outport_data_o  = {{15{out_sel[16]}}, out_sel};
// outport_idx = {k[2:0], k[5:3]} : the un-transpose permutation
assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

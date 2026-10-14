//-----------------------------------------------------------------
// jpeg_idct_2d - full 8x8 2-D IDCT of the baseline JPEG decoder
// (ultraembedded/core_jpeg, Apache 2.0): row pass (OUT_SHIFT=11) ->
// double-buffered transpose -> column pass (OUT_SHIFT=15), composed
// exactly as in the original jpeg_idct.v (downstream ready tied high,
// row pass accepts no backpressure).
//
// Pipeline structure, case schedules, FSMs and arithmetic are kept
// bit-exact with the original RTL.  Optimizations (all provably value
// preserving for any 16-bit signed input):
//   * multiplier operand registers narrowed (16/17 x 13 for the x pass,
//     21/22 x 13 for the y pass) instead of 32 x 32; products are still
//     evaluated in the 32-bit signed assignment context so the low-32
//     wrap-around behaviour is identical,
//   * cosine constants kept in 13-bit registers,
//   * transpose RAM storage narrowed to 21 bits
//     (|x-pass output| <= 709424128 >>> 11 = 346350 < 2^19),
//   * y-pass output buffer narrowed to 17 bits
//     (result = 32-bit sum >>> 15 lies in [-65536, 65535]).
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

localparam XW = 21;   // narrowed x-pass result width (|x| <= 346350 < 2^20)
localparam YW = 17;   // narrowed y-pass result width (|y| < 2^16)

//-----------------------------------------------------------------
// Row pass (jpeg_idct_x equivalent, OUT_SHIFT = 11)
//-----------------------------------------------------------------
wire          x_valid_w;
wire [XW-1:0] x_data_w;
wire [  5:0]  x_idx_w;

jpeg_idct_1d #(
     .OUT_SHIFT(11)
    ,.IN_W(16)
    ,.OUT_W(XW)
    ,.IDX_SWAP(1'b0)
) u_idct_x (
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

//-----------------------------------------------------------------
// Double-buffered corner turn
//-----------------------------------------------------------------
wire          t_valid_w;
wire [XW-1:0] t_data0_w;
wire [XW-1:0] t_data1_w;
wire [XW-1:0] t_data2_w;
wire [XW-1:0] t_data3_w;
wire [  2:0]  t_idx_w;

jpeg_idct_transpose #(
     .W(XW)
) u_transpose (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o()
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
    ,.outport_ready_i(1'b1)
);

//-----------------------------------------------------------------
// Column pass (jpeg_idct_y equivalent, OUT_SHIFT = 15)
//-----------------------------------------------------------------
wire          y_valid_w;
wire [YW-1:0] y_data_w;
wire [  5:0]  y_idx_w;

jpeg_idct_1d #(
     .OUT_SHIFT(15)
    ,.IN_W(XW)
    ,.OUT_W(YW)
    ,.IDX_SWAP(1'b1)
) u_idct_y (
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
    ,.outport_idx_o(y_idx_w)
);

assign outport_valid_o = y_valid_w;
assign outport_idx_o   = y_idx_w;
assign outport_data_o  = {{(32-YW){y_data_w[YW-1]}}, y_data_w};

endmodule

//-----------------------------------------------------------------
// One-dimensional IDCT pass (factorized AAN style, 8-beat schedule,
// 5 time-shared multipliers).  Shared by the row and column pass.
//   IDX_SWAP = 0 : output index = stream counter      (x pass)
//   IDX_SWAP = 1 : output index = {cnt[2:0],cnt[5:3]} (y pass)
//-----------------------------------------------------------------
module jpeg_idct_1d #(
     parameter OUT_SHIFT = 11
    ,parameter IN_W      = 16
    ,parameter OUT_W     = 21
    ,parameter IDX_SWAP  = 1'b0
)(
     input                 clk
    ,input                 rst
    ,input                 img_start_i
    ,input                 img_end_i
    ,input                 inport_valid_i
    ,input  [ IN_W-1:0]    inport_data0_i
    ,input  [ IN_W-1:0]    inport_data1_i
    ,input  [ IN_W-1:0]    inport_data2_i
    ,input  [ IN_W-1:0]    inport_data3_i
    ,input  [  2:0]        inport_idx_i
    ,output                outport_valid_o
    ,output [ OUT_W-1:0]   outport_data_o
    ,output [  5:0]        outport_idx_o
);

// cosine constants x4096 (all < 4096 => fit unsigned 12 bit, which is
// what makes the unsigned-multiply re formulation below exact)
localparam [11:0] C1 = 12'd4017; // cos( pi/16) x4096
localparam [11:0] C2 = 12'd3784; // cos(2pi/16) x4096
localparam [11:0] C3 = 12'd3406; // cos(3pi/16) x4096
localparam [11:0] C4 = 12'd2896; // cos(4pi/16) x4096
localparam [11:0] C5 = 12'd2276; // cos(5pi/16) x4096
localparam [11:0] C6 = 12'd1567; // cos(6pi/16) x4096
localparam [11:0] C7 = 12'd799;  // cos(7pi/16) x4096

// operand width for the multiplier input stage: holds the sign
// extended input data as well as i0 = d0 + d4 (IN_W+1 bits)
localparam MA_W = IN_W + 1;

wire signed [MA_W-1:0] block_in_0_1 = {{(MA_W-IN_W){inport_data0_i[IN_W-1]}}, inport_data0_i};
wire signed [MA_W-1:0] block_in_2_3 = {{(MA_W-IN_W){inport_data1_i[IN_W-1]}}, inport_data1_i};
wire signed [MA_W-1:0] block_in_4_5 = {{(MA_W-IN_W){inport_data2_i[IN_W-1]}}, inport_data2_i};
wire signed [MA_W-1:0] block_in_6_7 = {{(MA_W-IN_W){inport_data3_i[IN_W-1]}}, inport_data3_i};

//-----------------------------------------------------------------
// IDCT
//-----------------------------------------------------------------
reg signed [MA_W-1:0] i0;
reg signed [MA_W-1:0] mul0_a;
reg        [11:0]     mul0_b;
reg signed [MA_W-1:0] mul1_a;
reg        [11:0]     mul1_b;
reg signed [MA_W-1:0] mul2_a;
reg        [11:0]     mul2_b;
reg signed [MA_W-1:0] mul3_a;
reg        [11:0]     mul3_b;
reg signed [MA_W-1:0] mul4_a;
reg        [11:0]     mul4_b;

always @ (posedge clk )
if (rst)
begin
    i0     <= {MA_W{1'b0}};
    mul0_a <= {MA_W{1'b0}};
    mul0_b <= 12'b0;
    mul1_a <= {MA_W{1'b0}};
    mul1_b <= 12'b0;
    mul2_a <= {MA_W{1'b0}};
    mul2_b <= 12'b0;
    mul3_a <= {MA_W{1'b0}};
    mul3_b <= 12'b0;
    mul4_a <= {MA_W{1'b0}};
    mul4_b <= 12'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul0_b <= C2;
        mul1_a <= block_in_6_7;
        mul1_b <= C6;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C1;
        mul1_a <= block_in_6_7;
        mul1_b <= C7;
        mul2_a <= block_in_4_5;
        mul2_b <= C5;
        mul3_a <= block_in_2_3;
        mul3_b <= C3;
        mul4_a <= i0;
        mul4_b <= C4;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7;
        mul1_a <= block_in_6_7;
        mul1_b <= C1;
        mul2_a <= block_in_4_5;
        mul2_b <= C3;
        mul3_a <= block_in_2_3;
        mul3_b <= C5;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7;
        mul1_a <= block_in_6_7;
        mul1_b <= C1;
        mul2_a <= block_in_4_5;
        mul2_b <= C3;
        mul3_a <= block_in_2_3;
        mul3_b <= C5;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_b <= C6;
        mul1_a <= block_in_6_7;
        mul1_b <= C2;
        mul4_a <= i0;
        mul4_b <= C4;
    end
    default:
        ;
    endcase
end

// Products.
// The factor operand b is always one of the non-negative constants
// C1..C7 (12 bit unsigned, value < 4096), so a signed multiply of the
// sign-extended a with the zero-extended b is bit-identical to the
// original 32x32 signed wrap-around product.
reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk )
if (rst)
begin
    mul0_q <= 32'b0;
    mul1_q <= 32'b0;
    mul2_q <= 32'b0;
    mul3_q <= 32'b0;
    mul4_q <= 32'b0;
end
else
begin
    mul0_q <= $signed({{(32-MA_W){mul0_a[MA_W-1]}}, mul0_a}) * $signed({20'b0, mul0_b});
    mul1_q <= $signed({{(32-MA_W){mul1_a[MA_W-1]}}, mul1_a}) * $signed({20'b0, mul1_b});
    mul2_q <= $signed({{(32-MA_W){mul2_a[MA_W-1]}}, mul2_a}) * $signed({20'b0, mul2_b});
    mul3_q <= $signed({{(32-MA_W){mul3_a[MA_W-1]}}, mul3_a}) * $signed({20'b0, mul3_b});
    mul4_q <= $signed({{(32-MA_W){mul4_a[MA_W-1]}}, mul4_a}) * $signed({20'b0, mul4_b});
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk )
if (rst)
begin
    mul0 <= 32'b0;
    mul1 <= 32'b0;
    mul2 <= 32'b0;
    mul3 <= 32'b0;
    mul4 <= 32'b0;
end
else
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
end

reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
end
else
begin
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
end

reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
end
else
begin
    out_stg2_valid_q <= out_stg1_valid_q;
    out_stg2_idx_q   <= out_stg1_idx_q;
end

reg signed [31:0] o_s5;
reg signed [31:0] o_s6;
reg signed [31:0] o_s7;
reg signed [31:0] o_t0;
reg signed [31:0] o_t1;
reg signed [31:0] o_t2;
reg signed [31:0] o_t3;
reg signed [31:0] o_t4;
reg signed [31:0] o_t5;
reg signed [31:0] o_t6;
reg signed [31:0] o_t7;
reg signed [31:0] o_t6_5;
reg signed [31:0] o_t5_6;

// signed division of a 32-bit value by 256, truncating toward zero
function signed [31:0] trunc_div256;
    input signed [31:0] w;
    reg [23:0] s;
    begin
        s = w[31:8] + {23'b0, (w[31] & (|w[7:0]))};
        trunc_div256 = {{8{s[23]}}, s};
    end
endfunction

always @ (posedge clk )
if (rst)
begin
    o_s5   <= 32'b0;
    o_s6   <= 32'b0;
    o_s7   <= 32'b0;
    o_t0   <= 32'b0;
    o_t1   <= 32'b0;
    o_t2   <= 32'b0;
    o_t3   <= 32'b0;
    o_t4   <= 32'b0;
    o_t5   <= 32'b0;
    o_t6   <= 32'b0;
    o_t7   <= 32'b0;
    o_t6_5 <= 32'b0;
    o_t5_6 <= 32'b0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0:
    begin
        o_t3 <= mul0 + mul1; // s3
    end
    3'd1:
    begin
        o_s7 <= mul0 + mul1;
        o_s6 <= mul2 + mul3;
        o_t0 <= mul4;        // s0
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3; // t0
        o_t3 <= o_t0 - o_t3; // t3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3:
    begin
        o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    end
    3'd4:
    begin
        o_t0 <= mul0 - mul1; // s4
        o_s5 <= mul2 - mul3;
    end
    3'd5:
    begin
        o_t3 <= mul0 - mul1; // s2
        o_t4 <= mul4;        // s1
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5;
        o_t5_6 <= o_t5 + o_t6;
    end
    default:
    begin
        // (x * 181) / 256 : product wraps at 32 bits, signed division
        // truncating toward zero == arithmetic shift right 8 plus a
        // +1 correction for negative, non-multiple results
        o_s5 <= trunc_div256(o_t6_5 * 181); // 1/sqrt(2)
        o_s6 <= trunc_div256(o_t5_6 * 181); // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg3_valid_q <= 1'b0;
    out_stg3_idx_q   <= 3'b0;
end
else
begin
    out_stg3_valid_q <= out_stg2_valid_q;
    out_stg3_idx_q   <= out_stg2_idx_q;
end

// output rounding: arithmetic shift right, result guaranteed to fit
// in OUT_W bits (see header), truncate explicitly to keep lint clean
wire signed [31:0] r0 = (o_t0 + o_t7) >>> OUT_SHIFT;
wire signed [31:0] r1 = (o_t1 + o_s6) >>> OUT_SHIFT;
wire signed [31:0] r2 = (o_t2 + o_s5) >>> OUT_SHIFT;
wire signed [31:0] r3 = (o_t3 + o_t4) >>> OUT_SHIFT;
wire signed [31:0] r4 = (o_t3 - o_t4) >>> OUT_SHIFT;
wire signed [31:0] r5 = (o_t2 - o_s5) >>> OUT_SHIFT;
wire signed [31:0] r6 = (o_t1 - o_s6) >>> OUT_SHIFT;
wire signed [31:0] r7 = (o_t0 - o_t7) >>> OUT_SHIFT;

reg signed [OUT_W-1:0] block_out[0:7];
reg signed [OUT_W-1:0] block_out_tmp;

always @ (posedge clk )
if (rst)
begin
    block_out[0] <= {OUT_W{1'b0}};
    block_out[1] <= {OUT_W{1'b0}};
    block_out[2] <= {OUT_W{1'b0}};
    block_out[3] <= {OUT_W{1'b0}};
    block_out[4] <= {OUT_W{1'b0}};
    block_out[5] <= {OUT_W{1'b0}};
    block_out[6] <= {OUT_W{1'b0}};
    block_out[7] <= {OUT_W{1'b0}};
    block_out_tmp <= {OUT_W{1'b0}};
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= r0[OUT_W-1:0];
        block_out_tmp <= r7[OUT_W-1:0]; // block_out[7]
        block_out[3] <= r3[OUT_W-1:0];
        block_out[4] <= r4[OUT_W-1:0];
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= r2[OUT_W-1:0];
        block_out[5] <= r5[OUT_W-1:0];
        block_out[1] <= r1[OUT_W-1:0];
        block_out[6] <= r6[OUT_W-1:0];
    end
end

reg [7:0] valid_q;

always @ (posedge clk )
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk )
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = IDX_SWAP ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

endmodule

//-----------------------------------------------------------------
// Double-buffered transpose (corner turn).  Four banks, each 32
// entries x W bits; write side takes one value per accepted beat,
// read side walks the fixed column-pass schedule.
//-----------------------------------------------------------------
module jpeg_idct_transpose #(
     parameter W = 21
)(
     input               clk
    ,input               rst
    ,input               img_start_i
    ,input               img_end_i
    ,input               inport_valid_i
    ,input  [ W-1:0]     inport_data_i
    ,input  [  5:0]      inport_idx_i
    ,input               outport_ready_i
    ,output              inport_accept_o
    ,output              outport_valid_o
    ,output [ W-1:0]     outport_data0_o
    ,output [ W-1:0]     outport_data1_o
    ,output [ W-1:0]     outport_data2_o
    ,output [ W-1:0]     outport_data3_o
    ,output [  2:0]      outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

reg [W-1:0] ram0 [0:31];
reg [W-1:0] ram1 [0:31];
reg [W-1:0] ram2 [0:31];
reg [W-1:0] ram3 [0:31];

reg [W-1:0] ram_read0_q;
reg [W-1:0] ram_read1_q;
reg [W-1:0] ram_read2_q;
reg [W-1:0] ram_read3_q;

always @ (posedge clk)
begin
    if (wr0_w)
        ram0[wr_ptr_w] <= inport_data_i;
end

always @ (posedge clk)
begin
    if (wr1_w)
        ram1[wr_ptr_w] <= inport_data_i;
end

always @ (posedge clk)
begin
    if (wr2_w)
        ram2[wr_ptr_w] <= inport_data_i;
end

always @ (posedge clk)
begin
    if (wr3_w)
        ram3[wr_ptr_w] <= inport_data_i;
end

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

always @ (posedge clk )
if (rst)
begin
    block_ready_q     <= 2'b0;
    block_wr_q        <= 1'b0;
    block_rd_q        <= 1'b0;
end
else if (img_start_i)
begin
    block_ready_q     <= 2'b0;
    block_wr_q        <= 1'b0;
    block_rd_q        <= 1'b0;
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

//-----------------------------------------------------------------
// FSM
//-----------------------------------------------------------------
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
    begin
        if (block_ready_q[block_rd_q] && outport_ready_i)
            next_state_r = STATE_SETUP;
    end
    STATE_SETUP:
    begin
        next_state_r = STATE_ACTIVE;
    end
    STATE_ACTIVE:
    begin
        if (outport_valid_o && rd_idx_q == 6'd63)
            next_state_r = STATE_IDLE;
    end
    default: ;
    endcase

    if (img_start_i)
        next_state_r = STATE_IDLE;
end

always @ (posedge clk )
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @ (posedge clk )
if (rst)
    rd_idx_q <= 6'b0;
else if (img_start_i)
    rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE)
    rd_idx_q <= rd_idx_q + 6'd1;

always @ (posedge clk )
if (rst)
    rd_addr_q <= 4'b0;
else if (state_q == STATE_IDLE)
    rd_addr_q <= 4'b0;
else if (state_q == STATE_SETUP)
    rd_addr_q <= 4'd8;
else if (state_q == STATE_ACTIVE)
begin
    case (rd_idx_q[2:0])
    3'd0: rd_addr_q <= rd_addr_q - 4'd8;
    3'd1: rd_addr_q <= rd_addr_q + 4'd8;
    3'd2: ;
    3'd3: rd_addr_q <= rd_addr_q - 4'd8;
    3'd4: rd_addr_q <= rd_addr_q + 4'd8;
    3'd5: rd_addr_q <= rd_addr_q - 4'd8;
    3'd6: rd_addr_q <= rd_addr_q + 4'd1;
    3'd7: rd_addr_q <= rd_addr_q + 4'd8;
    endcase
end

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

//-----------------------------------------------------------------
// Read path: register a one-hot decode of the NEXT read address and
// use a low fan-out AND-OR read instead of a 32:1 mux driven by the
// address registers.  The next-address logic below mirrors the
// rd_addr_q / block_rd_q register updates exactly, so the one-cycle
// read latency (and read-first behaviour) is unchanged.
//-----------------------------------------------------------------
reg [3:0] rd_addr_nxt;

always @ *
begin
    rd_addr_nxt = rd_addr_q;
    if (rst || img_start_i)
        rd_addr_nxt = 4'd0;
    else if (state_q == STATE_IDLE)
        rd_addr_nxt = 4'd0;
    else if (state_q == STATE_SETUP)
        rd_addr_nxt = 4'd8;
    else if (state_q == STATE_ACTIVE)
    begin
        case (rd_idx_q[2:0])
        3'd0: rd_addr_nxt = rd_addr_q - 4'd8;
        3'd1: rd_addr_nxt = rd_addr_q + 4'd8;
        3'd2: rd_addr_nxt = rd_addr_q;
        3'd3: rd_addr_nxt = rd_addr_q - 4'd8;
        3'd4: rd_addr_nxt = rd_addr_q + 4'd8;
        3'd5: rd_addr_nxt = rd_addr_q - 4'd8;
        3'd6: rd_addr_nxt = rd_addr_q + 4'd1;
        3'd7: rd_addr_nxt = rd_addr_q + 4'd8;
        endcase
    end
end

wire rd_blk_nxt = (rst || img_start_i) ? 1'b0 :
                  (outport_valid_o && rd_idx_q == 6'd63) ? ~block_rd_q :
                  block_rd_q;

reg [31:0] rd_oh_q;

always @ (posedge clk)
    rd_oh_q <= 32'd1 << {rd_blk_nxt, rd_addr_nxt};

// AND-OR read over the registered one-hot select
reg [W-1:0] rd0_r;
reg [W-1:0] rd1_r;
reg [W-1:0] rd2_r;
reg [W-1:0] rd3_r;

integer bi;
always @ *
begin
    rd0_r = {W{1'b0}};
    rd1_r = {W{1'b0}};
    rd2_r = {W{1'b0}};
    rd3_r = {W{1'b0}};
    for (bi = 0; bi < 32; bi = bi + 1)
    begin
        rd0_r = rd0_r | ({W{rd_oh_q[bi]}} & ram0[bi]);
        rd1_r = rd1_r | ({W{rd_oh_q[bi]}} & ram1[bi]);
        rd2_r = rd2_r | ({W{rd_oh_q[bi]}} & ram2[bi]);
        rd3_r = rd3_r | ({W{rd_oh_q[bi]}} & ram3[bi]);
    end
end

always @ (posedge clk)
begin
    ram_read0_q <= rd0_r;
    ram_read1_q <= rd1_r;
    ram_read2_q <= rd2_r;
    ram_read3_q <= rd3_r;
end

assign outport_data0_o = ram_read0_q;
assign outport_data1_o = ram_read1_q;
assign outport_data2_o = ram_read2_q;
assign outport_data3_o = ram_read3_q;

endmodule

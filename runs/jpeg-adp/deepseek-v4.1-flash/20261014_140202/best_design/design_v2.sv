//-----------------------------------------------------------------
// jpeg_idct_2d — bit-exact optimized equivalent of the core_jpeg
// composition (jpeg_idct_x -> transpose -> jpeg_idct_y).
// Same schedule/arithmetic; multiplier operand widths narrowed to the
// true data ranges (product low-32 wrap preserved), transpose datapath
// 21 bits.
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

wire          x_valid;
wire [20:0]   x_data;
wire [ 5:0]   x_idx;

wire          t_valid;
wire [20:0]   t_data0, t_data1, t_data2, t_data3;
wire [ 2:0]   t_idx;

jpeg_idct_x u_idct_x
(
     .clk(clk), .rst(rst)
    ,.img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(inport_data0_i), .inport_data1_i(inport_data1_i)
    ,.inport_data2_i(inport_data2_i), .inport_data3_i(inport_data3_i)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(x_valid), .outport_data_o(x_data), .outport_idx_o(x_idx)
);

jpeg_idct_transpose u_transpose
(
     .clk(clk), .rst(rst)
    ,.img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(x_valid), .inport_data_i(x_data), .inport_idx_i(x_idx)
    ,.outport_valid_o(t_valid)
    ,.outport_data0_o(t_data0), .outport_data1_o(t_data1)
    ,.outport_data2_o(t_data2), .outport_data3_o(t_data3)
    ,.outport_idx_o(t_idx)
);

jpeg_idct_y u_idct_y
(
     .clk(clk), .rst(rst)
    ,.img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(t_valid)
    ,.inport_data0_i(t_data0), .inport_data1_i(t_data1)
    ,.inport_data2_i(t_data2), .inport_data3_i(t_data3)
    ,.inport_idx_i(t_idx)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//-----------------------------------------------------------------
// Row pass (OUT_SHIFT = 11)
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
    ,output [ 20:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam signed [12:0] C1_16 = 4017;
localparam signed [12:0] C2_16 = 3784;
localparam signed [12:0] C3_16 = 3406;
localparam signed [12:0] C4_16 = 2896;
localparam signed [12:0] C5_16 = 2276;
localparam signed [12:0] C6_16 = 1567;
localparam signed [12:0] C7_16 = 799;

wire signed [16:0] block_in_0_1 = $signed(inport_data0_i);
wire signed [16:0] block_in_2_3 = $signed(inport_data1_i);
wire signed [16:0] block_in_4_5 = $signed(inport_data2_i);
wire signed [16:0] block_in_6_7 = $signed(inport_data3_i);

reg signed [16:0] i0;
reg signed [16:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
reg signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 17'sd0;
    mul0_a <= 17'sd0; mul1_a <= 17'sd0; mul2_a <= 17'sd0; mul3_a <= 17'sd0; mul4_a <= 17'sd0;
    mul0_b <= 13'sd0; mul1_b <= 13'sd0; mul2_b <= 13'sd0; mul3_b <= 13'sd0; mul4_b <= 13'sd0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul0_b <= C2_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C1_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C7_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C5_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C3_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C2_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

wire signed [29:0] p0 = mul0_a * mul0_b;
wire signed [29:0] p1 = mul1_a * mul1_b;
wire signed [29:0] p2 = mul2_a * mul2_b;
wire signed [29:0] p3 = mul3_a * mul3_b;
wire signed [29:0] p4 = mul4_a * mul4_b;

reg signed [29:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 30'sd0; mul1_q <= 30'sd0; mul2_q <= 30'sd0; mul3_q <= 30'sd0; mul4_q <= 30'sd0;
end
else
begin
    mul0_q <= p0;
    mul1_q <= p1;
    mul2_q <= p2;
    mul3_q <= p3;
    mul4_q <= p4;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'sd0; mul1 <= 32'sd0; mul2 <= 32'sd0; mul3 <= 32'sd0; mul4 <= 32'sd0;
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

always @ (posedge clk)
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

always @ (posedge clk)
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

always @ (posedge clk)
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

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'sd0; o_s6 <= 32'sd0; o_s7 <= 32'sd0;
    o_t0 <= 32'sd0; o_t1 <= 32'sd0; o_t2 <= 32'sd0; o_t3 <= 32'sd0;
    o_t4 <= 32'sd0; o_t5 <= 32'sd0; o_t6 <= 32'sd0; o_t7 <= 32'sd0;
    o_t6_5 <= 32'sd0; o_t5_6 <= 32'sd0;
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
        o_s5 <= (o_t6_5 * 181) / 256; // 1/sqrt(2)
        o_s6 <= (o_t5_6 * 181) / 256; // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
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

reg signed [20:0] block_out[0:7];
reg signed [20:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 21'sd0; block_out[1] <= 21'sd0;
    block_out[2] <= 21'sd0; block_out[3] <= 21'sd0;
    block_out[4] <= 21'sd0; block_out[5] <= 21'sd0;
    block_out[6] <= 21'sd0; block_out[7] <= 21'sd0;
    block_out_tmp <= 21'sd0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 11);
        block_out_tmp <= ((o_t0 - o_t7) >>> 11);
        block_out[3] <= ((o_t3 + o_t4) >>> 11);
        block_out[4] <= ((o_t3 - o_t4) >>> 11);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> 11);
        block_out[5] <= ((o_t2 - o_s5) >>> 11);
        block_out[1] <= ((o_t1 + o_s6) >>> 11);
        block_out[6] <= ((o_t1 - o_s6) >>> 11);
    end
end

reg [7:0] valid_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = {{11{block_out[ptr_q[2:0]][20]}}, block_out[ptr_q[2:0]]};
assign outport_idx_o   = ptr_q;

endmodule

//-----------------------------------------------------------------
// Corner turn: double buffered transpose
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 20:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 20:0]  outport_data0_o
    ,output [ 20:0]  outport_data1_o
    ,output [ 20:0]  outport_data2_o
    ,output [ 20:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire [20:0] outport_data0_w;
wire [20:0] outport_data1_w;
wire [20:0] outport_data2_w;
wire [20:0] outport_data3_w;

wire inport_accept_o;

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram u_ram0
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr0_w)
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data0_w)
);
jpeg_idct_transpose_ram u_ram1
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr1_w)
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data1_w)
);
jpeg_idct_transpose_ram u_ram2
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr2_w)
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data2_w)
);
jpeg_idct_transpose_ram u_ram3
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr3_w)
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data3_w)
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

    if (outport_valid_o && rd_idx_q == 6'd63)
    begin
        block_ready_q[block_rd_q] <= 1'b0;
        block_rd_q <= ~block_rd_q;
    end
end

assign inport_accept_o = ~block_ready_q[block_wr_q];

localparam STATE_IDLE   = 2'd0;
localparam STATE_SETUP  = 2'd1;
localparam STATE_ACTIVE = 2'd2;

reg [1:0] state_q;
reg [1:0] next_state_r;

always @ *
begin
    next_state_r = state_q;

    case (state_q)
    STATE_IDLE:
    begin
        if (block_ready_q[block_rd_q])
            next_state_r = STATE_SETUP;
    end
    STATE_SETUP:
    begin
        next_state_r = STATE_ACTIVE;
    end
    STATE_ACTIVE:
    begin
        if (rd_idx_q == 6'd63)
            next_state_r = STATE_IDLE;
    end
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

always @ (posedge clk)
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
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

module jpeg_idct_transpose_ram
(
     input           clk
    ,input  [  4:0]  addr0_i
    ,input  [ 20:0]  data0_i
    ,input           wr0_i
    ,input  [  4:0]  addr1_i
    ,output [ 20:0]  data1_o
);

reg [20:0] ram [31:0];
reg [20:0] ram_read1_q;

always @ (posedge clk)
begin
    if (wr0_i)
        ram[addr0_i] <= data0_i;

    ram_read1_q <= ram[addr1_i];
end

assign data1_o = ram_read1_q;

endmodule

//-----------------------------------------------------------------
// Column pass (OUT_SHIFT = 15)
//-----------------------------------------------------------------
module jpeg_idct_y
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 20:0]  inport_data0_i
    ,input  [ 20:0]  inport_data1_i
    ,input  [ 20:0]  inport_data2_i
    ,input  [ 20:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam signed [12:0] C1_16 = 4017;
localparam signed [12:0] C2_16 = 3784;
localparam signed [12:0] C3_16 = 3406;
localparam signed [12:0] C4_16 = 2896;
localparam signed [12:0] C5_16 = 2276;
localparam signed [12:0] C6_16 = 1567;
localparam signed [12:0] C7_16 = 799;

wire signed [20:0] block_in_0_1 = inport_data0_i;
wire signed [20:0] block_in_2_3 = inport_data1_i;
wire signed [20:0] block_in_4_5 = inport_data2_i;
wire signed [20:0] block_in_6_7 = inport_data3_i;

reg signed [21:0] i0;
reg signed [21:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
reg signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 22'sd0;
    mul0_a <= 22'sd0; mul1_a <= 22'sd0; mul2_a <= 22'sd0; mul3_a <= 22'sd0; mul4_a <= 22'sd0;
    mul0_b <= 13'sd0; mul1_b <= 13'sd0; mul2_b <= 13'sd0; mul3_b <= 13'sd0; mul4_b <= 13'sd0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul0_b <= C2_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C1_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C7_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C5_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C3_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C2_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

// full 35-bit products; low 32 bits reproduce the original 32-bit wrap
wire signed [34:0] p0 = mul0_a * mul0_b;
wire signed [34:0] p1 = mul1_a * mul1_b;
wire signed [34:0] p2 = mul2_a * mul2_b;
wire signed [34:0] p3 = mul3_a * mul3_b;
wire signed [34:0] p4 = mul4_a * mul4_b;

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'sd0; mul1_q <= 32'sd0; mul2_q <= 32'sd0; mul3_q <= 32'sd0; mul4_q <= 32'sd0;
end
else
begin
    mul0_q <= p0[31:0];
    mul1_q <= p1[31:0];
    mul2_q <= p2[31:0];
    mul3_q <= p3[31:0];
    mul4_q <= p4[31:0];
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'sd0; mul1 <= 32'sd0; mul2 <= 32'sd0; mul3 <= 32'sd0; mul4 <= 32'sd0;
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

always @ (posedge clk)
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

always @ (posedge clk)
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

always @ (posedge clk)
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

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;
wire signed [31:0] sc5 = o_t6_5 * 181;
wire signed [31:0] sc6 = o_t5_6 * 181;

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'sd0; o_s6 <= 32'sd0; o_s7 <= 32'sd0;
    o_t0 <= 32'sd0; o_t1 <= 32'sd0; o_t2 <= 32'sd0; o_t3 <= 32'sd0;
    o_t4 <= 32'sd0; o_t5 <= 32'sd0; o_t6 <= 32'sd0; o_t7 <= 32'sd0;
    o_t6_5 <= 32'sd0; o_t5_6 <= 32'sd0;
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
        o_s5 <= sc5 / 256; // 1/sqrt(2)
        o_s6 <= sc6 / 256; // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
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

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 32'sd0; block_out[1] <= 32'sd0;
    block_out[2] <= 32'sd0; block_out[3] <= 32'sd0;
    block_out[4] <= 32'sd0; block_out[5] <= 32'sd0;
    block_out[6] <= 32'sd0; block_out[7] <= 32'sd0;
    block_out_tmp <= 32'sd0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 15);
        block_out_tmp <= ((o_t0 - o_t7) >>> 15);
        block_out[3] <= ((o_t3 + o_t4) >>> 15);
        block_out[4] <= ((o_t3 - o_t4) >>> 15);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> 15);
        block_out[5] <= ((o_t2 - o_s5) >>> 15);
        block_out[1] <= ((o_t1 + o_s6) >>> 15);
        block_out[6] <= ((o_t1 - o_s6) >>> 15);
    end
end

reg [7:0] valid_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

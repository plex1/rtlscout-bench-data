//-----------------------------------------------------------------
// jpeg_idct_2d : x-pass -> transpose -> y-pass
// v4: transpose read uses a registered one-hot select (kills the
//     huge fanout from the read-address register into the 32:1 mux).
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

jpeg_idct_x u_idct_x
(
     .clk(clk), .rst(rst),
     .img_start_i(img_start_i), .img_end_i(img_end_i),
     .inport_valid_i(inport_valid_i),
     .inport_data0_i(inport_data0_i), .inport_data1_i(inport_data1_i),
     .inport_data2_i(inport_data2_i), .inport_data3_i(inport_data3_i),
     .inport_idx_i(inport_idx_i),
     .outport_valid_o(idct_x_valid_w),
     .outport_data_o(idct_x_data_w),
     .outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [ 19:0]  transpose_data0_w;
wire [ 19:0]  transpose_data1_w;
wire [ 19:0]  transpose_data2_w;
wire [ 19:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose u_transpose
(
     .clk(clk), .rst(rst),
     .img_start_i(img_start_i),
     .inport_valid_i(idct_x_valid_w),
     .inport_data_i(idct_x_data_w),
     .inport_idx_i(idct_x_idx_w),
     .inport_accept_o(),
     .outport_valid_o(transpose_valid_w),
     .outport_data0_o(transpose_data0_w),
     .outport_data1_o(transpose_data1_w),
     .outport_data2_o(transpose_data2_w),
     .outport_data3_o(transpose_data3_w),
     .outport_idx_o(transpose_idx_w)
);

jpeg_idct_y u_idct_y
(
     .clk(clk), .rst(rst),
     .img_start_i(img_start_i), .img_end_i(img_end_i),
     .inport_valid_i(transpose_valid_w),
     .inport_data0_i(transpose_data0_w),
     .inport_data1_i(transpose_data1_w),
     .inport_data2_i(transpose_data2_w),
     .inport_data3_i(transpose_data3_w),
     .inport_idx_i(transpose_idx_w),
     .outport_valid_o(outport_valid_o),
     .outport_data_o(outport_data_o),
     .outport_idx_o(outport_idx_o)
);

endmodule

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

localparam signed [12:0] C1_16 = 4017;
localparam signed [12:0] C2_16 = 3784;
localparam signed [12:0] C3_16 = 3406;
localparam signed [12:0] C4_16 = 2896;
localparam signed [12:0] C5_16 = 2276;
localparam signed [12:0] C6_16 = 1567;
localparam signed [12:0] C7_16 = 799;

wire signed [15:0] b0 = inport_data0_i;
wire signed [15:0] b2 = inport_data1_i;
wire signed [15:0] b4 = inport_data2_i;
wire signed [15:0] b6 = inport_data3_i;

reg  signed [17:0] i0;
reg  signed [17:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
reg  signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 18'b0;
    mul0_a <= 18'b0; mul0_b <= 13'b0;
    mul1_a <= 18'b0; mul1_b <= 13'b0;
    mul2_a <= 18'b0; mul2_b <= 13'b0;
    mul3_a <= 18'b0; mul3_b <= 13'b0;
    mul4_a <= 18'b0; mul4_b <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= b0 + b4;
        mul0_a <= b2;  mul0_b <= C2_16;
        mul1_a <= b6;  mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= b0;  mul0_b <= C1_16;
        mul1_a <= b6;  mul1_b <= C7_16;
        mul2_a <= b4;  mul2_b <= C5_16;
        mul3_a <= b2;  mul3_b <= C3_16;
        mul4_a <= i0;  mul4_b <= C4_16;
    end
    3'd2: i0 <= b0 - b4;
    3'd3, 3'd4:
    begin
        mul0_a <= b0;  mul0_b <= C7_16;
        mul1_a <= b6;  mul1_b <= C1_16;
        mul2_a <= b4;  mul2_b <= C3_16;
        mul3_a <= b2;  mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= b2;  mul0_b <= C6_16;
        mul1_a <= b6;  mul1_b <= C2_16;
        mul4_a <= i0;  mul4_b <= C4_16;
    end
    default: ;
    endcase
end

reg signed [30:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 31'b0; mul1_q <= 31'b0; mul2_q <= 31'b0;
    mul3_q <= 31'b0; mul4_q <= 31'b0;
end
else
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [30:0] mul0, mul1, mul2, mul3, mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 31'b0; mul1 <= 31'b0; mul2 <= 31'b0; mul3 <= 31'b0; mul4 <= 31'b0;
end
else
begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q; mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;    out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;  out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q;  out_stg2_idx_q <= out_stg1_idx_q;
    out_stg3_valid_q <= out_stg2_valid_q;  out_stg3_idx_q <= out_stg2_idx_q;
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

wire signed [31:0] px_s5 = o_t6_5 * 181;
wire signed [31:0] px_s6 = o_t5_6 * 181;
wire signed [31:0] o_s5_n = (px_s5 + (px_s5[31] ? 32'sd255 : 32'sd0)) >>> 8;
wire signed [31:0] o_s6_n = (px_s6 + (px_s6[31] ? 32'sd255 : 32'sd0)) >>> 8;

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'b0; o_s6 <= 32'b0; o_s7 <= 32'b0;
    o_t0 <= 32'b0; o_t1 <= 32'b0; o_t2 <= 32'b0; o_t3 <= 32'b0;
    o_t4 <= 32'b0; o_t5 <= 32'b0; o_t6 <= 32'b0; o_t7 <= 32'b0;
    o_t6_5 <= 32'b0; o_t5_6 <= 32'b0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0: o_t3 <= mul0 + mul1;
    3'd1:
    begin
        o_s7 <= mul0 + mul1;
        o_s6 <= mul2 + mul3;
        o_t0 <= mul4;
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3;
        o_t3 <= o_t0 - o_t3;
        o_t7 <= o_s6 + o_s7;
    end
    3'd3: o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    3'd4:
    begin
        o_t0 <= mul0 - mul1;
        o_s5 <= mul2 - mul3;
    end
    3'd5:
    begin
        o_t3 <= mul0 - mul1;
        o_t4 <= mul4;
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
        o_s5 <= o_s5_n;
        o_s6 <= o_s6_n;
    end
    endcase
end

reg signed [19:0] block_out[0:7];
reg signed [19:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 20'b0; block_out[1] <= 20'b0; block_out[2] <= 20'b0;
    block_out[3] <= 20'b0; block_out[4] <= 20'b0; block_out[5] <= 20'b0;
    block_out[6] <= 20'b0; block_out[7] <= 20'b0;
    block_out_tmp <= 20'b0;
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
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = ptr_q;

endmodule

//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
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

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};
wire        wr_en_w  = inport_valid_i && inport_accept_o;

// Registered one-hot read select, computed from the NEXT address so that it
// is aligned with the address actually applied during the current cycle.
// This removes the huge fanout of rd_addr_q into the 32:1 read muxes.
reg  [3:0]  nxt_addr;
reg  [31:0] rd_sel_q;

always @ *
begin
    if (state_q == STATE_ACTIVE)
    begin
        case (rd_idx_q[2:0])
        3'd0:    nxt_addr = rd_addr_q - 4'd8;
        3'd1:    nxt_addr = rd_addr_q + 4'd8;
        3'd2:    nxt_addr = rd_addr_q;
        3'd3:    nxt_addr = rd_addr_q - 4'd8;
        3'd4:    nxt_addr = rd_addr_q + 4'd8;
        3'd5:    nxt_addr = rd_addr_q - 4'd8;
        3'd6:    nxt_addr = rd_addr_q + 4'd1;
        default: nxt_addr = rd_addr_q + 4'd8;
        endcase
    end
    else if (state_q == STATE_SETUP)
        nxt_addr = 4'd8;
    else
        nxt_addr = 4'd0;
end

always @ (posedge clk)
if (rst)
    rd_sel_q <= 32'd1;
else
    rd_sel_q <= 32'd1 << {block_rd_q, nxt_addr};

wire [19:0] outport_data0_w;
wire [19:0] outport_data1_w;
wire [19:0] outport_data2_w;
wire [19:0] outport_data3_w;

jpeg_idct_transpose_ram u_ram0
(
     .clk_i(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i),
     .wr0_i(wr_en_w && (inport_idx_i[5:4] == 2'd0)),
     .sel1_i(rd_sel_q), .data1_o(outport_data0_w)
);
jpeg_idct_transpose_ram u_ram1
(
     .clk_i(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i),
     .wr0_i(wr_en_w && (inport_idx_i[5:4] == 2'd1)),
     .sel1_i(rd_sel_q), .data1_o(outport_data1_w)
);
jpeg_idct_transpose_ram u_ram2
(
     .clk_i(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i),
     .wr0_i(wr_en_w && (inport_idx_i[5:4] == 2'd2)),
     .sel1_i(rd_sel_q), .data1_o(outport_data2_w)
);
jpeg_idct_transpose_ram u_ram3
(
     .clk_i(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i),
     .wr0_i(wr_en_w && (inport_idx_i[5:4] == 2'd3)),
     .sel1_i(rd_sel_q), .data1_o(outport_data3_w)
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
    STATE_IDLE:   if (block_ready_q[block_rd_q]) next_state_r = STATE_SETUP;
    STATE_SETUP:  next_state_r = STATE_ACTIVE;
    STATE_ACTIVE: if (outport_valid_o && rd_idx_q == 6'd63) next_state_r = STATE_IDLE;
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

// rd_addr_q holds the address used in the CURRENT read cycle; nxt_addr computes
// the following one combinationally so rd_sel_q can be registered.
always @ (posedge clk)
if (rst)
    rd_addr_q <= 4'b0;
else
    rd_addr_q <= nxt_addr;

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

module jpeg_idct_transpose_ram
(
     input           clk_i
    ,input  [  4:0]  addr0_i
    ,input  [ 19:0]  data0_i
    ,input           wr0_i
    ,input  [ 31:0]  sel1_i
    ,output [ 19:0]  data1_o
);

reg  [19:0] ram [0:31];
reg  [19:0] rd_q;
reg  [19:0] mux_r;

integer i;

always @ (posedge clk_i)
begin
    if (wr0_i)
        ram[addr0_i] <= data0_i;
end

always @ *
begin
    mux_r = 20'd0;
    for (i = 0; i < 32; i = i + 1)
        mux_r = mux_r | ({20{sel1_i[i]}} & ram[i]);
end

always @ (posedge clk_i)
    rd_q <= mux_r;

assign data1_o = rd_q;

endmodule

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

localparam signed [12:0] C1_16 = 4017;
localparam signed [12:0] C2_16 = 3784;
localparam signed [12:0] C3_16 = 3406;
localparam signed [12:0] C4_16 = 2896;
localparam signed [12:0] C5_16 = 2276;
localparam signed [12:0] C6_16 = 1567;
localparam signed [12:0] C7_16 = 799;

wire signed [19:0] b0 = inport_data0_i;
wire signed [19:0] b2 = inport_data1_i;
wire signed [19:0] b4 = inport_data2_i;
wire signed [19:0] b6 = inport_data3_i;

reg  signed [20:0] i0;
reg  signed [20:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
reg  signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 21'b0;
    mul0_a <= 21'b0; mul0_b <= 13'b0;
    mul1_a <= 21'b0; mul1_b <= 13'b0;
    mul2_a <= 21'b0; mul2_b <= 13'b0;
    mul3_a <= 21'b0; mul3_b <= 13'b0;
    mul4_a <= 21'b0; mul4_b <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= b0 + b4;
        mul0_a <= b2;  mul0_b <= C2_16;
        mul1_a <= b6;  mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= b0;  mul0_b <= C1_16;
        mul1_a <= b6;  mul1_b <= C7_16;
        mul2_a <= b4;  mul2_b <= C5_16;
        mul3_a <= b2;  mul3_b <= C3_16;
        mul4_a <= i0;  mul4_b <= C4_16;
    end
    3'd2: i0 <= b0 - b4;
    3'd3, 3'd4:
    begin
        mul0_a <= b0;  mul0_b <= C7_16;
        mul1_a <= b6;  mul1_b <= C1_16;
        mul2_a <= b4;  mul2_b <= C3_16;
        mul3_a <= b2;  mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= b2;  mul0_b <= C6_16;
        mul1_a <= b6;  mul1_b <= C2_16;
        mul4_a <= i0;  mul4_b <= C4_16;
    end
    default: ;
    endcase
end

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0; mul1_q <= 32'b0; mul2_q <= 32'b0;
    mul3_q <= 32'b0; mul4_q <= 32'b0;
end
else
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0; mul1 <= 32'b0; mul2 <= 32'b0; mul3 <= 32'b0; mul4 <= 32'b0;
end
else
begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q; mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;    out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;  out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q;  out_stg2_idx_q <= out_stg1_idx_q;
    out_stg3_valid_q <= out_stg2_valid_q;  out_stg3_idx_q <= out_stg2_idx_q;
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

wire signed [31:0] py_s5 = o_t6_5 * 181;
wire signed [31:0] py_s6 = o_t5_6 * 181;
wire signed [31:0] o_s5_n = (py_s5 + (py_s5[31] ? 32'sd255 : 32'sd0)) >>> 8;
wire signed [31:0] o_s6_n = (py_s6 + (py_s6[31] ? 32'sd255 : 32'sd0)) >>> 8;

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'b0; o_s6 <= 32'b0; o_s7 <= 32'b0;
    o_t0 <= 32'b0; o_t1 <= 32'b0; o_t2 <= 32'b0; o_t3 <= 32'b0;
    o_t4 <= 32'b0; o_t5 <= 32'b0; o_t6 <= 32'b0; o_t7 <= 32'b0;
    o_t6_5 <= 32'b0; o_t5_6 <= 32'b0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0: o_t3 <= mul0 + mul1;
    3'd1:
    begin
        o_s7 <= mul0 + mul1;
        o_s6 <= mul2 + mul3;
        o_t0 <= mul4;
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3;
        o_t3 <= o_t0 - o_t3;
        o_t7 <= o_s6 + o_s7;
    end
    3'd3: o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    3'd4:
    begin
        o_t0 <= mul0 - mul1;
        o_s5 <= mul2 - mul3;
    end
    3'd5:
    begin
        o_t3 <= mul0 - mul1;
        o_t4 <= mul4;
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
        o_s5 <= o_s5_n;
        o_s6 <= o_s6_n;
    end
    endcase
end

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 32'b0; block_out[1] <= 32'b0; block_out[2] <= 32'b0;
    block_out[3] <= 32'b0; block_out[4] <= 32'b0; block_out[5] <= 32'b0;
    block_out[6] <= 32'b0; block_out[7] <= 32'b0;
    block_out_tmp <= 32'b0;
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

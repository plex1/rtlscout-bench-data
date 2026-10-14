//-----------------------------------------------------------------
// jpeg_idct_2d - full 8x8 2-D IDCT of the baseline JPEG decoder
// (ultraembedded/core_jpeg). Optimized variant:
//   - narrowed datapath widths (the x-pass output provably fits in 22
//     bits, the cosine constants in 13 bits)
//   - the y-pass idx permutation replaced by an equivalent bit swap
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

localparam XW = 22;

wire          idct_x_valid_w;
wire [XW-1:0] idct_x_data_w;
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
wire [XW-1:0] transpose_data0_w;
wire [XW-1:0] transpose_data1_w;
wire [XW-1:0] transpose_data2_w;
wire [XW-1:0] transpose_data3_w;
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
// Row pass (1-D IDCT, OUT_SHIFT=11)
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
    ,output [ 21:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam OUT_SHIFT = 11;

localparam [12:0] C1_16 = 4017; // cos( pi/16) x4096
localparam [12:0] C2_16 = 3784; // cos(2pi/16) x4096
localparam [12:0] C3_16 = 3406; // cos(3pi/16) x4096
localparam [12:0] C4_16 = 2896; // cos(4pi/16) x4096
localparam [12:0] C5_16 = 2276; // cos(5pi/16) x4096
localparam [12:0] C6_16 = 1567; // cos(6pi/16) x4096
localparam [12:0] C7_16 = 799;  // cos(7pi/16) x4096

wire signed [16:0] block_in_0_1 = {{1{inport_data0_i[15]}}, inport_data0_i};
wire signed [16:0] block_in_2_3 = {{1{inport_data1_i[15]}}, inport_data1_i};
wire signed [16:0] block_in_4_5 = {{1{inport_data2_i[15]}}, inport_data2_i};
wire signed [16:0] block_in_6_7 = {{1{inport_data3_i[15]}}, inport_data3_i};

reg signed [17:0] i0;
reg signed [16:0] mul0_a;
reg signed [12:0] mul0_b;
reg signed [16:0] mul1_a;
reg signed [12:0] mul1_b;
reg signed [16:0] mul2_a;
reg signed [12:0] mul2_b;
reg signed [16:0] mul3_a;
reg signed [12:0] mul3_b;
reg signed [16:0] mul4_a;
reg signed [12:0] mul4_b;

always @ (posedge clk )
if (rst)
begin
    i0     <= 18'b0;
    mul0_a <= 17'b0;
    mul0_b <= 13'b0;
    mul1_a <= 17'b0;
    mul1_b <= 13'b0;
    mul2_a <= 17'b0;
    mul2_b <= 13'b0;
    mul3_a <= 17'b0;
    mul3_b <= 13'b0;
    mul4_a <= 17'b0;
    mul4_b <= 13'b0;
end
else
begin
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
        mul4_a <= i0[16:0];
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
        mul4_a <= i0[16:0];
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
end

reg signed [30:0] mul0_q;
reg signed [30:0] mul1_q;
reg signed [30:0] mul2_q;
reg signed [30:0] mul3_q;
reg signed [30:0] mul4_q;

always @ (posedge clk )
if (rst)
begin
    mul0_q <= 31'b0;
    mul1_q <= 31'b0;
    mul2_q <= 31'b0;
    mul3_q <= 31'b0;
    mul4_q <= 31'b0;
end
else
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [30:0] mul0;
reg signed [30:0] mul1;
reg signed [30:0] mul2;
reg signed [30:0] mul3;
reg signed [30:0] mul4;

always @ (posedge clk )
if (rst)
begin
    mul0 <= 31'b0;
    mul1 <= 31'b0;
    mul2 <= 31'b0;
    mul3 <= 31'b0;
    mul4 <= 31'b0;
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

reg signed [30:0] o_s5;
reg signed [30:0] o_s6;
reg signed [30:0] o_s7;
reg signed [30:0] o_t0;
reg signed [30:0] o_t1;
reg signed [30:0] o_t2;
reg signed [30:0] o_t3;
reg signed [30:0] o_t4;
reg signed [30:0] o_t5;
reg signed [30:0] o_t6;
reg signed [30:0] o_t7;
reg signed [30:0] o_t6_5;
reg signed [30:0] o_t5_6;

always @ (posedge clk )
if (rst)
begin
    o_s5   <= 31'b0;
    o_s6   <= 31'b0;
    o_s7   <= 31'b0;
    o_t0   <= 31'b0;
    o_t1   <= 31'b0;
    o_t2   <= 31'b0;
    o_t3   <= 31'b0;
    o_t4   <= 31'b0;
    o_t5   <= 31'b0;
    o_t6   <= 31'b0;
    o_t7   <= 31'b0;
    o_t6_5 <= 31'b0;
    o_t5_6 <= 31'b0;
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
        o_t4 <= mul4; // s1
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

reg signed [21:0] block_out[0:7];
reg signed [21:0] block_out_tmp;

always @ (posedge clk )
if (rst)
begin
    block_out[0] <= 22'b0;
    block_out[1] <= 22'b0;
    block_out[2] <= 22'b0;
    block_out[3] <= 22'b0;
    block_out[4] <= 22'b0;
    block_out[5] <= 22'b0;
    block_out[6] <= 22'b0;
    block_out[7] <= 22'b0;
    block_out_tmp <= 22'b0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= (o_t0 + o_t7) >>> OUT_SHIFT;
        block_out_tmp <= (o_t0 - o_t7) >>> OUT_SHIFT; // block_out[7]
        block_out[3] <= (o_t3 + o_t4) >>> OUT_SHIFT;
        block_out[4] <= (o_t3 - o_t4) >>> OUT_SHIFT;
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= (o_t2 + o_s5) >>> OUT_SHIFT;
        block_out[5] <= (o_t2 - o_s5) >>> OUT_SHIFT;
        block_out[1] <= (o_t1 + o_s6) >>> OUT_SHIFT;
        block_out[6] <= (o_t1 - o_s6) >>> OUT_SHIFT;
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
assign outport_idx_o   = ptr_q;

endmodule

//-----------------------------------------------------------------
// Corner turn (double buffered)
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 21:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          outport_valid_o
    ,output [ 21:0]  outport_data0_o
    ,output [ 21:0]  outport_data1_o
    ,output [ 21:0]  outport_data2_o
    ,output [ 21:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire [21:0] outport_data0_w;
wire [21:0] outport_data1_w;
wire [21:0] outport_data2_w;
wire [21:0] outport_data3_w;

wire wr0_w = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram u_ram0
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr0_w),
     .addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data0_w)
);

jpeg_idct_transpose_ram u_ram1
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr1_w),
     .addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data1_w)
);

jpeg_idct_transpose_ram u_ram2
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr2_w),
     .addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data2_w)
);

jpeg_idct_transpose_ram u_ram3
(
     .clk(clk), .addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr3_w),
     .addr1_i({block_rd_q, rd_addr_q}), .data1_o(outport_data3_w)
);

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
    if (inport_valid_i && inport_idx_i == 6'd63)
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

//-----------------------------------------------------------------
// FSM
//-----------------------------------------------------------------
localparam STATE_IDLE        = 2'd0;
localparam STATE_SETUP       = 2'd1;
localparam STATE_ACTIVE      = 2'd2;

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
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

//-----------------------------------------------------------------
// Transpose storage (single write port, single read port)
//-----------------------------------------------------------------
module jpeg_idct_transpose_ram
(
     input           clk
    ,input  [  4:0]  addr0_i
    ,input  [ 21:0]  data0_i
    ,input           wr0_i
    ,input  [  4:0]  addr1_i
    ,output [ 21:0]  data1_o
);

//-----------------------------------------------------------------
// Dual Port RAM
// Mode: Read First
//-----------------------------------------------------------------
/* verilator lint_off MULTIDRIVEN */
reg [21:0]   ram [31:0];
/* verilator lint_on MULTIDRIVEN */

reg [21:0] rdata0;
reg [21:0] rdata1;
reg [21:0] ram_read1_q;

// one-hot selects (expressed as decoders so the address flops only drive
// the decoder inputs, not hundreds of mux gates)
wire [31:0] sel_wr = 32'b1 << addr0_i;
wire [31:0] sel_rd = 32'b1 << addr1_i;

integer i;

always @ (posedge clk)
begin
    if (wr0_i)
    begin
        for (i = 0; i < 32; i = i + 1)
            if (sel_wr[i])
                ram[i] <= data0_i;
    end
end

always @ *
begin
    rdata0 = 22'b0;
    rdata1 = 22'b0;
    for (i = 0; i < 16; i = i + 1)
    begin
        rdata0 = rdata0 | ({22{sel_rd[i]}}       & ram[i]);
        rdata1 = rdata1 | ({22{sel_rd[i + 16]}} & ram[i + 16]);
    end
end

always @ (posedge clk)
begin
    ram_read1_q <= rdata0 | rdata1;
end

assign data1_o = ram_read1_q;

endmodule

//-----------------------------------------------------------------
// Column pass (1-D IDCT, OUT_SHIFT=15)
//-----------------------------------------------------------------
module jpeg_idct_y
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 21:0]  inport_data0_i
    ,input  [ 21:0]  inport_data1_i
    ,input  [ 21:0]  inport_data2_i
    ,input  [ 21:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam OUT_SHIFT = 15;

localparam [12:0] C1_16 = 4017; // cos( pi/16) x4096
localparam [12:0] C2_16 = 3784; // cos(2pi/16) x4096
localparam [12:0] C3_16 = 3406; // cos(3pi/16) x4096
localparam [12:0] C4_16 = 2896; // cos(4pi/16) x4096
localparam [12:0] C5_16 = 2276; // cos(5pi/16) x4096
localparam [12:0] C6_16 = 1567; // cos(6pi/16) x4096
localparam [12:0] C7_16 = 799;  // cos(7pi/16) x4096

wire signed [22:0] block_in_0_1 = {{1{inport_data0_i[21]}}, inport_data0_i};
wire signed [22:0] block_in_2_3 = {{1{inport_data1_i[21]}}, inport_data1_i};
wire signed [22:0] block_in_4_5 = {{1{inport_data2_i[21]}}, inport_data2_i};
wire signed [22:0] block_in_6_7 = {{1{inport_data3_i[21]}}, inport_data3_i};

reg signed [23:0] i0;
reg signed [22:0] mul0_a;
reg signed [12:0] mul0_b;
reg signed [22:0] mul1_a;
reg signed [12:0] mul1_b;
reg signed [22:0] mul2_a;
reg signed [12:0] mul2_b;
reg signed [22:0] mul3_a;
reg signed [12:0] mul3_b;
reg signed [22:0] mul4_a;
reg signed [12:0] mul4_b;

always @ (posedge clk )
if (rst)
begin
    i0     <= 24'b0;
    mul0_a <= 23'b0;
    mul0_b <= 13'b0;
    mul1_a <= 23'b0;
    mul1_b <= 13'b0;
    mul2_a <= 23'b0;
    mul2_b <= 13'b0;
    mul3_a <= 23'b0;
    mul3_b <= 13'b0;
    mul4_a <= 23'b0;
    mul4_b <= 13'b0;
end
else
begin
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
        mul4_a <= i0[22:0];
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
        mul4_a <= i0[22:0];
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
end

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
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
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
        o_t4 <= mul4; // s1
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

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk )
if (rst)
begin
    block_out[0] <= 32'b0;
    block_out[1] <= 32'b0;
    block_out[2] <= 32'b0;
    block_out[3] <= 32'b0;
    block_out[4] <= 32'b0;
    block_out[5] <= 32'b0;
    block_out[6] <= 32'b0;
    block_out[7] <= 32'b0;
    block_out_tmp <= 32'b0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= (o_t0 + o_t7) >>> OUT_SHIFT;
        block_out_tmp <= (o_t0 - o_t7) >>> OUT_SHIFT; // block_out[7]
        block_out[3] <= (o_t3 + o_t4) >>> OUT_SHIFT;
        block_out[4] <= (o_t3 - o_t4) >>> OUT_SHIFT;
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= (o_t2 + o_s5) >>> OUT_SHIFT;
        block_out[5] <= (o_t2 - o_s5) >>> OUT_SHIFT;
        block_out[1] <= (o_t1 + o_s6) >>> OUT_SHIFT;
        block_out[6] <= (o_t1 - o_s6) >>> OUT_SHIFT;
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

assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

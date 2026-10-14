//-----------------------------------------------------------------
// jpeg_idct_2d — Optimized 8x8 2-D IDCT (AAN-style)
// Based on ultraembedded/core_jpeg (Apache 2.0)
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
wire [ 31:0]  idct_x_data_w;
wire [  5:0]  idct_x_idx_w;
wire          idct_x_accept_w;

jpeg_idct_1d #(
    .OUT_SHIFT(11),
    .IS_Y(0)
) u_idct_x (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i({{16{inport_data0_i[15]}}, inport_data0_i})
    ,.inport_data1_i({{16{inport_data1_i[15]}}, inport_data1_i})
    ,.inport_data2_i({{16{inport_data2_i[15]}}, inport_data2_i})
    ,.inport_data3_i({{16{inport_data3_i[15]}}, inport_data3_i})
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(idct_x_valid_w)
    ,.outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;
wire          transpose_ready_w = 1'b1;

jpeg_idct_transpose u_transpose (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(idct_x_valid_w)
    ,.inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w)
    ,.inport_accept_o(idct_x_accept_w)
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
    ,.outport_ready_i(transpose_ready_w)
);

jpeg_idct_1d #(
    .OUT_SHIFT(15),
    .IS_Y(1)
) u_idct_y (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
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
// Parameterized 1-D IDCT (AAN-style)
//-----------------------------------------------------------------
module jpeg_idct_1d #(
    parameter OUT_SHIFT = 11,
    parameter IS_Y = 0
) (
    input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data0_i
    ,input  [ 31:0]  inport_data1_i
    ,input  [ 31:0]  inport_data2_i
    ,input  [ 31:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam signed [12:0] C1 = 13'sd4017;
localparam signed [12:0] C2 = 13'sd3784;
localparam signed [12:0] C3 = 13'sd3406;
localparam signed [12:0] C4 = 13'sd2896;
localparam signed [12:0] C5 = 13'sd2276;
localparam signed [12:0] C6 = 13'sd1567;
localparam signed [12:0] C7 = 13'sd799;

wire signed [31:0] block_in_0_1 = inport_data0_i;
wire signed [31:0] block_in_2_3 = inport_data1_i;
wire signed [31:0] block_in_4_5 = inport_data2_i;
wire signed [31:0] block_in_6_7 = inport_data3_i;

reg signed [31:0] i0;
reg signed [31:0] mul0_a; reg signed [12:0] mul0_b;
reg signed [31:0] mul1_a; reg signed [12:0] mul1_b;
reg signed [31:0] mul2_a; reg signed [12:0] mul2_b;
reg signed [31:0] mul3_a; reg signed [12:0] mul3_b;
reg signed [31:0] mul4_a; reg signed [12:0] mul4_b;

always @(posedge clk)
if (rst) begin
    i0 <= 32'b0;
    mul0_a <= 32'b0; mul0_b <= 13'b0;
    mul1_a <= 32'b0; mul1_b <= 13'b0;
    mul2_a <= 32'b0; mul2_b <= 13'b0;
    mul3_a <= 32'b0; mul3_b <= 13'b0;
    mul4_a <= 32'b0; mul4_b <= 13'b0;
end else begin
    case (inport_idx_i)
    3'd0: begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;  mul0_b <= C2;
        mul1_a <= block_in_6_7;   mul1_b <= C6;
    end
    3'd1: begin
        mul0_a <= block_in_0_1;  mul0_b <= C1;
        mul1_a <= block_in_6_7;   mul1_b <= C7;
        mul2_a <= block_in_4_5;  mul2_b <= C5;
        mul3_a <= block_in_2_3;  mul3_b <= C3;
        mul4_a <= i0;            mul4_b <= C4;
    end
    3'd2: begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3: begin
        mul0_a <= block_in_0_1;  mul0_b <= C7;
        mul1_a <= block_in_6_7;   mul1_b <= C1;
        mul2_a <= block_in_4_5;  mul2_b <= C3;
        mul3_a <= block_in_2_3;  mul3_b <= C5;
    end
    3'd4: begin
        mul0_a <= block_in_0_1;  mul0_b <= C7;
        mul1_a <= block_in_6_7;   mul1_b <= C1;
        mul2_a <= block_in_4_5;  mul2_b <= C3;
        mul3_a <= block_in_2_3;  mul3_b <= C5;
    end
    3'd5: begin
        mul0_a <= block_in_2_3;  mul0_b <= C6;
        mul1_a <= block_in_6_7;   mul1_b <= C2;
        mul4_a <= i0;            mul4_b <= C4;
    end
    default: ;
    endcase
end

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

always @(posedge clk)
if (rst) begin
    mul0_q <= 32'b0; mul1_q <= 32'b0; mul2_q <= 32'b0;
    mul3_q <= 32'b0; mul4_q <= 32'b0;
end else begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

always @(posedge clk)
if (rst) begin
    mul0 <= 32'b0; mul1 <= 32'b0; mul2 <= 32'b0;
    mul3 <= 32'b0; mul4 <= 32'b0;
end else begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q;
    mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q; reg [2:0] out_stg0_idx_q;
reg        out_stg1_valid_q; reg [2:0] out_stg1_idx_q;
reg        out_stg2_valid_q; reg [2:0] out_stg2_idx_q;

always @(posedge clk)
if (rst) begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end else begin
    out_stg0_valid_q <= inport_valid_i;  out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q;
end

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

wire signed [31:0] prod_6_5 = o_t6_5 * 32'sd181;
wire signed [31:0] prod_5_6 = o_t5_6 * 32'sd181;
wire signed [23:0] hi_6_5 = prod_6_5[31:8];
wire signed [23:0] hi_5_6 = prod_5_6[31:8];
wire corr_6_5 = prod_6_5[31] & |prod_6_5[7:0];
wire corr_5_6 = prod_5_6[31] & |prod_5_6[7:0];

always @(posedge clk)
if (rst) begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0;
    o_t0 <= 0; o_t1 <= 0; o_t2 <= 0; o_t3 <= 0; o_t4 <= 0;
    o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end else begin
    case (out_stg2_idx_q)
    3'd0: begin
        o_t3 <= mul0 + mul1;
    end
    3'd1: begin
        o_s7 <= mul0 + mul1;
        o_s6 <= mul2 + mul3;
        o_t0 <= mul4;
    end
    3'd2: begin
        o_t0 <= o_t0 + o_t3;
        o_t3 <= o_t0 - o_t3;
        o_t7 <= o_s6 + o_s7;
    end
    3'd3: begin
        o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    end
    3'd4: begin
        o_t0 <= mul0 - mul1;
        o_s5 <= mul2 - mul3;
    end
    3'd5: begin
        o_t3 <= mul0 - mul1;
        o_t4 <= mul4;
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6: begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5;
        o_t5_6 <= o_t5 + o_t6;
    end
    default: begin
        o_s5 <= {{8{hi_6_5[23]}}, hi_6_5} + corr_6_5;
        o_s6 <= {{8{hi_5_6[23]}}, hi_5_6} + corr_5_6;
    end
    endcase
end

reg        out_stg3_valid_q; reg [2:0] out_stg3_idx_q;

always @(posedge clk)
if (rst) begin
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
end else begin
    out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q;
end

reg signed [31:0] block_out [0:7];
reg signed [31:0] block_out_tmp;

always @(posedge clk)
if (rst) begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end else if (out_stg3_valid_q) begin
    if (out_stg3_idx_q == 3'd3) begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end
    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;
    if (out_stg3_idx_q == 3'd7) begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;
always @(posedge clk)
if (rst)        valid_q <= 8'b0;
else if (img_start_i) valid_q <= 8'b0;
else            valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;
always @(posedge clk)
if (rst)        ptr_q <= 6'd0;
else if (img_start_i) ptr_q <= 6'd0;
else if (outport_valid_o) ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];

generate
if (IS_Y) begin : y_idx
    assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
end else begin : x_idx
    assign outport_idx_o = ptr_q;
end
endgenerate

endmodule

//-----------------------------------------------------------------
// Transpose (corner-turn) with double-buffered RAMs
//-----------------------------------------------------------------
module jpeg_idct_transpose (
    input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data0_o
    ,output [ 31:0]  outport_data1_o
    ,output [ 31:0]  outport_data2_o
    ,output [ 31:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire [31:0] outport_data0_w;
wire [31:0] outport_data1_w;
wire [31:0] outport_data2_w;
wire [31:0] outport_data3_w;

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram u_ram0 (
     .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr0_w), .data0_o()
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data0_w)
);
jpeg_idct_transpose_ram u_ram1 (
     .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr1_w), .data0_o()
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data1_w)
);
jpeg_idct_transpose_ram u_ram2 (
     .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr2_w), .data0_o()
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data2_w)
);
jpeg_idct_transpose_ram u_ram3 (
     .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr3_w), .data0_o()
    ,.addr1_i({block_rd_q, rd_addr_q}), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data3_w)
);

reg [1:0] block_ready_q;

always @(posedge clk)
if (rst) begin
    block_ready_q <= 2'b0;
    block_wr_q     <= 1'b0;
    block_rd_q     <= 1'b0;
end else if (img_start_i) begin
    block_ready_q <= 2'b0;
    block_wr_q     <= 1'b0;
    block_rd_q     <= 1'b0;
end else begin
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o) begin
        block_ready_q[block_wr_q] <= 1'b1;
        block_wr_q                 <= ~block_wr_q;
    end
    if (outport_valid_o && rd_idx_q[5:0] == 6'd63) begin
        block_ready_q[block_rd_q] <= 1'b0;
        block_rd_q <= ~block_rd_q;
    end
end

assign inport_accept_o = ~block_ready_q[block_wr_q];

localparam STATE_IDLE = 2'd0;
localparam STATE_SETUP = 2'd1;
localparam STATE_ACTIVE = 2'd2;

reg [1:0] state_q;
reg [1:0] next_state_r;

always @* begin
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

always @(posedge clk)
if (rst) state_q <= STATE_IDLE;
else     state_q <= next_state_r;

always @(posedge clk)
if (rst)               rd_idx_q <= 6'b0;
else if (img_start_i)  rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE) rd_idx_q <= rd_idx_q + 6'd1;

always @(posedge clk)
if (rst)               rd_addr_q <= 4'b0;
else if (state_q == STATE_IDLE)  rd_addr_q <= 4'b0;
else if (state_q == STATE_SETUP) rd_addr_q <= 4'd8;
else if (state_q == STATE_ACTIVE) begin
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
// Dual-port RAM (register-based)
//-----------------------------------------------------------------
module jpeg_idct_transpose_ram (
    input           clk0_i
    ,input           rst0_i
    ,input  [  4:0]  addr0_i
    ,input  [ 31:0]  data0_i
    ,input           wr0_i
    ,input           clk1_i
    ,input           rst1_i
    ,input  [  4:0]  addr1_i
    ,input  [ 31:0]  data1_i
    ,input           wr1_i
    ,output [ 31:0]  data0_o
    ,output [ 31:0]  data1_o
);

reg [31:0] ram [0:31];
reg [31:0] ram_read0_q;
reg [31:0] ram_read1_q;

always @(posedge clk0_i) begin
    if (wr0_i)
        ram[addr0_i] <= data0_i;
    ram_read0_q <= ram[addr0_i];
end

always @(posedge clk1_i) begin
    if (wr1_i)
        ram[addr1_i] <= data1_i;
    ram_read1_q <= ram[addr1_i];
end

assign data0_o = ram_read0_q;
assign data1_o = ram_read1_q;

endmodule

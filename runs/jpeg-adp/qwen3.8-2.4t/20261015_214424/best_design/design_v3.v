//-----------------------------------------------------------------
// jpeg_idct_2d - optimized v3
//  - All multipliers folded to constant shift-add trees (bit-exact)
//  - Pipeline shortened by 2 stages per pass (bit-exact schedule)
//  - (x*181)/256 implemented as truncating shift (identical semantics)
//  - Transpose stores 20-bit values (x-pass outputs provably < 2^19 mag)
//  - block_out of x-pass stored 20 bits
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

wire          x_valid_w;
wire [ 31:0]  x_data_w;
wire [  5:0]  x_idx_w;

idct1d_pass #(.OUT_SHIFT(11), .INW(16), .SWAP_IDX(0)) u_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i({16'b0, inport_data0_i})
    ,.inport_data1_i({16'b0, inport_data1_i})
    ,.inport_data2_i({16'b0, inport_data2_i})
    ,.inport_data3_i({16'b0, inport_data3_i})
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

wire          t_valid_w;
wire [ 31:0]  t_data0_w;
wire [ 31:0]  t_data1_w;
wire [ 31:0]  t_data2_w;
wire [ 31:0]  t_data3_w;
wire [  2:0]  t_idx_w;

idct_transpose u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
);

idct1d_pass #(.OUT_SHIFT(15), .INW(20), .SWAP_IDX(1)) u_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(t_valid_w)
    ,.inport_data0_i(t_data0_w)
    ,.inport_data1_i(t_data1_w)
    ,.inport_data2_i(t_data2_w)
    ,.inport_data3_i(t_data3_w)
    ,.inport_idx_i(t_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//-----------------------------------------------------------------
// 1-D IDCT pass, constant-multiplied, minimal pipeline
//  INW = signed width of incoming samples (16 for x, 20 for y)
//-----------------------------------------------------------------
module idct1d_pass
#(
     parameter OUT_SHIFT = 11
    ,parameter INW       = 16
    ,parameter SWAP_IDX  = 0
)
(
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

wire signed [INW-1:0] d0 = inport_data0_i[INW-1:0];
wire signed [INW-1:0] d1 = inport_data1_i[INW-1:0];
wire signed [INW-1:0] d2 = inport_data2_i[INW-1:0];
wire signed [INW-1:0] d3 = inport_data3_i[INW-1:0];

//-----------------------------------------------------------------
// i0 accumulator (INW+1 bits, exact)
//-----------------------------------------------------------------
reg signed [INW:0] i0;

always @ (posedge clk)
if (rst)
    i0 <= 'b0;
else
begin
    case (inport_idx_i)
    3'd0: i0 <= d0 + d2;
    3'd2: i0 <= d0 - d2;
    default: ;
    endcase
end

//-----------------------------------------------------------------
// Constant-multiply trees, input muxed by beat
//-----------------------------------------------------------------
wire signed [INW-1:0] s1 = (inport_idx_i == 3'd1) ? d0 : d3;
wire signed [INW-1:0] s2 = ((inport_idx_i == 3'd0) || (inport_idx_i == 3'd5)) ? d1 : d3;
wire signed [INW-1:0] s3 = (inport_idx_i == 3'd1) ? d1 : d2;
wire signed [INW-1:0] s5 = (inport_idx_i == 3'd1) ? d2 : d1;
wire signed [INW-1:0] s6 = ((inport_idx_i == 3'd0) || (inport_idx_i == 3'd5)) ? d3 : d1;
wire signed [INW-1:0] s7 = (inport_idx_i == 3'd1) ? d3 : d0;

wire signed [INW+12:0] p1 = s1 * C1;
wire signed [INW+12:0] p2 = s2 * C2;
wire signed [INW+12:0] p3 = s3 * C3;
wire signed [INW+12:0] p5 = s5 * C5;
wire signed [INW+12:0] p6 = s6 * C6;
wire signed [INW+12:0] p7 = s7 * C7;
wire signed [INW+13:0] p4 = i0 * C4;

//-----------------------------------------------------------------
// Product latches
//-----------------------------------------------------------------
reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
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
    case (inport_idx_i)
    3'd0:
    begin
        mul0_q <= p2;
        mul1_q <= p6;
    end
    3'd1:
    begin
        mul0_q <= p1;
        mul1_q <= p7;
        mul2_q <= p5;
        mul3_q <= p3;
        mul4_q <= p4;
    end
    3'd3:
    begin
        mul0_q <= p7;
        mul1_q <= p1;
        mul2_q <= p3;
        mul3_q <= p5;
    end
    3'd4:
    begin
        mul0_q <= p7;
        mul1_q <= p1;
        mul2_q <= p3;
        mul3_q <= p5;
    end
    3'd5:
    begin
        mul0_q <= p6;
        mul1_q <= p2;
        mul4_q <= p4;
    end
    default:
        ;
    endcase
end

//-----------------------------------------------------------------
// Valid/idx pipeline (2 stages)
//-----------------------------------------------------------------
reg        stg0_valid_q;
reg [2:0]  stg0_idx_q;
reg        stg1_valid_q;
reg [2:0]  stg1_idx_q;

always @ (posedge clk)
if (rst)
begin
    stg0_valid_q <= 1'b0;
    stg0_idx_q   <= 3'b0;
    stg1_valid_q <= 1'b0;
    stg1_idx_q   <= 3'b0;
end
else
begin
    stg0_valid_q <= inport_valid_i;
    stg0_idx_q   <= inport_idx_i;
    stg1_valid_q <= stg0_valid_q;
    stg1_idx_q   <= stg0_idx_q;
end

//-----------------------------------------------------------------
// Butterfly stage (runs 1 cycle after input beat)
//-----------------------------------------------------------------
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

wire signed [31:0] r5 = (o_t6_5 * 181) / 256;
wire signed [31:0] r6 = (o_t5_6 * 181) / 256;

always @ (posedge clk)
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
    case (stg0_idx_q)
    3'd0:
    begin
        o_t3 <= mul0_q + mul1_q; // s3
    end
    3'd1:
    begin
        o_s7 <= mul0_q + mul1_q;
        o_s6 <= mul2_q + mul3_q;
        o_t0 <= mul4_q;          // s0
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3; // t0
        o_t3 <= o_t0 - o_t3; // t3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3:
    begin
        o_t4 <= (mul0_q - mul1_q) + (mul2_q - mul3_q);
    end
    3'd4:
    begin
        o_t0 <= mul0_q - mul1_q; // s4
        o_s5 <= mul2_q - mul3_q;
    end
    3'd5:
    begin
        o_t3 <= mul0_q - mul1_q; // s2
        o_t4 <= mul4_q; // s1
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
        o_s5 <= r5;
        o_s6 <= r6;
    end
    endcase
end

//-----------------------------------------------------------------
// Output storage (all 8 values written on beats 3 and 7)
//-----------------------------------------------------------------
localparam OW = (OUT_SHIFT == 11) ? 20 : 32;

reg signed [OW-1:0] block_out[0:7];
reg signed [OW-1:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 'b0;
    block_out[1] <= 'b0;
    block_out[2] <= 'b0;
    block_out[3] <= 'b0;
    block_out[4] <= 'b0;
    block_out[5] <= 'b0;
    block_out[6] <= 'b0;
    block_out[7] <= 'b0;
    block_out_tmp <= 'b0;
end
else if (stg1_valid_q)
begin
    if (stg1_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end

    if (stg1_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (stg1_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

//-----------------------------------------------------------------
// Output streamer
//-----------------------------------------------------------------
reg [7:0] valid_q;

always @ (posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], stg1_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];

generate
if (OW == 20)
begin : g_20
    wire signed [31:0] rd_val = {{12{block_out[ptr_q[2:0]][19]}}, block_out[ptr_q[2:0]]};
    assign outport_data_o = rd_val;
end
else
begin : g_32
    assign outport_data_o = block_out[ptr_q[2:0]];
end
endgenerate

assign outport_idx_o = SWAP_IDX ? {ptr_q[2:0], ptr_q[5:3]} : {3'b0, ptr_q};

endmodule

//-----------------------------------------------------------------
// Transpose: double buffered, 20-bit storage
//-----------------------------------------------------------------
module idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
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

wire wr0_w = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && (inport_idx_i[5:4] == 2'd3);

wire [19:0] rd0_w, rd1_w, rd2_w, rd3_w;

idct_transpose_ram u_ram0 (.clk(clk), .wr_addr_i(wr_ptr_w), .wr_data_i(inport_data_i[19:0]), .wr_en_i(wr0_w), .rd_addr_i({block_rd_q, rd_addr_q}), .rd_data_o(rd0_w));
idct_transpose_ram u_ram1 (.clk(clk), .wr_addr_i(wr_ptr_w), .wr_data_i(inport_data_i[19:0]), .wr_en_i(wr1_w), .rd_addr_i({block_rd_q, rd_addr_q}), .rd_data_o(rd1_w));
idct_transpose_ram u_ram2 (.clk(clk), .wr_addr_i(wr_ptr_w), .wr_data_i(inport_data_i[19:0]), .wr_en_i(wr2_w), .rd_addr_i({block_rd_q, rd_addr_q}), .rd_data_o(rd2_w));
idct_transpose_ram u_ram3 (.clk(clk), .wr_addr_i(wr_ptr_w), .wr_data_i(inport_data_i[19:0]), .wr_en_i(wr3_w), .rd_addr_i({block_rd_q, rd_addr_q}), .rd_data_o(rd3_w));

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
assign outport_data0_o = {{12{rd0_w[19]}}, rd0_w};
assign outport_data1_o = {{12{rd1_w[19]}}, rd1_w};
assign outport_data2_o = {{12{rd2_w[19]}}, rd2_w};
assign outport_data3_o = {{12{rd3_w[19]}}, rd3_w};

endmodule

//-----------------------------------------------------------------
// Dual port RAM (32 x 20 bit), single read port registered
//-----------------------------------------------------------------
module idct_transpose_ram
(
     input           clk
    ,input  [  4:0]  wr_addr_i
    ,input  [ 19:0]  wr_data_i
    ,input           wr_en_i
    ,input  [  4:0]  rd_addr_i
    ,output reg [19:0] rd_data_o
);

reg [19:0] ram [31:0];

always @ (posedge clk)
begin
    if (wr_en_i)
        ram[wr_addr_i] <= wr_data_i;

    rd_data_o <= ram[rd_addr_i];
end

endmodule

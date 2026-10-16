//-----------------------------------------------------------------
// jpeg_idct_2d - optimized v2
// Same architecture as core_jpeg (x-pass -> transpose -> y-pass),
// but all multipliers are constant-multipliers (shift-add trees),
// since in the original schedule one operand is always a fixed
// cosine constant. Bit-exact with the original.
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

jpeg_idct_pass_x
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
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
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(idct_x_valid_w)
    ,.inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w)
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
);

jpeg_idct_pass_y
u_idct_y
(
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
// Shared 1-D IDCT datapath (parameterized shift). All multiplies
// are by compile-time constants -> shift/add trees in synthesis.
//-----------------------------------------------------------------
module jpeg_idct_pass_core
#(
     parameter OUT_SHIFT = 11
    ,parameter INPUT_16B = 1
    ,parameter IDX_SWAP  = 0
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

localparam signed [31:0] C1 = 32'sd4017;
localparam signed [31:0] C2 = 32'sd3784;
localparam signed [31:0] C3 = 32'sd3406;
localparam signed [31:0] C4 = 32'sd2896;
localparam signed [31:0] C5 = 32'sd2276;
localparam signed [31:0] C6 = 32'sd1567;
localparam signed [31:0] C7 = 32'sd799;

wire signed [31:0] d0 = INPUT_16B ? { {16{inport_data0_i[15]}}, inport_data0_i[15:0]} : inport_data0_i;
wire signed [31:0] d1 = INPUT_16B ? { {16{inport_data1_i[15]}}, inport_data1_i[15:0]} : inport_data1_i;
wire signed [31:0] d2 = INPUT_16B ? { {16{inport_data2_i[15]}}, inport_data2_i[15:0]} : inport_data2_i;
wire signed [31:0] d3 = INPUT_16B ? { {16{inport_data3_i[15]}}, inport_data3_i[15:0]} : inport_data3_i;

//-----------------------------------------------------------------
// i0 accumulator
//-----------------------------------------------------------------
reg signed [31:0] i0;

always @ (posedge clk)
if (rst)
    i0 <= 32'b0;
else
begin
    case (inport_idx_i)
    3'd0: i0 <= d0 + d2;
    3'd2: i0 <= d0 - d2;
    default: ;
    endcase
end

//-----------------------------------------------------------------
// Constant multipliers (7 trees per pass), input-muxed by beat
//-----------------------------------------------------------------
wire signed [31:0] s1 = (inport_idx_i == 3'd1) ? d0 : d3;
wire signed [31:0] s2 = (inport_idx_i == 3'd0) ? d1 : d3;
wire signed [31:0] s3 = (inport_idx_i == 3'd1) ? d1 : d2;
wire signed [31:0] s5 = (inport_idx_i == 3'd1) ? d2 : d1;
wire signed [31:0] s6 = (inport_idx_i == 3'd0) ? d3 : d1;
wire signed [31:0] s7 = (inport_idx_i == 3'd1) ? d3 : d0;

/* verilator lint_off UNUSED */
wire signed [31:0] p1 = s1 * C1;
wire signed [31:0] p2 = s2 * C2;
wire signed [31:0] p3 = s3 * C3;
wire signed [31:0] p4 = i0 * C4;
wire signed [31:0] p5 = s5 * C5;
wire signed [31:0] p6 = s6 * C6;
wire signed [31:0] p7 = s7 * C7;
/* verilator lint_on UNUSED */

//-----------------------------------------------------------------
// Product latch (same values as original mul_a/mul_b + mul_q)
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
// Extra product pipeline stage (as original)
//-----------------------------------------------------------------
reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
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

//-----------------------------------------------------------------
// Valid/idx pipeline
//-----------------------------------------------------------------
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
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
    out_stg3_valid_q <= 1'b0;
    out_stg3_idx_q   <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q;
    out_stg2_idx_q   <= out_stg1_idx_q;
    out_stg3_valid_q <= out_stg2_valid_q;
    out_stg3_idx_q   <= out_stg2_idx_q;
end

//-----------------------------------------------------------------
// Butterfly accumulator stage
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

//-----------------------------------------------------------------
// Output storage
//-----------------------------------------------------------------
reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
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
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
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
assign outport_idx_o   = IDX_SWAP ? {ptr_q[2:0], ptr_q[5:3]} : {3'b0, ptr_q};

endmodule

//-----------------------------------------------------------------
// X pass wrapper (16-bit sign-extended inputs)
//-----------------------------------------------------------------
module jpeg_idct_pass_x
(
     input           clk
    ,input           rst
    ,input           img_start_i
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

jpeg_idct_pass_core
#(
    .OUT_SHIFT(11)
   ,.INPUT_16B(1)
   ,.IDX_SWAP(0)
)
u_core
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
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//-----------------------------------------------------------------
// Y pass wrapper (32-bit inputs, swapped output idx)
//-----------------------------------------------------------------
module jpeg_idct_pass_y
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

jpeg_idct_pass_core
#(
    .OUT_SHIFT(15)
   ,.INPUT_16B(0)
   ,.IDX_SWAP(1)
)
u_core
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(inport_data0_i)
    ,.inport_data1_i(inport_data1_i)
    ,.inport_data2_i(inport_data2_i)
    ,.inport_data3_i(inport_data3_i)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//-----------------------------------------------------------------
// Transpose (original structure)
//-----------------------------------------------------------------
module jpeg_idct_transpose
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

wire [31:0] outport_data0_w;
wire [31:0] outport_data1_w;
wire [31:0] outport_data2_w;
wire [31:0] outport_data3_w;

wire wr0_w = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram
u_ram0
(
     .clk_i(clk)
    ,.wr_addr_i(wr_ptr_w)
    ,.wr_data_i(inport_data_i)
    ,.wr_en_i(wr0_w)
    ,.rd_addr_i({block_rd_q, rd_addr_q})
    ,.rd_data_o(outport_data0_w)
);

jpeg_idct_transpose_ram
u_ram1
(
     .clk_i(clk)
    ,.wr_addr_i(wr_ptr_w)
    ,.wr_data_i(inport_data_i)
    ,.wr_en_i(wr1_w)
    ,.rd_addr_i({block_rd_q, rd_addr_q})
    ,.rd_data_o(outport_data1_w)
);

jpeg_idct_transpose_ram
u_ram2
(
     .clk_i(clk)
    ,.wr_addr_i(wr_ptr_w)
    ,.wr_data_i(inport_data_i)
    ,.wr_en_i(wr2_w)
    ,.rd_addr_i({block_rd_q, rd_addr_q})
    ,.rd_data_o(outport_data2_w)
);

jpeg_idct_transpose_ram
u_ram3
(
     .clk_i(clk)
    ,.wr_addr_i(wr_ptr_w)
    ,.wr_data_i(inport_data_i)
    ,.wr_en_i(wr3_w)
    ,.rd_addr_i({block_rd_q, rd_addr_q})
    ,.rd_data_o(outport_data3_w)
);

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;
wire      inport_accept_w = ~block_ready_q[block_wr_q];

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
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_w)
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
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

//-----------------------------------------------------------------
// Dual port RAM (write port + registered read port)
//-----------------------------------------------------------------
module jpeg_idct_transpose_ram
(
     input           clk_i
    ,input  [  4:0]  wr_addr_i
    ,input  [ 31:0]  wr_data_i
    ,input           wr_en_i
    ,input  [  4:0]  rd_addr_i
    ,output [ 31:0]  rd_data_o
);

reg [31:0] ram [31:0];
reg [31:0] ram_read_q;

always @ (posedge clk_i)
begin
    if (wr_en_i)
        ram[wr_addr_i] <= wr_data_i;

    ram_read_q <= ram[rd_addr_i];
end

assign rd_data_o = ram_read_q;

endmodule

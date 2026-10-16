//-----------------------------------------------------------------
// jpeg_idct_2d - v8: v6 + transpose read-path restructure.
// Each 32x20 RAM is split into 4 sub-banks of 8 entries
// ({block, addr[3]} select, addr[2:0] = group index). The read
// address is decoded combinationally from rd_idx_q (no rd_addr_q
// register, no +/-8 carry chain): addr(j,g) = g+8*b(j), except
// beat 7 which reads g+1 (b(j) = [1,0,1,1,0,1,0,0][j]).
// Read mux depth: 8:1 + 4:1 instead of 32:1.
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
wire [ 19:0]  transpose_data0_w;
wire [ 19:0]  transpose_data1_w;
wire [ 19:0]  transpose_data2_w;
wire [ 19:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose u_transpose
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

jpeg_idct_y u_idct_y
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

module const_mul_sel
(
    input  [2:0] idx,
    output reg signed [31:0] b0,
    output reg signed [31:0] b1,
    output reg signed [31:0] b2,
    output reg signed [31:0] b3
);
localparam signed [31:0] C1 = 32'sd4017;
localparam signed [31:0] C2 = 32'sd3784;
localparam signed [31:0] C3 = 32'sd3406;
localparam signed [31:0] C5 = 32'sd2276;
localparam signed [31:0] C6 = 32'sd1567;
localparam signed [31:0] C7 = 32'sd799;

always @ *
begin
    case (idx)
    3'd0: begin b0 = C2; b1 = C6; b2 = C3; b3 = C5; end
    3'd1: begin b0 = C1; b1 = C7; b2 = C5; b3 = C3; end
    3'd2: begin b0 = C1; b1 = C7; b2 = C5; b3 = C3; end
    3'd3: begin b0 = C7; b1 = C1; b2 = C3; b3 = C5; end
    3'd4: begin b0 = C7; b1 = C1; b2 = C3; b3 = C5; end
    3'd5: begin b0 = C6; b1 = C2; b2 = C3; b3 = C5; end
    3'd6: begin b0 = C6; b1 = C2; b2 = C3; b3 = C5; end
    default: begin b0 = C6; b1 = C2; b2 = C3; b3 = C5; end
    endcase
end
endmodule

//-----------------------------------------------------------------
// X pass
//-----------------------------------------------------------------
module jpeg_idct_x
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
    ,output [ 19:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam OUT_SHIFT = 11;
localparam signed [31:0] C4_32 = 32'sd2896;

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg signed [31:0] mul1_a;
reg signed [31:0] mul2_a;
reg signed [31:0] mul3_a;
reg signed [31:0] mul4_a;

always @ (posedge clk )
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0;
    mul1_a <= 32'b0;
    mul2_a <= 32'b0;
    mul3_a <= 32'b0;
    mul4_a <= 32'b0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul1_a <= block_in_6_7;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
        mul4_a <= i0;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul1_a <= block_in_6_7;
        mul4_a <= i0;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk )
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

wire signed [31:0] b0_w, b1_w, b2_w, b3_w;
const_mul_sel u_sel (.idx(out_stg0_idx_q), .b0(b0_w), .b1(b1_w), .b2(b2_w), .b3(b3_w));

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
    mul0_q <= mul0_a * b0_w;
    mul1_q <= mul1_a * b1_w;
    mul2_q <= mul2_a * b2_w;
    mul3_q <= mul3_a * b3_w;
    mul4_q <= mul4_a * C4_32;
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

// 1/sqrt(2) scaling, pipelined: pre5/pre6 (= t6-t5, t5+t6) are formed
// one cycle earlier so the x181 product is registered; the /256 (trunc
// toward zero) is then a short shift/correct path.
wire signed [31:0] pre5_w = (o_s7 - o_s6) - (o_t0 - o_s5);
wire signed [31:0] pre6_w = (o_t0 - o_s5) + (o_s7 - o_s6);

reg signed [31:0] p5_q;
reg signed [31:0] p6_q;

wire [23:0] adj5_w = p5_q[31:8] + {23'd0, (p5_q[31] & (p5_q[7:0] != 8'd0))};
wire [23:0] adj6_w = p6_q[31:8] + {23'd0, (p6_q[31] & (p6_q[7:0] != 8'd0))};
wire signed [31:0] r5_w = {{8{adj5_w[23]}}, adj5_w};
wire signed [31:0] r6_w = {{8{adj6_w[23]}}, adj6_w};

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
    p5_q   <= 32'b0;
    p6_q   <= 32'b0;
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
        p5_q <= pre5_w * 32'sd181; // (t6-t5)*181
        p6_q <= pre6_w * 32'sd181; // (t5+t6)*181
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_s5 <= r5_w; // 1/sqrt(2), signed trunc toward zero
        o_s6 <= r6_w;
    end
    default: ;
    endcase
end

reg signed [19:0] block_out[0:7];
reg signed [19:0] block_out_tmp;

always @ (posedge clk )
if (rst)
begin
    block_out[0] <= 20'b0;
    block_out[1] <= 20'b0;
    block_out[2] <= 20'b0;
    block_out[3] <= 20'b0;
    block_out[4] <= 20'b0;
    block_out[5] <= 20'b0;
    block_out[6] <= 20'b0;
    block_out[7] <= 20'b0;
    block_out_tmp <= 20'b0;
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
// Transpose with 4-way banked RAMs
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
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

//-----------------------------------------------------------------
// Read address decode. The registered RAM read latches the decode
// one cycle early, so during ACTIVE beat k the RAM input must carry
// the address of beat k+1: addr(beat) = sub=g, hi=[0,1,0,1,1,0,1,0]
//   -> live decode: hi(j)=[1,0,1,1,0,1,0,0][j], sub=g+(j==7)
// During SETUP the original presents the carried-over address
// {hi=0, sub=0} (needed for beat 0 of the new block).
//-----------------------------------------------------------------
wire        active_w = (state_q == STATE_ACTIVE);
wire [2:0]  j_w      = rd_idx_q[2:0];
wire [2:0]  g_w      = rd_idx_q[5:3];

wire        hi_w     = (j_w == 3'd0) || (j_w == 3'd2) ||
                       (j_w == 3'd3) || (j_w == 3'd5);
wire [2:0]  sub_w    = g_w + {2'd0, (j_w == 3'd7)};

wire [1:0]  bank_sel_w = {block_rd_q, active_w & hi_w};
wire [2:0]  sub_addr_w = active_w ? sub_w : 3'd0;

wire wr_r0 = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr_r1 = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr_r2 = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr_r3 = inport_valid_i && (inport_idx_i[5:4] == 2'd3);

idct_ram20_bank u_ram0 (.clk(clk), .wr_en_i(wr_r0), .wr_bank_i({block_wr_q, inport_idx_i[3]}), .wr_addr_i(inport_idx_i[2:0]), .wr_data_i(inport_data_i), .rd_bank_i(bank_sel_w), .rd_addr_i(sub_addr_w), .rd_data_o(outport_data0_o));
idct_ram20_bank u_ram1 (.clk(clk), .wr_en_i(wr_r1), .wr_bank_i({block_wr_q, inport_idx_i[3]}), .wr_addr_i(inport_idx_i[2:0]), .wr_data_i(inport_data_i), .rd_bank_i(bank_sel_w), .rd_addr_i(sub_addr_w), .rd_data_o(outport_data1_o));
idct_ram20_bank u_ram2 (.clk(clk), .wr_en_i(wr_r2), .wr_bank_i({block_wr_q, inport_idx_i[3]}), .wr_addr_i(inport_idx_i[2:0]), .wr_data_i(inport_data_i), .rd_bank_i(bank_sel_w), .rd_addr_i(sub_addr_w), .rd_data_o(outport_data2_o));
idct_ram20_bank u_ram3 (.clk(clk), .wr_en_i(wr_r3), .wr_bank_i({block_wr_q, inport_idx_i[3]}), .wr_addr_i(inport_idx_i[2:0]), .wr_data_i(inport_data_i), .rd_bank_i(bank_sel_w), .rd_addr_i(sub_addr_w), .rd_data_o(outport_data3_o));

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

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

endmodule

//-----------------------------------------------------------------
// Banked RAM: 4 banks x 8 entries x 20 bits, registered read
//-----------------------------------------------------------------
module idct_ram20_bank
(
     input           clk
    ,input           wr_en_i
    ,input  [  1:0]  wr_bank_i
    ,input  [  2:0]  wr_addr_i
    ,input  [ 19:0]  wr_data_i
    ,input  [  1:0]  rd_bank_i
    ,input  [  2:0]  rd_addr_i
    ,output reg [19:0] rd_data_o
);

reg [19:0] ram [3:0][7:0];

always @ (posedge clk)
begin
    if (wr_en_i)
        ram[wr_bank_i][wr_addr_i] <= wr_data_i;

    rd_data_o <= ram[rd_bank_i][rd_addr_i];
end

endmodule

//-----------------------------------------------------------------
// Y pass
//-----------------------------------------------------------------
module jpeg_idct_y
(
     input           clk
    ,input           rst
    ,input           img_start_i
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
localparam signed [31:0] C4_32 = 32'sd2896;

wire signed [31:0] block_in_0_1 = {{12{inport_data0_i[19]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{12{inport_data1_i[19]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{12{inport_data2_i[19]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{12{inport_data3_i[19]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg signed [31:0] mul1_a;
reg signed [31:0] mul2_a;
reg signed [31:0] mul3_a;
reg signed [31:0] mul4_a;

always @ (posedge clk )
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0;
    mul1_a <= 32'b0;
    mul2_a <= 32'b0;
    mul3_a <= 32'b0;
    mul4_a <= 32'b0;
end
else
begin
    /* verilator lint_off WIDTH */
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul1_a <= block_in_6_7;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
        mul4_a <= i0;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul1_a <= block_in_6_7;
        mul2_a <= block_in_4_5;
        mul3_a <= block_in_2_3;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul1_a <= block_in_6_7;
        mul4_a <= i0;
    end
    default:
        ;
    endcase
    /* verilator lint_on WIDTH */
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk )
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

wire signed [31:0] b0_w, b1_w, b2_w, b3_w;
const_mul_sel u_sel (.idx(out_stg0_idx_q), .b0(b0_w), .b1(b1_w), .b2(b2_w), .b3(b3_w));

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
    mul0_q <= mul0_a * b0_w;
    mul1_q <= mul1_a * b1_w;
    mul2_q <= mul2_a * b2_w;
    mul3_q <= mul3_a * b3_w;
    mul4_q <= mul4_a * C4_32;
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

// 1/sqrt(2) scaling, pipelined (see x-pass): products at beat 5,
// truncating /256 at beat 6 via shift + carry correction.
wire signed [31:0] pre5_w = (o_s7 - o_s6) - (o_t0 - o_s5);
wire signed [31:0] pre6_w = (o_t0 - o_s5) + (o_s7 - o_s6);

reg signed [31:0] p5_q;
reg signed [31:0] p6_q;

wire [23:0] adj5_w = p5_q[31:8] + {23'd0, (p5_q[31] & (p5_q[7:0] != 8'd0))};
wire [23:0] adj6_w = p6_q[31:8] + {23'd0, (p6_q[31] & (p6_q[7:0] != 8'd0))};
wire signed [31:0] r5_w = {{8{adj5_w[23]}}, adj5_w};
wire signed [31:0] r6_w = {{8{adj6_w[23]}}, adj6_w};

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
    p5_q   <= 32'b0;
    p6_q   <= 32'b0;
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
        p5_q <= pre5_w * 32'sd181; // (t6-t5)*181
        p6_q <= pre6_w * 32'sd181; // (t5+t6)*181
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_s5 <= r5_w; // 1/sqrt(2), signed trunc toward zero
        o_s6 <= r6_w;
    end
    default: ;
    endcase
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

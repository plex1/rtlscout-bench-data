//=====================================================================
// jpeg_idct_2d — full 8x8 2-D IDCT (core_jpeg algorithm, exact fixed
// point) — width/resource optimized version.
//
// Exactness arguments:
//  * x-pass outputs are (32-bit wrapped t) >>> 11, hence always in
//    [-2^20, 2^20)  -> stored/sign-extended with 21 bits, bit-exact.
//  * x-pass multiplier operands fit 17 bits (coeffs 16b, i0 = 17b),
//    constants fit 13 bits signed -> products fit 30 bits, no wrap,
//    sign-extension into the 32-bit accumulator context is exact.
//  * y-pass inputs fit 21 bits (i0 fits 22); products may exceed 32
//    bits, so the 35-bit true product is truncated to 32 bits (wrap),
//    matching the original 32x32 wrap multiply exactly.
//  * x*181 is computed as shift-adds keeping only the low 32 bits
//    (identical to the original 32-bit wrap multiply); /256 is signed
//    division truncating toward zero = arithmetic floor shift plus a
//    conditional increment for negative non-exact values.
//=====================================================================

module jpeg_idct_2d (
     input  wire        clk,
     input  wire        rst,
     input  wire        img_start_i,
     input  wire        img_end_i,
     input  wire        inport_valid_i,
     input  wire [15:0] inport_data0_i,
     input  wire [15:0] inport_data1_i,
     input  wire [15:0] inport_data2_i,
     input  wire [15:0] inport_data3_i,
     input  wire [ 2:0] inport_idx_i,
     output wire        outport_valid_o,
     output wire [31:0] outport_data_o,
     output wire [ 5:0] outport_idx_o
);

wire          x_valid_w;
wire [20:0]   x_data_w;
wire [ 5:0]   x_idx_w;
wire          x_accept_w;

wire          tr_valid_w;
wire [20:0]   tr_data0_w;
wire [20:0]   tr_data1_w;
wire [20:0]   tr_data2_w;
wire [20:0]   tr_data3_w;
wire [ 2:0]   tr_idx_w;

jpeg_idct_x u_idct_x (
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

jpeg_idct_transpose u_transpose (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o(x_accept_w)
    ,.outport_valid_o(tr_valid_w)
    ,.outport_data0_o(tr_data0_w)
    ,.outport_data1_o(tr_data1_w)
    ,.outport_data2_o(tr_data2_w)
    ,.outport_data3_o(tr_data3_w)
    ,.outport_idx_o(tr_idx_w)
    ,.outport_ready_i(1'b1)
);

jpeg_idct_y u_idct_y (
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(tr_valid_w)
    ,.inport_data0_i(tr_data0_w)
    ,.inport_data1_i(tr_data1_w)
    ,.inport_data2_i(tr_data2_w)
    ,.inport_data3_i(tr_data3_w)
    ,.inport_idx_i(tr_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

//=====================================================================
// 1-D IDCT row pass (OUT_SHIFT = 11)
//=====================================================================
module jpeg_idct_x (
     input  wire        clk,
     input  wire        rst,
     input  wire        img_start_i,
     input  wire        img_end_i,
     input  wire        inport_valid_i,
     input  wire [15:0] inport_data0_i,
     input  wire [15:0] inport_data1_i,
     input  wire [15:0] inport_data2_i,
     input  wire [15:0] inport_data3_i,
     input  wire [ 2:0] inport_idx_i,
     output wire        outport_valid_o,
     output wire [20:0] outport_data_o,
     output wire [ 5:0] outport_idx_o
);

localparam signed [12:0] C1_16 = 13'sd4017;
localparam signed [12:0] C2_16 = 13'sd3784;
localparam signed [12:0] C3_16 = 13'sd3406;
localparam signed [12:0] C4_16 = 13'sd2896;
localparam signed [12:0] C5_16 = 13'sd2276;
localparam signed [12:0] C6_16 = 13'sd1567;
localparam signed [12:0] C7_16 = 13'sd799;

wire signed [16:0] block_in_0_1 = {inport_data0_i[15], inport_data0_i};
wire signed [16:0] block_in_2_3 = {inport_data1_i[15], inport_data1_i};
wire signed [16:0] block_in_4_5 = {inport_data2_i[15], inport_data2_i};
wire signed [16:0] block_in_6_7 = {inport_data3_i[15], inport_data3_i};

reg signed [16:0] i0;
reg signed [16:0] mul0_a;
reg signed [16:0] mul1_a;
reg signed [16:0] mul2_a;
reg signed [16:0] mul3_a;
reg signed [16:0] mul4_a;
reg signed [12:0] mul0_b;
reg signed [12:0] mul1_b;
reg signed [12:0] mul2_b;
reg signed [12:0] mul3_b;
reg signed [12:0] mul4_b;

always @(posedge clk)
if (rst)
begin
    i0     <= 17'sd0;
    mul0_a <= 17'sd0; mul0_b <= 13'sd0;
    mul1_a <= 17'sd0; mul1_b <= 13'sd0;
    mul2_a <= 17'sd0; mul2_b <= 13'sd0;
    mul3_a <= 17'sd0; mul3_b <= 13'sd0;
    mul4_a <= 17'sd0; mul4_b <= 13'sd0;
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
        i0     <= block_in_0_1 - block_in_4_5;
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
    default:
        ;
    endcase
end

// products (17x13 -> 30 bits, never overflows)
reg signed [29:0] mul0;
reg signed [29:0] mul1;
reg signed [29:0] mul2;
reg signed [29:0] mul3;
reg signed [29:0] mul4;

always @(posedge clk)
if (rst)
begin
    mul0 <= 30'sd0;
    mul1 <= 30'sd0;
    mul2 <= 30'sd0;
    mul3 <= 30'sd0;
    mul4 <= 30'sd0;
end
else
begin
    mul0 <= mul0_a * mul0_b;
    mul1 <= mul1_a * mul1_b;
    mul2 <= mul2_a * mul2_b;
    mul3 <= mul3_a * mul3_b;
    mul4 <= mul4_a * mul4_b;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @(posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q;
    out_stg2_idx_q   <= out_stg1_idx_q;
end

// x*181 keeping low 32 bits (identical to 32-bit wrap multiply)
function [31:0] mul181_wrap;
    input signed [31:0] x;
    begin
        mul181_wrap = (x <<< 7) + (x <<< 5) + (x <<< 4) + (x <<< 2) + x;
    end
endfunction

// signed division by 256 truncating toward zero
function signed [31:0] div256_trunc;
    input signed [31:0] x;
    reg signed [31:0] f;
    begin
        f = x >>> 8;
        if (x[31] && (x[7:0] != 8'h00))
            f = f + 32'sd1;
        div256_trunc = f;
    end
endfunction

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

always @(posedge clk)
if (rst)
begin
    o_s5   <= 32'sd0;
    o_s6   <= 32'sd0;
    o_s7   <= 32'sd0;
    o_t0   <= 32'sd0;
    o_t1   <= 32'sd0;
    o_t2   <= 32'sd0;
    o_t3   <= 32'sd0;
    o_t4   <= 32'sd0;
    o_t5   <= 32'sd0;
    o_t6   <= 32'sd0;
    o_t7   <= 32'sd0;
    o_t6_5 <= 32'sd0;
    o_t5_6 <= 32'sd0;
end
else
begin
    case (out_stg1_idx_q)
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
        o_s5 <= div256_trunc(mul181_wrap(o_t6_5)); // 1/sqrt(2)
        o_s6 <= div256_trunc(mul181_wrap(o_t5_6)); // 1/sqrt(2)
    end
    endcase
end

// results of (32-bit >>> 11) always fit in 21 signed bits
reg signed [20:0] block_out [0:7];
reg signed [20:0] block_out_tmp;

always @(posedge clk)
if (rst)
begin
    block_out[0] <= 21'sd0;
    block_out[1] <= 21'sd0;
    block_out[2] <= 21'sd0;
    block_out[3] <= 21'sd0;
    block_out[4] <= 21'sd0;
    block_out[5] <= 21'sd0;
    block_out[6] <= 21'sd0;
    block_out[7] <= 21'sd0;
    block_out_tmp <= 21'sd0;
end
else if (out_stg2_valid_q)
begin
    if (out_stg2_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 11);
        block_out_tmp <= ((o_t0 - o_t7) >>> 11); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> 11);
        block_out[4] <= ((o_t3 - o_t4) >>> 11);
    end

    if (out_stg2_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg2_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> 11);
        block_out[5] <= ((o_t2 - o_s5) >>> 11);
        block_out[1] <= ((o_t1 + o_s6) >>> 11);
        block_out[6] <= ((o_t1 - o_s6) >>> 11);
    end
end

reg [7:0] valid_q;

always @(posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg2_valid_q};

reg [5:0] ptr_q;

always @(posedge clk)
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

//=====================================================================
// Double-buffered corner turn (21-bit words)
//=====================================================================
module jpeg_idct_transpose (
     input  wire        clk,
     input  wire        rst,
     input  wire        img_start_i,
     input  wire        img_end_i,
     input  wire        inport_valid_i,
     input  wire [20:0] inport_data_i,
     input  wire [ 5:0] inport_idx_i,
     input  wire        outport_ready_i,
     output wire        inport_accept_o,
     output wire        outport_valid_o,
     output wire [20:0] outport_data0_o,
     output wire [20:0] outport_data1_o,
     output wire [20:0] outport_data2_o,
     output wire [20:0] outport_data3_o,
     output wire [ 2:0] outport_idx_o
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

reg [20:0] ram0 [0:31];
reg [20:0] ram1 [0:31];
reg [20:0] ram2 [0:31];
reg [20:0] ram3 [0:31];

reg [20:0] ram0_q;
reg [20:0] ram1_q;
reg [20:0] ram2_q;
reg [20:0] ram3_q;

wire [4:0] rd_ptr_w = {block_rd_q, rd_addr_q};

always @(posedge clk)
begin
    if (wr0_w)
        ram0[wr_ptr_w] <= inport_data_i;
    ram0_q <= ram0[rd_ptr_w];
end

always @(posedge clk)
begin
    if (wr1_w)
        ram1[wr_ptr_w] <= inport_data_i;
    ram1_q <= ram1[rd_ptr_w];
end

always @(posedge clk)
begin
    if (wr2_w)
        ram2[wr_ptr_w] <= inport_data_i;
    ram2_q <= ram2[rd_ptr_w];
end

always @(posedge clk)
begin
    if (wr3_w)
        ram3[wr_ptr_w] <= inport_data_i;
    ram3_q <= ram3[rd_ptr_w];
end

assign outport_data0_o = ram0_q;
assign outport_data1_o = ram1_q;
assign outport_data2_o = ram2_q;
assign outport_data3_o = ram3_q;

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

always @(posedge clk )
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

always @(posedge clk )
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @(posedge clk )
if (rst)
    rd_idx_q <= 6'b0;
else if (img_start_i)
    rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE)
    rd_idx_q <= rd_idx_q + 6'd1;

always @(posedge clk )
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

endmodule

//=====================================================================
// 1-D IDCT column pass (OUT_SHIFT = 15)
//=====================================================================
module jpeg_idct_y (
     input  wire        clk,
     input  wire        rst,
     input  wire        img_start_i,
     input  wire        img_end_i,
     input  wire        inport_valid_i,
     input  wire [20:0] inport_data0_i,
     input  wire [20:0] inport_data1_i,
     input  wire [20:0] inport_data2_i,
     input  wire [20:0] inport_data3_i,
     input  wire [ 2:0] inport_idx_i,
     output wire        outport_valid_o,
     output wire [31:0] outport_data_o,
     output wire [ 5:0] outport_idx_o
);

localparam signed [12:0] C1_16 = 13'sd4017;
localparam signed [12:0] C2_16 = 13'sd3784;
localparam signed [12:0] C3_16 = 13'sd3406;
localparam signed [12:0] C4_16 = 13'sd2896;
localparam signed [12:0] C5_16 = 13'sd2276;
localparam signed [12:0] C6_16 = 13'sd1567;
localparam signed [12:0] C7_16 = 13'sd799;

// x-pass outputs always fit in 21 signed bits
wire signed [20:0] block_in_0_1 = inport_data0_i;
wire signed [20:0] block_in_2_3 = inport_data1_i;
wire signed [20:0] block_in_4_5 = inport_data2_i;
wire signed [20:0] block_in_6_7 = inport_data3_i;

reg signed [21:0] i0;
reg signed [21:0] mul0_a;
reg signed [21:0] mul1_a;
reg signed [21:0] mul2_a;
reg signed [21:0] mul3_a;
reg signed [21:0] mul4_a;
reg signed [12:0] mul0_b;
reg signed [12:0] mul1_b;
reg signed [12:0] mul2_b;
reg signed [12:0] mul3_b;
reg signed [12:0] mul4_b;

always @(posedge clk)
if (rst)
begin
    i0     <= 22'sd0;
    mul0_a <= 22'sd0; mul0_b <= 13'sd0;
    mul1_a <= 22'sd0; mul1_b <= 13'sd0;
    mul2_a <= 22'sd0; mul2_b <= 13'sd0;
    mul3_a <= 22'sd0; mul3_b <= 13'sd0;
    mul4_a <= 22'sd0; mul4_b <= 13'sd0;
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
        i0     <= block_in_0_1 - block_in_4_5;
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
    default:
        ;
    endcase
end

// full 22x13 products (35 bits); low 32 bits = 32-bit wrap multiply
reg signed [34:0] mul0_p;
reg signed [34:0] mul1_p;
reg signed [34:0] mul2_p;
reg signed [34:0] mul3_p;
reg signed [34:0] mul4_p;

wire signed [31:0] mul0 = mul0_p[31:0];
wire signed [31:0] mul1 = mul1_p[31:0];
wire signed [31:0] mul2 = mul2_p[31:0];
wire signed [31:0] mul3 = mul3_p[31:0];
wire signed [31:0] mul4 = mul4_p[31:0];

always @(posedge clk)
if (rst)
begin
    mul0_p <= 35'sd0;
    mul1_p <= 35'sd0;
    mul2_p <= 35'sd0;
    mul3_p <= 35'sd0;
    mul4_p <= 35'sd0;
end
else
begin
    mul0_p <= mul0_a * mul0_b;
    mul1_p <= mul1_a * mul1_b;
    mul2_p <= mul2_a * mul2_b;
    mul3_p <= mul3_a * mul3_b;
    mul4_p <= mul4_a * mul4_b;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @(posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q;
    out_stg2_idx_q   <= out_stg1_idx_q;
end

function [31:0] mul181_wrap;
    input signed [31:0] x;
    begin
        mul181_wrap = (x <<< 7) + (x <<< 5) + (x <<< 4) + (x <<< 2) + x;
    end
endfunction

function signed [31:0] div256_trunc;
    input signed [31:0] x;
    reg signed [31:0] f;
    begin
        f = x >>> 8;
        if (x[31] && (x[7:0] != 8'h00))
            f = f + 32'sd1;
        div256_trunc = f;
    end
endfunction

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

always @(posedge clk)
if (rst)
begin
    o_s5   <= 32'sd0;
    o_s6   <= 32'sd0;
    o_s7   <= 32'sd0;
    o_t0   <= 32'sd0;
    o_t1   <= 32'sd0;
    o_t2   <= 32'sd0;
    o_t3   <= 32'sd0;
    o_t4   <= 32'sd0;
    o_t5   <= 32'sd0;
    o_t6   <= 32'sd0;
    o_t7   <= 32'sd0;
    o_t6_5 <= 32'sd0;
    o_t5_6 <= 32'sd0;
end
else
begin
    case (out_stg1_idx_q)
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
        o_s5 <= div256_trunc(mul181_wrap(o_t6_5)); // 1/sqrt(2)
        o_s6 <= div256_trunc(mul181_wrap(o_t5_6)); // 1/sqrt(2)
    end
    endcase
end

// results of (32-bit >>> 15) always fit in 17 signed bits
reg signed [16:0] block_out [0:7];
reg signed [16:0] block_out_tmp;

always @(posedge clk)
if (rst)
begin
    block_out[0] <= 17'sd0;
    block_out[1] <= 17'sd0;
    block_out[2] <= 17'sd0;
    block_out[3] <= 17'sd0;
    block_out[4] <= 17'sd0;
    block_out[5] <= 17'sd0;
    block_out[6] <= 17'sd0;
    block_out[7] <= 17'sd0;
    block_out_tmp <= 17'sd0;
end
else if (out_stg2_valid_q)
begin
    if (out_stg2_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 15);
        block_out_tmp <= ((o_t0 - o_t7) >>> 15); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> 15);
        block_out[4] <= ((o_t3 - o_t4) >>> 15);
    end

    if (out_stg2_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg2_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> 15);
        block_out[5] <= ((o_t2 - o_s5) >>> 15);
        block_out[1] <= ((o_t1 + o_s6) >>> 15);
        block_out[6] <= ((o_t1 - o_s6) >>> 15);
    end
end

reg [7:0] valid_q;

always @(posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg2_valid_q};

reg [5:0] ptr_q;

always @(posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = {{15{block_out[ptr_q[2:0]][16]}}, block_out[ptr_q[2:0]]};
assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

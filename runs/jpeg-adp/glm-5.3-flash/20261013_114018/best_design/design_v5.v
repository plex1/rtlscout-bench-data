//=====================================================================
// jpeg_idct_2d : full 8x8 2-D IDCT (core_jpeg behavior, bit-exact)
//
// v5: 1-D pass re-scheduled to use only FOUR multipliers (the DC term
// i0*C4 is computed as c0*C4 + c4*C4 on the beat-2 even set and folded
// into the t0/t3 butterfly), constant operands narrowed to 13 bits, and
// explicit-width product wires (exact low-32-bit products).
//=====================================================================

//---------------------------------------------------------------------
// Two-stage pipelined signed multiplier: AW x 13b -> low 32 bits.
// Stage A: partial products + carry-save tree -> 2 vectors (registered)
// Stage B: final 32-bit carry-propagate add (registered)
// Bit-exact with a single-cycle multiply truncated to 32 bits.
//---------------------------------------------------------------------
module jpeg_idct_mul2 #(
    parameter AW = 23
)(
     input            clk
    ,input            rst
    ,input  [AW-1:0]  a_i
    ,input  [12:0]    b_i
    ,output reg [31:0] p_o
);

wire [31:0] a_sext = {{(32-AW){a_i[AW-1]}}, a_i};
wire [31:0] a_not  = ~a_sext;

wire [31:0] pp0  = b_i[0]  ? a_sext        : 32'b0;
wire [31:0] pp1  = b_i[1]  ? (a_sext << 1)  : 32'b0;
wire [31:0] pp2  = b_i[2]  ? (a_sext << 2)  : 32'b0;
wire [31:0] pp3  = b_i[3]  ? (a_sext << 3)  : 32'b0;
wire [31:0] pp4  = b_i[4]  ? (a_sext << 4)  : 32'b0;
wire [31:0] pp5  = b_i[5]  ? (a_sext << 5)  : 32'b0;
wire [31:0] pp6  = b_i[6]  ? (a_sext << 6)  : 32'b0;
wire [31:0] pp7  = b_i[7]  ? (a_sext << 7)  : 32'b0;
wire [31:0] pp8  = b_i[8]  ? (a_sext << 8)  : 32'b0;
wire [31:0] pp9  = b_i[9]  ? (a_sext << 9)  : 32'b0;
wire [31:0] pp10 = b_i[10] ? (a_sext << 10) : 32'b0;
wire [31:0] pp11 = b_i[11] ? (a_sext << 11) : 32'b0;
wire [31:0] pp12 = b_i[12] ? (a_not  << 12) : 32'b0;
wire [31:0] pp13 = b_i[12] ? 32'h1000       : 32'b0;

// level 1
wire [31:0] maj1_0 = (pp0 & pp1) | (pp0 & pp2) | (pp1 & pp2);
wire [31:0] maj1_1 = (pp3 & pp4) | (pp3 & pp5) | (pp4 & pp5);
wire [31:0] maj1_2 = (pp6 & pp7) | (pp6 & pp8) | (pp7 & pp8);
wire [31:0] maj1_3 = (pp9 & pp10) | (pp9 & pp11) | (pp10 & pp11);
wire [31:0] s1_0 = pp0 ^ pp1 ^ pp2;
wire [31:0] c1_0 = {maj1_0[30:0], 1'b0};
wire [31:0] s1_1 = pp3 ^ pp4 ^ pp5;
wire [31:0] c1_1 = {maj1_1[30:0], 1'b0};
wire [31:0] s1_2 = pp6 ^ pp7 ^ pp8;
wire [31:0] c1_2 = {maj1_2[30:0], 1'b0};
wire [31:0] s1_3 = pp9 ^ pp10 ^ pp11;
wire [31:0] c1_3 = {maj1_3[30:0], 1'b0};

// level 2 (leftover pp12)
wire [31:0] maj2_0 = (s1_0 & s1_1) | (s1_0 & s1_2) | (s1_1 & s1_2);
wire [31:0] maj2_1 = (s1_3 & c1_0) | (s1_3 & c1_1) | (c1_0 & c1_1);
wire [31:0] maj2_2 = (c1_2 & c1_3) | (c1_2 & pp12) | (c1_3 & pp12);
wire [31:0] s2_0 = s1_0 ^ s1_1 ^ s1_2;
wire [31:0] c2_0 = {maj2_0[30:0], 1'b0};
wire [31:0] s2_1 = s1_3 ^ c1_0 ^ c1_1;
wire [31:0] c2_1 = {maj2_1[30:0], 1'b0};
wire [31:0] s2_2 = c1_2 ^ c1_3 ^ pp12;
wire [31:0] c2_2 = {maj2_2[30:0], 1'b0};

// level 3 (leftover pp13)
wire [31:0] maj3_0 = (s2_0 & s2_1) | (s2_0 & s2_2) | (s2_1 & s2_2);
wire [31:0] maj3_1 = (c2_0 & c2_1) | (c2_0 & c2_2) | (c2_1 & c2_2);
wire [31:0] s3_0 = s2_0 ^ s2_1 ^ s2_2;
wire [31:0] c3_0 = {maj3_0[30:0], 1'b0};
wire [31:0] s3_1 = c2_0 ^ c2_1 ^ c2_2;
wire [31:0] c3_1 = {maj3_1[30:0], 1'b0};

// level 4 (leftovers c3_0, c3_1)
wire [31:0] maj4_0 = (s3_0 & s3_1) | (s3_0 & pp13) | (s3_1 & pp13);
wire [31:0] s4_0 = s3_0 ^ s3_1 ^ pp13;
wire [31:0] c4_0 = {maj4_0[30:0], 1'b0};

// level 5 (leftover c4_0)
wire [31:0] maj5_0 = (s4_0 & c3_0) | (s4_0 & c3_1) | (c3_0 & c3_1);
wire [31:0] s5_0 = s4_0 ^ c3_0 ^ c3_1;
wire [31:0] c5_0 = {maj5_0[30:0], 1'b0};

// level 6
wire [31:0] maj6_0 = (s5_0 & c5_0) | (s5_0 & c4_0) | (c5_0 & c4_0);
wire [31:0] s6_0 = s5_0 ^ c5_0 ^ c4_0;
wire [31:0] c6_0 = {maj6_0[30:0], 1'b0};

reg [31:0] sa_q;
reg [31:0] sc_q;

always @ (posedge clk)
if (rst)
begin
    sa_q <= 32'b0;
    sc_q <= 32'b0;
    p_o  <= 32'b0;
end
else
begin
    sa_q <= s6_0;
    sc_q <= c6_0;
    p_o  <= sa_q + sc_q;
end

endmodule

//---------------------------------------------------------------------
// 1-D IDCT pass (factorized AAN algorithm, 4 time-shared multipliers)
//---------------------------------------------------------------------
module jpeg_idct_1d #(
     parameter IN_W      = 16
    ,parameter OUT_W     = 32
    ,parameter OUT_SHIFT = 11
    ,parameter CONV_IDX  = 0
)(
     input            clk
    ,input            rst
    ,input            img_start_i
    ,input            img_end_i
    ,input            inport_valid_i
    ,input  [IN_W-1:0] inport_data0_i
    ,input  [IN_W-1:0] inport_data1_i
    ,input  [IN_W-1:0] inport_data2_i
    ,input  [IN_W-1:0] inport_data3_i
    ,input  [ 2:0]    inport_idx_i
    ,output           outport_valid_o
    ,output [OUT_W-1:0] outport_data_o
    ,output [ 5:0]    outport_idx_o
);

// cosine constants scaled by 4096 (all < 4096 -> 13-bit signed)
localparam signed [12:0] C1_16 = 4017; // cos( pi/16) x4096
localparam signed [12:0] C2_16 = 3784; // cos(2pi/16) x4096
localparam signed [12:0] C3_16 = 3406; // cos(3pi/16) x4096
localparam signed [12:0] C4_16 = 2896; // cos(4pi/16) x4096
localparam signed [12:0] C5_16 = 2276; // cos(5pi/16) x4096
localparam signed [12:0] C6_16 = 1567; // cos(6pi/16) x4096
localparam signed [12:0] C7_16 = 799;  // cos(7pi/16) x4096

localparam PW = IN_W + 14;             // max product width ((IN_W+1)+13)

wire signed [IN_W-1:0] block_in_0_1 = inport_data0_i;
wire signed [IN_W-1:0] block_in_2_3 = inport_data1_i;
wire signed [IN_W-1:0] block_in_4_5 = inport_data2_i;
wire signed [IN_W-1:0] block_in_6_7 = inport_data3_i;

reg signed [IN_W:0]   i1;             // c0 - c4
reg signed [22:0]     mul0_a;
reg signed [12:0]     mul0_b;
reg signed [22:0]     mul1_a;
reg signed [12:0]     mul1_b;
reg signed [22:0]     mul2_a;
reg signed [12:0]     mul2_b;
reg signed [22:0]     mul3_a;
reg signed [12:0]     mul3_b;

always @ (posedge clk)
if (rst)
begin
    i1     <= {(IN_W+1){1'b0}};
    mul0_a <= 23'b0;
    mul0_b <= 13'b0;
    mul1_a <= 23'b0;
    mul1_b <= 13'b0;
    mul2_a <= 23'b0;
    mul2_b <= 13'b0;
    mul3_a <= 23'b0;
    mul3_b <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
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
    end
    3'd2:
    begin
        i1     <= block_in_0_1 - block_in_4_5;
        mul0_a <= block_in_0_1;
        mul0_b <= C4_16;
        mul1_a <= block_in_4_5;
        mul1_b <= C4_16;
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
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C2_16;
        mul3_a <= i1;
        mul3_b <= C4_16;
    end
    default:
        ;
    endcase
end

// products via 2-stage pipelined multipliers (23b x 13b -> low 32 bits)
wire [31:0] prod0_w;
wire [31:0] prod1_w;
wire [31:0] prod2_w;
wire [31:0] prod3_w;

jpeg_idct_mul2 u_mul0 (.clk(clk), .rst(rst), .a_i(mul0_a), .b_i(mul0_b), .p_o(prod0_w));
jpeg_idct_mul2 u_mul1 (.clk(clk), .rst(rst), .a_i(mul1_a), .b_i(mul1_b), .p_o(prod1_w));
jpeg_idct_mul2 u_mul2 (.clk(clk), .rst(rst), .a_i(mul2_a), .b_i(mul2_b), .p_o(prod2_w));
jpeg_idct_mul2 u_mul3 (.clk(clk), .rst(rst), .a_i(mul3_a), .b_i(mul3_b), .p_o(prod3_w));

// the mul2 modules already carry two register stages -> same latency as
// the original (operand latch -> use) pipeline
wire signed [31:0] mul0 = prod0_w;
wire signed [31:0] mul1 = prod1_w;
wire signed [31:0] mul2 = prod2_w;
wire signed [31:0] mul3 = prod3_w;

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
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

// signed division by 256 truncating toward zero (== Verilog '/' semantics)
function signed [31:0] div256;
    input signed [31:0] x;
begin
    div256 = (x[31]) ? ((x + 256'sd255) >>> 8) : (x >>> 8);
end
endfunction

wire signed [31:0] s0_w = mul0 + mul1;   // (c0*C4 + c4*C4) == i0*C4

// *181/256 split across two beats: multiply on beat 6, divide on beat 7
reg signed [31:0] m5_q;
reg signed [31:0] m6_q;

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
    m5_q   <= 32'b0;
    m6_q   <= 32'b0;
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
    end
    3'd2:
    begin
        o_t0 <= s0_w + o_t3; // t0 = s0 + s3
        o_t3 <= s0_w - o_t3; // t3 = s0 - s3
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
        o_t4 <= mul3; // s1
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        m5_q <= o_t5 * 181; // 1/sqrt(2) scaling, stage 1
        m6_q <= o_t6 * 181;
    end
    default:
    begin
        o_s5 <= div256(m6_q - m5_q); // 1/sqrt(2), stage 2
        o_s6 <= div256(m5_q + m6_q);
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

reg signed [OUT_W-1:0] block_out [0:7];
reg signed [OUT_W-1:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= {(OUT_W){1'b0}};
    block_out[1] <= {(OUT_W){1'b0}};
    block_out[2] <= {(OUT_W){1'b0}};
    block_out[3] <= {(OUT_W){1'b0}};
    block_out[4] <= {(OUT_W){1'b0}};
    block_out[5] <= {(OUT_W){1'b0}};
    block_out[6] <= {(OUT_W){1'b0}};
    block_out[7] <= {(OUT_W){1'b0}};
    block_out_tmp <= {(OUT_W){1'b0}};
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

function [5:0] ptr_conv;
    input [5:0] idx;
    reg [5:0] out_idx;
begin
    case (idx)
    6'd0:  out_idx = 6'd0;
    6'd1:  out_idx = 6'd8;
    6'd2:  out_idx = 6'd16;
    6'd3:  out_idx = 6'd24;
    6'd4:  out_idx = 6'd32;
    6'd5:  out_idx = 6'd40;
    6'd6:  out_idx = 6'd48;
    6'd7:  out_idx = 6'd56;
    6'd8:  out_idx = 6'd1;
    6'd9:  out_idx = 6'd9;
    6'd10:  out_idx = 6'd17;
    6'd11:  out_idx = 6'd25;
    6'd12:  out_idx = 6'd33;
    6'd13:  out_idx = 6'd41;
    6'd14:  out_idx = 6'd49;
    6'd15:  out_idx = 6'd57;
    6'd16:  out_idx = 6'd2;
    6'd17:  out_idx = 6'd10;
    6'd18:  out_idx = 6'd18;
    6'd19:  out_idx = 6'd26;
    6'd20:  out_idx = 6'd34;
    6'd21:  out_idx = 6'd42;
    6'd22:  out_idx = 6'd50;
    6'd23:  out_idx = 6'd58;
    6'd24:  out_idx = 6'd3;
    6'd25:  out_idx = 6'd11;
    6'd26:  out_idx = 6'd19;
    6'd27:  out_idx = 6'd27;
    6'd28:  out_idx = 6'd35;
    6'd29:  out_idx = 6'd43;
    6'd30:  out_idx = 6'd51;
    6'd31:  out_idx = 6'd59;
    6'd32:  out_idx = 6'd4;
    6'd33:  out_idx = 6'd12;
    6'd34:  out_idx = 6'd20;
    6'd35:  out_idx = 6'd28;
    6'd36:  out_idx = 6'd36;
    6'd37:  out_idx = 6'd44;
    6'd38:  out_idx = 6'd52;
    6'd39:  out_idx = 6'd60;
    6'd40:  out_idx = 6'd5;
    6'd41:  out_idx = 6'd13;
    6'd42:  out_idx = 6'd21;
    6'd43:  out_idx = 6'd29;
    6'd44:  out_idx = 6'd37;
    6'd45:  out_idx = 6'd45;
    6'd46:  out_idx = 6'd53;
    6'd47:  out_idx = 6'd61;
    6'd48:  out_idx = 6'd6;
    6'd49:  out_idx = 6'd14;
    6'd50:  out_idx = 6'd22;
    6'd51:  out_idx = 6'd30;
    6'd52:  out_idx = 6'd38;
    6'd53:  out_idx = 6'd46;
    6'd54:  out_idx = 6'd54;
    6'd55:  out_idx = 6'd62;
    6'd56:  out_idx = 6'd7;
    6'd57:  out_idx = 6'd15;
    6'd58:  out_idx = 6'd23;
    6'd59:  out_idx = 6'd31;
    6'd60:  out_idx = 6'd39;
    6'd61:  out_idx = 6'd47;
    6'd62:  out_idx = 6'd55;
    default:  out_idx = 6'd63;
    endcase

    ptr_conv = out_idx;
end
endfunction

assign outport_idx_o   = (CONV_IDX != 0) ? ptr_conv(ptr_q) : ptr_q;

endmodule

//---------------------------------------------------------------------
// Corner turn RAM: 16 x XW, one write port + one-hot read select
//---------------------------------------------------------------------
module jpeg_idct_transpose_ram #(
    parameter XW = 20
)(
     input            clk0_i
    ,input  [ 3:0]    addr0_i
    ,input  [XW-1:0]  data0_i
    ,input            wr0_i
    ,input            clk1_i
    ,input  [15:0]    addr1_hot_i
    ,output [XW-1:0]  data1_o
);

reg [XW-1:0] ram [15:0];
reg [XW-1:0] ram_read1_q;

integer j;
reg [XW-1:0] acc;

always @ (posedge clk0_i)
begin
    if (wr0_i)
        ram[addr0_i] <= data0_i;
end

always @ (posedge clk1_i)
begin
    acc = {XW{1'b0}};
    for (j = 0; j < 16; j = j + 1)
        acc = acc | ({XW{addr1_hot_i[j]}} & ram[j]);
    ram_read1_q <= acc;
end

assign data1_o = ram_read1_q;

endmodule

//---------------------------------------------------------------------
// Corner turn (double-buffered transpose)
//---------------------------------------------------------------------
module jpeg_idct_transpose #(
    parameter XW = 20
)(
    // Inputs
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [XW-1:0] inport_data_i
    ,input  [ 5:0]   inport_idx_i
    ,input           outport_ready_i

    // Outputs
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [XW-1:0] outport_data0_o
    ,output [XW-1:0] outport_data1_o
    ,output [XW-1:0] outport_data2_o
    ,output [XW-1:0] outport_data3_o
    ,output [ 2:0]   outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [15:0]  rd_hot_q;      // one-hot slot select for the read port

wire [3:0]  wr_ptr_w = inport_idx_i[3:0];

wire wr_bank0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr_bank1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr_bank2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr_bank3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

wire [XW-1:0] outport_data0_a_w;
wire [XW-1:0] outport_data0_b_w;
wire [XW-1:0] outport_data1_a_w;
wire [XW-1:0] outport_data1_b_w;
wire [XW-1:0] outport_data2_a_w;
wire [XW-1:0] outport_data2_b_w;
wire [XW-1:0] outport_data3_a_w;
wire [XW-1:0] outport_data3_b_w;

// lane 0: two banks
jpeg_idct_transpose_ram #(.XW(XW))
u_ram0a
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank0_w & ~block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data0_a_w)
);

jpeg_idct_transpose_ram #(.XW(XW))
u_ram0b
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank0_w & block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data0_b_w)
);

// lane 1
jpeg_idct_transpose_ram #(.XW(XW))
u_ram1a
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank1_w & ~block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data1_a_w)
);

jpeg_idct_transpose_ram #(.XW(XW))
u_ram1b
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank1_w & block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data1_b_w)
);

// lane 2
jpeg_idct_transpose_ram #(.XW(XW))
u_ram2a
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank2_w & ~block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data2_a_w)
);

jpeg_idct_transpose_ram #(.XW(XW))
u_ram2b
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank2_w & block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data2_b_w)
);

// lane 3
jpeg_idct_transpose_ram #(.XW(XW))
u_ram3a
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank3_w & ~block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data3_a_w)
);

jpeg_idct_transpose_ram #(.XW(XW))
u_ram3b
(
     .clk0_i(clk)
    ,.addr0_i(wr_ptr_w)
    ,.data0_i(inport_data_i)
    ,.wr0_i(wr_bank3_w & block_wr_q)
    ,.clk1_i(clk)
    ,.addr1_hot_i(rd_hot_q)
    ,.data1_o(outport_data3_b_w)
);

// bank select after the RAM output registers
assign outport_data0_o = block_rd_q ? outport_data0_b_w : outport_data0_a_w;
assign outport_data1_o = block_rd_q ? outport_data1_b_w : outport_data1_a_w;
assign outport_data2_o = block_rd_q ? outport_data2_b_w : outport_data2_a_w;
assign outport_data3_o = block_rd_q ? outport_data3_b_w : outport_data3_a_w;

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

// one-hot read address: replicates rd_addr_q sequence
//  +-8  -> swap the two 8-slot halves
//  +1   -> rotate the 16-bit one-hot left by one
wire [15:0] rd_hot_swap_w   = {rd_hot_q[7:0], rd_hot_q[15:8]};
wire [15:0] rd_hot_rot_w    = {rd_hot_q[14:0], rd_hot_q[15]};

always @ (posedge clk )
if (rst)
    rd_hot_q <= 16'h0001;
else if (img_start_i)
    rd_hot_q <= 16'h0001;
else if (state_q == STATE_IDLE)
    rd_hot_q <= 16'h0001;
else if (state_q == STATE_SETUP)
    rd_hot_q <= 16'h0100;
else if (state_q == STATE_ACTIVE)
begin
    case (rd_idx_q[2:0])
    3'd0: rd_hot_q <= rd_hot_swap_w;
    3'd1: rd_hot_q <= rd_hot_swap_w;
    3'd2: rd_hot_q <= rd_hot_q;
    3'd3: rd_hot_q <= rd_hot_swap_w;
    3'd4: rd_hot_q <= rd_hot_swap_w;
    3'd5: rd_hot_q <= rd_hot_swap_w;
    3'd6: rd_hot_q <= rd_hot_rot_w;
    3'd7: rd_hot_q <= rd_hot_swap_w;
    endcase
end

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

endmodule

//=====================================================================
// Top level
//=====================================================================
module jpeg_idct_2d #(
    parameter XW = 20
)
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
wire [XW-1:0] idct_x_data_w;
wire [ 5:0]   idct_x_idx_w;
wire          idct_x_accept_w;

jpeg_idct_1d #(.IN_W(16), .OUT_W(XW), .OUT_SHIFT(11), .CONV_IDX(0))
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
wire          transpose_ready_w = 1'b1;

jpeg_idct_transpose #(.XW(XW))
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
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

// pipeline register between the corner turn and the y pass: splits the
// bank-select mux from the y-pass input logic (uniform 1-cycle delay)
reg         y_in_valid_q;
reg [ 2:0]  y_in_idx_q;
reg [XW-1:0] y_in0_q;
reg [XW-1:0] y_in1_q;
reg [XW-1:0] y_in2_q;
reg [XW-1:0] y_in3_q;

always @ (posedge clk)
if (rst)
begin
    y_in_valid_q <= 1'b0;
    y_in_idx_q   <= 3'b0;
    y_in0_q      <= {XW{1'b0}};
    y_in1_q      <= {XW{1'b0}};
    y_in2_q      <= {XW{1'b0}};
    y_in3_q      <= {XW{1'b0}};
end
else if (img_start_i)
begin
    y_in_valid_q <= 1'b0;
    y_in_idx_q   <= 3'b0;
end
else
begin
    y_in_valid_q <= transpose_valid_w;
    y_in_idx_q   <= transpose_idx_w;
    y_in0_q      <= transpose_data0_w;
    y_in1_q      <= transpose_data1_w;
    y_in2_q      <= transpose_data2_w;
    y_in3_q      <= transpose_data3_w;
end

jpeg_idct_1d #(.IN_W(XW), .OUT_W(32), .OUT_SHIFT(15), .CONV_IDX(1))
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(y_in_valid_q)
    ,.inport_data0_i(y_in0_q)
    ,.inport_data1_i(y_in1_q)
    ,.inport_data2_i(y_in2_q)
    ,.inport_data3_i(y_in3_q)
    ,.inport_idx_i(y_in_idx_q)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

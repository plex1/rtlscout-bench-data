//=====================================================================
// Optimized full 8x8 2-D IDCT (core_jpeg bit-exact behavior).
//
//   jpeg_idct_x (row pass, OUT_SHIFT=11)
//     -> double-buffered corner turn (transpose)
//     -> jpeg_idct_y (column pass, OUT_SHIFT=15)
//
// Bit-exactness notes:
//  * Constant multipliers (C1..C7 = 4017,3784,3406,2896,2276,1567,799,
//    i.e. cos*4096) are expanded into balanced shift/add trees.
//    - Row pass: operands are 16-bit coefficients or (c0+c4) [17-bit];
//      every product |a*C| < 2^28, so a 29-bit exact tree is used.
//    - Column pass: operands are the row-pass outputs (provably
//      bounded to 20 bits signed, so 21-bit operand registers) and
//      products are computed with 32-bit wrap arithmetic, exactly
//      equal to (a*C) mod 2^32 as produced by the original 32x32
//      signed multiply.
//  * (x*181)/256: p = (x*181) mod 2^32 (as the original), then
//    truncating signed division by 256:
//        p >= 0 : p >>> 8 ;  p < 0 : (p + 255) >>> 8.
//  * Row-pass outputs are bounded by 2^19 in magnitude (all sums stay
//    below 2^30), so the corner-turn storage is 20-bit signed.
//  * Column-pass outputs are (32-bit value >>> 15), i.e. 17-bit signed.
//=====================================================================

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

wire        x_valid_w;
wire [19:0] x_data_w;
wire [ 5:0] x_idx_w;

jpeg_idct_1d
#(
     .OUT_SHIFT (11),
     .IW        (16),
     .OPW       (17),
     .MUL_WRAP  (0),
     .OW        (31),
     .BOW       (20),
     .PERM_IDX  (0)
)
u_idct_x
(
     .clk             (clk)
    ,.rst             (rst)
    ,.img_start_i     (img_start_i)
    ,.inport_valid_i  (inport_valid_i)
    ,.inport_data0_i  (inport_data0_i)
    ,.inport_data1_i  (inport_data1_i)
    ,.inport_data2_i  (inport_data2_i)
    ,.inport_data3_i  (inport_data3_i)
    ,.inport_idx_i    (inport_idx_i)
    ,.outport_valid_o (x_valid_w)
    ,.outport_data_o  (x_data_w)
    ,.outport_idx_o   (x_idx_w)
);

wire        t_valid_w;
wire [31:0] t_data0_w;
wire [31:0] t_data1_w;
wire [31:0] t_data2_w;
wire [31:0] t_data3_w;
wire [ 2:0] t_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk             (clk)
    ,.rst             (rst)
    ,.img_start_i     (img_start_i)
    ,.inport_valid_i  (x_valid_w)
    ,.inport_data_i   (x_data_w)
    ,.inport_idx_i    (x_idx_w)
    ,.outport_ready_i (1'b1)
    ,.outport_valid_o (t_valid_w)
    ,.outport_data0_o (t_data0_w)
    ,.outport_data1_o (t_data1_w)
    ,.outport_data2_o (t_data2_w)
    ,.outport_data3_o (t_data3_w)
    ,.outport_idx_o   (t_idx_w)
);

jpeg_idct_1d
#(
     .OUT_SHIFT (15),
     .IW        (32),
     .OPW       (21),
     .MUL_WRAP  (1),
     .OW        (32),
     .BOW       (17),
     .PERM_IDX  (1)
)
u_idct_y
(
     .clk             (clk)
    ,.rst             (rst)
    ,.img_start_i     (img_start_i)
    ,.inport_valid_i  (t_valid_w)
    ,.inport_data0_i  (t_data0_w)
    ,.inport_data1_i  (t_data1_w)
    ,.inport_data2_i  (t_data2_w)
    ,.inport_data3_i  (t_data3_w)
    ,.inport_idx_i    (t_idx_w)
    ,.outport_valid_o (outport_valid_o)
    ,.outport_data_o  (outport_data_o)
    ,.outport_idx_o   (outport_idx_o)
);

endmodule

//=====================================================================
// Shared 1-D factorized (AAN-style) IDCT pass.
//   IW       : input lane width (16 = row pass, 32 = column pass)
//   OPW      : multiplier operand register width (17 = row, 21 = col)
//   MUL_WRAP : 0 = exact 29-bit product tree (row pass)
//              1 = 32-bit wrap product tree (column pass)
//   OW       : combining-stage register width
//   BOW      : block_out register width
//=====================================================================
module jpeg_idct_1d
#(
     parameter OUT_SHIFT = 11,
     parameter IW        = 16,
     parameter OPW       = 17,
     parameter MUL_WRAP  = 0,
     parameter OW        = 31,
     parameter BOW       = 20,
     parameter PERM_IDX  = 0
)
(
     input            clk
    ,input            rst
    ,input            img_start_i
    ,input            inport_valid_i
    ,input  [ IW-1:0] inport_data0_i
    ,input  [ IW-1:0] inport_data1_i
    ,input  [ IW-1:0] inport_data2_i
    ,input  [ IW-1:0] inport_data3_i
    ,input  [  2:0]   inport_idx_i
    ,output           outport_valid_o
    ,output [ BOW-1:0] outport_data_o
    ,output [  5:0]   outport_idx_o
);

localparam MW = (MUL_WRAP != 0) ? 32 : 29;

// Operand narrowing: row pass inputs are 16-bit signed; column pass
// inputs are 20-bit signed values (sign-extended in 32 bits).
wire signed [OPW-1:0] a0 = (IW == 16) ? {inport_data0_i[15], inport_data0_i[15:0]}
                                      : {inport_data0_i[19], inport_data0_i[19:0]};
wire signed [OPW-1:0] a2 = (IW == 16) ? {inport_data1_i[15], inport_data1_i[15:0]}
                                      : {inport_data1_i[19], inport_data1_i[19:0]};
wire signed [OPW-1:0] a4 = (IW == 16) ? {inport_data2_i[15], inport_data2_i[15:0]}
                                      : {inport_data2_i[19], inport_data2_i[19:0]};
wire signed [OPW-1:0] a6 = (IW == 16) ? {inport_data3_i[15], inport_data3_i[15:0]}
                                      : {inport_data3_i[19], inport_data3_i[19:0]};

//---------------------------------------------------------------------
// Constant multiplier: b in {1..7} selects C1..C7.
//   C1 = 4017 = 4096 - 64 - 15
//   C2 = 3784 = 4096 - 256 - 56
//   C3 = 3406 = 4096 - 512 - 160 - 18
//   C4 = 2896 = 4096 - 1024 - 160 - 16
//   C5 = 2276 = 2048 + 128 + 100
//   C6 = 1567 = 2048 - 512 + 31
//   C7 =  799 = 1024 - 256 + 31
//---------------------------------------------------------------------
function signed [MW-1:0] cmul;
    input signed [OPW-1:0] a;
    input         [ 2:0] b;
    reg signed [28:0] n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12;
    reg signed [28:0] q1, q2, q3, acc29;
    reg signed [31:0] w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12;
    reg signed [31:0] r1, r2, r3, acc32;
    begin
        if (MUL_WRAP == 0)
        begin
            n0  = a;
            n1  = n0 <<< 1;
            n2  = n0 <<< 2;
            n3  = n0 <<< 3;
            n4  = n0 <<< 4;
            n5  = n0 <<< 5;
            n6  = n0 <<< 6;
            n7  = n0 <<< 7;
            n8  = n0 <<< 8;
            n9  = n0 <<< 9;
            n10 = n0 <<< 10;
            n11 = n0 <<< 11;
            n12 = n0 <<< 12;
            case (b)
            3'd1: begin q1 = n12 - n6;  q2 = n0 - n4;  acc29 = q1 + q2; end
            3'd2: begin q1 = n12 - n8;  q2 = n3 - n6;  acc29 = q1 + q2; end
            3'd3: begin q1 = n12 - n9;  q2 = n7 + n5;  q3 = n4 + n1;
                        acc29 = (q1 - q2) - q3; end
            3'd4: begin q1 = n12 - n10; q2 = n7 + n5;  acc29 = (q1 - q2) - n4; end
            3'd5: begin q1 = n11 + n7;  q2 = n6 + n5;  acc29 = (q1 + q2) + n2; end
            3'd6: begin q1 = n11 - n9;  q2 = n5 - n0;  acc29 = q1 + q2; end
            default: begin q1 = n10 - n8; q2 = n5 - n0; acc29 = q1 + q2; end
            endcase
            cmul = acc29;
        end
        else
        begin
            w1  = a <<< 1;
            w2  = a <<< 2;
            w3  = a <<< 3;
            w4  = a <<< 4;
            w5  = a <<< 5;
            w6  = a <<< 6;
            w7  = a <<< 7;
            w8  = a <<< 8;
            w9  = a <<< 9;
            w10 = a <<< 10;
            w11 = a <<< 11;
            w12 = a <<< 12;
            case (b)
            3'd1: begin r1 = w12 - w6;  r2 = a  - w4;  acc32 = r1 + r2; end
            3'd2: begin r1 = w12 - w8;  r2 = w3 - w6;  acc32 = r1 + r2; end
            3'd3: begin r1 = w12 - w9;  r2 = w7 + w5;  r3 = w4 + w1;
                        acc32 = (r1 - r2) - r3; end
            3'd4: begin r1 = w12 - w10; r2 = w7 + w5;  acc32 = (r1 - r2) - w4; end
            3'd5: begin r1 = w11 + w7;  r2 = w6 + w5;  acc32 = (r1 + r2) + w2; end
            3'd6: begin r1 = w11 - w9;  r2 = w5 - a;   acc32 = r1 + r2; end
            default: begin r1 = w10 - w8; r2 = w5 - a; acc32 = r1 + r2; end
            endcase
            cmul = acc32;
        end
    end
endfunction

// (x*181)/256 with signed division truncating toward zero.
// p = (x*181) mod 2^32 exactly as the original; then
//   p >= 0 : p >>> 8 ;  p < 0 : (p + 255) >>> 8.
function signed [31:0] m181_div256;
    input signed [OW-1:0] x;
    reg signed [31:0] p1, p2, p;
    begin
        p1 = (x <<< 7) + (x <<< 5);   // x*160
        p2 = (x <<< 4) + (x <<< 2);   // x*20
        p  = (p1 + p2) + x;           // x*181 (mod 2^32)
        m181_div256 = (p[31] ? (p + 32'sd255) : p) >>> 8;
    end
endfunction

//---------------------------------------------------------------------
// Operand/constant registers (5 time-shared multipliers)
//---------------------------------------------------------------------
reg signed [OPW-1:0] i0;
reg signed [OPW-1:0] mul0_a;
reg         [ 2:0]   mul0_b;
reg signed [OPW-1:0] mul1_a;
reg         [ 2:0]   mul1_b;
reg signed [OPW-1:0] mul2_a;
reg         [ 2:0]   mul2_b;
reg signed [OPW-1:0] mul3_a;
reg         [ 2:0]   mul3_b;
reg signed [OPW-1:0] mul4_a;
reg         [ 2:0]   mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= {OPW{1'b0}};
    mul0_a <= {OPW{1'b0}}; mul0_b <= 3'b0;
    mul1_a <= {OPW{1'b0}}; mul1_b <= 3'b0;
    mul2_a <= {OPW{1'b0}}; mul2_b <= 3'b0;
    mul3_a <= {OPW{1'b0}}; mul3_b <= 3'b0;
    mul4_a <= {OPW{1'b0}}; mul4_b <= 3'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= a0 + a4;
        mul0_a <= a2; mul0_b <= 3'd2;
        mul1_a <= a6; mul1_b <= 3'd6;
    end
    3'd1:
    begin
        mul0_a <= a0; mul0_b <= 3'd1;
        mul1_a <= a6; mul1_b <= 3'd7;
        mul2_a <= a4; mul2_b <= 3'd5;
        mul3_a <= a2; mul3_b <= 3'd3;
        mul4_a <= i0; mul4_b <= 3'd4;
    end
    3'd2:
    begin
        i0     <= a0 - a4;
    end
    3'd3:
    begin
        mul0_a <= a0; mul0_b <= 3'd7;
        mul1_a <= a6; mul1_b <= 3'd1;
        mul2_a <= a4; mul2_b <= 3'd3;
        mul3_a <= a2; mul3_b <= 3'd5;
    end
    3'd4:
    begin
        mul0_a <= a0; mul0_b <= 3'd7;
        mul1_a <= a6; mul1_b <= 3'd1;
        mul2_a <= a4; mul2_b <= 3'd3;
        mul3_a <= a2; mul3_b <= 3'd5;
    end
    3'd5:
    begin
        mul0_a <= a2; mul0_b <= 3'd6;
        mul1_a <= a6; mul1_b <= 3'd2;
        mul4_a <= i0; mul4_b <= 3'd4;
    end
    default:
        ;
    endcase
end

reg signed [MW-1:0] mul0_q;
reg signed [MW-1:0] mul1_q;
reg signed [MW-1:0] mul2_q;
reg signed [MW-1:0] mul3_q;
reg signed [MW-1:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= {MW{1'b0}}; mul1_q <= {MW{1'b0}}; mul2_q <= {MW{1'b0}};
    mul3_q <= {MW{1'b0}}; mul4_q <= {MW{1'b0}};
end
else
begin
    mul0_q <= cmul(mul0_a, mul0_b);
    mul1_q <= cmul(mul1_a, mul1_b);
    mul2_q <= cmul(mul2_a, mul2_b);
    mul3_q <= cmul(mul3_a, mul3_b);
    mul4_q <= cmul(mul4_a, mul4_b);
end

reg signed [MW-1:0] mul0;
reg signed [MW-1:0] mul1;
reg signed [MW-1:0] mul2;
reg signed [MW-1:0] mul3;
reg signed [MW-1:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= {MW{1'b0}}; mul1 <= {MW{1'b0}}; mul2 <= {MW{1'b0}};
    mul3 <= {MW{1'b0}}; mul4 <= {MW{1'b0}};
end
else
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

// Widen to 32 bits for the combining-stage expressions.
wire signed [31:0] mul0_w = {{(32-MW){mul0[MW-1]}}, mul0};
wire signed [31:0] mul1_w = {{(32-MW){mul1[MW-1]}}, mul1};
wire signed [31:0] mul2_w = {{(32-MW){mul2[MW-1]}}, mul2};
wire signed [31:0] mul3_w = {{(32-MW){mul3[MW-1]}}, mul3};
wire signed [31:0] mul4_w = {{(32-MW){mul4[MW-1]}}, mul4};

//---------------------------------------------------------------------
// Pipeline valid/idx tracking
//---------------------------------------------------------------------
reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
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

//---------------------------------------------------------------------
// Combining stage (beat-scheduled, same as the original)
//---------------------------------------------------------------------
reg signed [OW-1:0] o_s5;
reg signed [OW-1:0] o_s6;
reg signed [OW-1:0] o_s7;
reg signed [OW-1:0] o_t0;
reg signed [OW-1:0] o_t1;
reg signed [OW-1:0] o_t2;
reg signed [OW-1:0] o_t3;
reg signed [OW-1:0] o_t4;
reg signed [OW-1:0] o_t5;
reg signed [OW-1:0] o_t6;
reg signed [OW-1:0] o_t7;
reg signed [OW-1:0] o_t6_5;
reg signed [OW-1:0] o_t5_6;

always @ (posedge clk)
if (rst)
begin
    o_s5 <= {OW{1'b0}}; o_s6 <= {OW{1'b0}}; o_s7 <= {OW{1'b0}};
    o_t0 <= {OW{1'b0}}; o_t1 <= {OW{1'b0}}; o_t2 <= {OW{1'b0}};
    o_t3 <= {OW{1'b0}}; o_t4 <= {OW{1'b0}}; o_t5 <= {OW{1'b0}};
    o_t6 <= {OW{1'b0}}; o_t7 <= {OW{1'b0}};
    o_t6_5 <= {OW{1'b0}}; o_t5_6 <= {OW{1'b0}};
end
else
begin
    case (out_stg2_idx_q)
    3'd0:
    begin
        o_t3 <= mul0_w + mul1_w; // s3
    end
    3'd1:
    begin
        o_s7 <= mul0_w + mul1_w;
        o_s6 <= mul2_w + mul3_w;
        o_t0 <= mul4_w;          // s0
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3;     // t0
        o_t3 <= o_t0 - o_t3;     // t3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3:
    begin
        o_t4 <= (mul0_w - mul1_w) + (mul2_w - mul3_w);
    end
    3'd4:
    begin
        o_t0 <= mul0_w - mul1_w; // s4
        o_s5 <= mul2_w - mul3_w;
    end
    3'd5:
    begin
        o_t3 <= mul0_w - mul1_w; // s2
        o_t4 <= mul4_w;          // s1
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6:
    begin
        o_t1   <= o_t4 + o_t3;
        o_t2   <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5;
        o_t5_6 <= o_t5 + o_t6;
    end
    default:
    begin
        o_s5 <= m181_div256(o_t6_5); // 1/sqrt(2)
        o_s6 <= m181_div256(o_t5_6);
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

//---------------------------------------------------------------------
// Output register file + streamer
//---------------------------------------------------------------------
reg signed [BOW-1:0] block_out [0:7];
reg signed [BOW-1:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= {BOW{1'b0}}; block_out[1] <= {BOW{1'b0}};
    block_out[2] <= {BOW{1'b0}}; block_out[3] <= {BOW{1'b0}};
    block_out[4] <= {BOW{1'b0}}; block_out[5] <= {BOW{1'b0}};
    block_out[6] <= {BOW{1'b0}}; block_out[7] <= {BOW{1'b0}};
    block_out_tmp <= {BOW{1'b0}};
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0]   <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp  <= ((o_t0 - o_t7) >>> OUT_SHIFT); // block_out[7]
        block_out[3]   <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4]   <= ((o_t3 - o_t4) >>> OUT_SHIFT);
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
reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

wire [BOW-1:0] out_sel_u = block_out[ptr_q[2:0]];

assign outport_valid_o = valid_q[6];
assign outport_data_o  = {{(32-BOW){out_sel_u[BOW-1]}}, out_sel_u};
// Column pass emits column-major (slot = {column, row}); the streamed
// idx must label the row-major pixel position {slot[2:0], slot[5:3]}.
assign outport_idx_o   = (PERM_IDX != 0) ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

endmodule

//=====================================================================
// Double-buffered corner turn.
//
// The row pass emits slot {row[2:0], j[2:0]}; storage word address is
// {buffer, lane, offset} with lane = row[2:1] and offset = {row[0], j}
// — exactly the original's four 32x32 RAMs (lane selects the RAM,
// offset the word).  The column pass reads, per beat, all four lanes
// at one offset, so each read lane is a 32:1 mux of 20-bit words.
//=====================================================================
module jpeg_idct_transpose
(
     input            clk
    ,input            rst
    ,input            img_start_i
    ,input            inport_valid_i
    ,input  [ 19:0]   inport_data_i
    ,input  [  5:0]   inport_idx_i
    ,input            outport_ready_i
    ,output           outport_valid_o
    ,output [ 31:0]   outport_data0_o
    ,output [ 31:0]   outport_data1_o
    ,output [ 31:0]   outport_data2_o
    ,output [ 31:0]   outport_data3_o
    ,output [  2:0]   outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire wr0_w = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && (inport_idx_i[5:4] == 2'd3);
wire wr_any_w = wr0_w | wr1_w | wr2_w | wr3_w;

reg [19:0] bufs [0:127];
reg [19:0] rd0_q;
reg [19:0] rd1_q;
reg [19:0] rd2_q;
reg [19:0] rd3_q;

wire [6:0] wr_addr_w = {block_wr_q, inport_idx_i[5:4], inport_idx_i[3:0]};

always @ (posedge clk)
begin
    if (wr_any_w)
        bufs[wr_addr_w] <= inport_data_i;

    rd0_q <= bufs[{block_rd_q, 2'd0, rd_addr_q}];
    rd1_q <= bufs[{block_rd_q, 2'd1, rd_addr_q}];
    rd2_q <= bufs[{block_rd_q, 2'd2, rd_addr_q}];
    rd3_q <= bufs[{block_rd_q, 2'd3, rd_addr_q}];
end

assign outport_data0_o = {{12{rd0_q[19]}}, rd0_q};
assign outport_data1_o = {{12{rd1_q[19]}}, rd1_q};
assign outport_data2_o = {{12{rd2_q[19]}}, rd2_q};
assign outport_data3_o = {{12{rd3_q[19]}}, rd3_q};

//-----------------------------------------------------------------
// Buffer tracking
//-----------------------------------------------------------------
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
    if (wr_any_w && inport_idx_i == 6'd63)
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
// Read FSM (identical schedule to the original)
//-----------------------------------------------------------------
localparam STATE_W      = 2;
localparam STATE_IDLE   = 2'd0;
localparam STATE_SETUP  = 2'd1;
localparam STATE_ACTIVE = 2'd2;

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

endmodule

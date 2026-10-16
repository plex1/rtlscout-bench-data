//-----------------------------------------------------------------
// jpeg_idct_2d - full 8x8 2-D IDCT of the baseline JPEG decoder
// (ultraembedded/core_jpeg algorithm, AAN-style factorized IDCT).
//
// x-pass (OUT_SHIFT=11) -> double-buffered corner-turn -> y-pass
// (OUT_SHIFT=15).  Bit-exact with the original RTL: all arithmetic
// evaluated in 32-bit signed context with wrap-around; outputs of
// each pass arithmetic-shifted right; the 1/sqrt(2) odd scaling is
// (x*181)/256 with wrapping multiply and signed division truncating
// toward zero.
//
// Optimizations vs. the original:
//  - 7 fixed-constant multiplier units per pass (every beat of the
//    schedule uses distinct constants) with 2:1 data muxes, instead
//    of 5 time-shared 32x32 multipliers; constant CSD shift-add
//    trees split over the two product pipeline stages.
//  - 4 multipliers' worth of products re-scheduled onto free beats
//    (operand register persistence where the original re-wrote
//    identical values).
//  - (x*181)/256 as carry-save shift/add with truncating-divide
//    fixup carried into the output adder as a carry-in bit.
//  - Transpose buffers hold only the significant 20 bits of the
//    x-pass outputs; registered one-hot word select for the read.
//  - y-pass output index is the wiring permutation {p[2:0],p[5:3]}.
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

//-----------------------------------------------------------------
// Row pass (1-D IDCT on rows of the coefficient block)
//-----------------------------------------------------------------
wire        x_valid_w;
wire [19:0] x_data_w;
wire [ 5:0] x_idx_w;

jpeg_idct_1d #(.OUT_SHIFT(11), .IN_W(16), .OUT_W(20), .IDX_PERM(0))
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
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

//-----------------------------------------------------------------
// Corner turn (double buffered)
//-----------------------------------------------------------------
wire        t_valid_w;
wire [19:0] t_data0_w;
wire [19:0] t_data1_w;
wire [19:0] t_data2_w;
wire [19:0] t_data3_w;
wire [ 2:0] t_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o()
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
    ,.outport_ready_i(1'b1)
);

//-----------------------------------------------------------------
// Column pass (1-D IDCT on columns of the row-transformed block)
//-----------------------------------------------------------------
wire        y_valid_w;
wire [16:0] y_data_w;
wire [ 5:0] y_idx_w;

jpeg_idct_1d #(.OUT_SHIFT(15), .IN_W(20), .OUT_W(17), .IDX_PERM(1))
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(t_valid_w)
    ,.inport_data0_i(t_data0_w)
    ,.inport_data1_i(t_data1_w)
    ,.inport_data2_i(t_data2_w)
    ,.inport_data3_i(t_data3_w)
    ,.inport_idx_i(t_idx_w)
    ,.outport_valid_o(y_valid_w)
    ,.outport_data_o(y_data_w)
    ,.outport_idx_o(y_idx_w)
);

assign outport_valid_o = y_valid_w;
assign outport_data_o  = {{(32-17){y_data_w[16]}}, y_data_w};
assign outport_idx_o   = y_idx_w;

endmodule

//-----------------------------------------------------------------
// jpeg_idct_1d - one factorized (AAN-style) 1-D IDCT pass.
// Constant multiplier units U1..U7 compute data*C1..data*C7; the
// fixed 8-beat schedule always uses distinct constants per beat.
//-----------------------------------------------------------------
module jpeg_idct_1d #(
     parameter OUT_SHIFT = 11
    ,parameter IN_W      = 16
    ,parameter OUT_W     = 20
    ,parameter IDX_PERM  = 0
)
(
    // Inputs
     input               clk
    ,input               rst
    ,input               img_start_i
    ,input               img_end_i
    ,input               inport_valid_i
    ,input  [IN_W-1:0]   inport_data0_i
    ,input  [IN_W-1:0]   inport_data1_i
    ,input  [IN_W-1:0]   inport_data2_i
    ,input  [IN_W-1:0]   inport_data3_i
    ,input  [ 2:0]       inport_idx_i

    // Outputs
    ,output              outport_valid_o
    ,output [OUT_W-1:0]  outport_data_o
    ,output [ 5:0]       outport_idx_o
);

localparam AW = IN_W + 1;

wire signed [AW-1:0] block_in_0_1 = {inport_data0_i[IN_W-1], inport_data0_i};
wire signed [AW-1:0] block_in_2_3 = {inport_data1_i[IN_W-1], inport_data1_i};
wire signed [AW-1:0] block_in_4_5 = {inport_data2_i[IN_W-1], inport_data2_i};
wire signed [AW-1:0] block_in_6_7 = {inport_data3_i[IN_W-1], inport_data3_i};

//-----------------------------------------------------------------
// 3:2 carry-save add of three 32-bit operands (+ carry-in), result
// modulo 2^32:  one XOR/MAJ level followed by one carry-propagate
// adder (faster than two chained adders).
//-----------------------------------------------------------------
function automatic signed [31:0] add3ci;
    input signed [31:0] a;
    input signed [31:0] b;
    input signed [31:0] c;
    input               ci;
    reg [31:0] sm;
    reg [31:0] cy;
    begin
        sm = a ^ b ^ c;
        cy = {((a & b) | (b & c) | (a & c)), 1'b0};
        add3ci = $signed(sm) + $signed(cy | {31'b0, ci});
    end
endfunction

//-----------------------------------------------------------------
// Sum of four 32-bit terms with per-term negation (negmask[i]=1
// subtracts t_i), modulo 2^32:  two XOR/MAJ levels followed by one
// carry-propagate adder.  The compensation +1 for every negated
// term must be folded into t0 by the caller (t0 has low zero bits
// in all uses).
//-----------------------------------------------------------------
function automatic signed [31:0] csa4;
    input signed [31:0] t0;
    input signed [31:0] t1;
    input signed [31:0] t2;
    input signed [31:0] t3;
    input [3:0] negmask;
    reg [31:0] n0, n1, n2, n3;
    reg [31:0] s1, c1, s2, c2;
    begin
        n0 = negmask[0] ? ~t0 : t0;
        n1 = negmask[1] ? ~t1 : t1;
        n2 = negmask[2] ? ~t2 : t2;
        n3 = negmask[3] ? ~t3 : t3;
        s1 = n0 ^ n1 ^ n2;
        c1 = {((n0 & n1) | (n1 & n2) | (n0 & n2)), 1'b0};
        s2 = s1 ^ c1 ^ n3;
        c2 = {((s1 & c1) | (c1 & n3) | (s1 & n3)), 1'b0};
        csa4 = $signed(s2) + $signed(c2);
    end
endfunction

//-----------------------------------------------------------------
// Compress six 32-bit terms into two rows (sum, carry<<1), four
// XOR/MAJ levels; the rows are registered and the final carry
// propagate adder lives in the next pipeline stage.
//-----------------------------------------------------------------
function automatic [63:0] csa6to2;
    input [31:0] t0;
    input [31:0] t1;
    input [31:0] t2;
    input [31:0] t3;
    input [31:0] t4;
    input [31:0] t5;
    reg [31:0] s1, c1, s2, c2, s3, c3, s4, c4;
    begin
        s1 = t0 ^ t1 ^ t2;
        c1 = {((t0 & t1) | (t1 & t2) | (t0 & t2)), 1'b0};
        s2 = s1 ^ c1 ^ t3;
        c2 = {((s1 & c1) | (c1 & t3) | (s1 & t3)), 1'b0};
        s3 = s2 ^ c2 ^ t4;
        c3 = {((s2 & c2) | (c2 & t4) | (s2 & t4)), 1'b0};
        s4 = s3 ^ c3 ^ t5;
        c4 = {((s3 & c3) | (c3 & t5) | (s3 & t5)), 1'b0};
        csa6to2 = {s4, c4};
    end
endfunction

//-----------------------------------------------------------------
// Operand capture stage (per beat of the fixed input schedule).
// Each constant unit gets its data from a 2:1 lane mux; U4 (C4)
// takes i0 = c0+-c4.
//-----------------------------------------------------------------
reg signed [AW-1:0] i0s;
reg signed [AW-1:0] i0d;
reg signed [AW-1:0] u1_a;
reg signed [AW-1:0] u2_a;
reg signed [AW-1:0] u3_a;
reg signed [AW-1:0] u4_a;
reg signed [AW-1:0] u5_a;
reg signed [AW-1:0] u6_a;
reg signed [AW-1:0] u7_a;

always @ (posedge clk )
begin
    if (inport_idx_i == 3'd0)
        i0s <= block_in_0_1 + block_in_4_5;
    if (inport_idx_i == 3'd2)
        i0d <= block_in_0_1 - block_in_4_5;

    case (inport_idx_i)
    3'd0:
    begin
        u2_a <= block_in_2_3;     // c2*C2
        u6_a <= block_in_6_7;     // c6*C6
    end
    3'd1:
    begin
        u1_a <= block_in_0_1;     // c1*C1
        u7_a <= block_in_6_7;     // c7*C7
        u5_a <= block_in_4_5;     // c5*C5
        u3_a <= block_in_2_3;     // c3*C3
    end
    3'd2:
    begin
        u4_a <= i0s;              // (c0+c4)*C4
    end
    3'd3:
    begin
        u7_a <= block_in_0_1;     // c1*C7
        u1_a <= block_in_6_7;     // c7*C1
        u3_a <= block_in_4_5;     // c5*C3
        u5_a <= block_in_2_3;     // c3*C5
    end
    3'd5:
    begin
        u6_a <= block_in_2_3;     // c2*C6
        u2_a <= block_in_6_7;     // c6*C2
        u4_a <= i0d;              // (c0-c4)*C4
    end
    default:
        ;
    endcase
end

//-----------------------------------------------------------------
// Constant multiplier units (cosine constants x4096):
//   C1=4017, C2=3784, C3=3406, C4=2896, C5=2276, C6=1567, C7=799
// Stage 1 compresses the CSD terms to two rows, stage 2 is the
// final carry-propagate adder (products land two cycles after the
// operand capture, consumed at out idx = beat).
//-----------------------------------------------------------------
wire signed [31:0] u1_a32 = {{(32-AW){u1_a[AW-1]}}, u1_a};
wire signed [31:0] u2_a32 = {{(32-AW){u2_a[AW-1]}}, u2_a};
wire signed [31:0] u3_a32 = {{(32-AW){u3_a[AW-1]}}, u3_a};
wire signed [31:0] u4_a32 = {{(32-AW){u4_a[AW-1]}}, u4_a};
wire signed [31:0] u5_a32 = {{(32-AW){u5_a[AW-1]}}, u5_a};
wire signed [31:0] u6_a32 = {{(32-AW){u6_a[AW-1]}}, u6_a};
wire signed [31:0] u7_a32 = {{(32-AW){u7_a[AW-1]}}, u7_a};

// negated terms pre-complemented; the +1 per negation folded into a
// high shift term of the same unit (low bits zero there)
wire [63:0] u1_rows_w = csa6to2((u1_a32 <<< 12) | 32'd2, ~(u1_a32 <<< 6), ~(u1_a32 <<< 4), u1_a32, 32'd0, 32'd0);
wire [63:0] u2_rows_w = csa6to2((u2_a32 <<< 12) | 32'd2, ~(u2_a32 <<< 8), ~(u2_a32 <<< 6), (u2_a32 <<< 3), 32'd0, 32'd0);
wire [63:0] u3_rows_w = csa6to2((u3_a32 <<< 12) | 32'd2, ~(u3_a32 <<< 10), (u3_a32 <<< 8), (u3_a32 <<< 6), (u3_a32 <<< 4), ~(u3_a32 <<< 1));
wire [63:0] u4_rows_w = csa6to2((u4_a32 <<< 12) | 32'd2, ~(u4_a32 <<< 10), ~(u4_a32 <<< 8), (u4_a32 <<< 6), (u4_a32 <<< 4), 32'd0);
wire [63:0] u5_rows_w = csa6to2((u5_a32 <<< 11) | 32'd1, (u5_a32 <<< 8), ~(u5_a32 <<< 5), (u5_a32 <<< 2), 32'd0, 32'd0);
wire [63:0] u6_rows_w = csa6to2((u6_a32 <<< 11) | 32'd2, ~(u6_a32 <<< 9), (u6_a32 <<< 5), ~u6_a32, 32'd0, 32'd0);
wire [63:0] u7_rows_w = csa6to2((u7_a32 <<< 10) | 32'd2, ~(u7_a32 <<< 8), (u7_a32 <<< 5), ~u7_a32, 32'd0, 32'd0);

reg [63:0] u1_r_q;
reg [63:0] u2_r_q;
reg [63:0] u3_r_q;
reg [63:0] u4_r_q;
reg [63:0] u5_r_q;
reg [63:0] u6_r_q;
reg [63:0] u7_r_q;

always @ (posedge clk )
begin
    u1_r_q <= u1_rows_w;
    u2_r_q <= u2_rows_w;
    u3_r_q <= u3_rows_w;
    u4_r_q <= u4_rows_w;
    u5_r_q <= u5_rows_w;
    u6_r_q <= u6_rows_w;
    u7_r_q <= u7_rows_w;
end

reg signed [31:0] u1_p;
reg signed [31:0] u2_p;
reg signed [31:0] u3_p;
reg signed [31:0] u4_p;
reg signed [31:0] u5_p;
reg signed [31:0] u6_p;
reg signed [31:0] u7_p;

always @ (posedge clk )
begin
    u1_p <= $signed(u1_r_q[63:32]) + $signed(u1_r_q[31:0]);
    u2_p <= $signed(u2_r_q[63:32]) + $signed(u2_r_q[31:0]);
    u3_p <= $signed(u3_r_q[63:32]) + $signed(u3_r_q[31:0]);
    u4_p <= $signed(u4_r_q[63:32]) + $signed(u4_r_q[31:0]);
    u5_p <= $signed(u5_r_q[63:32]) + $signed(u5_r_q[31:0]);
    u6_p <= $signed(u6_r_q[63:32]) + $signed(u6_r_q[31:0]);
    u7_p <= $signed(u7_r_q[63:32]) + $signed(u7_r_q[31:0]);
end

//-----------------------------------------------------------------
// Valid / index pipelines
//-----------------------------------------------------------------
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

//-----------------------------------------------------------------
// Butterfly stage
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
reg               o_c5;
reg               o_c6;
reg signed [31:0] o_r65s;
reg signed [31:0] o_r65b;
reg signed [31:0] o_r56s;
reg signed [31:0] o_r56b;

// (x * 181) / 256 with wrapping multiply and signed division toward
// zero, spread over the pipeline:
//   181 = 256-64-16+4+1
//   step 1 (idx 6): s = 5x ;  b = 80x   (x = o_t6 -+ o_t5)
//   step 2 (idx 7): p = 256x + s - b = 181x (all mod 2^32)
//                   q0 = p >>> 8 (floor);  c = (p<0 && p[7:0]!=0)
//                   trunc(p/256) = q0 + c
//   step 3 (out):   block_out = (t +/- q0 +/- c) >>> OUT_SHIFT
wire signed [31:0] s65s_w = csa4((o_t6 <<< 2) | 32'd2, o_t6, (o_t5 <<< 2), o_t5, 4'b1100);
wire signed [31:0] s65b_w = csa4((o_t6 <<< 6) | 32'd2, (o_t6 <<< 4), (o_t5 <<< 6), (o_t5 <<< 4), 4'b1100);
wire signed [31:0] s56s_w = csa4((o_t5 <<< 2), o_t5, (o_t6 <<< 2), o_t6, 4'b0000);
wire signed [31:0] s56b_w = csa4((o_t5 <<< 6), (o_t5 <<< 4), (o_t6 <<< 6), (o_t6 <<< 4), 4'b0000);

wire signed [31:0] p65_w = csa4((o_t6 <<< 8) | 32'd2, o_r65s, (o_t5 <<< 8), o_r65b, 4'b1100);
wire signed [31:0] p56_w = csa4((o_t5 <<< 8) | 32'd1, (o_t6 <<< 8), o_r56s, o_r56b, 4'b1000);
wire               c65_w = p65_w[31] & (|p65_w[7:0]);
wire               c56_w = p56_w[31] & (|p56_w[7:0]);
wire signed [31:0] c5_ext = {31'b0, o_c5};
wire signed [31:0] c6_ext = {31'b0, o_c6};

always @ (posedge clk )
begin
    case (out_stg2_idx_q)
    3'd0:
    begin
        o_t3 <= u2_p + u6_p; // s3 = c2*C2 + c6*C6
    end
    3'd1:
    begin
        o_s7 <= u1_p + u7_p; // c1*C1 + c7*C7
        o_s6 <= u5_p + u3_p; // c5*C5 + c3*C3
    end
    3'd2:
    begin
        o_t0 <= u4_p + o_t3; // t0 = (c0+c4)*C4 + s3
        o_t3 <= u4_p - o_t3; // t3 = (c0+c4)*C4 - s3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3:
    begin
        o_t4 <= (u7_p - u1_p) + (u3_p - u5_p);
    end
    3'd4:
    begin
        o_t0 <= u7_p - u1_p; // s4 (products persist from beat 3)
        o_s5 <= u3_p - u5_p;
    end
    3'd5:
    begin
        o_t3 <= u6_p - u2_p; // s2
        o_t4 <= u4_p;        // s1 = (c0-c4)*C4
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_r65s <= s65s_w;
        o_r65b <= s65b_w;
        o_r56s <= s56s_w;
        o_r56b <= s56b_w;
    end
    default:
    begin
        o_s5 <= p65_w >>> 8; // 1/sqrt(2), fixup carry in o_c5
        o_s6 <= p56_w >>> 8; // 1/sqrt(2)
        o_c5 <= c65_w;
        o_c6 <= c56_w;
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

//-----------------------------------------------------------------
// Output block registers
//-----------------------------------------------------------------
wire signed [31:0] sh_r0 = (o_t0 + o_t7) >>> OUT_SHIFT;
wire signed [31:0] sh_r7 = (o_t0 - o_t7) >>> OUT_SHIFT;
wire signed [31:0] sh_r3 = (o_t3 + o_t4) >>> OUT_SHIFT;
wire signed [31:0] sh_r4 = (o_t3 - o_t4) >>> OUT_SHIFT;
wire signed [31:0] sh_r2 = add3ci(o_t2, o_s5, c5_ext, 1'b0) >>> OUT_SHIFT;
wire signed [31:0] sh_r5 = add3ci(o_t2, ~o_s5, {31'b0, ~o_c5}, 1'b0) >>> OUT_SHIFT;
wire signed [31:0] sh_r1 = add3ci(o_t1, o_s6, c6_ext, 1'b0) >>> OUT_SHIFT;
wire signed [31:0] sh_r6 = add3ci(o_t1, ~o_s6, {31'b0, ~o_c6}, 1'b0) >>> OUT_SHIFT;

reg signed [OUT_W-1:0] block_out[0:7];
reg signed [OUT_W-1:0] block_out_tmp;

always @ (posedge clk )
if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0]  <= sh_r0[OUT_W-1:0];
        block_out_tmp <= sh_r7[OUT_W-1:0]; // block_out[7]
        block_out[3]  <= sh_r3[OUT_W-1:0];
        block_out[4]  <= sh_r4[OUT_W-1:0];
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= sh_r2[OUT_W-1:0];
        block_out[5] <= sh_r5[OUT_W-1:0];
        block_out[1] <= sh_r1[OUT_W-1:0];
        block_out[6] <= sh_r6[OUT_W-1:0];
    end
end

//-----------------------------------------------------------------
// Output stream
//-----------------------------------------------------------------
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

wire [2:0] ptr_lo = ptr_q[2:0];
wire [2:0] ptr_hi = ptr_q[5:3];

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_lo];
assign outport_idx_o   = (IDX_PERM != 0) ? {ptr_lo, ptr_hi} : ptr_q;

endmodule

//-----------------------------------------------------------------
// jpeg_idct_transpose - double-buffered corner turn.
// Four 32 x 20 banks (written by lane group), read back in
// transposed order using a registered one-hot word select.
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
    // Inputs
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i

    // Outputs
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 19:0]  outport_data0_o
    ,output [ 19:0]  outport_data1_o
    ,output [ 19:0]  outport_data2_o
    ,output [ 19:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

localparam STATE_W      = 2;
localparam STATE_IDLE   = 2'd0;
localparam STATE_SETUP  = 2'd1;
localparam STATE_ACTIVE = 2'd2;

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [1:0]   block_ready_q;
reg [STATE_W-1:0] state_q;

//-----------------------------------------------------------------
// Storage: 4 banks x (2 blocks x 16 words) x 20 bits
//-----------------------------------------------------------------
wire [4:0] wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire wr_rdy_w = inport_valid_i && inport_accept_o;
wire wr0_w = wr_rdy_w && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = wr_rdy_w && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = wr_rdy_w && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = wr_rdy_w && (inport_idx_i[5:4] == 2'd3);

reg [19:0] ram0 [0:31];
reg [19:0] ram1 [0:31];
reg [19:0] ram2 [0:31];
reg [19:0] ram3 [0:31];

always @ (posedge clk)
begin
    if (wr0_w) ram0[wr_ptr_w] <= inport_data_i;
    if (wr1_w) ram1[wr_ptr_w] <= inport_data_i;
    if (wr2_w) ram2[wr_ptr_w] <= inport_data_i;
    if (wr3_w) ram3[wr_ptr_w] <= inport_data_i;
end

//-----------------------------------------------------------------
// Read address sequencing (identical to the original rd_addr_q):
//   IDLE   -> 0
//   SETUP  -> 8
//   ACTIVE -> step pattern = f(rd_idx_q + 1) (next cycle's value)
// The one-hot word select is registered with the NEXT address value
// so the mux select lines are flip-flop outputs (low fanout path).
//-----------------------------------------------------------------
wire [3:0] qm_w = {1'b0, rd_idx_q[5:3]};

reg [3:0] rd_addr_next_r;
always @ *
begin
    case (rd_idx_q[2:0])
    3'd0: rd_addr_next_r = qm_w;         // next j=1 -> {0,m}
    3'd1: rd_addr_next_r = qm_w + 4'd8;  // next j=2 -> {1,m}
    3'd2: rd_addr_next_r = qm_w + 4'd8;  // next j=3 -> {1,m}
    3'd3: rd_addr_next_r = qm_w;         // next j=4 -> {0,m}
    3'd4: rd_addr_next_r = qm_w + 4'd8;  // next j=5 -> {1,m}
    3'd5: rd_addr_next_r = qm_w;         // next j=6 -> {0,m}
    3'd6: rd_addr_next_r = qm_w + 4'd1;  // next j=7 -> m+1
    default: rd_addr_next_r = qm_w + 4'd9; // next j=0 -> {1,m+1}
    endcase
end

wire [3:0] addr_sel_w = (state_q == STATE_IDLE)  ? 4'd0 :
                        (state_q == STATE_SETUP) ? 4'd8 :
                                                   rd_addr_next_r;

reg [15:0] sel_q;
reg        blk_sel_q;

always @ (posedge clk)
if (rst)
begin
    sel_q     <= 16'd1;
    blk_sel_q <= 1'b0;
end
else if (img_start_i)
begin
    sel_q     <= 16'd1;
    blk_sel_q <= 1'b0;
end
else
begin
    sel_q     <= (16'd1 << addr_sel_w);
    blk_sel_q <= block_rd_q;
end

//-----------------------------------------------------------------
// Read mux (one-hot AND-OR per block half + registered block 2:1)
// + output registers
//-----------------------------------------------------------------
reg [19:0] rd0_lo;
reg [19:0] rd0_hi;
reg [19:0] rd1_lo;
reg [19:0] rd1_hi;
reg [19:0] rd2_lo;
reg [19:0] rd2_hi;
reg [19:0] rd3_lo;
reg [19:0] rd3_hi;

integer i;
always @ *
begin
    rd0_lo = 20'b0; rd0_hi = 20'b0;
    rd1_lo = 20'b0; rd1_hi = 20'b0;
    rd2_lo = 20'b0; rd2_hi = 20'b0;
    rd3_lo = 20'b0; rd3_hi = 20'b0;
    for (i = 0; i < 16; i = i + 1)
    begin
        rd0_lo = rd0_lo | (ram0[i]     & {20{sel_q[i]}});
        rd0_hi = rd0_hi | (ram0[16+i] & {20{sel_q[i]}});
        rd1_lo = rd1_lo | (ram1[i]     & {20{sel_q[i]}});
        rd1_hi = rd1_hi | (ram1[16+i] & {20{sel_q[i]}});
        rd2_lo = rd2_lo | (ram2[i]     & {20{sel_q[i]}});
        rd2_hi = rd2_hi | (ram2[16+i] & {20{sel_q[i]}});
        rd3_lo = rd3_lo | (ram3[i]     & {20{sel_q[i]}});
        rd3_hi = rd3_hi | (ram3[16+i] & {20{sel_q[i]}});
    end
end

reg [19:0] data0_q;
reg [19:0] data1_q;
reg [19:0] data2_q;
reg [19:0] data3_q;

always @ (posedge clk)
begin
    data0_q <= blk_sel_q ? rd0_hi : rd0_lo;
    data1_q <= blk_sel_q ? rd1_hi : rd1_lo;
    data2_q <= blk_sel_q ? rd2_hi : rd2_lo;
    data3_q <= blk_sel_q ? rd3_hi : rd3_lo;
end

//-----------------------------------------------------------------
// Control
//-----------------------------------------------------------------
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

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];
assign outport_data0_o = data0_q;
assign outport_data1_o = data1_q;
assign outport_data2_o = data2_q;
assign outport_data3_o = data3_q;

endmodule

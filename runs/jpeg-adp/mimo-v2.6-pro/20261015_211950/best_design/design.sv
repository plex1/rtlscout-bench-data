//-----------------------------------------------------------------
// jpeg_idct_2d - full 8x8 2-D IDCT (baseline JPEG decoder)
// x-pass (OUT_SHIFT=11) -> double-buffered transpose -> y-pass
// (OUT_SHIFT=15).  Behaviorally identical to the original core_jpeg
// composition (ultraembedded/core_jpeg, Apache 2.0):
//  * factorized AAN-style 1-D pass, C1..C7 scaled by 4096,
//  * 5 time-shared multipliers over the fixed 8-beat schedule,
//  * all arithmetic in 32-bit signed context with wrap-around,
//  * outputs arithmetic-shifted right (11 then 15),
//  * odd-part scaling (x * 181) / 256 with signed truncating
//    division toward zero.
// Optimizations vs. the original RTL (same behavior):
//  * multiplies by constants expanded to CSD shift/add trees,
//  * (x*181)/256 computed as shift/add + truncating shift,
//  * shared parameterized 1-D pass for row and column passes.
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
// Row pass (sign-extends 16-bit inputs)
//-----------------------------------------------------------------
wire signed [31:0] x_in0 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] x_in1 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] x_in2 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] x_in3 = {{16{inport_data3_i[15]}}, inport_data3_i};

wire          idct_x_valid_w;
wire [ 31:0]  idct_x_data_w;
wire [  5:0]  idct_x_idx_w;

idct1d #(.OUT_SHIFT(11), .FLIP(0))
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(x_in0)
    ,.inport_data1_i(x_in1)
    ,.inport_data2_i(x_in2)
    ,.inport_data3_i(x_in3)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(idct_x_valid_w)
    ,.outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

//-----------------------------------------------------------------
// Double-buffered corner-turn
//-----------------------------------------------------------------
wire          transpose_valid_w;
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

transpose
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

//-----------------------------------------------------------------
// Column pass
//-----------------------------------------------------------------
idct1d #(.OUT_SHIFT(15), .FLIP(1))
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
// idct1d - one factorized (AAN) 8-point IDCT pass.
// 4 input lanes carry (c0,c2,c4,c6) on beats idx in {0,2,5,7} and
// (c1,c3,c5,c7) on beats idx in {1,3,4,6}.  The whole datapath runs
// every clock (no backpressure); valid/idx are pipelined alongside.
//-----------------------------------------------------------------
module idct1d
#(
     parameter OUT_SHIFT = 11
    ,parameter FLIP      = 0
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

wire signed [31:0] block_in_0_1 = inport_data0_i;
wire signed [31:0] block_in_2_3 = inport_data1_i;
wire signed [31:0] block_in_4_5 = inport_data2_i;
wire signed [31:0] block_in_6_7 = inport_data3_i;

//-----------------------------------------------------------------
// Constant multiplies as CSD shift/add (result wraps at 32 bits,
// identical to a 32-bit signed multiply truncated to 32 bits).
// Each product is split so that no cycle contains more than two
// 32-bit adder levels: for a two-constant multiplier one constant is
// registered as a complete product (2 levels) and the other as two
// halves (2 levels) that are summed and selected one stage later.
// m01x/m01y (sel[0]): halves of the sel[1]=0 product (3784 / 4017)
// m01b (sel[0]):      full  sel[1]=1 product ( 799 / 1567)
// m23a:               full  s=0 product (2276)
// m23x/m23y:          halves of the s=1 product (3406)
// m4x/m4y:            halves of the only product (2896)
function [31:0] m01x;
    input [31:0] a;
    input        s;
    begin
        if (s)
            m01x = (a << 12) - (a << 6);                         // 4032
        else
            m01x = (a << 12) - (a << 8);                         // 3840
    end
endfunction

function [31:0] m01y;
    input [31:0] a;
    input        s;
    begin
        if (s)
            m01y = a - (a << 4);                                 // -15
        else
            m01y = (a << 3) - (a << 6);                          // -56
    end
endfunction

function [31:0] m01bx;
    input [31:0] a;
    input        s;
    begin
        if (s)
            m01bx = (a << 11) - (a << 9);                        // 1536
        else
            m01bx = (a << 10) - (a << 8);                        //  768
    end
endfunction

function [31:0] m01by;
    input [31:0] a;
    begin
        m01by = (a << 5) - a;                                    //   31
    end
endfunction

// Multi-operand signed sums use carry-save (3:2) compression:
// s = a^b^c, k = maj(a,b,c)<<1 with s + k = a+b+c (mod 2^32).  The
// +1 of a two's-complement negation ~z+1 folds into the carry wiring
// ({k[30:0],1'b1}) when one operand is negated; with two negated
// operands one input is ~z+2 instead.  Result: four-operand sums in
// two bitwise CSA levels plus one adder instead of three adders.
function [31:0] m23a;
    input [31:0] a;
    begin
        m23a = (a << 11) + (a << 8) + (a << 2) - (a << 5);       // 2276
    end
endfunction

function [31:0] m23x;
    input [31:0] a;
    begin
        m23x = (a << 12) - (a << 10);                            // 3072
    end
endfunction

function [31:0] m23y;
    input [31:0] a;
    begin
        m23y = (a << 8) + (a << 6) + (a << 4) - (a << 1);        //  334
    end
endfunction

function [31:0] m4x;
    input [31:0] a;
    begin
        m4x = (a << 12) - (a << 10);                             // 3072
    end
endfunction

function [31:0] m4y;
    input [31:0] a;
    begin
        m4y = (a << 6) + (a << 4) - (a << 8);                    // -176
    end
endfunction

// 181*x split in halves (wraps at 32 bits); the /256 with signed
// division truncating toward zero is applied one beat later.
function [31:0] p181h;
    input [31:0] x;
    begin
        p181h = (x << 8) - (x << 6);                             // 192*x
    end
endfunction

function [31:0] p181l;
    input [31:0] x;
    begin
        p181l = (x << 2) + x - (x << 4);                          // -11*x
    end
endfunction

function [31:0] trunc256;
    input [31:0] p;
    reg signed [31:0] q;
    begin
        q = $signed(p) >>> 8;
        if (p[31] && (p[7:0] != 8'd0))
            q = q + 32'sd1;
        trunc256 = q;
    end
endfunction

//-----------------------------------------------------------------
// IDCT
//-----------------------------------------------------------------
reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg        [1:0]  mul0_s;
reg signed [31:0] mul1_a;
reg        [1:0]  mul1_s;
reg signed [31:0] mul2_a;
reg               mul2_s;
reg signed [31:0] mul3_a;
reg               mul3_s;
reg signed [31:0] mul4_a;

always @ (posedge clk )
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0;
    mul0_s <= 2'b0;
    mul1_a <= 32'b0;
    mul1_s <= 2'b0;
    mul2_a <= 32'b0;
    mul2_s <= 1'b0;
    mul3_a <= 32'b0;
    mul3_s <= 1'b0;
    mul4_a <= 32'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul0_s <= 2'd0; // C2
        mul1_a <= block_in_6_7;
        mul1_s <= 2'd3; // C6
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul0_s <= 2'd1; // C1
        mul1_a <= block_in_6_7;
        mul1_s <= 2'd2; // C7
        mul2_a <= block_in_4_5;
        mul2_s <= 1'd0; // C5
        mul3_a <= block_in_2_3;
        mul3_s <= 1'd1; // C3
        mul4_a <= i0;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul0_s <= 2'd2; // C7
        mul1_a <= block_in_6_7;
        mul1_s <= 2'd1; // C1
        mul2_a <= block_in_4_5;
        mul2_s <= 1'd1; // C3
        mul3_a <= block_in_2_3;
        mul3_s <= 1'd0; // C5
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul0_s <= 2'd2; // C7
        mul1_a <= block_in_6_7;
        mul1_s <= 2'd1; // C1
        mul2_a <= block_in_4_5;
        mul2_s <= 1'd1; // C3
        mul3_a <= block_in_2_3;
        mul3_s <= 1'd0; // C5
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_s <= 2'd3; // C6
        mul1_a <= block_in_6_7;
        mul1_s <= 2'd0; // C2
        mul4_a <= i0;
    end
    default:
        ;
    endcase
end

reg signed [31:0] mul0_qx;
reg signed [31:0] mul0_qy;
reg signed [31:0] mul0_qbx;
reg signed [31:0] mul0_qby;
reg               mul0_sq;
reg signed [31:0] mul1_qx;
reg signed [31:0] mul1_qy;
reg signed [31:0] mul1_qbx;
reg signed [31:0] mul1_qby;
reg               mul1_sq;
reg signed [31:0] mul2_qa;
reg signed [31:0] mul2_qx;
reg signed [31:0] mul2_qy;
reg               mul2_sq;
reg signed [31:0] mul3_qa;
reg signed [31:0] mul3_qx;
reg signed [31:0] mul3_qy;
reg               mul3_sq;
reg signed [31:0] mul4_qx;
reg signed [31:0] mul4_qy;

always @ (posedge clk )
if (rst)
begin
    mul0_qx <= 32'b0;
    mul0_qy <= 32'b0;
    mul0_qbx <= 32'b0;
    mul0_qby <= 32'b0;
    mul0_sq <= 1'b0;
    mul1_qx <= 32'b0;
    mul1_qy <= 32'b0;
    mul1_qbx <= 32'b0;
    mul1_qby <= 32'b0;
    mul1_sq <= 1'b0;
    mul2_qa <= 32'b0;
    mul2_qx <= 32'b0;
    mul2_qy <= 32'b0;
    mul2_sq <= 1'b0;
    mul3_qa <= 32'b0;
    mul3_qx <= 32'b0;
    mul3_qy <= 32'b0;
    mul3_sq <= 1'b0;
    mul4_qx <= 32'b0;
    mul4_qy <= 32'b0;
end
else
begin
    mul0_qx <= m01x(mul0_a, mul0_s[0]);
    mul0_qy <= m01y(mul0_a, mul0_s[0]);
    mul0_qbx <= m01bx(mul0_a, mul0_s[0]);
    mul0_qby <= m01by(mul0_a);
    mul0_sq <= mul0_s[1];
    mul1_qx <= m01x(mul1_a, mul1_s[0]);
    mul1_qy <= m01y(mul1_a, mul1_s[0]);
    mul1_qbx <= m01bx(mul1_a, mul1_s[0]);
    mul1_qby <= m01by(mul1_a);
    mul1_sq <= mul1_s[1];
    mul2_qa <= m23a(mul2_a);
    mul2_qx <= m23x(mul2_a);
    mul2_qy <= m23y(mul2_a);
    mul2_sq <= mul2_s;
    mul3_qa <= m23a(mul3_a);
    mul3_qx <= m23x(mul3_a);
    mul3_qy <= m23y(mul3_a);
    mul3_sq <= mul3_s;
    mul4_qx <= m4x(mul4_a);
    mul4_qy <= m4y(mul4_a);
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
    mul0 <= mul0_sq ? (mul0_qbx + mul0_qby) : (mul0_qx + mul0_qy);
    mul1 <= mul1_sq ? (mul1_qbx + mul1_qby) : (mul1_qx + mul1_qy);
    mul2 <= mul2_sq ? (mul2_qx + mul2_qy) : mul2_qa;
    mul3 <= mul3_sq ? (mul3_qx + mul3_qy) : mul3_qa;
    mul4 <= mul4_qx + mul4_qy;
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
reg signed [31:0] o_t7;
reg signed [31:0] o_t6_5;
reg signed [31:0] o_t5_6;
reg signed [31:0] o_p5h;
reg signed [31:0] o_p5l;
reg signed [31:0] o_p6h;
reg signed [31:0] o_p6l;

// 4-operand signed sums, written flat so the adder chains can be
// merged into a single multi-operand adder network.
wire signed [31:0] t4_w  = mul0 + mul2 - mul1 - mul3;
wire signed [31:0] t56_w = o_t0 + o_s7 - o_s5 - o_s6;
wire signed [31:0] t65_w = o_s7 + o_s5 - o_s6 - o_t0;

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
    o_t7   <= 32'b0;
    o_t6_5 <= 32'b0;
    o_t5_6 <= 32'b0;
    o_p5h  <= 32'b0;
    o_p5l  <= 32'b0;
    o_p6h  <= 32'b0;
    o_p6l  <= 32'b0;
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
        o_t4 <= t4_w;
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
        o_t6_5 <= t65_w;
        o_t5_6 <= t56_w;
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_p5h <= p181h(o_t6_5);
        o_p5l <= p181l(o_t6_5);
        o_p6h <= p181h(o_t5_6);
        o_p6l <= p181l(o_t5_6);
    end
    default:
    begin
        o_s5 <= trunc256(o_p5h + o_p5l); // 1/sqrt(2)
        o_s6 <= trunc256(o_p6h + o_p6l); // 1/sqrt(2)
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

reg signed [31:0] block_out0;
reg signed [31:0] block_out1;
reg signed [31:0] block_out2;
reg signed [31:0] block_out3;
reg signed [31:0] block_out4;
reg signed [31:0] block_out5;
reg signed [31:0] block_out6;
reg signed [31:0] block_out7;
reg signed [31:0] block_out_tmp;

always @ (posedge clk )
if (rst)
begin
    block_out0 <= 32'b0;
    block_out1 <= 32'b0;
    block_out2 <= 32'b0;
    block_out3 <= 32'b0;
    block_out4 <= 32'b0;
    block_out5 <= 32'b0;
    block_out6 <= 32'b0;
    block_out7 <= 32'b0;
    block_out_tmp <= 32'b0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out0 <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT); // block_out7
        block_out3 <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out4 <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out7 <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out2 <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out5 <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out1 <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out6 <= ((o_t1 - o_s6) >>> OUT_SHIFT);
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

reg [31:0] block_out_mux;
always @ *
begin
    case (ptr_q[2:0])
    3'd0: block_out_mux = block_out0;
    3'd1: block_out_mux = block_out1;
    3'd2: block_out_mux = block_out2;
    3'd3: block_out_mux = block_out3;
    3'd4: block_out_mux = block_out4;
    3'd5: block_out_mux = block_out5;
    3'd6: block_out_mux = block_out6;
    default: block_out_mux = block_out7;
    endcase
end

assign outport_data_o  = block_out_mux;
assign outport_idx_o   = FLIP ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

endmodule

//-----------------------------------------------------------------
// transpose - double-buffered corner turn with four 32x32 dual-port
// RAMs (register arrays).  Identical control to the original.
//-----------------------------------------------------------------
module transpose
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
reg [1:0]   block_ready_q;
reg [1:0]   state_q;
reg [1:0]   next_state_r;
reg         out_valid_q;
reg [2:0]   out_idx_q;

wire        inport_accept_w = ~block_ready_q[block_wr_q];
wire        wr_en_w = inport_valid_i && inport_accept_w;

wire        wr0_w = wr_en_w && (inport_idx_i[5:4] == 2'd0);
wire        wr1_w = wr_en_w && (inport_idx_i[5:4] == 2'd1);
wire        wr2_w = wr_en_w && (inport_idx_i[5:4] == 2'd2);
wire        wr3_w = wr_en_w && (inport_idx_i[5:4] == 2'd3);

// Last read beat of the current block and the periodic re-alignment shift
wire        rd_last_w = (state_q == STATE_ACTIVE) && (rd_idx_q == 6'd63);
wire        rd_sh_w   = (state_q == STATE_ACTIVE) && (rd_idx_q[2:0] == 3'd7);

// Odd/even input beat pattern of the 1-D pass (idx 1,3,4,6 = odd set)
wire        tap_sel_w = (rd_idx_q[2:0] == 3'd1) || (rd_idx_q[2:0] == 3'd3) ||
                        (rd_idx_q[2:0] == 3'd4) || (rd_idx_q[2:0] == 3'd6);

wire [31:0] outport_data0_w;
wire [31:0] outport_data1_w;
wire [31:0] outport_data2_w;
wire [31:0] outport_data3_w;

tr_bank u_bank0
(
     .clk(clk)
    ,.wr_en(wr0_w)
    ,.wr_half(block_wr_q)
    ,.wr_data(inport_data_i)
    ,.rd_sh(rd_sh_w)
    ,.rd_half(block_rd_q)
    ,.tap_sel(tap_sel_w)
    ,.rd_data(outport_data0_w)
);

tr_bank u_bank1
(
     .clk(clk)
    ,.wr_en(wr1_w)
    ,.wr_half(block_wr_q)
    ,.wr_data(inport_data_i)
    ,.rd_sh(rd_sh_w)
    ,.rd_half(block_rd_q)
    ,.tap_sel(tap_sel_w)
    ,.rd_data(outport_data1_w)
);

tr_bank u_bank2
(
     .clk(clk)
    ,.wr_en(wr2_w)
    ,.wr_half(block_wr_q)
    ,.wr_data(inport_data_i)
    ,.rd_sh(rd_sh_w)
    ,.rd_half(block_rd_q)
    ,.tap_sel(tap_sel_w)
    ,.rd_data(outport_data2_w)
);

tr_bank u_bank3
(
     .clk(clk)
    ,.wr_en(wr3_w)
    ,.wr_half(block_wr_q)
    ,.wr_data(inport_data_i)
    ,.rd_sh(rd_sh_w)
    ,.rd_half(block_rd_q)
    ,.tap_sel(tap_sel_w)
    ,.rd_data(outport_data3_w)
);

//-----------------------------------------------------------------
// Input Buffer
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
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_w)
    begin
        block_ready_q[block_wr_q] <= 1'b1;
        block_wr_q                <= ~block_wr_q;
    end

    if (rd_last_w)
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
        if (rd_last_w)
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

// Registered outputs keep the (data, idx) pairing of the beat; the
// stream is therefore one cycle later than the original composition,
// which is invisible at the interface.
always @ (posedge clk )
if (rst)
begin
    out_valid_q <= 1'b0;
    out_idx_q   <= 3'b0;
end
else
begin
    out_valid_q <= (state_q == STATE_ACTIVE);
    out_idx_q   <= rd_idx_q[2:0];
end

assign outport_valid_o = out_valid_q;
assign outport_idx_o   = out_idx_q;
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

//-----------------------------------------------------------------
// tr_bank - one bank of the corner-turn: two 16-word halves of 32
// bits, implemented as shift registers.  Writes arrive in strictly
// increasing word order (x-pass idx[3:0] = 0..15), so a write is a
// shift-in at position 0: word a ends at position 15-a.  The fixed
// read pattern visits word pairs (8p+g, g=0..7) in order of g, so a
// one-position shift after every 8 read beats keeps the pair of the
// active column at the two fixed taps 15 (p=0) and 7 (p=1).  Read
// data is registered.
//-----------------------------------------------------------------
module tr_bank
(
     input           clk
    ,input           wr_en
    ,input           wr_half
    ,input  [ 31:0]  wr_data
    ,input           rd_sh
    ,input           rd_half
    ,input           tap_sel
    ,output [ 31:0]  rd_data
);

reg [31:0] h0 [0:15];
reg [31:0] h1 [0:15];
reg [31:0] rd_q;

wire sh0_w = (wr_en & ~wr_half) | (rd_sh & ~rd_half);
wire sh1_w = (wr_en &  wr_half) | (rd_sh &  rd_half);

wire [31:0] din0_w = (wr_en & ~wr_half) ? wr_data : 32'b0;
wire [31:0] din1_w = (wr_en &  wr_half) ? wr_data : 32'b0;

integer i;
always @ (posedge clk )
begin
    if (sh0_w)
    begin
        h0[0] <= din0_w;
        for (i = 1; i < 16; i = i + 1)
            h0[i] <= h0[i-1];
    end

    if (sh1_w)
    begin
        h1[0] <= din1_w;
        for (i = 1; i < 16; i = i + 1)
            h1[i] <= h1[i-1];
    end

    case ({rd_half, tap_sel})
    2'b00:   rd_q <= h0[15];
    2'b01:   rd_q <= h0[7];
    2'b10:   rd_q <= h1[15];
    default: rd_q <= h1[7];
    endcase
end

assign rd_data = rd_q;

endmodule

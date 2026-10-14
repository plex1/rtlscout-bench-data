//=====================================================================
// jpeg_idct_2d — full 8x8 2-D IDCT (core_jpeg algorithm, exact fixed
// point) — width + structure optimized, carry-save multipliers and a
// 2-cycle *181/256 path.
//
// Exactness arguments:
//  * x-pass outputs are (32-bit wrapped t) >>> 11.  All x-pass
//    intermediates are bounded by 2^30 in magnitude (max |t| =
//    2^15*(4017+3784+3406+2896) + 2^15*(4017+799+2276+1567) < 2^30),
//    so no wrap occurs and the results fit in 20 signed bits
//    [-2^19, 2^19).  Stored/extended with 20 bits — bit-exact.
//  * x-pass multiplier operands fit 17 bits, constants fit 13 bits.
//  * y-pass inputs fit 20 bits (i0 fits 21).  Products may exceed 32
//    bits, so all product arithmetic is done mod 2^32, matching the
//    original 32x32 wrap multiply exactly (the carry-save tree and the
//    final 32-bit add both operate mod 2^32; dropped carry-outs ARE
//    the wrap).
//  * The multiply a*C (C > 0, registered) uses the biased-operand
//    identity a*C = a_u*C - 2^(AW-1)*C with a_u = a + 2^(AW-1)
//    (unsigned).  Rows: a_u gated by C[k], shifted by k, plus one
//    constant correction row -(C<<(AW-1)) mod 2^32.  Reduced by a 3:2
//    CSA tree to (sum, carry), registered, then summed mod 2^32.
//  * x*181 mod 2^32 = CSA reduction of {x, x<<2, x<<4, x<<5, x<<7}
//    plus one 32-bit wrap add (identical to the original expression).
//    /256 truncating toward zero = arithmetic shift by 8 plus a
//    conditional increment for negative non-exact values; the
//    increment is folded into the carry-in of the following adder.
//  * Schedule note: the scaled values are registered one beat before
//    the block_out write that consumes them, so the *181 add and the
//    final truncating add live in different cycles.
//
// Corner turn restructure:
//  The y-pass at output row r, beat j consumes exactly the x-pass
//  outputs x[r + 8*(2m+o)] on lane m, where o = j[0]^j[2].  So: 8
//  banks of 16 entries (2 block parities x 8 slots); a row latch
//  grabs the whole bank at the start of each y row and the per-beat
//  lane muxes are just 2:1s — no 32:1 address decode anywhere.
//=====================================================================

// 2-stage pipelined multiply by a positive 13-bit constant (bit-exact
// 32-bit wrap semantics):
//   stage 1: biased-operand partial rows -> CSA tree -> (sum,carry) reg
//   stage 2: 32-bit wrap add of the two vectors
module idct_mul13 #(
    parameter AW = 21
) (
    input  wire                   clk,
    input  wire                   rst,
    input  wire signed [AW-1:0]   a_i,
    input  wire [12:0]            b_i,   // positive constant < 4096
    output reg  [31:0]            p_o
);
    // biased unsigned operand: a_u = a + 2^(AW-1)  (just flips the MSB)
    wire [AW-1:0] a_u = {~a_i[AW-1], a_i[AW-2:0]};

    // partial rows, all mod 2^32 (32-bit wrap == original semantics)
    wire [31:0] r [0:12];
    genvar k;
    generate
        for (k = 0; k < 12; k = k + 1) begin : ROWS
            assign r[k] = (a_u & {AW{b_i[k]}}) << k;
        end
    endgenerate
    // correction row: -(C<<(AW-1)) mod 2^32  (a positive 32-bit value:
    // the two's complement), making the total == a*C mod 2^32
    wire [31:0] c32 = {19'd0, b_i};
    assign r[12] = -(c32 << (AW-1));

    // 3:2 carry-save reducers (all mod 2^32)
    function [31:0] csa_s;
        input [31:0] x, y, z;
        begin
            csa_s = x ^ y ^ z;
        end
    endfunction
    function [31:0] csa_c;
        input [31:0] x, y, z;
        begin
            csa_c = ((x & y) | (x & z) | (y & z)) << 1;
        end
    endfunction

    // 13 rows -> 2 vectors, 5 CSA levels
    wire [31:0] a0 = csa_s(r[0],  r[1],  r[2]);
    wire [31:0] b0 = csa_c(r[0],  r[1],  r[2]);
    wire [31:0] a1 = csa_s(r[3],  r[4],  r[5]);
    wire [31:0] b1 = csa_c(r[3],  r[4],  r[5]);
    wire [31:0] a2 = csa_s(r[6],  r[7],  r[8]);
    wire [31:0] b2 = csa_c(r[6],  r[7],  r[8]);
    wire [31:0] a3 = csa_s(r[9],  r[10], r[11]);
    wire [31:0] b3 = csa_c(r[9],  r[10], r[11]);

    wire [31:0] c0 = csa_s(a0, b0, a1);
    wire [31:0] d0 = csa_c(a0, b0, a1);
    wire [31:0] c1 = csa_s(b1, a2, b2);
    wire [31:0] d1 = csa_c(b1, a2, b2);
    wire [31:0] c2 = csa_s(a3, b3, r[12]);
    wire [31:0] d2 = csa_c(a3, b3, r[12]);

    wire [31:0] e0 = csa_s(c0, d0, c1);
    wire [31:0] f0 = csa_c(c0, d0, c1);
    wire [31:0] e1 = csa_s(d1, c2, d2);
    wire [31:0] f1 = csa_c(d1, c2, d2);

    wire [31:0] g0 = csa_s(e0, f0, e1);
    wire [31:0] h0 = csa_c(e0, f0, e1);

    wire [31:0] s_w = csa_s(g0, h0, f1);
    wire [31:0] c_w = csa_c(g0, h0, f1);

    reg [31:0] s_q;
    reg [31:0] c_q;

    always @(posedge clk)
    if (rst)
    begin
        s_q <= 32'b0;
        c_q <= 32'b0;
    end
    else
    begin
        s_q <= s_w;
        c_q <= c_w;
    end

    always @(posedge clk)
    if (rst)
        p_o <= 32'b0;
    else
        p_o <= s_q + c_q;   // c_q already carries the 2x weight
endmodule

// combinational: m = 181*x mod 2^32 (wrap, as in the original), plus
// the /256 truncate-toward-zero correction bit
module idct_scale181 (
    input  wire signed [31:0] x,
    output wire        [31:0] m,
    output wire               incr
);
    wire [31:0] r0 = x;
    wire [31:0] r1 = x << 2;
    wire [31:0] r2 = x << 4;
    wire [31:0] r3 = x << 5;
    wire [31:0] r4 = x << 7;

    wire [31:0] s1 = r0 ^ r1 ^ r2;
    wire [31:0] c1 = ((r0 & r1) | (r0 & r2) | (r1 & r2)) << 1;
    wire [31:0] s2 = r3 ^ r4;
    wire [31:0] c2 = (r3 & r4) << 1;

    wire [31:0] s3 = s1 ^ c1 ^ s2;
    wire [31:0] c3 = ((s1 & c1) | (s1 & s2) | (c1 & s2)) << 1;

    wire [31:0] s4 = s3 ^ c3 ^ c2;
    wire [31:0] c4 = ((s3 & c3) | (s3 & c2) | (c3 & c2)) << 1;

    assign m = s4 + c4;
    // signed division by 256 truncating toward zero:
    // (m>>>8) + (m<0 && (m & 255) != 0)
    assign incr = m[31] & (m[7:0] != 8'b0);
endmodule

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
wire [19:0]   x_data_w;
wire [ 5:0]   x_idx_w;
wire          x_accept_w;

wire          tr_valid_w;
wire [19:0]   tr_data0_w;
wire [19:0]   tr_data1_w;
wire [19:0]   tr_data2_w;
wire [19:0]   tr_data3_w;
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
     output wire [19:0] outport_data_o,
     output wire [ 5:0] outport_idx_o
);

localparam [12:0] C1_16 = 13'd4017;
localparam [12:0] C2_16 = 13'd3784;
localparam [12:0] C3_16 = 13'd3406;
localparam [12:0] C4_16 = 13'd2896;
localparam [12:0] C5_16 = 13'd2276;
localparam [12:0] C6_16 = 13'd1567;
localparam [12:0] C7_16 = 13'd799;

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
reg [12:0] mul0_b;
reg [12:0] mul1_b;
reg [12:0] mul2_b;
reg [12:0] mul3_b;
reg [12:0] mul4_b;

always @(posedge clk)
if (rst)
begin
    i0     <= 17'sd0;
    mul0_a <= 17'sd0; mul0_b <= 13'd0;
    mul1_a <= 17'sd0; mul1_b <= 13'd0;
    mul2_a <= 17'sd0; mul2_b <= 13'd0;
    mul3_a <= 17'sd0; mul3_b <= 13'd0;
    mul4_a <= 17'sd0; mul4_b <= 13'd0;
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

wire [31:0] mul0;
wire [31:0] mul1;
wire [31:0] mul2;
wire [31:0] mul3;
wire [31:0] mul4;

idct_mul13 #(.AW(17)) u_mul0 (.clk(clk), .rst(rst), .a_i(mul0_a), .b_i(mul0_b), .p_o(mul0));
// raw beat-4 s5 / beat-1 s6 temporaries (time-multiplexed, reset not needed)
reg signed [31:0] o_s5;
reg signed [31:0] o_s6;
idct_mul13 #(.AW(17)) u_mul1 (.clk(clk), .rst(rst), .a_i(mul1_a), .b_i(mul1_b), .p_o(mul1));
idct_mul13 #(.AW(17)) u_mul2 (.clk(clk), .rst(rst), .a_i(mul2_a), .b_i(mul2_b), .p_o(mul2));
idct_mul13 #(.AW(17)) u_mul3 (.clk(clk), .rst(rst), .a_i(mul3_a), .b_i(mul3_b), .p_o(mul3));
idct_mul13 #(.AW(17)) u_mul4 (.clk(clk), .rst(rst), .a_i(mul4_a), .b_i(mul4_b), .p_o(mul4));

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @(posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
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

// *181/256 scaling, split over two beats (see idct_scale181)
wire [31:0] m5_w;
wire        i5_w;
wire [31:0] m6_w;
wire        i6_w;

idct_scale181 u_sc5 (.x(o_t6_5), .m(m5_w), .incr(i5_w));
idct_scale181 u_sc6 (.x(o_t5_6), .m(m6_w), .incr(i6_w));

reg [31:0] m5_q;
reg        i5_q;
reg [31:0] m6_q;
reg        i6_q;

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
    m5_q   <= 32'b0;
    i5_q   <= 1'b0;
    m6_q   <= 32'b0;
    i6_q   <= 1'b0;
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
        m5_q <= m5_w;
        i5_q <= i5_w;
        m6_q <= m6_w;
        i6_q <= i6_w;
    end
    endcase
end

// beat B: fold the >>>8 / truncation into the block_out adders
wire signed [31:0] a5 = $signed(m5_q) >>> 8;
wire signed [31:0] a6 = $signed(m6_q) >>> 8;

wire [31:0] sum_p2 = o_t2 + a5 + {31'b0, i5_q};      // o_t2 + s5
wire [31:0] sum_n5 = o_t2 + ~a5 + {31'b0, ~i5_q};    // o_t2 - s5
wire [31:0] sum_p1 = o_t1 + a6 + {31'b0, i6_q};      // o_t1 + s6
wire [31:0] sum_n6 = o_t1 + ~a6 + {31'b0, ~i6_q};    // o_t1 - s6

// results of (32-bit >>> 11) always fit in 20 signed bits
reg signed [19:0] block_out [0:7];
reg signed [19:0] block_out_tmp;

always @(posedge clk)
if (rst)
begin
    block_out[0] <= 20'sd0;
    block_out[1] <= 20'sd0;
    block_out[2] <= 20'sd0;
    block_out[3] <= 20'sd0;
    block_out[4] <= 20'sd0;
    block_out[5] <= 20'sd0;
    block_out[6] <= 20'sd0;
    block_out[7] <= 20'sd0;
    block_out_tmp <= 20'sd0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 11);
        block_out_tmp <= ((o_t0 - o_t7) >>> 11); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> 11);
        block_out[4] <= ((o_t3 - o_t4) >>> 11);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= sum_p2[31:11];
        block_out[5] <= sum_n5[31:11];
        block_out[1] <= sum_p1[31:11];
        block_out[6] <= sum_n6[31:11];
    end
end

reg [7:0] valid_q;

always @(posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

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
// Double-buffered corner turn, banked (see header note)
//=====================================================================
module jpeg_idct_transpose (
     input  wire        clk,
     input  wire        rst,
     input  wire        img_start_i,
     input  wire        img_end_i,
     input  wire        inport_valid_i,
     input  wire [19:0] inport_data_i,
     input  wire [ 5:0] inport_idx_i,
     input  wire        outport_ready_i,
     output wire        inport_accept_o,
     output wire        outport_valid_o,
     output wire [19:0] outport_data0_o,
     output wire [19:0] outport_data1_o,
     output wire [19:0] outport_data2_o,
     output wire [19:0] outport_data3_o,
     output wire [ 2:0] outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [1:0]   block_ready_q;

//-----------------------------------------------------------------
// Storage: 8 banks x (2 parities x 8 slots), 20-bit entries.
// x-out beat k writes bank k[2:0], slot k[5:3], parity block_wr_q.
// y row r reads bank r (slots 0..7) of parity block_rd_q.
//-----------------------------------------------------------------
reg [19:0] stor [0:127];    // index {bank[2:0], parity, slot[2:0]}

wire       wen    = inport_valid_i && inport_accept_o;
wire [6:0] waddr  = {inport_idx_i[2:0], block_wr_q, inport_idx_i[5:3]};

always @(posedge clk)
    if (wen)
        stor[waddr] <= inport_data_i;

//-----------------------------------------------------------------
// Row latch read path
//-----------------------------------------------------------------
wire [2:0] r_r    = rd_idx_q[5:3];
wire [2:0] j_r    = rd_idx_q[2:0];
wire       cap    = (rd_idx_q != 6'd63) && (j_r == 3'd7);
// at idx 63: block_rd_q toggles, so load row 0 of the parity AFTER toggle;
// at idx 0 (first ever row): load row 0 of current parity only if a block
// has actually been received (block_ready_q != 0).
wire       cap0   = (rd_idx_q == 6'd0) && (block_ready_q != 2'd0);
wire [2:0] cap_r  = cap ? (r_r + 3'd1) : 3'd0;

reg [19:0] row_latch [0:7];

always @(posedge clk)
if (cap || cap0)
begin
    row_latch[0] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd0}];
    row_latch[1] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd1}];
    row_latch[2] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd2}];
    row_latch[3] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd3}];
    row_latch[4] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd4}];
    row_latch[5] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd5}];
    row_latch[6] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd6}];
    row_latch[7] <= stor[{cap_r, block_rd_q ^ (rd_idx_q == 6'd63), 3'd7}];
end

// o = j[0] ^ j[2]  (1 for j in {1,3,4,6})
wire       o_w    = j_r[0] ^ j_r[2];

assign outport_data0_o = o_w ? row_latch[1] : row_latch[0];
assign outport_data1_o = o_w ? row_latch[3] : row_latch[2];
assign outport_data2_o = o_w ? row_latch[5] : row_latch[4];
assign outport_data3_o = o_w ? row_latch[7] : row_latch[6];

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
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
     input  wire [19:0] inport_data0_i,
     input  wire [19:0] inport_data1_i,
     input  wire [19:0] inport_data2_i,
     input  wire [19:0] inport_data3_i,
     input  wire [ 2:0] inport_idx_i,
     output wire        outport_valid_o,
     output wire [31:0] outport_data_o,
     output wire [ 5:0] outport_idx_o
);

localparam [12:0] C1_16 = 13'd4017;
localparam [12:0] C2_16 = 13'd3784;
localparam [12:0] C3_16 = 13'd3406;
localparam [12:0] C4_16 = 13'd2896;
localparam [12:0] C5_16 = 13'd2276;
localparam [12:0] C6_16 = 13'd1567;
localparam [12:0] C7_16 = 13'd799;

// x-pass outputs always fit in 20 signed bits
wire signed [19:0] block_in_0_1 = inport_data0_i;
wire signed [19:0] block_in_2_3 = inport_data1_i;
wire signed [19:0] block_in_4_5 = inport_data2_i;
wire signed [19:0] block_in_6_7 = inport_data3_i;

reg signed [20:0] i0;
reg signed [20:0] mul0_a;
reg signed [20:0] mul1_a;
reg signed [20:0] mul2_a;
reg signed [20:0] mul3_a;
reg signed [20:0] mul4_a;
reg [12:0] mul0_b;
reg [12:0] mul1_b;
reg [12:0] mul2_b;
reg [12:0] mul3_b;
reg [12:0] mul4_b;

always @(posedge clk)
if (rst)
begin
    i0     <= 21'sd0;
    mul0_a <= 21'sd0; mul0_b <= 13'd0;
    mul1_a <= 21'sd0; mul1_b <= 13'd0;
    mul2_a <= 21'sd0; mul2_b <= 13'd0;
    mul3_a <= 21'sd0; mul3_b <= 13'd0;
    mul4_a <= 21'sd0; mul4_b <= 13'd0;
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

wire [31:0] mul0;
wire [31:0] mul1;
wire [31:0] mul2;
wire [31:0] mul3;
wire [31:0] mul4;

idct_mul13 #(.AW(21)) u_mul0 (.clk(clk), .rst(rst), .a_i(mul0_a), .b_i(mul0_b), .p_o(mul0));
// raw beat-4 s5 / beat-1 s6 temporaries (time-multiplexed, reset not needed)
reg signed [31:0] o_s5;
reg signed [31:0] o_s6;
idct_mul13 #(.AW(21)) u_mul1 (.clk(clk), .rst(rst), .a_i(mul1_a), .b_i(mul1_b), .p_o(mul1));
idct_mul13 #(.AW(21)) u_mul2 (.clk(clk), .rst(rst), .a_i(mul2_a), .b_i(mul2_b), .p_o(mul2));
idct_mul13 #(.AW(21)) u_mul3 (.clk(clk), .rst(rst), .a_i(mul3_a), .b_i(mul3_b), .p_o(mul3));
idct_mul13 #(.AW(21)) u_mul4 (.clk(clk), .rst(rst), .a_i(mul4_a), .b_i(mul4_b), .p_o(mul4));

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @(posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
    out_stg3_valid_q <= 1'b0; out_stg3_idx_q <= 3'b0;
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

// *181/256 scaling, split over two beats (see idct_scale181)
wire [31:0] m5_w;
wire        i5_w;
wire [31:0] m6_w;
wire        i6_w;

idct_scale181 u_sc5 (.x(o_t6_5), .m(m5_w), .incr(i5_w));
idct_scale181 u_sc6 (.x(o_t5_6), .m(m6_w), .incr(i6_w));

reg [31:0] m5_q;
reg        i5_q;
reg [31:0] m6_q;
reg        i6_q;

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
    m5_q   <= 32'b0;
    i5_q   <= 1'b0;
    m6_q   <= 32'b0;
    i6_q   <= 1'b0;
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
        m5_q <= m5_w;
        i5_q <= i5_w;
        m6_q <= m6_w;
        i6_q <= i6_w;
    end
    endcase
end

// beat B: fold the >>>8 / truncation into the block_out adders
wire signed [31:0] a5 = $signed(m5_q) >>> 8;
wire signed [31:0] a6 = $signed(m6_q) >>> 8;

wire [31:0] sum_p2 = o_t2 + a5 + {31'b0, i5_q};      // o_t2 + s5
wire [31:0] sum_n5 = o_t2 + ~a5 + {31'b0, ~i5_q};    // o_t2 - s5
wire [31:0] sum_p1 = o_t1 + a6 + {31'b0, i6_q};      // o_t1 + s6
wire [31:0] sum_n6 = o_t1 + ~a6 + {31'b0, ~i6_q};    // o_t1 - s6

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
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> 15);
        block_out_tmp <= ((o_t0 - o_t7) >>> 15); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> 15);
        block_out[4] <= ((o_t3 - o_t4) >>> 15);
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= sum_p2[31:15];
        block_out[5] <= sum_n5[31:15];
        block_out[1] <= sum_p1[31:15];
        block_out[6] <= sum_n6[31:15];
    end
end

reg [7:0] valid_q;

always @(posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

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

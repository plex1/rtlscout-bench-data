// Optimized full 8x8 2-D IDCT (core_jpeg bit-exact).
// x-pass (OUT_SHIFT=11) -> double-buffered corner turn -> y-pass (OUT_SHIFT=15).
//
// Savings vs. original:
//  * Constant multipliers (C1..C7) as balanced shift/add trees (32-bit wrap
//    context == (a*C) mod 2^32 of the original 32x32 signed multiply).
//  * Row-pass operand registers 17-bit (coeffs are 16-bit signed;
//    (c0+c4)/(c0-c4) fit in [-2^16,2^16]).
//  * Row-pass outputs bounded by 2^19 -> 20-bit corner-turn storage
//    (all internal sums stay below 2^30, so no wrap).
//  * Column-pass operand registers 21-bit (20-bit stored values).
//  * (x*181)/256: p = (x*181) mod 2^32; truncating signed div by 256:
//    p>=0 : p>>>8 ; p<0 : (p+255)>>>8.
//  * Corner turn read path: the original read address sequence is
//    addr(i) = g(i) + 8*b(k(i)) with g = i[5:3] (stable for 8 beats) and
//    b = k[0]^k[2]; each lane therefore only ever needs the two words
//    {0,g} and {1,g} of its group.  Both are prefetched through 8:1 muxes
//    whose select is the next group id, and the per-beat value is a 2:1
//    mux -> very short logic depth.

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
     .IW_16     (1),
     .OPW       (17),
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
wire [19:0] t_data0_w;
wire [19:0] t_data1_w;
wire [19:0] t_data2_w;
wire [19:0] t_data3_w;
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

wire [19:0] y_data_w;

jpeg_idct_1d
#(
     .OUT_SHIFT (15),
     .IW_16     (0),
     .OPW       (21),
     .BOW       (20),
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
    ,.outport_data_o  (y_data_w)
    ,.outport_idx_o   (outport_idx_o)
);

assign outport_data_o = {{12{y_data_w[19]}}, y_data_w};

endmodule

//=====================================================================
// Shared 1-D factorized (AAN-style) IDCT pass.
//   IW_16 = 1 : 16-bit signed input lanes (row pass)
//   IW_16 = 0 : 20-bit signed input lanes (column pass)
//=====================================================================
module jpeg_idct_1d
#(
     parameter OUT_SHIFT = 11,
     parameter IW_16     = 1,
     parameter OPW       = 17,
     parameter BOW       = 20,
     parameter PERM_IDX  = 0
)
(
     input            clk
    ,input            rst
    ,input            img_start_i
    ,input            inport_valid_i
    ,input  [ 19:0]   inport_data0_i
    ,input  [ 19:0]   inport_data1_i
    ,input  [ 19:0]   inport_data2_i
    ,input  [ 19:0]   inport_data3_i
    ,input  [  2:0]   inport_idx_i
    ,output           outport_valid_o
    ,output [ 19:0]   outport_data_o
    ,output [  5:0]   outport_idx_o
);

// Operand registers: sign-extended copies of the inputs.
wire [OPW-1:0] a0 = (IW_16 != 0) ? {inport_data0_i[15], inport_data0_i[15:0]}
                                 : {inport_data0_i[19], inport_data0_i[19:0]};
wire [OPW-1:0] a2 = (IW_16 != 0) ? {inport_data1_i[15], inport_data1_i[15:0]}
                                 : {inport_data1_i[19], inport_data1_i[19:0]};
wire [OPW-1:0] a4 = (IW_16 != 0) ? {inport_data2_i[15], inport_data2_i[15:0]}
                                 : {inport_data2_i[19], inport_data2_i[19:0]};
wire [OPW-1:0] a6 = (IW_16 != 0) ? {inport_data3_i[15], inport_data3_i[15:0]}
                                 : {inport_data3_i[19], inport_data3_i[19:0]};

//---------------------------------------------------------------------
// Constant multiplier: b in {1..7} selects C1..C7.
//   C1 = 4017 = 4096 - 64 - 15  = (s12 - s6) - (s4 - a)
//   C2 = 3784 = 4096 - 256 - 56 = (s12 - s8) - (s6 - s3)
//   C3 = 3406 = 4096 - 512 - 160 - 18
//   C4 = 2896 = 4096 - 1024 - 160 - 16
//   C5 = 2276 = 2048 + 128 + 100
//   C6 = 1567 = 2048 - 512 + 31
//   C7 =  799 = 1024 - 256 + 31
// Evaluated in 32-bit wrap-around context == (a*C) mod 2^32.
//---------------------------------------------------------------------
function signed [31:0] cmul;
    input signed [OPW-1:0] a;
    input         [ 2:0] b;
    reg signed [31:0] a32;
    reg signed [31:0] s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12;
    reg signed [31:0] r1, r2, r3;
    begin
        a32 = a;
        s1  = a32 <<< 1;
        s2  = a32 <<< 2;
        s3  = a32 <<< 3;
        s4  = a32 <<< 4;
        s5  = a32 <<< 5;
        s6  = a32 <<< 6;
        s7  = a32 <<< 7;
        s8  = a32 <<< 8;
        s9  = a32 <<< 9;
        s10 = a32 <<< 10;
        s11 = a32 <<< 11;
        s12 = a32 <<< 12;
        case (b)
        3'd1: begin r1 = s12 - s6;  r2 = s4 - a32; cmul = r1 - r2; end
        3'd2: begin r1 = s12 - s8;  r2 = s6 - s3;  cmul = r1 - r2; end
        3'd3: begin r1 = s12 - s9;  r2 = s7 + s5;  r3 = s4 + s1;
                    cmul = (r1 - r2) - r3; end
        3'd4: begin r1 = s12 - s10; r2 = s7 + s5;  cmul = (r1 - r2) - s4; end
        3'd5: begin r1 = s11 + s7;  r2 = s6 + s5;  cmul = (r1 + r2) + s2; end
        3'd6: begin r1 = s11 - s9;  r2 = s5 - a32; cmul = r1 + r2; end
        default: begin r1 = s10 - s8; r2 = s5 - a32; cmul = r1 + r2; end
        endcase
    end
endfunction

// (x*181)/256 with signed division truncating toward zero.
// p = (x*181) mod 2^32 (as the original); then
//   p >= 0 : p >>> 8 (plain truncating division, no bias)
//   p <  0 : (p + 255) >>> 8
function signed [31:0] m181_div256;
    input signed [31:0] x;
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

reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0; mul1_q <= 32'b0; mul2_q <= 32'b0;
    mul3_q <= 32'b0; mul4_q <= 32'b0;
end
else
begin
    mul0_q <= cmul(mul0_a, mul0_b);
    mul1_q <= cmul(mul1_a, mul1_b);
    mul2_q <= cmul(mul2_a, mul2_b);
    mul3_q <= cmul(mul3_a, mul3_b);
    mul4_q <= cmul(mul4_a, mul4_b);
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0; mul1 <= 32'b0; mul2 <= 32'b0; mul3 <= 32'b0; mul4 <= 32'b0;
end
else
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

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
// Combining stage (beat-scheduled, 32-bit signed wrap-around)
//---------------------------------------------------------------------
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
    o_s5 <= 32'b0; o_s6 <= 32'b0; o_s7 <= 32'b0;
    o_t0 <= 32'b0; o_t1 <= 32'b0; o_t2 <= 32'b0; o_t3 <= 32'b0;
    o_t4 <= 32'b0; o_t5 <= 32'b0; o_t6 <= 32'b0; o_t7 <= 32'b0;
    o_t6_5 <= 32'b0; o_t5_6 <= 32'b0;
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
assign outport_data_o  = out_sel_u;
// Column pass emits column-major (slot = {column, row}); the streamed
// idx must label the row-major pixel position {slot[2:0], slot[5:3]}.
assign outport_idx_o   = (PERM_IDX != 0) ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

endmodule

//=====================================================================
// Double-buffered corner turn.
//
// x-pass output slot {row[2:0], j[2:0]}: lane = row[2:1], word = {row[0], j}.
// Storage: 8 arrays (buffer x lane) x 16 words x 20b.
//
// Read: the original address sequence is addr(i) = g(i) + 8*b(k(i)) with
// g = i[5:3] and b = k[0]^k[2] (verified against the original recurrence),
// so per lane and group only the two words {0,g} and {1,g} are needed.
// Both are prefetched (8:1 muxes, select = next group id, stable 8 cycles)
// and the per-beat word is a final 2:1 by b.
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
    ,output [ 19:0]   outport_data0_o
    ,output [ 19:0]   outport_data1_o
    ,output [ 19:0]   outport_data2_o
    ,output [ 19:0]   outport_data3_o
    ,output [  2:0]   outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;

wire wr0_w = inport_valid_i && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && (inport_idx_i[5:4] == 2'd3);

// storage: mem{buf*4 + lane}[word], word = {row[0], j[2:0]}
reg [19:0] mem0 [0:15];
reg [19:0] mem1 [0:15];
reg [19:0] mem2 [0:15];
reg [19:0] mem3 [0:15];
reg [19:0] mem4 [0:15];
reg [19:0] mem5 [0:15];
reg [19:0] mem6 [0:15];
reg [19:0] mem7 [0:15];

always @ (posedge clk)
begin
    if (wr0_w && !block_wr_q) mem0[inport_idx_i[3:0]] <= inport_data_i;
    if (wr1_w && !block_wr_q) mem1[inport_idx_i[3:0]] <= inport_data_i;
    if (wr2_w && !block_wr_q) mem2[inport_idx_i[3:0]] <= inport_data_i;
    if (wr3_w && !block_wr_q) mem3[inport_idx_i[3:0]] <= inport_data_i;
    if (wr0_w &&  block_wr_q) mem4[inport_idx_i[3:0]] <= inport_data_i;
    if (wr1_w &&  block_wr_q) mem5[inport_idx_i[3:0]] <= inport_data_i;
    if (wr2_w &&  block_wr_q) mem6[inport_idx_i[3:0]] <= inport_data_i;
    if (wr3_w &&  block_wr_q) mem7[inport_idx_i[3:0]] <= inport_data_i;
end

// prefetch: pf[buf][lane][b] <= word {b, g_next} of that buffer/lane
reg [19:0] pf [0:1][0:3][0:1];

wire [5:0] idx1_w   = rd_idx_q + 6'd1;
wire [2:0] g_next_w = idx1_w[5:3];
wire [3:0] glo_w    = {1'b0, g_next_w};
wire [3:0] ghi_w    = {1'b1, g_next_w};

always @ (posedge clk)
begin
    pf[0][0][0] <= mem0[glo_w];  pf[0][0][1] <= mem0[ghi_w];
    pf[0][1][0] <= mem1[glo_w];  pf[0][1][1] <= mem1[ghi_w];
    pf[0][2][0] <= mem2[glo_w];  pf[0][2][1] <= mem2[ghi_w];
    pf[0][3][0] <= mem3[glo_w];  pf[0][3][1] <= mem3[ghi_w];
    pf[1][0][0] <= mem4[glo_w];  pf[1][0][1] <= mem4[ghi_w];
    pf[1][1][0] <= mem5[glo_w];  pf[1][1][1] <= mem5[ghi_w];
    pf[1][2][0] <= mem6[glo_w];  pf[1][2][1] <= mem6[ghi_w];
    pf[1][3][0] <= mem7[glo_w];  pf[1][3][1] <= mem7[ghi_w];
end

// per-beat word select: b = k[0] ^ k[2] (k = beat index within the group)
wire b_w = rd_idx_q[0] ^ rd_idx_q[2];

assign outport_data0_o = block_rd_q ? (b_w ? pf[1][0][1] : pf[1][0][0])
                                    : (b_w ? pf[0][0][1] : pf[0][0][0]);
assign outport_data1_o = block_rd_q ? (b_w ? pf[1][1][1] : pf[1][1][0])
                                    : (b_w ? pf[0][1][1] : pf[0][1][0]);
assign outport_data2_o = block_rd_q ? (b_w ? pf[1][2][1] : pf[1][2][0])
                                    : (b_w ? pf[0][2][1] : pf[0][2][0]);
assign outport_data3_o = block_rd_q ? (b_w ? pf[1][3][1] : pf[1][3][0])
                                    : (b_w ? pf[0][3][1] : pf[0][3][0]);

//-----------------------------------------------------------------
// Buffer tracking
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

wire wr_any_w = wr0_w | wr1_w | wr2_w | wr3_w;

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

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

endmodule

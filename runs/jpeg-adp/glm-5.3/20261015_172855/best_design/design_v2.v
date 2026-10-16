//-----------------------------------------------------------------
// jpeg_idct_2d — optimized 8x8 2-D IDCT (row pass -> transpose -> col pass)
// Structure and per-cycle schedule identical to ultraembedded/core_jpeg.
// Optimizations (all value-preserving):
//  * transpose RAM narrowed to 20 bits (|x-pass out| <= 347K < 2^19, proven)
//  * transpose RAM split into per-buffer 16-entry arrays (shallower muxes)
//  * multiplier operands narrowed (a: 18b row / 21b col, b: 13b constant);
//    products computed exactly then truncated to 32 bits == original wrap
//  * /256 (signed, trunc toward 0) implemented explicitly
//  * dead logic removed (unused RAM ports, LUT idx permutation)
//-----------------------------------------------------------------
module jpeg_idct_2d
(
     input  wire        clk
    ,input  wire        rst
    ,input  wire        img_start_i
    ,input  wire        img_end_i
    ,input  wire        inport_valid_i
    ,input  wire [15:0] inport_data0_i
    ,input  wire [15:0] inport_data1_i
    ,input  wire [15:0] inport_data2_i
    ,input  wire [15:0] inport_data3_i
    ,input  wire [2:0]  inport_idx_i
    ,output wire        outport_valid_o
    ,output wire [31:0] outport_data_o
    ,output wire [5:0]  outport_idx_o
);

    wire        x_valid_w;
    wire [19:0] x_data_w;
    wire [5:0]  x_idx_w;

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
        ,.outport_valid_o(x_valid_w)
        ,.outport_data_o(x_data_w)
        ,.outport_idx_o(x_idx_w)
    );

    wire        t_valid_w;
    wire [19:0] t_data0_w;
    wire [19:0] t_data1_w;
    wire [19:0] t_data2_w;
    wire [19:0] t_data3_w;
    wire [2:0]  t_idx_w;

    jpeg_idct_transpose u_transpose
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

    jpeg_idct_y u_idct_y
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
// Row pass (OUT_SHIFT = 11)
//-----------------------------------------------------------------
module jpeg_idct_x
(
     input  wire        clk
    ,input  wire        rst
    ,input  wire        img_start_i
    ,input  wire        inport_valid_i
    ,input  wire signed [15:0] inport_data0_i
    ,input  wire signed [15:0] inport_data1_i
    ,input  wire signed [15:0] inport_data2_i
    ,input  wire signed [15:0] inport_data3_i
    ,input  wire [2:0]  inport_idx_i
    ,output wire        outport_valid_o
    ,output wire [19:0] outport_data_o
    ,output wire [5:0]  outport_idx_o
);

    localparam signed [12:0] C1_16 = 13'sd4017;
    localparam signed [12:0] C2_16 = 13'sd3784;
    localparam signed [12:0] C3_16 = 13'sd3406;
    localparam signed [12:0] C4_16 = 13'sd2896;
    localparam signed [12:0] C5_16 = 13'sd2276;
    localparam signed [12:0] C6_16 = 13'sd1567;
    localparam signed [12:0] C7_16 = 13'sd799;

    reg signed [17:0] i0;
    reg signed [17:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
    reg signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

    always @(posedge clk)
    if (rst)
    begin
        i0     <= 18'sd0;
        mul0_a <= 18'sd0; mul0_b <= 13'sd0;
        mul1_a <= 18'sd0; mul1_b <= 13'sd0;
        mul2_a <= 18'sd0; mul2_b <= 13'sd0;
        mul3_a <= 18'sd0; mul3_b <= 13'sd0;
        mul4_a <= 18'sd0; mul4_b <= 13'sd0;
    end
    else
    begin
        case (inport_idx_i)
        3'd0:
        begin
            i0     <= inport_data0_i + inport_data2_i;
            mul0_a <= inport_data1_i; mul0_b <= C2_16;
            mul1_a <= inport_data3_i; mul1_b <= C6_16;
        end
        3'd1:
        begin
            mul0_a <= inport_data0_i; mul0_b <= C1_16;
            mul1_a <= inport_data3_i; mul1_b <= C7_16;
            mul2_a <= inport_data2_i; mul2_b <= C5_16;
            mul3_a <= inport_data1_i; mul3_b <= C3_16;
            mul4_a <= i0;            mul4_b <= C4_16;
        end
        3'd2:
        begin
            i0 <= inport_data0_i - inport_data2_i;
        end
        3'd3, 3'd4:
        begin
            mul0_a <= inport_data0_i; mul0_b <= C7_16;
            mul1_a <= inport_data3_i; mul1_b <= C1_16;
            mul2_a <= inport_data2_i; mul2_b <= C3_16;
            mul3_a <= inport_data1_i; mul3_b <= C5_16;
        end
        3'd5:
        begin
            mul0_a <= inport_data1_i; mul0_b <= C6_16;
            mul1_a <= inport_data3_i; mul1_b <= C2_16;
            mul4_a <= i0;            mul4_b <= C4_16;
        end
        default:
            ;
        endcase
    end

    // 18 x 13 signed -> exact 31-bit product (|a*b| <= 65536*4017 < 2^31)
    wire signed [30:0] p0 = mul0_a * mul0_b;
    wire signed [30:0] p1 = mul1_a * mul1_b;
    wire signed [30:0] p2 = mul2_a * mul2_b;
    wire signed [30:0] p3 = mul3_a * mul3_b;
    wire signed [30:0] p4 = mul4_a * mul4_b;

    reg signed [30:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

    always @(posedge clk)
    if (rst)
    begin
        mul0_q <= 31'sd0; mul1_q <= 31'sd0; mul2_q <= 31'sd0;
        mul3_q <= 31'sd0; mul4_q <= 31'sd0;
    end
    else
    begin
        mul0_q <= p0; mul1_q <= p1; mul2_q <= p2; mul3_q <= p3; mul4_q <= p4;
    end

    reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

    always @(posedge clk)
    if (rst)
    begin
        mul0 <= 32'sd0; mul1 <= 32'sd0; mul2 <= 32'sd0;
        mul3 <= 32'sd0; mul4 <= 32'sd0;
    end
    else
    begin
        mul0 <= {{1{mul0_q[30]}}, mul0_q};
        mul1 <= {{1{mul1_q[30]}}, mul1_q};
        mul2 <= {{1{mul2_q[30]}}, mul2_q};
        mul3 <= {{1{mul3_q[30]}}, mul3_q};
        mul4 <= {{1{mul4_q[30]}}, mul4_q};
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

    reg signed [31:0] o_s5, o_s6, o_s7;
    reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
    reg signed [31:0] o_t6_5, o_t5_6;

    // (x * 181) with 32-bit wrap, then signed /256 truncating toward zero
    wire signed [39:0] pw5 = o_t6_5 * 40'sd181;
    wire signed [39:0] pw6 = o_t5_6 * 40'sd181;
    wire signed [31:0] pw5_32 = pw5[31:0];
    wire signed [31:0] pw6_32 = pw6[31:0];
    wire signed [31:0] pw5_sh = pw5_32 >>> 8;
    wire signed [31:0] pw6_sh = pw6_32 >>> 8;
    wire pw5_c = pw5_32[31] & (|pw5_32[7:0]);
    wire pw6_c = pw6_32[31] & (|pw6_32[7:0]);

    always @(posedge clk)
    if (rst)
    begin
        o_s5 <= 32'sd0; o_s6 <= 32'sd0; o_s7 <= 32'sd0;
        o_t0 <= 32'sd0; o_t1 <= 32'sd0; o_t2 <= 32'sd0; o_t3 <= 32'sd0;
        o_t4 <= 32'sd0; o_t5 <= 32'sd0; o_t6 <= 32'sd0; o_t7 <= 32'sd0;
        o_t6_5 <= 32'sd0; o_t5_6 <= 32'sd0;
    end
    else
    begin
        case (out_stg2_idx_q)
        3'd0:
        begin
            o_t3 <= mul0 + mul1;
        end
        3'd1:
        begin
            o_s7 <= mul0 + mul1;
            o_s6 <= mul2 + mul3;
            o_t0 <= mul4;
        end
        3'd2:
        begin
            o_t0 <= o_t0 + o_t3;
            o_t3 <= o_t0 - o_t3;
            o_t7 <= o_s6 + o_s7;
        end
        3'd3:
        begin
            o_t4 <= (mul0 - mul1) + (mul2 - mul3);
        end
        3'd4:
        begin
            o_t0 <= mul0 - mul1;
            o_s5 <= mul2 - mul3;
        end
        3'd5:
        begin
            o_t3 <= mul0 - mul1;
            o_t4 <= mul4;
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
            o_s5 <= pw5_c ? (pw5_sh + 32'sd1) : pw5_sh;
            o_s6 <= pw6_c ? (pw6_sh + 32'sd1) : pw6_sh;
        end
        endcase
    end

    reg        out_stg3_valid_q;
    reg [2:0]  out_stg3_idx_q;

    always @(posedge clk)
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

    reg signed [19:0] block_out[0:7];
    reg signed [19:0] block_out_tmp;

    always @(posedge clk)
    if (rst)
    begin
        block_out[0] <= 20'sd0; block_out[1] <= 20'sd0;
        block_out[2] <= 20'sd0; block_out[3] <= 20'sd0;
        block_out[4] <= 20'sd0; block_out[5] <= 20'sd0;
        block_out[6] <= 20'sd0; block_out[7] <= 20'sd0;
        block_out_tmp <= 20'sd0;
    end
    else if (out_stg3_valid_q)
    begin
        if (out_stg3_idx_q == 3'd3)
        begin
            block_out[0] <= (o_t0 + o_t7) >>> 11;
            block_out_tmp <= (o_t0 - o_t7) >>> 11;
            block_out[3] <= (o_t3 + o_t4) >>> 11;
            block_out[4] <= (o_t3 - o_t4) >>> 11;
        end

        if (out_stg3_idx_q == 3'd6)
            block_out[7] <= block_out_tmp;

        if (out_stg3_idx_q == 3'd7)
        begin
            block_out[2] <= (o_t2 + o_s5) >>> 11;
            block_out[5] <= (o_t2 - o_s5) >>> 11;
            block_out[1] <= (o_t1 + o_s6) >>> 11;
            block_out[6] <= (o_t1 - o_s6) >>> 11;
        end
    end

    reg [6:0] valid_q;

    always @(posedge clk)
    if (rst)
        valid_q <= 7'b0;
    else if (img_start_i)
        valid_q <= 7'b0;
    else
        valid_q <= {valid_q[5:0], out_stg3_valid_q};

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

//-----------------------------------------------------------------
// Column pass (OUT_SHIFT = 15)
//-----------------------------------------------------------------
module jpeg_idct_y
(
     input  wire        clk
    ,input  wire        rst
    ,input  wire        img_start_i
    ,input  wire        inport_valid_i
    ,input  wire signed [19:0] inport_data0_i
    ,input  wire signed [19:0] inport_data1_i
    ,input  wire signed [19:0] inport_data2_i
    ,input  wire signed [19:0] inport_data3_i
    ,input  wire [2:0]  inport_idx_i
    ,output wire        outport_valid_o
    ,output wire [31:0] outport_data_o
    ,output wire [5:0]  outport_idx_o
);

    localparam signed [12:0] C1_16 = 13'sd4017;
    localparam signed [12:0] C2_16 = 13'sd3784;
    localparam signed [12:0] C3_16 = 13'sd3406;
    localparam signed [12:0] C4_16 = 13'sd2896;
    localparam signed [12:0] C5_16 = 13'sd2276;
    localparam signed [12:0] C6_16 = 13'sd1567;
    localparam signed [12:0] C7_16 = 13'sd799;

    reg signed [20:0] i0;
    reg signed [20:0] mul0_a, mul1_a, mul2_a, mul3_a, mul4_a;
    reg signed [12:0] mul0_b, mul1_b, mul2_b, mul3_b, mul4_b;

    always @(posedge clk)
    if (rst)
    begin
        i0     <= 21'sd0;
        mul0_a <= 21'sd0; mul0_b <= 13'sd0;
        mul1_a <= 21'sd0; mul1_b <= 13'sd0;
        mul2_a <= 21'sd0; mul2_b <= 13'sd0;
        mul3_a <= 21'sd0; mul3_b <= 13'sd0;
        mul4_a <= 21'sd0; mul4_b <= 13'sd0;
    end
    else
    begin
        case (inport_idx_i)
        3'd0:
        begin
            i0     <= inport_data0_i + inport_data2_i;
            mul0_a <= inport_data1_i; mul0_b <= C2_16;
            mul1_a <= inport_data3_i; mul1_b <= C6_16;
        end
        3'd1:
        begin
            mul0_a <= inport_data0_i; mul0_b <= C1_16;
            mul1_a <= inport_data3_i; mul1_b <= C7_16;
            mul2_a <= inport_data2_i; mul2_b <= C5_16;
            mul3_a <= inport_data1_i; mul3_b <= C3_16;
            mul4_a <= i0;            mul4_b <= C4_16;
        end
        3'd2:
        begin
            i0 <= inport_data0_i - inport_data2_i;
        end
        3'd3, 3'd4:
        begin
            mul0_a <= inport_data0_i; mul0_b <= C7_16;
            mul1_a <= inport_data3_i; mul1_b <= C1_16;
            mul2_a <= inport_data2_i; mul2_b <= C3_16;
            mul3_a <= inport_data1_i; mul3_b <= C5_16;
        end
        3'd5:
        begin
            mul0_a <= inport_data1_i; mul0_b <= C6_16;
            mul1_a <= inport_data3_i; mul1_b <= C2_16;
            mul4_a <= i0;            mul4_b <= C4_16;
        end
        default:
            ;
        endcase
    end

    // 21 x 13 signed -> exact 34-bit product, truncated to 32 (== wrap)
    wire signed [33:0] p0 = mul0_a * mul0_b;
    wire signed [33:0] p1 = mul1_a * mul1_b;
    wire signed [33:0] p2 = mul2_a * mul2_b;
    wire signed [33:0] p3 = mul3_a * mul3_b;
    wire signed [33:0] p4 = mul4_a * mul4_b;

    reg signed [33:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;

    always @(posedge clk)
    if (rst)
    begin
        mul0_q <= 34'sd0; mul1_q <= 34'sd0; mul2_q <= 34'sd0;
        mul3_q <= 34'sd0; mul4_q <= 34'sd0;
    end
    else
    begin
        mul0_q <= p0; mul1_q <= p1; mul2_q <= p2; mul3_q <= p3; mul4_q <= p4;
    end

    reg signed [31:0] mul0, mul1, mul2, mul3, mul4;

    always @(posedge clk)
    if (rst)
    begin
        mul0 <= 32'sd0; mul1 <= 32'sd0; mul2 <= 32'sd0;
        mul3 <= 32'sd0; mul4 <= 32'sd0;
    end
    else
    begin
        mul0 <= mul0_q[31:0];
        mul1 <= mul1_q[31:0];
        mul2 <= mul2_q[31:0];
        mul3 <= mul3_q[31:0];
        mul4 <= mul4_q[31:0];
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

    reg signed [31:0] o_s5, o_s6, o_s7;
    reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
    reg signed [31:0] o_t6_5, o_t5_6;

    wire signed [39:0] pw5 = o_t6_5 * 40'sd181;
    wire signed [39:0] pw6 = o_t5_6 * 40'sd181;
    wire signed [31:0] pw5_32 = pw5[31:0];
    wire signed [31:0] pw6_32 = pw6[31:0];
    wire signed [31:0] pw5_sh = pw5_32 >>> 8;
    wire signed [31:0] pw6_sh = pw6_32 >>> 8;
    wire pw5_c = pw5_32[31] & (|pw5_32[7:0]);
    wire pw6_c = pw6_32[31] & (|pw6_32[7:0]);

    always @(posedge clk)
    if (rst)
    begin
        o_s5 <= 32'sd0; o_s6 <= 32'sd0; o_s7 <= 32'sd0;
        o_t0 <= 32'sd0; o_t1 <= 32'sd0; o_t2 <= 32'sd0; o_t3 <= 32'sd0;
        o_t4 <= 32'sd0; o_t5 <= 32'sd0; o_t6 <= 32'sd0; o_t7 <= 32'sd0;
        o_t6_5 <= 32'sd0; o_t5_6 <= 32'sd0;
    end
    else
    begin
        case (out_stg2_idx_q)
        3'd0:
        begin
            o_t3 <= mul0 + mul1;
        end
        3'd1:
        begin
            o_s7 <= mul0 + mul1;
            o_s6 <= mul2 + mul3;
            o_t0 <= mul4;
        end
        3'd2:
        begin
            o_t0 <= o_t0 + o_t3;
            o_t3 <= o_t0 - o_t3;
            o_t7 <= o_s6 + o_s7;
        end
        3'd3:
        begin
            o_t4 <= (mul0 - mul1) + (mul2 - mul3);
        end
        3'd4:
        begin
            o_t0 <= mul0 - mul1;
            o_s5 <= mul2 - mul3;
        end
        3'd5:
        begin
            o_t3 <= mul0 - mul1;
            o_t4 <= mul4;
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
            o_s5 <= pw5_c ? (pw5_sh + 32'sd1) : pw5_sh;
            o_s6 <= pw6_c ? (pw6_sh + 32'sd1) : pw6_sh;
        end
        endcase
    end

    reg        out_stg3_valid_q;
    reg [2:0]  out_stg3_idx_q;

    always @(posedge clk)
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

    reg signed [31:0] block_out[0:7];
    reg signed [31:0] block_out_tmp;

    always @(posedge clk)
    if (rst)
    begin
        block_out[0] <= 32'sd0; block_out[1] <= 32'sd0;
        block_out[2] <= 32'sd0; block_out[3] <= 32'sd0;
        block_out[4] <= 32'sd0; block_out[5] <= 32'sd0;
        block_out[6] <= 32'sd0; block_out[7] <= 32'sd0;
        block_out_tmp <= 32'sd0;
    end
    else if (out_stg3_valid_q)
    begin
        if (out_stg3_idx_q == 3'd3)
        begin
            block_out[0] <= (o_t0 + o_t7) >>> 15;
            block_out_tmp <= (o_t0 - o_t7) >>> 15;
            block_out[3] <= (o_t3 + o_t4) >>> 15;
            block_out[4] <= (o_t3 - o_t4) >>> 15;
        end

        if (out_stg3_idx_q == 3'd6)
            block_out[7] <= block_out_tmp;

        if (out_stg3_idx_q == 3'd7)
        begin
            block_out[2] <= (o_t2 + o_s5) >>> 15;
            block_out[5] <= (o_t2 - o_s5) >>> 15;
            block_out[1] <= (o_t1 + o_s6) >>> 15;
            block_out[6] <= (o_t1 - o_s6) >>> 15;
        end
    end

    reg [6:0] valid_q;

    always @(posedge clk)
    if (rst)
        valid_q <= 7'b0;
    else if (img_start_i)
        valid_q <= 7'b0;
    else
        valid_q <= {valid_q[5:0], out_stg3_valid_q};

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
    assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

//-----------------------------------------------------------------
// Corner-turn double buffer: 4 banks x 2 buffers x 16 x 20b
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input  wire        clk
    ,input  wire        rst
    ,input  wire        img_start_i
    ,input  wire        inport_valid_i
    ,input  wire [19:0] inport_data_i
    ,input  wire [5:0]  inport_idx_i
    ,output wire        inport_accept_o
    ,output wire        outport_valid_o
    ,output wire [19:0] outport_data0_o
    ,output wire [19:0] outport_data1_o
    ,output wire [19:0] outport_data2_o
    ,output wire [19:0] outport_data3_o
    ,output wire [2:0]  outport_idx_o
);

    reg        block_wr_q;
    reg        block_rd_q;
    reg [1:0]  block_ready_q;
    reg [5:0]  rd_idx_q;
    reg [3:0]  rd_addr_q;

    assign inport_accept_o = ~block_ready_q[block_wr_q];

    wire wr_en_w = inport_valid_i && inport_accept_o;

    wire [19:0] rd_data_w [0:3];

    genvar gb;
    generate
    for (gb = 0; gb < 4; gb = gb + 1)
    begin : g_bank
        reg [19:0] mem0 [0:15];
        reg [19:0] mem1 [0:15];

        wire bank_sel_w = (inport_idx_i[5:4] == gb[1:0]);

        always @(posedge clk)
        begin
            if (wr_en_w && bank_sel_w)
            begin
                if (block_wr_q)
                    mem1[inport_idx_i[3:0]] <= inport_data_i;
                else
                    mem0[inport_idx_i[3:0]] <= inport_data_i;
            end
        end

        reg [19:0] rd0_q, rd1_q;

        always @(posedge clk)
        begin
            rd0_q <= mem0[rd_addr_q];
            rd1_q <= mem1[rd_addr_q];
        end

        assign rd_data_w[gb] = block_rd_q ? rd1_q : rd0_q;
    end
    endgenerate

    assign outport_data0_o = rd_data_w[0];
    assign outport_data1_o = rd_data_w[1];
    assign outport_data2_o = rd_data_w[2];
    assign outport_data3_o = rd_data_w[3];

    //-----------------------------------------------------------------
    // Input buffer bookkeeping
    //-----------------------------------------------------------------
    always @(posedge clk)
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
        if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o)
        begin
            block_ready_q[block_wr_q] <= 1'b1;
            block_wr_q                <= ~block_wr_q;
        end

        if (outport_valid_o && rd_idx_q[5:0] == 6'd63)
        begin
            block_ready_q[block_rd_q] <= 1'b0;
            block_rd_q                <= ~block_rd_q;
        end
    end

    //-----------------------------------------------------------------
    // FSM
    //-----------------------------------------------------------------
    localparam STATE_IDLE   = 2'd0;
    localparam STATE_SETUP  = 2'd1;
    localparam STATE_ACTIVE = 2'd2;

    reg [1:0] state_q;
    reg [1:0] next_state_r;

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

    always @(posedge clk)
    if (rst)
        state_q <= STATE_IDLE;
    else
        state_q <= next_state_r;

    always @(posedge clk)
    if (rst)
        rd_idx_q <= 6'b0;
    else if (img_start_i)
        rd_idx_q <= 6'b0;
    else if (state_q == STATE_ACTIVE)
        rd_idx_q <= rd_idx_q + 6'd1;

    always @(posedge clk)
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

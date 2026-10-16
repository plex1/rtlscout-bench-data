// DIAGNOSTIC ONLY: v2 with the transpose replaced by a stub (area breakdown probe)
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

    assign t_valid_w = 1'b0;
    assign t_data0_w = 20'b0;
    assign t_data1_w = 20'b0;
    assign t_data2_w = 20'b0;
    assign t_data3_w = 20'b0;
    assign t_idx_w   = 3'b0;

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
            o_s5 <= pw5_c ? (pw5_sh + 32'sd1) : pw6_sh;
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

// PROBE: x-pass only (original core_jpeg jpeg_idct_x logic verbatim)
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

localparam [15:0] C1_16 = 4017;
localparam [15:0] C2_16 = 3784;
localparam [15:0] C3_16 = 3406;
localparam [15:0] C4_16 = 2896;
localparam [15:0] C5_16 = 2276;
localparam [15:0] C6_16 = 1567;
localparam [15:0] C7_16 = 799;

localparam OUT_SHIFT = 11;

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a, mul0_b, mul1_a, mul1_b, mul2_a, mul2_b;
reg signed [31:0] mul3_a, mul3_b, mul4_a, mul4_b;

always @ (posedge clk)
if (rst) begin
    i0 <= 32'b0; mul0_a <= 0; mul0_b <= 0; mul1_a <= 0; mul1_b <= 0;
    mul2_a <= 0; mul2_b <= 0; mul3_a <= 0; mul3_b <= 0; mul4_a <= 0; mul4_b <= 0;
end else begin
    case (inport_idx_i)
    3'd0: begin
        i0 <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= C2_16;
        mul1_a <= block_in_6_7; mul1_b <= C6_16;
    end
    3'd1: begin
        mul0_a <= block_in_0_1; mul0_b <= C1_16;
        mul1_a <= block_in_6_7; mul1_b <= C7_16;
        mul2_a <= block_in_4_5; mul2_b <= C5_16;
        mul3_a <= block_in_2_3; mul3_b <= C3_16;
        mul4_a <= i0;           mul4_b <= C4_16;
    end
    3'd2: i0 <= block_in_0_1 - block_in_4_5;
    3'd3: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd4: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd5: begin
        mul0_a <= block_in_2_3; mul0_b <= C6_16;
        mul1_a <= block_in_6_7; mul1_b <= C2_16;
        mul4_a <= i0;           mul4_b <= C4_16;
    end
    default: ;
    endcase
end

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;
always @ (posedge clk)
if (rst) begin
    mul0_q <= 0; mul1_q <= 0; mul2_q <= 0; mul3_q <= 0; mul4_q <= 0;
end else begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;
always @ (posedge clk)
if (rst) begin
    mul0 <= 0; mul1 <= 0; mul2 <= 0; mul3 <= 0; mul4 <= 0;
end else begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q; mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q; reg [2:0] out_stg0_idx_q;
reg        out_stg1_valid_q; reg [2:0] out_stg1_idx_q;
reg        out_stg2_valid_q; reg [2:0] out_stg2_idx_q;

always @ (posedge clk)
if (rst) begin
    out_stg0_valid_q <= 0; out_stg0_idx_q <= 0;
    out_stg1_valid_q <= 0; out_stg1_idx_q <= 0;
    out_stg2_valid_q <= 0; out_stg2_idx_q <= 0;
end else begin
    out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q;
    out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q;
end

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

always @ (posedge clk)
if (rst) begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0; o_t0 <= 0; o_t1 <= 0; o_t2 <= 0;
    o_t3 <= 0; o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end else begin
    case (out_stg2_idx_q)
    3'd0: o_t3 <= mul0 + mul1;
    3'd1: begin
        o_s7 <= mul0 + mul1;
        o_s6 <= mul2 + mul3;
        o_t0 <= mul4;
    end
    3'd2: begin
        o_t0 <= o_t0 + o_t3;
        o_t3 <= o_t0 - o_t3;
        o_t7 <= o_s6 + o_s7;
    end
    3'd3: o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    3'd4: begin
        o_t0 <= mul0 - mul1;
        o_s5 <= mul2 - mul3;
    end
    3'd5: begin
        o_t3 <= mul0 - mul1;
        o_t4 <= mul4;
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6: begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5;
        o_t5_6 <= o_t5 + o_t6;
    end
    default: begin
        o_s5 <= (o_t6_5 * 181) / 256;
        o_s6 <= (o_t5_6 * 181) / 256;
    end
    endcase
end

reg        out_stg3_valid_q; reg [2:0] out_stg3_idx_q;
always @ (posedge clk)
if (rst) begin
    out_stg3_valid_q <= 0; out_stg3_idx_q <= 0;
end else begin
    out_stg3_valid_q <= out_stg2_valid_q;
    out_stg3_idx_q   <= out_stg2_idx_q;
end

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
if (rst) begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end else if (out_stg3_valid_q) begin
    if (out_stg3_idx_q == 3'd3) begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end
    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;
    if (out_stg3_idx_q == 3'd7) begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;
always @ (posedge clk)
if (rst) valid_q <= 8'b0;
else if (img_start_i) valid_q <= 8'b0;
else valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;
always @ (posedge clk)
if (rst) ptr_q <= 6'd0;
else if (img_start_i) ptr_q <= 6'd0;
else if (outport_valid_o) ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = ptr_q;

endmodule

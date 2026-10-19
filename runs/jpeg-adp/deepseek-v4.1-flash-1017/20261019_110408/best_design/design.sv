//-----------------------------------------------------------------
// jpeg_idct_2d  -- full 8x8 2-D IDCT (baseline JPEG)
// Row pass (OUT_SHIFT=11) -> double-buffered corner turn -> column
// pass (OUT_SHIFT=15).  Bit-exact reproduction of the fixed point
// behaviour of the original core_jpeg implementation.
//-----------------------------------------------------------------

//-----------------------------------------------------------------
// 1-D IDCT (AAN style), time shared 5 multipliers over an 8 beat
// schedule.  Input lanes:  even beats {0,2,5,7} carry (c0,c2,c4,c6)
//                          odd  beats {1,3,4,6} carry (c1,c3,c5,c7)
//-----------------------------------------------------------------
module idct1d #(
    parameter integer IN_W      = 16,
    parameter integer OUT_SHIFT = 11,
    parameter integer OUT_W     = 32
)(
    input  wire                       clk,
    input  wire                       rst,
    input  wire                       valid_i,
    input  wire [2:0]                 idx_i,
    input  wire signed [IN_W-1:0]     d0_i,
    input  wire signed [IN_W-1:0]     d1_i,
    input  wire signed [IN_W-1:0]     d2_i,
    input  wire signed [IN_W-1:0]     d3_i,
    output wire                       valid_o,
    output wire [5:0]                 cnt_o,
    output wire signed [OUT_W-1:0]    data_o
);

localparam signed [15:0] C1 = 16'sd4017;
localparam signed [15:0] C2 = 16'sd3784;
localparam signed [15:0] C3 = 16'sd3406;
localparam signed [15:0] C4 = 16'sd2896;
localparam signed [15:0] C5 = 16'sd2276;
localparam signed [15:0] C6 = 16'sd1567;
localparam signed [15:0] C7 = 16'sd799;

//-----------------------------------------------------------------
// multiplier operand select
//-----------------------------------------------------------------
wire signed [IN_W:0] i0_sum = d0_i + d2_i;   // c0 + c4 (even beat)
wire signed [IN_W:0] i0_dif = d0_i - d2_i;   // c0 - c4 (even beat)

reg signed [IN_W:0] i0q;

always @(posedge clk)
if (rst)
    i0q <= 0;
else if (idx_i == 3'd0)
    i0q <= i0_sum;
else if (idx_i == 3'd2)
    i0q <= i0_dif;

reg signed [IN_W-1:0] a0,a1,a2,a3;
reg signed [IN_W:0]   a4;
reg signed [15:0]     b0,b1,b2,b3,b4;

always @* begin
    a0 = d1_i; b0 = C2;
    a1 = d3_i; b1 = C6;
    a2 = d2_i; b2 = C5;
    a3 = d1_i; b3 = C3;
    a4 = i0q; b4 = C4;
    case (idx_i)
    3'd1: begin // odd set: c1,c3,c5,c7
        a0 = d0_i; b0 = C1;
        a1 = d3_i; b1 = C7;
        a2 = d2_i; b2 = C5;
        a3 = d1_i; b3 = C3;
    end
    3'd3, 3'd4: begin
        a0 = d0_i; b0 = C7;
        a1 = d3_i; b1 = C1;
        a2 = d2_i; b2 = C3;
        a3 = d1_i; b3 = C5;
    end
    3'd5: begin // even set
        a0 = d1_i; b0 = C6;
        a1 = d3_i; b1 = C2;
    end
    default: ;
    endcase
end

reg signed [31:0] p0,p1,p2,p3,p4;

always @(posedge clk)
if (rst) begin
    p0 <= 32'sd0; p1 <= 32'sd0; p2 <= 32'sd0; p3 <= 32'sd0; p4 <= 32'sd0;
end
else begin
    p0 <= a0 * b0;
    p1 <= a1 * b1;
    p2 <= a2 * b2;
    p3 <= a3 * b3;
    p4 <= a4 * b4;
end

//-----------------------------------------------------------------
// pipeline index / valid tracking (one extra stage for the adder tree)
//-----------------------------------------------------------------
reg [2:0] idx_d1, idx_d2;
reg       vld_d1, vld_d2;

always @(posedge clk)
if (rst) begin
    idx_d1 <= 3'd0; idx_d2 <= 3'd0; vld_d1 <= 1'b0; vld_d2 <= 1'b0;
end
else begin
    idx_d1 <= idx_i;  vld_d1 <= valid_i;
    idx_d2 <= idx_d1; vld_d2 <= vld_d1;
end

//-----------------------------------------------------------------
// butterfly accumulators
//-----------------------------------------------------------------
reg signed [31:0] o_s5,o_s6,o_s7,o_t0,o_t1,o_t2,o_t3,o_t4,o_t5,o_t6,o_t7;
reg signed [31:0] o_t6_5,o_t5_6;

always @(posedge clk)
if (rst) begin
    o_s5<=0; o_s6<=0; o_s7<=0; o_t0<=0; o_t1<=0; o_t2<=0; o_t3<=0;
    o_t4<=0; o_t5<=0; o_t6<=0; o_t7<=0; o_t6_5<=0; o_t5_6<=0;
end
else begin
    case (idx_d1)
    3'd0: begin
        o_t3 <= p0 + p1;                    // s3
    end
    3'd1: begin
        o_s7 <= p0 + p1;
        o_s6 <= p2 + p3;
        o_t0 <= p4;                         // s0
    end
    3'd2: begin
        o_t0 <= o_t0 + o_t3;                // t0
        o_t3 <= o_t0 - o_t3;                // t3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3: begin
        o_t4 <= (p0 - p1) + (p2 - p3);
    end
    3'd4: begin
        o_t0 <= p0 - p1;                    // s4
        o_s5 <= p2 - p3;
    end
    3'd5: begin
        o_t3 <= p0 - p1;                    // s2
        o_t4 <= p4;                         // s1
        o_t5 <= o_t0 - o_s5;
        o_t6 <= o_s7 - o_s6;
    end
    3'd6: begin
        o_t1   <= o_t4 + o_t3;
        o_t2   <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5;
        o_t5_6 <= o_t5 + o_t6;
    end
    default: begin
        o_s5 <= (o_t6_5 * 181) / 256;       // 1/sqrt(2)
        o_s6 <= (o_t5_6 * 181) / 256;
    end
    endcase
end

//-----------------------------------------------------------------
// row output buffer (double banked) / output streamer
//-----------------------------------------------------------------
reg signed [OUT_W-1:0] blk [0:7];
reg signed [31:0]      blk_tmp;

always @(posedge clk)
if (rst) begin
    blk[0]<=0; blk[1]<=0; blk[2]<=0; blk[3]<=0;
    blk[4]<=0; blk[5]<=0; blk[6]<=0; blk[7]<=0;
    blk_tmp <= 0;
end
else if (vld_d2) begin
    if (idx_d2 == 3'd3) begin
        blk[0]  <= (o_t0 + o_t7) >>> OUT_SHIFT;
        blk_tmp <= (o_t0 - o_t7) >>> OUT_SHIFT;
        blk[3]  <= (o_t3 + o_t4) >>> OUT_SHIFT;
        blk[4]  <= (o_t3 - o_t4) >>> OUT_SHIFT;
    end

    if (idx_d2 == 3'd6)
        blk[7] <= blk_tmp;

    if (idx_d2 == 3'd7) begin
        blk[2] <= (o_t2 + o_s5) >>> OUT_SHIFT;
        blk[5] <= (o_t2 - o_s5) >>> OUT_SHIFT;
        blk[1] <= (o_t1 + o_s6) >>> OUT_SHIFT;
        blk[6] <= (o_t1 - o_s6) >>> OUT_SHIFT;
    end
end

reg [7:0] vq;

always @(posedge clk)
if (rst)
    vq <= 8'd0;
else
    vq <= {vq[6:0], vld_d2};

reg [5:0] ptr_q;

always @(posedge clk)
if (rst) begin
    ptr_q <= 6'd0;
end
else begin
    if (valid_o) begin
        ptr_q <= ptr_q + 6'd1;
    end
end

assign valid_o = vq[6];
assign cnt_o   = ptr_q;
assign data_o  = blk[ptr_q[2:0]];

endmodule

//-----------------------------------------------------------------
// top level
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
// row pass
//-----------------------------------------------------------------
wire                x_valid;
wire [5:0]          x_cnt;
wire signed [31:0]  x_data;

idct1d #(.IN_W(16), .OUT_SHIFT(11), .OUT_W(32)) u_x
(
     .clk(clk)
    ,.rst(rst)
    ,.valid_i(inport_valid_i)
    ,.idx_i(inport_idx_i)
    ,.d0_i(inport_data0_i)
    ,.d1_i(inport_data1_i)
    ,.d2_i(inport_data2_i)
    ,.d3_i(inport_data3_i)
    ,.valid_o(x_valid)
    ,.cnt_o(x_cnt)
    ,.data_o(x_data)
);

//-----------------------------------------------------------------
// corner turn memory:  mem[blk][row][col]
//-----------------------------------------------------------------
reg signed [31:0] mem [0:1][0:7][0:7];

reg        wr_blk;
reg [1:0]  blk_rdy;
reg        rd_blk;
reg        rd_act;
reg [5:0]  rd_cnt;

wire [2:0] wr_row = x_cnt[5:3];
wire [2:0] wr_col = x_cnt[2:0];

always @(posedge clk)
if (x_valid)
    mem[wr_blk][wr_row][wr_col] <= x_data;

always @(posedge clk)
if (rst) begin
    wr_blk  <= 1'b0;
    blk_rdy <= 2'b0;
    rd_blk  <= 1'b0;
    rd_act  <= 1'b0;
    rd_cnt  <= 6'd0;
end
else begin
    if (x_valid && x_cnt == 6'd63) begin
        wr_blk          <= ~wr_blk;
        blk_rdy[wr_blk] <= 1'b1;
    end

    if (rd_act) begin
        if (rd_cnt == 6'd63) begin
            rd_act          <= 1'b0;
            blk_rdy[rd_blk] <= 1'b0;
            rd_blk          <= ~rd_blk;
            rd_cnt          <= 6'd0;
        end
        else
            rd_cnt <= rd_cnt + 6'd1;
    end
    else if (blk_rdy[rd_blk]) begin
        rd_act <= 1'b1;
        rd_cnt <= 6'd0;
    end
end

wire        rd_valid = rd_act;
wire [2:0]  rd_beat  = rd_cnt[2:0];
wire [2:0]  rd_col   = rd_cnt[5:3];

// parity: odd beats {1,3,4,6} take the odd rows
wire rd_par = ~(rd_beat == 3'd0) & ~(rd_beat == 3'd2) &
              ~(rd_beat == 3'd5) & ~(rd_beat == 3'd7);

// lane data:  lane j = mem[rd_blk][2j+parity][rd_col]
wire signed [31:0] q0 = mem[rd_blk][{2'd0,rd_par}][rd_col];
wire signed [31:0] q1 = mem[rd_blk][{2'd1,rd_par}][rd_col];
wire signed [31:0] q2 = mem[rd_blk][{2'd2,rd_par}][rd_col];
wire signed [31:0] q3 = mem[rd_blk][{2'd3,rd_par}][rd_col];

//-----------------------------------------------------------------
// column pass
//-----------------------------------------------------------------
wire                y_valid;
wire [5:0]          y_cnt;
wire signed [31:0]  y_data;

idct1d #(.IN_W(32), .OUT_SHIFT(15), .OUT_W(32)) u_y
(
     .clk(clk)
    ,.rst(rst)
    ,.valid_i(rd_valid)
    ,.idx_i(rd_beat)
    ,.d0_i(q0)
    ,.d1_i(q1)
    ,.d2_i(q2)
    ,.d3_i(q3)
    ,.valid_o(y_valid)
    ,.cnt_o(y_cnt)
    ,.data_o(y_data)
);

assign outport_valid_o = y_valid;
assign outport_data_o  = y_data;
assign outport_idx_o   = {y_cnt[2:0], y_cnt[5:3]};

endmodule

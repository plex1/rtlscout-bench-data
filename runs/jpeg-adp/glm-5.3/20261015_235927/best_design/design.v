//-----------------------------------------------------------------
// 1-D IDCT stage (AAN factorized, 5 multipliers, 8-beat schedule)
// Streams 8 results per row starting 7 cycles after the row's idx-3 beat.
//-----------------------------------------------------------------
module idct1d #(parameter SH = 11)
(
    input  clk, rst, img_start,
    input  in_valid,
    input  [2:0] in_idx,
    input  signed [31:0] d0, d1, d2, d3,
    output reg        out_valid,
    output reg [31:0] out_data,
    output reg [5:0]  out_cnt
);
    localparam signed [15:0] C1 = 4017, C2 = 3784, C3 = 3406,
                             C4 = 2896, C5 = 2276, C6 = 1567, C7 = 799;

    // stage 0: operand setup
    reg signed [31:0] i0, ma0, mb0, ma1, mb1, ma2, mb2, ma3, mb3, ma4, mb4;
    always @(posedge clk) begin
        if (rst) begin
            i0<=0; ma0<=0; mb0<=0; ma1<=0; mb1<=0; ma2<=0; mb2<=0;
            ma3<=0; mb3<=0; ma4<=0; mb4<=0;
        end else begin
            case (in_idx)
            3'd0: begin i0 <= d0 + d2; ma0 <= d1; mb0 <= C2;
                         ma1 <= d3; mb1 <= C6; end
            3'd1: begin ma0 <= d0; mb0 <= C1; ma1 <= d3; mb1 <= C7;
                         ma2 <= d2; mb2 <= C5; ma3 <= d1; mb3 <= C3;
                         ma4 <= i0; mb4 <= C4; end
            3'd2: begin i0 <= d0 - d2; end
            3'd3: begin ma0 <= d0; mb0 <= C7; ma1 <= d3; mb1 <= C1;
                         ma2 <= d2; mb2 <= C3; ma3 <= d1; mb3 <= C5; end
            3'd4: begin ma0 <= d0; mb0 <= C7; ma1 <= d3; mb1 <= C1;
                         ma2 <= d2; mb2 <= C3; ma3 <= d1; mb3 <= C5; end
            3'd5: begin ma0 <= d1; mb0 <= C6; ma1 <= d3; mb1 <= C2;
                         ma4 <= i0; mb4 <= C4; end
            default: ;
            endcase
        end
    end

    // stage 1: multiply
    reg signed [31:0] q0, q1, q2, q3, q4;
    always @(posedge clk) begin
        if (rst) begin q0<=0; q1<=0; q2<=0; q3<=0; q4<=0; end
        else begin
            q0 <= ma0 * mb0;
            q1 <= ma1 * mb1;
            q2 <= ma2 * mb2;
            q3 <= ma3 * mb3;
            q4 <= ma4 * mb4;
        end
    end

    // stage 2: pipeline
    reg signed [31:0] m0, m1, m2, m3, m4;
    always @(posedge clk) begin
        if (rst) begin m0<=0; m1<=0; m2<=0; m3<=0; m4<=0; end
        else begin m0<=q0; m1<=q1; m2<=q2; m3<=q3; m4<=q4; end
    end

    reg       d1_v, d2_v, d3_v, d4_v;
    reg [2:0] d1_i, d2_i, d3_i, d4_i;
    always @(posedge clk) begin
        if (rst) begin
            d1_v<=0; d2_v<=0; d3_v<=0; d4_v<=0;
            d1_i<=0; d2_i<=0; d3_i<=0; d4_i<=0;
        end else begin
            d1_v <= in_valid; d1_i <= in_idx;
            d2_v <= d1_v;     d2_i <= d1_i;
            d3_v <= d2_v;     d3_i <= d2_i;
            d4_v <= d3_v;     d4_i <= d3_i;
        end
    end

    // stage 3: butterfly network
    reg signed [31:0] s5, s6, s7, t0, t1, t2, t3, t4, t5, t6, t7, t65, t56;
    always @(posedge clk) begin
        if (rst) begin
            s5<=0;s6<=0;s7<=0;t0<=0;t1<=0;t2<=0;t3<=0;t4<=0;
            t5<=0;t6<=0;t7<=0;t65<=0;t56<=0;
        end else begin
            case (d3_i)
            3'd0: t3 <= m0 + m1;
            3'd1: begin s7 <= m0 + m1; s6 <= m2 + m3; t0 <= m4; end
            3'd2: begin t0 <= t0 + t3; t3 <= t0 - t3; t7 <= s6 + s7; end
            3'd3: t4 <= (m0 - m1) + (m2 - m3);
            3'd4: begin t0 <= m0 - m1; s5 <= m2 - m3; end
            3'd5: begin t3 <= m0 - m1; t4 <= m4; t5 <= t0 - s5; t6 <= s7 - s6; end
            3'd6: begin t1 <= t4 + t3; t2 <= t4 - t3;
                         t65 <= t6 - t5; t56 <= t5 + t6; end
            default: begin s5 <= (t65 * 181) / 256; s6 <= (t56 * 181) / 256; end
            endcase
        end
    end

    // stage 4: result registers
    reg signed [31:0] bo [0:7];
    reg signed [31:0] bo_tmp;
    always @(posedge clk) begin
        if (rst) begin
            bo[0]<=0;bo[1]<=0;bo[2]<=0;bo[3]<=0;
            bo[4]<=0;bo[5]<=0;bo[6]<=0;bo[7]<=0;bo_tmp<=0;
        end else if (d4_v) begin
            if (d4_i == 3'd3) begin
                bo[0] <= (t0 + t7) >>> SH;
                bo_tmp <= (t0 - t7) >>> SH;
                bo[3] <= (t3 + t4) >>> SH;
                bo[4] <= (t3 - t4) >>> SH;
            end
            if (d4_i == 3'd6)
                bo[7] <= bo_tmp;
            if (d4_i == 3'd7) begin
                bo[2] <= (t2 + s5) >>> SH;
                bo[5] <= (t2 - s5) >>> SH;
                bo[1] <= (t1 + s6) >>> SH;
                bo[6] <= (t1 - s6) >>> SH;
            end
        end
    end

    // stage 5: output streamer (7-cycle delayed valid, 8:1 mux)
    reg [7:0] vld_q;
    always @(posedge clk) begin
        if (rst)      vld_q <= 8'b0;
        else if (img_start) vld_q <= 8'b0;
        else          vld_q <= {vld_q[6:0], d4_v};
    end

    always @(posedge clk) begin
        if (rst)            out_cnt <= 6'd0;
        else if (img_start) out_cnt <= 6'd0;
        else if (vld_q[6])  out_cnt <= out_cnt + 6'd1;
    end

    always @* begin
        out_valid = vld_q[6];
        case (out_cnt[2:0])
        3'd0: out_data = bo[0];
        3'd1: out_data = bo[1];
        3'd2: out_data = bo[2];
        3'd3: out_data = bo[3];
        3'd4: out_data = bo[4];
        3'd5: out_data = bo[5];
        3'd6: out_data = bo[6];
        default: out_data = bo[7];
        endcase
    end
endmodule

//-----------------------------------------------------------------
// Optimized 2-D IDCT (jpeg_idct_2d)
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

    wire        rx_v;
    wire [31:0] rx_d;
    wire [5:0]  rx_cnt;

    idct1d #(.SH(11)) u_x (
        .clk(clk), .rst(rst), .img_start(img_start_i),
        .in_valid(inport_valid_i), .in_idx(inport_idx_i),
        .d0({{16{inport_data0_i[15]}}, inport_data0_i}),
        .d1({{16{inport_data1_i[15]}}, inport_data1_i}),
        .d2({{16{inport_data2_i[15]}}, inport_data2_i}),
        .d3({{16{inport_data3_i[15]}}, inport_data3_i}),
        .out_valid(rx_v), .out_data(rx_d), .out_cnt(rx_cnt)
    );

    wire        tp_accept;
    wire        tp_valid;
    wire [31:0] tpd0, tpd1, tpd2, tpd3;
    wire [2:0]  tp_idx;

    jpeg_idct_transpose u_tp (
        .clk(clk), .rst(rst), .img_start_i(img_start_i), .img_end_i(img_end_i),
        .inport_valid_i(rx_v), .inport_data_i(rx_d),
        .inport_idx_i(rx_cnt), .inport_accept_o(tp_accept),
        .outport_valid_o(tp_valid),
        .outport_data0_o(tpd0), .outport_data1_o(tpd1),
        .outport_data2_o(tpd2), .outport_data3_o(tpd3),
        .outport_idx_o(tp_idx), .outport_ready_i(1'b1)
    );

    wire        ry_v;
    wire [31:0] ry_d;
    wire [5:0]  ry_cnt;

    idct1d #(.SH(15)) u_y (
        .clk(clk), .rst(rst), .img_start(img_start_i),
        .in_valid(tp_valid), .in_idx(tp_idx),
        .d0(tpd0), .d1(tpd1), .d2(tpd2), .d3(tpd3),
        .out_valid(ry_v), .out_data(ry_d), .out_cnt(ry_cnt)
    );

    assign outport_valid_o = ry_v;
    assign outport_data_o  = ry_d;
    assign outport_idx_o   = {ry_cnt[2:0], ry_cnt[5:3]};

endmodule

//-----------------------------------------------------------------
// Transpose: 4 x 32x32 dual-port RAM, double buffered
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data0_o
    ,output [ 31:0]  outport_data1_o
    ,output [ 31:0]  outport_data2_o
    ,output [ 31:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);
    reg block_wr_q, block_rd_q;
    reg [1:0] block_ready_q;
    reg [5:0] rd_idx_q;
    reg [3:0] rd_addr_q;

    wire [4:0] wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

    wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
    wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
    wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
    wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

    wire [31:0] rd0_w, rd1_w, rd2_w, rd3_w;

    ram32x32 r0 (.clk(clk), .addr0(wr_ptr_w), .data0(inport_data_i), .wr0(wr0_w),
                 .addr1({block_rd_q, rd_addr_q}), .q1(rd0_w));
    ram32x32 r1 (.clk(clk), .addr0(wr_ptr_w), .data0(inport_data_i), .wr0(wr1_w),
                 .addr1({block_rd_q, rd_addr_q}), .q1(rd1_w));
    ram32x32 r2 (.clk(clk), .addr0(wr_ptr_w), .data0(inport_data_i), .wr0(wr2_w),
                 .addr1({block_rd_q, rd_addr_q}), .q1(rd2_w));
    ram32x32 r3 (.clk(clk), .addr0(wr_ptr_w), .data0(inport_data_i), .wr0(wr3_w),
                 .addr1({block_rd_q, rd_addr_q}), .q1(rd3_w));

    always @(posedge clk) begin
        if (rst) begin
            block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
        end else if (img_start_i) begin
            block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
        end else begin
            if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o) begin
                block_ready_q[block_wr_q] <= 1'b1;
                block_wr_q <= ~block_wr_q;
            end
            if (outport_valid_o && rd_idx_q[5:0] == 6'd63) begin
                block_ready_q[block_rd_q] <= 1'b0;
                block_rd_q <= ~block_rd_q;
            end
        end
    end

    assign inport_accept_o = ~block_ready_q[block_wr_q];

    localparam STATE_IDLE = 2'd0, STATE_SETUP = 2'd1, STATE_ACTIVE = 2'd2;
    reg [1:0] state_q;
    reg [1:0] next_state_r;

    always @* begin
        next_state_r = state_q;
        case (state_q)
        STATE_IDLE:  if (block_ready_q[block_rd_q] && outport_ready_i) next_state_r = STATE_SETUP;
        STATE_SETUP: next_state_r = STATE_ACTIVE;
        STATE_ACTIVE: if (outport_valid_o && rd_idx_q == 6'd63) next_state_r = STATE_IDLE;
        default: ;
        endcase
        if (img_start_i) next_state_r = STATE_IDLE;
    end

    always @(posedge clk) begin
        if (rst) state_q <= STATE_IDLE;
        else     state_q <= next_state_r;
    end

    always @(posedge clk) begin
        if (rst)                 rd_idx_q <= 6'd0;
        else if (img_start_i)    rd_idx_q <= 6'd0;
        else if (state_q == STATE_ACTIVE) rd_idx_q <= rd_idx_q + 6'd1;
    end

    always @(posedge clk) begin
        if (rst)                rd_addr_q <= 4'b0;
        else if (state_q == STATE_IDLE)  rd_addr_q <= 4'b0;
        else if (state_q == STATE_SETUP) rd_addr_q <= 4'd8;
        else if (state_q == STATE_ACTIVE) begin
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
    end

    assign outport_valid_o = (state_q == STATE_ACTIVE);
    assign outport_idx_o   = rd_idx_q[2:0];
    assign outport_data0_o = rd0_w;
    assign outport_data1_o = rd1_w;
    assign outport_data2_o = rd2_w;
    assign outport_data3_o = rd3_w;
endmodule

//-----------------------------------------------------------------
// 32x32 dual-port RAM (write port 0, read port 1)
//-----------------------------------------------------------------
module ram32x32
(
    input            clk,
    input  [4:0]     addr0,
    input  [31:0]    data0,
    input            wr0,
    input  [4:0]     addr1,
    output reg [31:0] q1
);
    /* verilator lint_off MULTIDRIVEN */
    reg [31:0] ram [0:31];
    /* verilator lint_on MULTIDRIVEN */

    always @(posedge clk) begin
        if (wr0) ram[addr0] <= data0;
    end
    always @(posedge clk) begin
        q1 <= ram[addr1];
    end
endmodule

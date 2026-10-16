//-----------------------------------------------------------------
// 1-D IDCT stage (AAN factorized, 5 multipliers, 8-beat schedule)
//-----------------------------------------------------------------
module idct1d #(parameter SH = 11, parameter IW = 17, parameter OW = 20)
(
    input  clk, rst, img_start,
    input  in_valid,
    input  [2:0] in_idx,
    input  signed [IW-1:0] d0, d1, d2, d3,
    output reg        out_valid,
    output     [31:0] out_data,
    output reg [5:0]  out_cnt
);
    localparam signed [12:0] C1 = 13'sd4017, C2 = 13'sd3784, C3 = 13'sd3406,
                             C4 = 13'sd2896, C5 = 13'sd2276, C6 = 13'sd1567,
                             C7 = 13'sd799;

    // stage 0: operand setup
    reg signed [IW:0]   i0;
    reg signed [IW-1:0] ma0, ma1, ma2, ma3;
    reg signed [IW:0]   ma4;
    reg signed [12:0]   mb0, mb1, mb2, mb3, mb4;
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
    wire signed [IW+12:0] p0 = ma0 * mb0;
    wire signed [IW+12:0] p1 = ma1 * mb1;
    wire signed [IW+12:0] p2 = ma2 * mb2;
    wire signed [IW+12:0] p3 = ma3 * mb3;
    wire signed [IW+13:0] p4 = ma4 * mb4;

    reg signed [31:0] q0, q1, q2, q3, q4;
    always @(posedge clk) begin
        if (rst) begin q0<=0; q1<=0; q2<=0; q3<=0; q4<=0; end
        else begin
            q0 <= p0; q1 <= p1; q2 <= p2; q3 <= p3; q4 <= p4;
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
    wire signed [31:0] e65 = t65 * 181;
    wire signed [31:0] e56 = t56 * 181;
    wire signed [31:0] o65 = (e65 >>> 8) + ((e65[31] & (|e65[7:0])) ? 32'sd1 : 32'sd0);
    wire signed [31:0] o56 = (e56 >>> 8) + ((e56[31] & (|e56[7:0])) ? 32'sd1 : 32'sd0);

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
            default: begin s5 <= o65; s6 <= o56; end
            endcase
        end
    end

    // stage 4: result registers
    wire signed [31:0] v0 = (t0 + t7) >>> SH;
    wire signed [31:0] v1 = (t1 + s6) >>> SH;
    wire signed [31:0] v2 = (t2 + s5) >>> SH;
    wire signed [31:0] v3 = (t3 + t4) >>> SH;
    wire signed [31:0] v4 = (t3 - t4) >>> SH;
    wire signed [31:0] v5 = (t2 - s5) >>> SH;
    wire signed [31:0] v6 = (t1 - s6) >>> SH;
    wire signed [31:0] v7 = (t0 - t7) >>> SH;

    reg signed [OW-1:0] bo0, bo1, bo2, bo3, bo4, bo5, bo6, bo7, bo_tmp;
    always @(posedge clk) begin
        if (rst) begin
            bo0<=0;bo1<=0;bo2<=0;bo3<=0;bo4<=0;bo5<=0;bo6<=0;bo7<=0;bo_tmp<=0;
        end else if (d4_v) begin
            if (d4_i == 3'd3) begin
                bo0 <= v0[OW-1:0];
                bo_tmp <= v7[OW-1:0];
                bo3 <= v3[OW-1:0];
                bo4 <= v4[OW-1:0];
            end
            if (d4_i == 3'd6)
                bo7 <= bo_tmp;
            if (d4_i == 3'd7) begin
                bo2 <= v2[OW-1:0];
                bo5 <= v5[OW-1:0];
                bo1 <= v1[OW-1:0];
                bo6 <= v6[OW-1:0];
            end
        end
    end

    // stage 5: output streamer
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

    reg signed [OW-1:0] bos;
    always @* begin
        case (out_cnt[2:0])
        3'd0: bos = bo0;
        3'd1: bos = bo1;
        3'd2: bos = bo2;
        3'd3: bos = bo3;
        3'd4: bos = bo4;
        3'd5: bos = bo5;
        3'd6: bos = bo6;
        default: bos = bo7;
        endcase
    end

    always @* out_valid = vld_q[6];
    assign out_data = {{(32-OW){bos[OW-1]}}, bos};
endmodule

//-----------------------------------------------------------------
// Optimized 2-D IDCT
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

    idct1d #(.SH(11), .IW(17), .OW(20)) u_x (
        .clk(clk), .rst(rst), .img_start(img_start_i),
        .in_valid(inport_valid_i), .in_idx(inport_idx_i),
        .d0({inport_data0_i[15], inport_data0_i}),
        .d1({inport_data1_i[15], inport_data1_i}),
        .d2({inport_data2_i[15], inport_data2_i}),
        .d3({inport_data3_i[15], inport_data3_i}),
        .out_valid(rx_v), .out_data(rx_d), .out_cnt(rx_cnt)
    );

    wire        tp_accept;
    wire        tp_valid;
    wire [19:0] tpd0, tpd1, tpd2, tpd3;
    wire [2:0]  tp_idx;

    jpeg_idct_transpose u_tp (
        .clk(clk), .rst(rst), .img_start_i(img_start_i), .img_end_i(img_end_i),
        .inport_valid_i(rx_v), .inport_data_i(rx_d[19:0]),
        .inport_idx_i(rx_cnt), .inport_accept_o(tp_accept),
        .outport_valid_o(tp_valid),
        .outport_data0_o(tpd0), .outport_data1_o(tpd1),
        .outport_data2_o(tpd2), .outport_data3_o(tpd3),
        .outport_idx_o(tp_idx), .outport_ready_i(1'b1)
    );

    wire        ry_v;
    wire [31:0] ry_d;
    wire [5:0]  ry_cnt;

    idct1d #(.SH(15), .IW(20), .OW(17)) u_y (
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
// Transpose: 8 banks (one per source row) x 16 words (2 blocks x 8
// columns) x 20 bits. Read: all banks at word {rd_blk, group}; the
// four output lanes pick bank 2b / 2b+1 by beat parity.
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 19:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 19:0]  outport_data0_o
    ,output [ 19:0]  outport_data1_o
    ,output [ 19:0]  outport_data2_o
    ,output [ 19:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);
    reg block_wr_q, block_rd_q;
    reg [1:0] block_ready_q;
    reg [5:0] rd_idx_q;

    localparam STATE_IDLE = 2'd0, STATE_SETUP = 2'd1, STATE_ACTIVE = 2'd2;
    reg [1:0] state_q;
    reg [1:0] next_state_r;

    assign outport_valid_o = (state_q == STATE_ACTIVE);

    // ---------------- write side ----------------
    wire wr_en = inport_valid_i && inport_accept_o;
    wire [2:0] wbank = inport_idx_i[5:3];
    wire [3:0] wword = {block_wr_q, inport_idx_i[2:0]};

    wire we0 = wr_en && (wbank == 3'd0);
    wire we1 = wr_en && (wbank == 3'd1);
    wire we2 = wr_en && (wbank == 3'd2);
    wire we3 = wr_en && (wbank == 3'd3);
    wire we4 = wr_en && (wbank == 3'd4);
    wire we5 = wr_en && (wbank == 3'd5);
    wire we6 = wr_en && (wbank == 3'd6);
    wire we7 = wr_en && (wbank == 3'd7);

    // ---------------- read address ----------------
    // RA = {rd_blk, group}; must lead the beat it serves by one cycle
    reg [3:0] RA;
    wire [5:0] nxt2 = rd_idx_q + 6'd2;
    always @(posedge clk) begin
        if (rst)               RA <= 4'd0;
        else if (img_start_i)  RA <= 4'd0;
        else if (state_q == STATE_IDLE)  RA <= {block_rd_q, 3'd0};
        else if (state_q == STATE_ACTIVE) RA <= {block_rd_q, nxt2[5:3]};
    end

    reg [19:0] rdr0, rdr1, rdr2, rdr3, rdr4, rdr5, rdr6, rdr7;

    /* verilator lint_off MULTIDRIVEN */
    reg [19:0] mem0 [0:15];
    reg [19:0] mem1 [0:15];
    reg [19:0] mem2 [0:15];
    reg [19:0] mem3 [0:15];
    reg [19:0] mem4 [0:15];
    reg [19:0] mem5 [0:15];
    reg [19:0] mem6 [0:15];
    reg [19:0] mem7 [0:15];
    /* verilator lint_on MULTIDRIVEN */

    always @(posedge clk) begin
        if (we0) mem0[wword] <= inport_data_i;
        rdr0 <= mem0[RA];
    end
    always @(posedge clk) begin
        if (we1) mem1[wword] <= inport_data_i;
        rdr1 <= mem1[RA];
    end
    always @(posedge clk) begin
        if (we2) mem2[wword] <= inport_data_i;
        rdr2 <= mem2[RA];
    end
    always @(posedge clk) begin
        if (we3) mem3[wword] <= inport_data_i;
        rdr3 <= mem3[RA];
    end
    always @(posedge clk) begin
        if (we4) mem4[wword] <= inport_data_i;
        rdr4 <= mem4[RA];
    end
    always @(posedge clk) begin
        if (we5) mem5[wword] <= inport_data_i;
        rdr5 <= mem5[RA];
    end
    always @(posedge clk) begin
        if (we6) mem6[wword] <= inport_data_i;
        rdr6 <= mem6[RA];
    end
    always @(posedge clk) begin
        if (we7) mem7[wword] <= inport_data_i;
        rdr7 <= mem7[RA];
    end

    // lane select: p = 0 on beats {0,2,5,7}, 1 on {1,3,4,6}
    wire psel = rd_idx_q[0] ^ rd_idx_q[2];

    assign outport_data0_o = psel ? rdr1 : rdr0;
    assign outport_data1_o = psel ? rdr3 : rdr2;
    assign outport_data2_o = psel ? rdr5 : rdr4;
    assign outport_data3_o = psel ? rdr7 : rdr6;
    assign outport_idx_o   = rd_idx_q[2:0];

    // ---------------- control ----------------
    always @(posedge clk) begin
        if (rst) begin
            block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
        end else if (img_start_i) begin
            block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
        end else begin
            if (wr_en && inport_idx_i == 6'd63) begin
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

endmodule

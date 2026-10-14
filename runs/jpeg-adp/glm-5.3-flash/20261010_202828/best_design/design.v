//-----------------------------------------------------------------
// jpeg_idct_2d — full 8x8 2-D IDCT of the baseline JPEG decoder
// (ultraembedded/core_jpeg fixed-point behavior, exactly reproduced).
//
// Optimized implementation:
//  - 1-D passes (AAN-style factorization, constants x4096): 4 time-shared
//    constant-per-beat multiplier lanes (the original uses 5: the s0
//    product (x0+x4)*C4 is produced by lane 2 loaded at beat 0 and
//    consumed by accumulator case 0; the (x0-x4)*C4 product is produced
//    by lane 2 at beat 5; the beat-4 operand reload of the original is
//    redundant since beat 3 and beat 4 deliver identical data in both
//    passes). Datapath widths narrowed to the provable value ranges so
//    every 32-bit wrap point of the original is preserved exactly.
//  - single product register stage (case m of the accumulator consumes
//    the products loaded at beat m, so the second product register of
//    the original is redundant).
//  - accumulator written as one dedicated register + enable per value
//    (no per-value case muxes).
//  - (x*181)/256 with signed truncation toward zero: the *181 products
//    are computed in accumulator case 6 (from the t5/t6 datapath) and
//    only the sign fixup + arithmetic shift remain in case 7 — no
//    divider, 32-bit wrap on x*181 kept.
//  - Corner turn: 8 row-banks x 2 block-copies of 8-deep shift registers.
//    Writes push in row order; each bank pops exactly once per column
//    and the read tap holds the value across the column's beats,
//    reproducing the original's re-read pattern (even row set at beats
//    0/2/5, odd set at 1/3/4). No address decode, no read muxes beyond
//    a 2:1 copy select and an even/odd lane select.
//-----------------------------------------------------------------

module jpeg_idct_1d #(
    parameter integer OUT_SHIFT = 11,
    parameter integer IN_W      = 16,
    parameter integer A_W       = 18,
    parameter integer PROD_W    = 29,
    parameter integer T_W       = 31,
    parameter integer OUT_W     = 20,
    parameter integer IDX_SWAP  = 0
)(
    input  wire             clk,
    input  wire             rst,
    input  wire             img_start_i,
    input  wire             inport_valid_i,
    input  wire [IN_W-1:0]  inport_data0_i,
    input  wire [IN_W-1:0]  inport_data1_i,
    input  wire [IN_W-1:0]  inport_data2_i,
    input  wire [IN_W-1:0]  inport_data3_i,
    input  wire [2:0]       inport_idx_i,
    output wire             outport_valid_o,
    output wire [OUT_W-1:0] outport_data_o,
    output wire [5:0]       outport_idx_o
);

    localparam signed [12:0] C1_16 = 13'sd4017; // cos( pi/16) x4096
    localparam signed [12:0] C2_16 = 13'sd3784; // cos(2pi/16) x4096
    localparam signed [12:0] C3_16 = 13'sd3406; // cos(3pi/16) x4096
    localparam signed [12:0] C4_16 = 13'sd2896; // cos(4pi/16) x4096
    localparam signed [12:0] C5_16 = 13'sd2276; // cos(5pi/16) x4096
    localparam signed [12:0] C6_16 = 13'sd1567; // cos(6pi/16) x4096
    localparam signed [12:0] C7_16 = 13'sd799;  // cos(7pi/16) x4096

    wire signed [A_W-1:0] bin0 = $signed(inport_data0_i);
    wire signed [A_W-1:0] bin1 = $signed(inport_data1_i);
    wire signed [A_W-1:0] bin2 = $signed(inport_data2_i);
    wire signed [A_W-1:0] bin3 = $signed(inport_data3_i);

    //-----------------------------------------------------------------
    // Multiplier operand load (8-beat schedule).
    // Input lane groups: beats {0,2,5,7} carry the even coefficient set
    // (x0,x2,x4,x6) on (lane0..lane3); beats {1,3,4,6} the odd set
    // (x1,x3,x5,x7).
    //-----------------------------------------------------------------
    reg signed [A_W-1:0] i0_r;
    reg signed [A_W-1:0] m0a, m1a, m2a, m3a;
    reg signed [12:0]    m0b, m1b, m2b, m3b;

    always @(posedge clk) begin
        if (rst) begin
            i0_r <= {A_W{1'b0}};
            m0a  <= {A_W{1'b0}};
            m1a  <= {A_W{1'b0}};
            m2a  <= {A_W{1'b0}};
            m3a  <= {A_W{1'b0}};
            m0b  <= 13'sd0;
            m1b  <= 13'sd0;
            m2b  <= 13'sd0;
            m3b  <= 13'sd0;
        end else begin
            case (inport_idx_i)
            3'd0: begin
                i0_r <= bin0 + bin2;        // x0+x4
                m0a  <= bin1; m0b <= C2_16; // x2*C2
                m1a  <= bin3; m1b <= C6_16; // x6*C6
                m2a  <= bin0 + bin2; m2b <= C4_16; // (x0+x4)*C4 = s0
            end
            3'd1: begin
                m0a  <= bin0; m0b <= C1_16; // x1*C1
                m1a  <= bin3; m1b <= C7_16; // x7*C7
                m2a  <= bin2; m2b <= C5_16; // x5*C5
                m3a  <= bin1; m3b <= C3_16; // x3*C3
            end
            3'd2: begin
                i0_r <= bin0 - bin2;        // x0-x4
            end
            3'd3: begin
                m0a  <= bin0; m0b <= C7_16; // x1*C7
                m1a  <= bin3; m1b <= C1_16; // x7*C1
                m2a  <= bin2; m2b <= C3_16; // x5*C3
                m3a  <= bin1; m3b <= C5_16; // x3*C5
            end
            // 3'd4: no reload — beat-3 products persist (identical data)
            3'd5: begin
                m0a  <= bin1; m0b <= C6_16; // x2*C6
                m1a  <= bin3; m1b <= C2_16; // x6*C2
                m2a  <= i0_r; m2b <= C4_16; // (x0-x4)*C4 = s1
            end
            default: ;
            endcase
        end
    end

    //-----------------------------------------------------------------
    // Multipliers (single product register stage)
    //-----------------------------------------------------------------
    reg signed [PROD_W-1:0] p0, p1, p2, p3;

    always @(posedge clk) begin
        if (rst) begin
            p0 <= {PROD_W{1'b0}};
            p1 <= {PROD_W{1'b0}};
            p2 <= {PROD_W{1'b0}};
            p3 <= {PROD_W{1'b0}};
        end else begin
            p0 <= m0a * m0b;
            p1 <= m1a * m1b;
            p2 <= m2a * m2b;
            p3 <= m3a * m3b;
        end
    end

    //-----------------------------------------------------------------
    // Pipeline valid/idx (case m of the accumulator fires at beat m+2
    // and consumes the products loaded at beat m; the output stage
    // fires at beat m+3)
    //-----------------------------------------------------------------
    reg       v1, v2, v3;
    reg [2:0] d0, d1, d2;

    always @(posedge clk) begin
        if (rst) begin
            v1 <= 1'b0; v2 <= 1'b0; v3 <= 1'b0;
            d0 <= 3'd0; d1 <= 3'd0; d2 <= 3'd0;
        end else begin
            v1 <= inport_valid_i; d0 <= inport_idx_i;
            v2 <= v1;             d1 <= d0;
            v3 <= v2;             d2 <= d1;
        end
    end

    //-----------------------------------------------------------------
    // Accumulation stage: one dedicated register + enable per value.
    //  r0=s0 r1=s3 r2=s7 r3=s6 r4=t0 r5=t3 r6=t7 r7=t4
    //  r8=s4 r9=s5 r10=s2 r11=s1 r12=t1 r13=t2 r14=t5 r15=t6
    //  r16=t65 r17=t56 r18=s5o r19=s6o
    //-----------------------------------------------------------------
    reg signed [T_W-1:0] r0, r1, r2, r3, r4, r5, r6, r7;
    reg signed [T_W-1:0] r8, r9, r10, r11, r12, r13, r14, r15;
    reg signed [T_W-1:0] r16, r17, r18, r19;

    // (x*181)/256 with signed truncation toward zero, 32-bit wrap on x*181.
    // The *181 products are registered in case 6 (from the t5/t6 datapath
    // being written there); case 7 applies the sign fixup + shift.
    reg signed [31:0] m5, m6;
    wire signed [31:0] div_in_5 = m5 + (m5[31] ? 32'sd255 : 32'sd0);
    wire signed [31:0] div_in_6 = m6 + (m6[31] ? 32'sd255 : 32'sd0);
    wire signed [31:0] div_q_5  = div_in_5 >>> 8;
    wire signed [31:0] div_q_6  = div_in_6 >>> 8;

    always @(posedge clk) begin
        if (rst) begin
            r0  <= {T_W{1'b0}}; r1  <= {T_W{1'b0}};
            r2  <= {T_W{1'b0}}; r3  <= {T_W{1'b0}};
            r4  <= {T_W{1'b0}}; r5  <= {T_W{1'b0}};
            r6  <= {T_W{1'b0}}; r7  <= {T_W{1'b0}};
            r8  <= {T_W{1'b0}}; r9  <= {T_W{1'b0}};
            r10 <= {T_W{1'b0}}; r11 <= {T_W{1'b0}};
            r12 <= {T_W{1'b0}}; r13 <= {T_W{1'b0}};
            r14 <= {T_W{1'b0}}; r15 <= {T_W{1'b0}};
            r16 <= {T_W{1'b0}}; r17 <= {T_W{1'b0}};
            r18 <= {T_W{1'b0}}; r19 <= {T_W{1'b0}};
            m5  <= 32'd0;
            m6  <= 32'd0;
        end else begin
            case (d1)
            3'd0: begin
                r1 <= p0 + p1;  // s3
                r0 <= p2;       // s0
            end
            3'd1: begin
                r2 <= p0 + p1;  // s7
                r3 <= p2 + p3;  // s6
            end
            3'd2: begin
                r4 <= r0 + r1;  // t0 = s0 + s3
                r5 <= r0 - r1;  // t3 = s0 - s3
                r6 <= r3 + r2;  // t7
            end
            3'd3: begin
                r7 <= (p0 - p1) + (p2 - p3); // t4
            end
            3'd4: begin
                r8 <= p0 - p1;  // s4
                r9 <= p2 - p3;  // s5
            end
            3'd5: begin
                r10 <= p0 - p1; // s2
                r11 <= p2;      // s1
                r14 <= r8 - r9; // t5 = s4 - s5
                r15 <= r2 - r3; // t6 = s7 - s6
            end
            3'd6: begin
                r12 <= r11 + r10; // t1 = s1 + s2
                r13 <= r11 - r10; // t2 = s1 - s2
                r16 <= r15 - r14;
                r17 <= r14 + r15;
                m5  <= (r15 - r14) * 32'sd181; // t65*181 (32-bit wrap)
                m6  <= (r14 + r15) * 32'sd181; // t56*181
            end
            default: begin
                r18 <= $signed(div_q_5[24:0]); // 1/sqrt(2)
                r19 <= $signed(div_q_6[24:0]); // 1/sqrt(2)
            end
            endcase
        end
    end

    //-----------------------------------------------------------------
    // Output stage (fires at beat m+3)
    //-----------------------------------------------------------------
    wire signed [31:0] sum0 = r4 + r6;
    wire signed [31:0] sum7 = r4 - r6;
    wire signed [31:0] sum3 = r5 + r7;
    wire signed [31:0] sum4 = r5 - r7;
    wire signed [31:0] sum2 = r13 + r18;
    wire signed [31:0] sum5 = r13 - r18;
    wire signed [31:0] sum1 = r12 + r19;
    wire signed [31:0] sum6 = r12 - r19;

    wire signed [31:0] sh0 = sum0 >>> OUT_SHIFT;
    wire signed [31:0] sh7 = sum7 >>> OUT_SHIFT;
    wire signed [31:0] sh3 = sum3 >>> OUT_SHIFT;
    wire signed [31:0] sh4 = sum4 >>> OUT_SHIFT;
    wire signed [31:0] sh2 = sum2 >>> OUT_SHIFT;
    wire signed [31:0] sh5 = sum5 >>> OUT_SHIFT;
    wire signed [31:0] sh1 = sum1 >>> OUT_SHIFT;
    wire signed [31:0] sh6 = sum6 >>> OUT_SHIFT;

    reg signed [OUT_W-1:0] bo0, bo1, bo2, bo3, bo4, bo5, bo6, bo7;
    reg signed [OUT_W-1:0] bo_tmp;

    always @(posedge clk) begin
        if (rst) begin
            bo0 <= {OUT_W{1'b0}};
            bo1 <= {OUT_W{1'b0}};
            bo2 <= {OUT_W{1'b0}};
            bo3 <= {OUT_W{1'b0}};
            bo4 <= {OUT_W{1'b0}};
            bo5 <= {OUT_W{1'b0}};
            bo6 <= {OUT_W{1'b0}};
            bo7 <= {OUT_W{1'b0}};
            bo_tmp <= {OUT_W{1'b0}};
        end else if (v3) begin
            if (d2 == 3'd3) begin
                bo0    <= sh0[OUT_W-1:0];
                bo_tmp <= sh7[OUT_W-1:0]; // bo7
                bo3    <= sh3[OUT_W-1:0];
                bo4    <= sh4[OUT_W-1:0];
            end
            if (d2 == 3'd6)
                bo7 <= bo_tmp;
            if (d2 == 3'd7) begin
                bo2 <= sh2[OUT_W-1:0];
                bo5 <= sh5[OUT_W-1:0];
                bo1 <= sh1[OUT_W-1:0];
                bo6 <= sh6[OUT_W-1:0];
            end
        end
    end

    //-----------------------------------------------------------------
    // Output streamer
    //-----------------------------------------------------------------
    reg [7:0] valid_q;
    reg [5:0] ptr_q;

    always @(posedge clk) begin
        if (rst || img_start_i)
            valid_q <= 8'd0;
        else
            valid_q <= {valid_q[6:0], v3};
    end

    always @(posedge clk) begin
        if (rst || img_start_i)
            ptr_q <= 6'd0;
        else if (outport_valid_o)
            ptr_q <= ptr_q + 6'd1;
    end

    assign outport_valid_o = valid_q[6];
    assign outport_data_o  = (ptr_q[2:0] == 3'd0) ? bo0 :
                             (ptr_q[2:0] == 3'd1) ? bo1 :
                             (ptr_q[2:0] == 3'd2) ? bo2 :
                             (ptr_q[2:0] == 3'd3) ? bo3 :
                             (ptr_q[2:0] == 3'd4) ? bo4 :
                             (ptr_q[2:0] == 3'd5) ? bo5 :
                             (ptr_q[2:0] == 3'd6) ? bo6 : bo7;

    generate
        if (IDX_SWAP != 0) begin : g_swap
            assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
        end else begin : g_noswap
            assign outport_idx_o = ptr_q;
        end
    endgenerate

endmodule

//-----------------------------------------------------------------
// Corner turn: 8 row-banks x 2 block-copies, each an 8-deep shift reg.
// Write phase: the x-pass streams a block row-major (idx = {row, col});
// bank r receives its 8 values (columns 0..7) sequentially and pushes.
// Read phase: 64 consecutive beats (column c = cnt[5:3], beat j = cnt[2:0]).
// Odd banks (rows 1,3,5,7) pop at j==5, even banks (rows 0,2,4,6) at j==6;
// between pops the read tap holds, reproducing the original delivery
// pattern: even row set at beats {0,2,5}, odd row set at {1,3,4}.
//-----------------------------------------------------------------
module jpeg_idct_tr #(
    parameter integer DW = 20
)(
    input  wire          clk,
    input  wire          rst,
    input  wire          img_start_i,
    input  wire          inport_valid_i,
    input  wire [DW-1:0] inport_data_i,
    input  wire [5:0]    inport_idx_i,
    output wire          outport_valid_o,
    output wire [DW-1:0] outport_data0_o,
    output wire [DW-1:0] outport_data1_o,
    output wire [DW-1:0] outport_data2_o,
    output wire [DW-1:0] outport_data3_o,
    output wire [2:0]    outport_idx_o
);

    reg        wr_copy_q;   // copy being written
    reg        rd_copy_q;   // copy being read
    reg        rd_active_q;
    reg [5:0]  rd_cnt_q;

    wire wr_done_w = inport_valid_i & (inport_idx_i == 6'd63);

    always @(posedge clk) begin
        if (rst || img_start_i) begin
            wr_copy_q   <= 1'b0;
            rd_copy_q   <= 1'b0;
            rd_active_q <= 1'b0;
            rd_cnt_q    <= 6'd0;
        end else if (wr_done_w) begin
            wr_copy_q   <= ~wr_copy_q;
            rd_copy_q   <= wr_copy_q;    // read the block just completed
            rd_active_q <= 1'b1;
            rd_cnt_q    <= 6'd0;
        end else if (rd_active_q) begin
            if (rd_cnt_q == 6'd63)
                rd_active_q <= 1'b0;
            rd_cnt_q <= rd_cnt_q + 6'd1;
        end
    end

    wire [2:0] j_w = rd_cnt_q[2:0];
    wire j_eq5_w = j_w[2] & ~j_w[1] &  j_w[0];
    wire j_eq6_w = j_w[2] &  j_w[1] & ~j_w[0];
    // odd row set on beats {1,3,4}, even row set on {0,2,5}
    wire sel_O_w = (~j_w[2] & j_w[0]) | (j_w[2] & ~j_w[1] & ~j_w[0]);

    wire [7:0] wr_sel_w = inport_valid_i ? (8'h01 << inport_idx_i[5:3]) : 8'h00;

    // [copy][bank][stage]; push at stage 0, read tap at stage 7
    reg [DW-1:0] mem0 [0:63];
    reg [DW-1:0] mem1 [0:63];

    wire wr0 [0:7];
    wire wr1 [0:7];
    wire rd0 [0:7];
    wire rd1 [0:7];

    genvar b;
    generate
        for (b = 0; b < 8; b = b + 1) begin : g_en
            assign wr0[b] = wr_sel_w[b] & ~wr_copy_q;
            assign wr1[b] = wr_sel_w[b] &  wr_copy_q;
            assign rd0[b] = rd_active_q & ~rd_copy_q & (b[0] ? j_eq5_w : j_eq6_w);
            assign rd1[b] = rd_active_q &  rd_copy_q & (b[0] ? j_eq5_w : j_eq6_w);
        end
    endgenerate

    integer bb;
    always @(posedge clk) begin
        for (bb = 0; bb < 8; bb = bb + 1) begin
            if (wr0[bb] | rd0[bb]) begin
                mem0[bb*8+0] <= wr0[bb] ? inport_data_i : {DW{1'b0}};
                mem0[bb*8+1] <= mem0[bb*8+0];
                mem0[bb*8+2] <= mem0[bb*8+1];
                mem0[bb*8+3] <= mem0[bb*8+2];
                mem0[bb*8+4] <= mem0[bb*8+3];
                mem0[bb*8+5] <= mem0[bb*8+4];
                mem0[bb*8+6] <= mem0[bb*8+5];
                mem0[bb*8+7] <= mem0[bb*8+6];
            end
            if (wr1[bb] | rd1[bb]) begin
                mem1[bb*8+0] <= wr1[bb] ? inport_data_i : {DW{1'b0}};
                mem1[bb*8+1] <= mem1[bb*8+0];
                mem1[bb*8+2] <= mem1[bb*8+1];
                mem1[bb*8+3] <= mem1[bb*8+2];
                mem1[bb*8+4] <= mem1[bb*8+3];
                mem1[bb*8+5] <= mem1[bb*8+4];
                mem1[bb*8+6] <= mem1[bb*8+5];
                mem1[bb*8+7] <= mem1[bb*8+6];
            end
        end
    end

    // read taps (oldest entry of each bank of the copy being read)
    wire [DW-1:0] tap [0:7];
    generate
        for (b = 0; b < 8; b = b + 1) begin : g_tap
            assign tap[b] = rd_copy_q ? mem1[b*8+7] : mem0[b*8+7];
        end
    endgenerate

    assign outport_data0_o = sel_O_w ? tap[1] : tap[0];
    assign outport_data1_o = sel_O_w ? tap[3] : tap[2];
    assign outport_data2_o = sel_O_w ? tap[5] : tap[4];
    assign outport_data3_o = sel_O_w ? tap[7] : tap[6];

    assign outport_idx_o   = j_w;
    assign outport_valid_o = rd_active_q;

endmodule

//-----------------------------------------------------------------
// Top level
//-----------------------------------------------------------------
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
    input  wire [2:0]  inport_idx_i,
    output wire        outport_valid_o,
    output wire [31:0] outport_data_o,
    output wire [5:0]  outport_idx_o
);

    wire        x_valid;
    wire [19:0] x_data;
    wire [5:0]  x_idx;

    jpeg_idct_1d #(
        .OUT_SHIFT(11), .IN_W(16), .A_W(18),
        .PROD_W(29), .T_W(31), .OUT_W(20), .IDX_SWAP(0)
    ) u_idct_x (
        .clk            (clk),
        .rst            (rst),
        .img_start_i    (img_start_i),
        .inport_valid_i (inport_valid_i),
        .inport_data0_i (inport_data0_i),
        .inport_data1_i (inport_data1_i),
        .inport_data2_i (inport_data2_i),
        .inport_data3_i (inport_data3_i),
        .inport_idx_i   (inport_idx_i),
        .outport_valid_o(x_valid),
        .outport_data_o (x_data),
        .outport_idx_o  (x_idx)
    );

    wire        t_valid;
    wire [19:0] t_data0, t_data1, t_data2, t_data3;
    wire [2:0]  t_idx;

    jpeg_idct_tr #(.DW(20)) u_transpose (
        .clk            (clk),
        .rst            (rst),
        .img_start_i    (img_start_i),
        .inport_valid_i (x_valid),
        .inport_data_i  (x_data),
        .inport_idx_i   (x_idx),
        .outport_valid_o(t_valid),
        .outport_data0_o(t_data0),
        .outport_data1_o(t_data1),
        .outport_data2_o(t_data2),
        .outport_data3_o(t_data3),
        .outport_idx_o  (t_idx)
    );

    wire        y_valid;
    wire [17:0] y_data;
    wire [5:0]  y_idx;

    jpeg_idct_1d #(
        .OUT_SHIFT(15), .IN_W(20), .A_W(21),
        .PROD_W(32), .T_W(32), .OUT_W(18), .IDX_SWAP(1)
    ) u_idct_y (
        .clk            (clk),
        .rst            (rst),
        .img_start_i    (img_start_i),
        .inport_valid_i (t_valid),
        .inport_data0_i (t_data0),
        .inport_data1_i (t_data1),
        .inport_data2_i (t_data2),
        .inport_data3_i (t_data3),
        .inport_idx_i   (t_idx),
        .outport_valid_o(y_valid),
        .outport_data_o (y_data),
        .outport_idx_o  (y_idx)
    );

    assign outport_valid_o = y_valid;
    assign outport_data_o  = {{14{y_data[17]}}, y_data};
    assign outport_idx_o   = y_idx;

endmodule

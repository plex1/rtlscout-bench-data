//-----------------------------------------------------------------
// jpeg_idct_2d — optimized 8x8 2-D IDCT (row pass -> transpose -> col pass)
// Per-cycle schedule identical to ultraembedded/core_jpeg.
// Optimizations (all value-preserving):
//  * Each cosine constant C1..C7 gets a DEDICATED shift-add (CSD) constant
//    multiplier instead of a shared general multiplier. The 8-beat schedule
//    needs at most one product per constant per beat, so 7 fixed-constant
//    units replace 5 variable multipliers (exact same products).
//  * transpose RAM narrowed to 20 bits (|x-pass out| <= 374,555 provably)
//  * y-pass output narrowed to 17 bits ((v>>>15) of a 32-bit value <= 2^16)
//  * /256 (signed, trunc toward 0) implemented explicitly
//  * dead logic removed
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
    wire signed [19:0] x_data_w;
    wire [5:0]  x_idx_w;

    idct_pass #(.IN_W(16), .I0_W(18), .OUT_SHIFT(11), .OUT_W(20), .COL(0))
    u_idct_x
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
    wire signed [19:0] t_data0_w;
    wire signed [19:0] t_data1_w;
    wire signed [19:0] t_data2_w;
    wire signed [19:0] t_data3_w;
    wire [2:0]  t_idx_w;

    jpeg_idct_transpose u_transpose
    (
         .clk(clk)
        ,.rst(rst)
        ,.img_start_i(img_start_i)
        ,.inport_valid_i(x_valid_w)
        ,.inport_data_i(x_data_w)        ,.inport_idx_i(x_idx_w)
        ,.outport_valid_o(t_valid_w)
        ,.outport_data0_o(t_data0_w)
        ,.outport_data1_o(t_data1_w)
        ,.outport_data2_o(t_data2_w)
        ,.outport_data3_o(t_data3_w)
        ,.outport_idx_o(t_idx_w)
    );

    idct_pass #(.IN_W(20), .I0_W(21), .OUT_SHIFT(15), .OUT_W(17), .COL(1))
    u_idct_y
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
// One 1-D IDCT pass (AAN-style, 8-beat schedule, constant multipliers)
//-----------------------------------------------------------------
module idct_pass
#(
     parameter IN_W      = 16
    ,parameter I0_W      = 18
    ,parameter OUT_SHIFT = 11
    ,parameter OUT_W     = 20
    ,parameter COL       = 0
)
(
     input  wire                  clk
    ,input  wire                  rst
    ,input  wire                  img_start_i
    ,input  wire                  inport_valid_i
    ,input  wire signed [IN_W-1:0] inport_data0_i
    ,input  wire signed [IN_W-1:0] inport_data1_i
    ,input  wire signed [IN_W-1:0] inport_data2_i
    ,input  wire signed [IN_W-1:0] inport_data3_i
    ,input  wire [2:0]            inport_idx_i
    ,output wire                  outport_valid_o
    ,output wire [31:0]           outport_data_o
    ,output wire [5:0]            outport_idx_o
);

    //-----------------------------------------------------------------
    // Stage 0: operand registers (one per constant multiplier)
    //-----------------------------------------------------------------
    reg signed [IN_W-1:0] a1, a2, a3, a5, a6, a7;
    reg signed [I0_W-1:0] a4;
    reg signed [I0_W-1:0] i0;

    always @(posedge clk)
    if (rst)
    begin
        a1 <= 0; a2 <= 0; a3 <= 0; a5 <= 0; a6 <= 0; a7 <= 0;
        a4 <= 0; i0 <= 0;
    end
    else
    begin
        case (inport_idx_i)
        3'd0:
        begin
            i0 <= inport_data0_i + inport_data2_i;
            a2 <= inport_data1_i;   // c2 * C2
            a6 <= inport_data3_i;   // c6 * C6
        end
        3'd1:
        begin
            a1 <= inport_data0_i;   // c0 * C1
            a7 <= inport_data3_i;   // c6 * C7
            a5 <= inport_data2_i;   // c4 * C5
            a3 <= inport_data1_i;   // c2 * C3
            a4 <= i0;               // i0 * C4
        end
        3'd2:
        begin
            i0 <= inport_data0_i - inport_data2_i;
        end
        3'd3, 3'd4:
        begin
            a7 <= inport_data0_i;   // c0 * C7
            a1 <= inport_data3_i;   // c6 * C1
            a3 <= inport_data2_i;   // c4 * C3
            a5 <= inport_data1_i;   // c2 * C5
        end
        3'd5:
        begin
            a6 <= inport_data1_i;   // c2 * C6
            a2 <= inport_data3_i;   // c6 * C2
            a4 <= i0;               // i0 * C4
        end
        default:
            ;
        endcase
    end

    //-----------------------------------------------------------------
    // Constant multipliers (exact, 32-bit wrap-around context)
    //   C1=4017 C2=3784 C3=3406 C4=2896 C5=2276 C6=1567 C7=799
    //-----------------------------------------------------------------
    wire signed [31:0] e1 = {{(32-IN_W){a1[IN_W-1]}}, a1};
    wire signed [31:0] e2 = {{(32-IN_W){a2[IN_W-1]}}, a2};
    wire signed [31:0] e3 = {{(32-IN_W){a3[IN_W-1]}}, a3};
    wire signed [31:0] e5 = {{(32-IN_W){a5[IN_W-1]}}, a5};
    wire signed [31:0] e6 = {{(32-IN_W){a6[IN_W-1]}}, a6};
    wire signed [31:0] e7 = {{(32-IN_W){a7[IN_W-1]}}, a7};
    wire signed [31:0] e4 = {{(32-I0_W){a4[I0_W-1]}}, a4};

    wire signed [31:0] k1 = (e1 << 12) - (e1 << 6) - (e1 << 3) - (e1 << 2) - (e1 << 1) - e1;
    wire signed [31:0] k2 = (e2 << 12) - (e2 << 8) - (e2 << 5) - (e2 << 4) - (e2 << 3);
    wire signed [31:0] k3 = (e3 << 12) - (e3 << 9) - (e3 << 7) - (e3 << 5) - (e3 << 4) - (e3 << 1);
    wire signed [31:0] k4 = (e4 << 11) + (e4 << 9) + (e4 << 8) + (e4 << 6) + (e4 << 4);
    wire signed [31:0] k5 = (e5 << 11) + (e5 << 7) + (e5 << 6) + (e5 << 5) + (e5 << 2);
    wire signed [31:0] k6 = (e6 << 10) + (e6 << 9) + (e6 << 5) - e6;
    wire signed [31:0] k7 = (e7 << 9) + (e7 << 8) + (e7 << 5) - e7;

    reg signed [31:0] k1_q, k2_q, k3_q, k4_q, k5_q, k6_q, k7_q;

    always @(posedge clk)
    if (rst)
    begin
        k1_q <= 0; k2_q <= 0; k3_q <= 0;
        k4_q <= 0; k5_q <= 0; k6_q <= 0; k7_q <= 0;
    end
    else
    begin
        k1_q <= k1; k2_q <= k2; k3_q <= k3;
        k4_q <= k4; k5_q <= k5; k6_q <= k6; k7_q <= k7;
    end

    reg signed [31:0] k1_r, k2_r, k3_r, k4_r, k5_r, k6_r, k7_r;

    always @(posedge clk)
    if (rst)
    begin
        k1_r <= 0; k2_r <= 0; k3_r <= 0;
        k4_r <= 0; k5_r <= 0; k6_r <= 0; k7_r <= 0;
    end
    else
    begin
        k1_r <= k1_q; k2_r <= k2_q; k3_r <= k3_q;
        k4_r <= k4_q; k5_r <= k5_q; k6_r <= k6_q; k7_r <= k7_q;
    end

    //-----------------------------------------------------------------
    // Pipeline valid/idx
    //-----------------------------------------------------------------
    reg        out_stg0_valid_q, out_stg1_valid_q, out_stg2_valid_q, out_stg3_valid_q;
    reg [2:0]  out_stg0_idx_q, out_stg1_idx_q, out_stg2_idx_q, out_stg3_idx_q;

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

    //-----------------------------------------------------------------
    // Stage 3: butterfly accumulation (keyed on out_stg2_idx_q)
    //-----------------------------------------------------------------
    reg signed [31:0] o_s5, o_s6, o_s7;
    reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
    reg signed [31:0] o_t6_5, o_t5_6;

    wire signed [31:0] d17 = k7_r - k1_r;
    wire signed [31:0] d35 = k3_r - k5_r;

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
        o_s5 <= 0; o_s6 <= 0; o_s7 <= 0;
        o_t0 <= 0; o_t1 <= 0; o_t2 <= 0; o_t3 <= 0;
        o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
        o_t6_5 <= 0; o_t5_6 <= 0;
    end
    else
    begin
        case (out_stg2_idx_q)
        3'd0:
        begin
            o_t3 <= k2_r + k6_r;
        end
        3'd1:
        begin
            o_s7 <= k1_r + k7_r;
            o_s6 <= k5_r + k3_r;
            o_t0 <= k4_r;
        end
        3'd2:
        begin
            o_t0 <= o_t0 + o_t3;
            o_t3 <= o_t0 - o_t3;
            o_t7 <= o_s6 + o_s7;
        end
        3'd3:
        begin
            o_t4 <= d17 + d35;
        end
        3'd4:
        begin
            o_t0 <= d17;
            o_s5 <= d35;
        end
        3'd5:
        begin
            o_t3 <= k6_r - k2_r;
            o_t4 <= k4_r;
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

    //-----------------------------------------------------------------
    // Stage 4: output registers (keyed on out_stg3_idx_q)
    //-----------------------------------------------------------------
    reg signed [OUT_W-1:0] block_out[0:7];
    reg signed [OUT_W-1:0] block_out_tmp;

    always @(posedge clk)
    if (rst)
    begin
        block_out[0] <= 0; block_out[1] <= 0;
        block_out[2] <= 0; block_out[3] <= 0;
        block_out[4] <= 0; block_out[5] <= 0;
        block_out[6] <= 0; block_out[7] <= 0;
        block_out_tmp <= 0;
    end
    else if (out_stg3_valid_q)
    begin
        if (out_stg3_idx_q == 3'd3)
        begin
            block_out[0] <= (o_t0 + o_t7) >>> OUT_SHIFT;
            block_out_tmp <= (o_t0 - o_t7) >>> OUT_SHIFT;
            block_out[3] <= (o_t3 + o_t4) >>> OUT_SHIFT;
            block_out[4] <= (o_t3 - o_t4) >>> OUT_SHIFT;
        end

        if (out_stg3_idx_q == 3'd6)
            block_out[7] <= block_out_tmp;

        if (out_stg3_idx_q == 3'd7)
        begin
            block_out[2] <= (o_t2 + o_s5) >>> OUT_SHIFT;
            block_out[5] <= (o_t2 - o_s5) >>> OUT_SHIFT;
            block_out[1] <= (o_t1 + o_s6) >>> OUT_SHIFT;
            block_out[6] <= (o_t1 - o_s6) >>> OUT_SHIFT;
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

    wire signed [OUT_W-1:0] out_sel_w = block_out[ptr_q[2:0]];

    assign outport_valid_o = valid_q[6];
    assign outport_data_o  = {{(32-OUT_W){out_sel_w[OUT_W-1]}}, out_sel_w};
    assign outport_idx_o   = COL ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

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
    ,input  wire signed [19:0] inport_data_i
    ,input  wire [5:0]  inport_idx_i
    ,output wire        inport_accept_o
    ,output wire        outport_valid_o
    ,output wire signed [19:0] outport_data0_o
    ,output wire signed [19:0] outport_data1_o
    ,output wire signed [19:0] outport_data2_o
    ,output wire signed [19:0] outport_data3_o
    ,output wire [2:0]  outport_idx_o
);

    reg        block_wr_q;
    reg        block_rd_q;
    reg [1:0]  block_ready_q;
    reg [5:0]  rd_idx_q;
    reg [3:0]  rd_addr_q;

    assign inport_accept_o = ~block_ready_q[block_wr_q];

    wire wr_en_w = inport_valid_i && inport_accept_o;

    wire signed [19:0] rd_data_w [0:3];

    genvar gb;
    generate
    for (gb = 0; gb < 4; gb = gb + 1)
    begin : g_bank
        reg signed [19:0] mem0 [0:15];
        reg signed [19:0] mem1 [0:15];

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

        reg signed [19:0] rd0_q, rd1_q;

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

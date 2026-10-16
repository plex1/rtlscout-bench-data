//-----------------------------------------------------------------
// jpeg_idct_2d -- optimized full 8x8 2-D IDCT (core_jpeg compatible)
// Row pass (OUT_SHIFT=11) -> double-buffered corner turn -> column
// pass (OUT_SHIFT=15).  Bit-exact with the original core_jpeg RTL
// including 32-bit wrap-around arithmetic semantics.
//
// Optimizations vs. original:
//  * 5 generic multipliers per pass -> 7 shared constant shift-add
//    multipliers (one per cosine constant), temporally shared over
//    the 8-beat schedule; all arithmetic mod 2^32 as in original.
//    s1 (capture) -> s2 (hold) product pipeline preserves the exact
//    combine timing (same-unit captures >=2 beats apart).
//  * (x*181)/256 as shift-add + trunc-toward-zero adjust.
//  * Transpose storage 20 bits (x-pass outputs bounded by 414k, fit
//    20-bit signed), y block_out 17 bits (>>>15 result = sext of 17
//    bits).  o_t6_5/o_t5_6 registers folded into the scaling stage.
//  * No sync resets on datapath registers (all consumed values are
//    freshly computed within a row; outputs gated by reset valid_q).
//-----------------------------------------------------------------

module jpeg_idct_1d
#(
     parameter integer AW        = 17
    ,parameter integer OUT_SHIFT = 11
    ,parameter integer OUT_W     = 20
    ,parameter integer SWAP_IDX  = 0
)
(
     input                       clk
    ,input                       rst
    ,input                       img_start_i
    ,input                       inport_valid_i
    ,input  signed [AW-1:0]      in0_i
    ,input  signed [AW-1:0]      in1_i
    ,input  signed [AW-1:0]      in2_i
    ,input  signed [AW-1:0]      in3_i
    ,input  [  2:0]              inport_idx_i
    ,output                      outport_valid_o
    ,output [ 31:0]              outport_data_o
    ,output [  5:0]              outport_idx_o
);

//-----------------------------------------------------------------
// Lane operands, sign-extended to 32 bits (exact)
//-----------------------------------------------------------------
wire signed [31:0] e0 = {{(32-AW){in0_i[AW-1]}}, in0_i};
wire signed [31:0] e1 = {{(32-AW){in1_i[AW-1]}}, in1_i};
wire signed [31:0] e2 = {{(32-AW){in2_i[AW-1]}}, in2_i};
wire signed [31:0] e3 = {{(32-AW){in3_i[AW-1]}}, in3_i};

// beats 1,3 use group A units; beats 0,5 use group B units
wire selA = inport_idx_i[1];
wire selB = inport_idx_i[2];
wire enA  = (inport_idx_i == 3'd1) || (inport_idx_i == 3'd3);
wire enB  = (inport_idx_i == 3'd0) || (inport_idx_i == 3'd5);

//-----------------------------------------------------------------
// Constant multiply operands (7 units, temporally shared)
//-----------------------------------------------------------------
wire signed [31:0] u1_a = selA ? e3 : e0;               // C1: c1@1, c7@3
wire signed [31:0] u2_a = selB ? e3 : e1;               // C2: c2@0, c6@5
wire signed [31:0] u3_a = selA ? e2 : e1;               // C3: c3@1, c5@3
wire signed [31:0] u4_a = selB ? (e0 - e2) : (e0 + e2); // C4: c0+c4@0, c0-c4@5
wire signed [31:0] u5_a = selA ? e1 : e2;               // C5: c5@1, c3@3
wire signed [31:0] u6_a = selB ? e1 : e3;               // C6: c6@0, c2@5
wire signed [31:0] u7_a = selA ? e0 : e3;               // C7: c7@1, c1@3

// C1=4017 = 4096-64-16+1
wire signed [31:0] u1_p = ((u1_a <<< 12) - (u1_a <<< 6)) - ((u1_a <<< 4) - u1_a);
// C2=3784 = 4096-256-64+8
wire signed [31:0] u2_p = ((u2_a <<< 12) - (u2_a <<< 8)) - ((u2_a <<< 6) - (u2_a <<< 3));
// C3=3406 = 4096-512-128-32-16-2
wire signed [31:0] u3_p = ((u3_a <<< 12) - (u3_a <<< 9))
                         - (((u3_a <<< 7) + (u3_a <<< 5))
                         + ((u3_a <<< 4) + (u3_a <<< 1)));
// C4=2896 = 4096-1024-128-32-16
wire signed [31:0] u4_p = ((u4_a <<< 12) - (u4_a <<< 10))
                         - (((u4_a <<< 7) + (u4_a <<< 5)) + (u4_a <<< 4));
// C5=2276 = 2048+256-32+4
wire signed [31:0] u5_p = ((u5_a <<< 11) + (u5_a <<< 8)) - ((u5_a <<< 5) - (u5_a <<< 2));
// C6=1567 = 2048-512+32-1
wire signed [31:0] u6_p = ((u6_a <<< 11) - (u6_a <<< 9)) + ((u6_a <<< 5) - u6_a);
// C7=799 = 512+256+32-1
wire signed [31:0] u7_p = ((u7_a <<< 9) + (u7_a <<< 8)) + ((u7_a <<< 5) - u7_a);

//-----------------------------------------------------------------
// Product pipeline: s1 (capture) -> s2 (hold)
//-----------------------------------------------------------------
reg signed [31:0] s1_1, s1_2, s1_3, s1_4, s1_5, s1_6, s1_7;
reg signed [31:0] s2_1, s2_2, s2_3, s2_4, s2_5, s2_6, s2_7;

always @ (posedge clk)
begin
    if (enA)
    begin
        s1_1 <= u1_p;
        s1_3 <= u3_p;
        s1_5 <= u5_p;
        s1_7 <= u7_p;
    end
    if (enB)
    begin
        s1_2 <= u2_p;
        s1_4 <= u4_p;
        s1_6 <= u6_p;
    end
end

always @ (posedge clk)
begin
    s2_1 <= s1_1;
    s2_2 <= s1_2;
    s2_3 <= s1_3;
    s2_4 <= s1_4;
    s2_5 <= s1_5;
    s2_6 <= s1_6;
    s2_7 <= s1_7;
end

//-----------------------------------------------------------------
// valid / idx pipelines (fixed latency, no backpressure)
//-----------------------------------------------------------------
reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;
reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
    out_stg3_valid_q <= 1'b0;
    out_stg3_idx_q   <= 3'b0;
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
// Combine stage -- every register has a single load point
//-----------------------------------------------------------------
reg signed [31:0] o_s0;    // idx0: s0 = (c0+c4)*C4
reg signed [31:0] o_s3;    // idx0: s3 = c2*C2 + c6*C6
reg signed [31:0] o_s7;    // idx1: c1*C1 + c7*C7
reg signed [31:0] o_s6a;   // idx1: c5*C5 + c3*C3
reg signed [31:0] o_t0r;   // idx2: t0 = s0 + s3
reg signed [31:0] o_t3r;   // idx2: t3 = s0 - s3
reg signed [31:0] o_t7;    // idx2: t7 = s6a + s7
reg signed [31:0] o_t4r;   // idx3: t4
reg signed [31:0] o_s4r;   // idx4: s4
reg signed [31:0] o_s5a;   // idx4: s5
reg signed [31:0] o_s2r;   // idx5: s2
reg signed [31:0] o_s1r;   // idx5: s1
reg signed [31:0] o_d5;    // idx5: t6 - t5
reg signed [31:0] o_d6;    // idx5: t5 + t6
reg signed [31:0] o_t1;    // idx6: t1 = s1 + s2
reg signed [31:0] o_t2;    // idx6: t2 = s1 - s2
reg signed [31:0] o_p5;    // idx6: (t6-t5)*181
reg signed [31:0] o_p6;    // idx6: (t5+t6)*181
reg signed [31:0] o_s5b;   // idx7: o_p5 / 256
reg signed [31:0] o_s6b;   // idx7: o_p6 / 256

//-----------------------------------------------------------------
// 1/sqrt(2) odd scaling: (x * 181) / 256, 32-bit wrap, signed
// division truncating toward zero.  Split across two cycles:
// idx6 registers (x*181), idx7 registers the /256 result.
//-----------------------------------------------------------------
wire signed [31:0] p5_w = ((o_d5 <<< 7) + (o_d5 <<< 5))
                        + (((o_d5 <<< 4) + (o_d5 <<< 2)) + o_d5);
wire signed [31:0] p6_w = ((o_d6 <<< 7) + (o_d6 <<< 5))
                        + (((o_d6 <<< 4) + (o_d6 <<< 2)) + o_d6);
wire signed [31:0] q5_w = (o_p5 >>> 8) +
                          ((o_p5[31] && (o_p5[7:0] != 8'd0)) ? 32'sd1 : 32'sd0);
wire signed [31:0] q6_w = (o_p6 >>> 8) +
                          ((o_p6[31] && (o_p6[7:0] != 8'd0)) ? 32'sd1 : 32'sd0);

always @ (posedge clk)
begin
    if (out_stg2_idx_q == 3'd0)
    begin
        o_s3 <= s2_2 + s2_6;
        o_s0 <= s2_4;
    end
    if (out_stg2_idx_q == 3'd1)
    begin
        o_s7  <= s2_1 + s2_7;
        o_s6a <= s2_5 + s2_3;
    end
    if (out_stg2_idx_q == 3'd2)
    begin
        o_t0r <= o_s0 + o_s3;
        o_t3r <= o_s0 - o_s3;
        o_t7  <= o_s6a + o_s7;
    end
    if (out_stg2_idx_q == 3'd3)
        o_t4r <= (s2_7 - s2_1) + (s2_3 - s2_5);
    if (out_stg2_idx_q == 3'd4)
    begin
        o_s4r <= s2_7 - s2_1;
        o_s5a <= s2_3 - s2_5;
    end
    if (out_stg2_idx_q == 3'd5)
    begin
        o_s2r <= s2_6 - s2_2;
        o_s1r <= s2_4;
        o_d5  <= (o_s7 - o_s6a) - (o_s4r - o_s5a);
        o_d6  <= (o_s7 - o_s6a) + (o_s4r - o_s5a);
    end
    if (out_stg2_idx_q == 3'd6)
    begin
        o_t1 <= o_s1r + o_s2r;
        o_t2 <= o_s1r - o_s2r;
        o_p5 <= p5_w;
        o_p6 <= p6_w;
    end
    if (out_stg2_idx_q == 3'd7)
    begin
        o_s5b <= q5_w;
        o_s6b <= q6_w;
    end
end

//-----------------------------------------------------------------
// Block output (exact original timing)
//-----------------------------------------------------------------
wire signed [31:0] sh0_w = (o_t0r + o_t7) >>> OUT_SHIFT;
wire signed [31:0] sh1_w = (o_t0r - o_t7) >>> OUT_SHIFT;
wire signed [31:0] sh2_w = (o_t3r + o_t4r) >>> OUT_SHIFT;
wire signed [31:0] sh3_w = (o_t3r - o_t4r) >>> OUT_SHIFT;
wire signed [31:0] sh4_w = (o_t2 + o_s5b) >>> OUT_SHIFT;
wire signed [31:0] sh5_w = (o_t2 - o_s5b) >>> OUT_SHIFT;
wire signed [31:0] sh6_w = (o_t1 + o_s6b) >>> OUT_SHIFT;
wire signed [31:0] sh7_w = (o_t1 - o_s6b) >>> OUT_SHIFT;

reg signed [OUT_W-1:0] block_out0;
reg signed [OUT_W-1:0] block_out1;
reg signed [OUT_W-1:0] block_out2;
reg signed [OUT_W-1:0] block_out3;
reg signed [OUT_W-1:0] block_out4;
reg signed [OUT_W-1:0] block_out5;
reg signed [OUT_W-1:0] block_out6;
reg signed [OUT_W-1:0] block_out7;
reg signed [OUT_W-1:0] block_out_tmp;

always @ (posedge clk)
if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out0   <= sh0_w[OUT_W-1:0];
        block_out_tmp <= sh1_w[OUT_W-1:0]; // block_out7
        block_out3   <= sh2_w[OUT_W-1:0];
        block_out4   <= sh3_w[OUT_W-1:0];
    end

    if (out_stg3_idx_q == 3'd6)
        block_out7 <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out2 <= sh4_w[OUT_W-1:0];
        block_out5 <= sh5_w[OUT_W-1:0];
        block_out1 <= sh6_w[OUT_W-1:0];
        block_out6 <= sh7_w[OUT_W-1:0];
    end
end

//-----------------------------------------------------------------
// Output streamer
//-----------------------------------------------------------------
reg [7:0] valid_q;

always @ (posedge clk)
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];

reg signed [OUT_W-1:0] rd_val;
always @ *
begin
    case (ptr_q[2:0])
    3'd0: rd_val = block_out0;
    3'd1: rd_val = block_out1;
    3'd2: rd_val = block_out2;
    3'd3: rd_val = block_out3;
    3'd4: rd_val = block_out4;
    3'd5: rd_val = block_out5;
    3'd6: rd_val = block_out6;
    default: rd_val = block_out7;
    endcase
end

assign outport_data_o = {{(32-OUT_W){rd_val[OUT_W-1]}}, rd_val};

generate
if (SWAP_IDX != 0)
    assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
else
    assign outport_idx_o = ptr_q;
endgenerate

endmodule

//-----------------------------------------------------------------
// Double-buffered corner-turn, DW-bit words, 4 read lanes
//-----------------------------------------------------------------
module jpeg_idct_transpose
#(
     parameter integer DW = 20
)
(
     input                  clk
    ,input                  rst
    ,input                  img_start_i
    ,input                  inport_valid_i
    ,input  [DW-1:0]        inport_data_i
    ,input  [  5:0]         inport_idx_i
    ,input                  outport_ready_i
    ,output                 inport_accept_o
    ,output                 outport_valid_o
    ,output [DW-1:0]        outport_data0_o
    ,output [DW-1:0]        outport_data1_o
    ,output [DW-1:0]        outport_data2_o
    ,output [DW-1:0]        outport_data3_o
    ,output [  2:0]         outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [1:0]   block_ready_q;

localparam STATE_W           = 2;
localparam STATE_IDLE        = 2'd0;
localparam STATE_SETUP       = 2'd1;
localparam STATE_ACTIVE      = 2'd2;

reg [STATE_W-1:0] state_q;

assign inport_accept_o = ~block_ready_q[block_wr_q];

wire wr_en_w = inport_valid_i && inport_accept_o;

wire [3:0] wr_bank_oh;
assign wr_bank_oh[0] = (inport_idx_i[5:4] == 2'd0);
assign wr_bank_oh[1] = (inport_idx_i[5:4] == 2'd1);
assign wr_bank_oh[2] = (inport_idx_i[5:4] == 2'd2);
assign wr_bank_oh[3] = (inport_idx_i[5:4] == 2'd3);

wire [1:0] wr_pb  = {block_wr_q, inport_idx_i[3]};
wire [2:0] wr_col = inport_idx_i[2:0];

wire [3:0] wr_pb_oh;
assign wr_pb_oh[0] = (wr_pb == 2'd0);
assign wr_pb_oh[1] = (wr_pb == 2'd1);
assign wr_pb_oh[2] = (wr_pb == 2'd2);
assign wr_pb_oh[3] = (wr_pb == 2'd3);

wire [DW*4-1:0] rdq_flat;

genvar gb, gp;
generate
for (gb = 0; gb < 4; gb = gb + 1) begin : g_bank
    reg [3:0] addr_q;
    always @ (posedge clk )
    if (rst)
        addr_q <= (gb == 0) ? 4'd0 : (gb == 1) ? 4'd1 : (gb == 2) ? 4'd2 : 4'd3;
    else if (img_start_i)
        addr_q <= 4'd0;
    else if (state_q == STATE_IDLE)
        addr_q <= 4'd0;
    else if (state_q == STATE_SETUP)
        addr_q <= 4'd8;
    else if (state_q == STATE_ACTIVE)
    begin
        case (rd_idx_q[2:0])
        3'd0: addr_q <= addr_q - 4'd8;
        3'd1: addr_q <= addr_q + 4'd8;
        3'd2: ;
        3'd3: addr_q <= addr_q - 4'd8;
        3'd4: addr_q <= addr_q + 4'd8;
        3'd5: addr_q <= addr_q - 4'd8;
        3'd6: addr_q <= addr_q + 4'd1;
        3'd7: addr_q <= addr_q + 4'd8;
        endcase
    end

    wire [1:0] rd_pb  = {block_rd_q, addr_q[3]};
    wire [2:0] rd_col = addr_q[2:0];

    wire [DW*4-1:0] rdcol;
    for (gp = 0; gp < 4; gp = gp + 1) begin : g_pb
        reg [DW-1:0] mem [0:7];
        always @ (posedge clk)
        begin
            if (wr_en_w && wr_bank_oh[gb] && wr_pb_oh[gp])
                mem[wr_col] <= inport_data_i;
        end
        wire [DW-1:0] rd_w;
        assign rd_w = rd_col[2] ?
                      (rd_col[1] ? (rd_col[0] ? mem[7] : mem[6])
                                 : (rd_col[0] ? mem[5] : mem[4]))
                    : (rd_col[1] ? (rd_col[0] ? mem[3] : mem[2])
                                 : (rd_col[0] ? mem[1] : mem[0]));
        assign rdcol[(gp*DW) +: DW] = rd_w;
    end

    wire [DW-1:0] m0 = rdcol[0*DW +: DW];
    wire [DW-1:0] m1 = rdcol[1*DW +: DW];
    wire [DW-1:0] m2 = rdcol[2*DW +: DW];
    wire [DW-1:0] m3 = rdcol[3*DW +: DW];
    wire [DW-1:0] p_sel = rd_pb[0] ? (rd_pb[1] ? m3 : m1)
                                   : (rd_pb[1] ? m2 : m0);
    reg  [DW-1:0] rd_q;
    always @ (posedge clk)
        rd_q <= p_sel;
    assign rdq_flat[(gb*DW) +: DW] = rd_q;
end
endgenerate

assign outport_data0_o = rdq_flat[0*DW +: DW];
assign outport_data1_o = rdq_flat[1*DW +: DW];
assign outport_data2_o = rdq_flat[2*DW +: DW];
assign outport_data3_o = rdq_flat[3*DW +: DW];

//-----------------------------------------------------------------
// Input buffer status / FSM (verbatim behavior)
//-----------------------------------------------------------------
always @ (posedge clk )
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

always @ (posedge clk )
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @ (posedge clk )
if (rst)
    rd_idx_q <= 6'b0;
else if (img_start_i)
    rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE)
    rd_idx_q <= rd_idx_q + 6'd1;

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

endmodule

//-----------------------------------------------------------------
// Top: x pass -> transpose -> y pass
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

wire signed [16:0] x_in0_w = {inport_data0_i[15], inport_data0_i};
wire signed [16:0] x_in1_w = {inport_data1_i[15], inport_data1_i};
wire signed [16:0] x_in2_w = {inport_data2_i[15], inport_data2_i};
wire signed [16:0] x_in3_w = {inport_data3_i[15], inport_data3_i};

wire        x_valid_w;
wire [31:0] x_data_w;
wire [ 5:0] x_idx_w;

jpeg_idct_1d
#(
     .AW(17)
    ,.OUT_SHIFT(11)
    ,.OUT_W(20)
    ,.SWAP_IDX(0)
)
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.in0_i(x_in0_w)
    ,.in1_i(x_in1_w)
    ,.in2_i(x_in2_w)
    ,.in3_i(x_in3_w)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

wire          transpose_valid_w;
wire [ 19:0]  transpose_data0_w;
wire [ 19:0]  transpose_data1_w;
wire [ 19:0]  transpose_data2_w;
wire [ 19:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose
#(
     .DW(20)
)
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w[19:0])
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o()
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
    ,.outport_ready_i(1'b1)
);

wire signed [20:0] y_in0_w = {transpose_data0_w[19], transpose_data0_w};
wire signed [20:0] y_in1_w = {transpose_data1_w[19], transpose_data1_w};
wire signed [20:0] y_in2_w = {transpose_data2_w[19], transpose_data2_w};
wire signed [20:0] y_in3_w = {transpose_data3_w[19], transpose_data3_w};

jpeg_idct_1d
#(
     .AW(21)
    ,.OUT_SHIFT(15)
    ,.OUT_W(17)
    ,.SWAP_IDX(1)
)
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(transpose_valid_w)
    ,.in0_i(y_in0_w)
    ,.in1_i(y_in1_w)
    ,.in2_i(y_in2_w)
    ,.in3_i(y_in3_w)
    ,.inport_idx_i(transpose_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

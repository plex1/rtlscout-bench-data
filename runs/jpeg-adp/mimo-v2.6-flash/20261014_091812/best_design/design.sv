//-----------------------------------------------------------------
// jpeg_idct_2d — optimized, bit-exact re-implementation of the
// core_jpeg (ultraembedded, Apache 2.0) 8x8 2-D IDCT composition:
//   row pass (OUT_SHIFT=11) -> double-buffered transpose -> col pass (15)
//
// Bit-exactness notes:
//  * Arithmetic that can wrap in the original (32-bit wrap-around) is
//    kept in a 32-bit context.  Registers are narrowed ONLY where the
//    original values are provably inside range for ANY input:
//      - row pass: 16-bit inputs => intermediates <= 688,217,344 < 2^30
//        (except the deliberate (x*181) 32-bit wrap); products <= 1.9e8.
//      - row outputs after >>>11 fit in 20 signed bits (max 346,424),
//        so the transpose stores 20 bits (sign-extended on read).
//      - col pass lanes are those row outputs => mul lanes fit 20 bits
//        (i0 fits 21); col pass sums CAN wrap => kept 32-bit; col
//        outputs after >>>15 fit 18 signed bits.
//  * (x*181)/256: multiply kept at 32-bit (wraps), division replaced
//    by an exact truncation-toward-zero shift network.
//  * Pipeline structure, latency, FSMs and streaming schedule are
//    identical to the original.
//-----------------------------------------------------------------

module jpeg_idct_1d
#(
     parameter LANE_W    = 17
    ,parameter MULQ_W    = 29
    ,parameter O_W       = 31
    ,parameter OUT_W     = 20
    ,parameter OUT_SHIFT = 11
    ,parameter IDX_SWAP  = 1'b0
)
(
     input             clk
    ,input             rst
    ,input             img_start_i
    ,input             inport_valid_i
    ,input      [31:0] inport_data0_i
    ,input      [31:0] inport_data1_i
    ,input      [31:0] inport_data2_i
    ,input      [31:0] inport_data3_i
    ,input      [ 2:0] inport_idx_i
    ,output            outport_valid_o
    ,output     [31:0] outport_data_o
    ,output     [ 5:0] outport_idx_o
);
/* verilator lint_off WIDTH */

localparam [12:0] C1_16 = 13'd4017; // cos( pi/16) x4096
localparam [12:0] C2_16 = 13'd3784; // cos(2pi/16) x4096
localparam [12:0] C3_16 = 13'd3406; // cos(3pi/16) x4096
localparam [12:0] C4_16 = 13'd2896; // cos(4pi/16) x4096
localparam [12:0] C5_16 = 13'd2276; // cos(5pi/16) x4096
localparam [12:0] C6_16 = 13'd1567; // cos(6pi/16) x4096
localparam [12:0] C7_16 = 13'd799;  // cos(7pi/16) x4096

wire signed [LANE_W-1:0] block_in_0_1 = inport_data0_i[LANE_W-1:0];
wire signed [LANE_W-1:0] block_in_2_3 = inport_data1_i[LANE_W-1:0];
wire signed [LANE_W-1:0] block_in_4_5 = inport_data2_i[LANE_W-1:0];
wire signed [LANE_W-1:0] block_in_6_7 = inport_data3_i[LANE_W-1:0];

//-----------------------------------------------------------------
// IDCT input stage (schedule: even beats {0,2,5,7}, odd {1,3,4,6})
//-----------------------------------------------------------------
reg signed [LANE_W-1:0] i0;
reg signed [LANE_W-1:0] mul0_a;
reg signed [12:0]       mul0_b;
reg signed [LANE_W-1:0] mul1_a;
reg signed [12:0]       mul1_b;
reg signed [LANE_W-1:0] mul2_a;
reg signed [12:0]       mul2_b;
reg signed [LANE_W-1:0] mul3_a;
reg signed [12:0]       mul3_b;
reg signed [LANE_W-1:0] mul4_a;
reg signed [12:0]       mul4_b;

always @ (posedge clk)
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3;
        mul0_b <= C2_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C1_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C7_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C5_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C3_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1;
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5;
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3;
        mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3;
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7;
        mul1_b <= C2_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
end

reg signed [MULQ_W-1:0] mul0_q;
reg signed [MULQ_W-1:0] mul1_q;
reg signed [MULQ_W-1:0] mul2_q;
reg signed [MULQ_W-1:0] mul3_q;
reg signed [MULQ_W-1:0] mul4_q;

always @ (posedge clk)
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [MULQ_W-1:0] mul0;
reg signed [MULQ_W-1:0] mul1;
reg signed [MULQ_W-1:0] mul2;
reg signed [MULQ_W-1:0] mul3;
reg signed [MULQ_W-1:0] mul4;

always @ (posedge clk)
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
end

reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
end
else
begin
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
end

reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk )
if (rst)
begin
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
end
else
begin
    out_stg2_valid_q <= out_stg1_valid_q;
    out_stg2_idx_q   <= out_stg1_idx_q;
end

//-----------------------------------------------------------------
// (x*181)/256 : multiply wraps at 32 bits (as in the original),
// divide truncates toward zero (as Verilog signed "/" does).
//-----------------------------------------------------------------
reg signed [O_W-1:0] o_s5;
reg signed [O_W-1:0] o_s6;
reg signed [O_W-1:0] o_s7;
reg signed [O_W-1:0] o_t0;
reg signed [O_W-1:0] o_t1;
reg signed [O_W-1:0] o_t2;
reg signed [O_W-1:0] o_t3;
reg signed [O_W-1:0] o_t4;
reg signed [O_W-1:0] o_t5;
reg signed [O_W-1:0] o_t6;
reg signed [O_W-1:0] o_t7;
reg signed [O_W-1:0] o_t6_5;
reg signed [O_W-1:0] o_t5_6;

//-----------------------------------------------------------------
// (x*181)/256 : multiply wraps at 32 bits (as in the original),
// divide truncates toward zero (as Verilog signed "/" does).
//-----------------------------------------------------------------
wire signed [31:0] t6_5_ext = o_t6_5;
wire signed [31:0] t5_6_ext = o_t5_6;
wire signed [31:0] prod65   = t6_5_ext * 32'sd181;
wire signed [31:0] prod56   = t5_6_ext * 32'sd181;
wire signed [31:0] div65    = prod65[31] ? ((prod65 + 32'sd255) >>> 8) : (prod65 >>> 8);
wire signed [31:0] div56    = prod56[31] ? ((prod56 + 32'sd255) >>> 8) : (prod56 >>> 8);

always @ (posedge clk)
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
        o_t4 <= mul4;        // s1
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
        o_s5 <= div65; // 1/sqrt(2)
        o_s6 <= div56; // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk )
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

reg signed [OUT_W-1:0] block_out[0:7];
reg signed [OUT_W-1:0] block_out_tmp;

always @ (posedge clk)
if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT); // block_out[7]
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
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

always @ (posedge clk )
if (rst)
    valid_q  <= 8'b0;
else if (img_start_i)
    valid_q  <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;

always @ (posedge clk )
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];

wire signed [OUT_W-1:0] block_out_r = block_out[ptr_q[2:0]];
assign outport_data_o  = {{(32-OUT_W){block_out_r[OUT_W-1]}}, block_out_r};

// row pass: idx = ptr;  col pass: idx = {ptr[2:0], ptr[5:3]}
assign outport_idx_o   = IDX_SWAP ? {ptr_q[2:0], ptr_q[5:3]} : ptr_q;

/* verilator lint_on WIDTH */

endmodule

//-----------------------------------------------------------------
// Double-buffered corner-turn (transpose).  Behaviour identical to
// the original jpeg_idct_transpose; data path narrowed to DW bits
// (row-pass outputs provably fit 20 signed bits).
//-----------------------------------------------------------------
module jpeg_idct_transpose
#(
     parameter DW = 20
)
(
     input             clk
    ,input             rst
    ,input             img_start_i
    ,input             img_end_i
    ,input             inport_valid_i
    ,input      [DW-1:0] inport_data_i
    ,input      [ 5:0] inport_idx_i
    ,input             outport_ready_i
    ,output            inport_accept_o
    ,output            outport_valid_o
    ,output     [DW-1:0] outport_data0_o
    ,output     [DW-1:0] outport_data1_o
    ,output     [DW-1:0] outport_data2_o
    ,output     [DW-1:0] outport_data3_o
    ,output     [ 2:0] outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};
wire [4:0]  rd_ptr_w = {block_rd_q, rd_addr_q};

wire        wr_en_w = inport_valid_i && inport_accept_o;
wire [1:0]  wr_bank_w = inport_idx_i[5:4];

// Four banks x 32 words (bank selected by idx[5:4], word by
// {block_wr, idx[3:0]}); all banks read at {block_rd, rd_addr_q}.
reg [DW-1:0] mem0 [0:31];
reg [DW-1:0] mem1 [0:31];
reg [DW-1:0] mem2 [0:31];
reg [DW-1:0] mem3 [0:31];

reg [DW-1:0] ram_read0_q;
reg [DW-1:0] ram_read1_q;
reg [DW-1:0] ram_read2_q;
reg [DW-1:0] ram_read3_q;

// One-hot read decode: keeps fanout of the read-address nets tiny
// (16 loads per address bit instead of a shared 32:1 mux tree).
wire [31:0] rsel_w;
genvar gi;
generate
    for (gi = 0; gi < 32; gi = gi + 1) begin : g_rsel
        assign rsel_w[gi] = (rd_ptr_w == gi);
    end
endgenerate

reg [DW-1:0] rdc0, rdc1, rdc2, rdc3;
integer ri;
always @* begin
    rdc0 = {DW{1'b0}};
    rdc1 = {DW{1'b0}};
    rdc2 = {DW{1'b0}};
    rdc3 = {DW{1'b0}};
    for (ri = 0; ri < 32; ri = ri + 1) begin
        rdc0 = rdc0 | ({DW{rsel_w[ri]}} & mem0[ri]);
        rdc1 = rdc1 | ({DW{rsel_w[ri]}} & mem1[ri]);
        rdc2 = rdc2 | ({DW{rsel_w[ri]}} & mem2[ri]);
        rdc3 = rdc3 | ({DW{rsel_w[ri]}} & mem3[ri]);
    end
end

always @ (posedge clk)
begin
    if (wr_en_w && wr_bank_w == 2'd0)
        mem0[wr_ptr_w] <= inport_data_i;
    if (wr_en_w && wr_bank_w == 2'd1)
        mem1[wr_ptr_w] <= inport_data_i;
    if (wr_en_w && wr_bank_w == 2'd2)
        mem2[wr_ptr_w] <= inport_data_i;
    if (wr_en_w && wr_bank_w == 2'd3)
        mem3[wr_ptr_w] <= inport_data_i;

    ram_read0_q <= rdc0;
    ram_read1_q <= rdc1;
    ram_read2_q <= rdc2;
    ram_read3_q <= rdc3;
end

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

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

assign inport_accept_o = ~block_ready_q[block_wr_q];

//-----------------------------------------------------------------
// FSM
//-----------------------------------------------------------------
localparam STATE_W           = 2;
localparam STATE_IDLE        = 2'd0;
localparam STATE_SETUP       = 2'd1;
localparam STATE_ACTIVE      = 2'd2;

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

always @ (posedge clk )
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
assign outport_data0_o = ram_read0_q;
assign outport_data1_o = ram_read1_q;
assign outport_data2_o = ram_read2_q;
assign outport_data3_o = ram_read3_q;

endmodule

//-----------------------------------------------------------------
// jpeg_idct_2d — row pass -> transpose -> column pass
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

// ---- row pass inputs: sign-extend 16 -> 32 ------------------------
wire [31:0] xd0 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire [31:0] xd1 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire [31:0] xd2 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire [31:0] xd3 = {{16{inport_data3_i[15]}}, inport_data3_i};

wire        x_valid_w;
wire [31:0] x_data_w;
wire [ 5:0] x_idx_w;
wire        x_accept_w;

jpeg_idct_1d
#(
     .LANE_W(17)
    ,.MULQ_W(29)
    ,.O_W(31)
    ,.OUT_W(20)
    ,.OUT_SHIFT(11)
    ,.IDX_SWAP(1'b0)
)
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(xd0)
    ,.inport_data1_i(xd1)
    ,.inport_data2_i(xd2)
    ,.inport_data3_i(xd3)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

// ---- transpose (20-bit data path) --------------------------------
wire          t_valid_w;
wire [19:0]   t_data0_w;
wire [19:0]   t_data1_w;
wire [19:0]   t_data2_w;
wire [19:0]   t_data3_w;
wire [ 2:0]   t_idx_w;

jpeg_idct_transpose
#(
     .DW(20)
)
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w[19:0])
    ,.inport_idx_i(x_idx_w)
    ,.inport_accept_o(x_accept_w)
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
    ,.outport_ready_i(1'b1)
);

// ---- column pass inputs: sign-extend 20 -> 32 ---------------------
wire [31:0] yd0 = {{12{t_data0_w[19]}}, t_data0_w};
wire [31:0] yd1 = {{12{t_data1_w[19]}}, t_data1_w};
wire [31:0] yd2 = {{12{t_data2_w[19]}}, t_data2_w};
wire [31:0] yd3 = {{12{t_data3_w[19]}}, t_data3_w};

jpeg_idct_1d
#(
     .LANE_W(21)
    ,.MULQ_W(32)
    ,.O_W(32)
    ,.OUT_W(18)
    ,.OUT_SHIFT(15)
    ,.IDX_SWAP(1'b1)
)
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(t_valid_w)
    ,.inport_data0_i(yd0)
    ,.inport_data1_i(yd1)
    ,.inport_data2_i(yd2)
    ,.inport_data3_i(yd3)
    ,.inport_idx_i(t_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

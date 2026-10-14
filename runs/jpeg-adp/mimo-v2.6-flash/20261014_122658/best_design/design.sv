//-----------------------------------------------------------------
// jpeg_idct_2d — full 8x8 2-D IDCT of the baseline JPEG decoder
// (ultraembedded/core_jpeg, Apache 2.0): jpeg_idct_x -> jpeg_idct_transpose
// -> jpeg_idct_y. Re-structured: the two identical 1-D passes share one
// parameterized module. Behavior is bit-exact with the original RTL.
//-----------------------------------------------------------------
module jpeg_idct_pass
#(
     parameter OUT_SHIFT = 11
    ,parameter IDX_CONV  = 0
    ,parameter A_W       = 22   // width of multiplier operand a
)(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data0_i
    ,input  [ 31:0]  inport_data1_i
    ,input  [ 31:0]  inport_data2_i
    ,input  [ 31:0]  inport_data3_i
    ,input  [  2:0]  inport_idx_i
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam [15:0] C1_16 = 4017; // cos( pi/16) x4096
localparam [15:0] C2_16 = 3784; // cos(2pi/16) x4096
localparam [15:0] C3_16 = 3406; // cos(3pi/16) x4096
localparam [15:0] C4_16 = 2896; // cos(4pi/16) x4096
localparam [15:0] C5_16 = 2276; // cos(5pi/16) x4096
localparam [15:0] C6_16 = 1567; // cos(6pi/16) x4096
localparam [15:0] C7_16 = 799;  // cos(7pi/16) x4096

wire signed [A_W-1:0] block_in_0_1 = inport_data0_i[A_W-1:0];
wire signed [A_W-1:0] block_in_2_3 = inport_data1_i[A_W-1:0];
wire signed [A_W-1:0] block_in_4_5 = inport_data2_i[A_W-1:0];
wire signed [A_W-1:0] block_in_6_7 = inport_data3_i[A_W-1:0];

//-----------------------------------------------------------------
// IDCT
//-----------------------------------------------------------------
// Data-path operands are bounded (16-bit inputs / 16-bit inter-pass values),
// so the multiplier operands are narrowed: a is <= 18 bits signed, b is one
// of the cosine constants (<= 4017, 13 bits signed). Products are identical
// to the original 32-bit computation (no wrap in range).
reg signed [A_W:0] i0;
reg signed [A_W-1:0] mul_a [0:4];
reg signed [12:0]    mul_b [0:4];

always @ (posedge clk )
if (rst)
begin
    i0 <= '0;
    mul_a[0] <= '0; mul_b[0] <= 13'b0;
    mul_a[1] <= '0; mul_b[1] <= 13'b0;
    mul_a[2] <= '0; mul_b[2] <= 13'b0;
    mul_a[3] <= '0; mul_b[3] <= 13'b0;
    mul_a[4] <= '0; mul_b[4] <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0      <= block_in_0_1 + block_in_4_5;
        mul_a[0] <= block_in_2_3;
        mul_b[0] <= C2_16;
        mul_a[1] <= block_in_6_7;
        mul_b[1] <= C6_16;
    end
    3'd1:
    begin
        mul_a[0] <= block_in_0_1;
        mul_b[0] <= C1_16;
        mul_a[1] <= block_in_6_7;
        mul_b[1] <= C7_16;
        mul_a[2] <= block_in_4_5;
        mul_b[2] <= C5_16;
        mul_a[3] <= block_in_2_3;
        mul_b[3] <= C3_16;
        mul_a[4] <= i0;
        mul_b[4] <= C4_16;
    end
    3'd2:
    begin
        i0      <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul_a[0] <= block_in_0_1;
        mul_b[0] <= C7_16;
        mul_a[1] <= block_in_6_7;
        mul_b[1] <= C1_16;
        mul_a[2] <= block_in_4_5;
        mul_b[2] <= C3_16;
        mul_a[3] <= block_in_2_3;
        mul_b[3] <= C5_16;
    end
    3'd4:
    begin
        mul_a[0] <= block_in_0_1;
        mul_b[0] <= C7_16;
        mul_a[1] <= block_in_6_7;
        mul_b[1] <= C1_16;
        mul_a[2] <= block_in_4_5;
        mul_b[2] <= C3_16;
        mul_a[3] <= block_in_2_3;
        mul_b[3] <= C5_16;
    end
    3'd5:
    begin
        mul_a[0] <= block_in_2_3;
        mul_b[0] <= C6_16;
        mul_a[1] <= block_in_6_7;
        mul_b[1] <= C2_16;
        mul_a[4] <= i0;
        mul_b[4] <= C4_16;
    end
    default:
        ;
    endcase
end

// Split multiplier: b is split into three parts (b[4:0], b[9:5], b[12:10])
// so the 13-row partial-product array is replaced by three small arrays
// (5/5/3 rows) registered separately; the final carry-propagate add occupies
// the next pipeline register. Bit-exact with a full 32-bit-wrapping product.
wire signed [12:0] mb_lo  [0:4];
wire signed [12:0] mb_mid [0:4];
wire signed [12:0] mb_hi  [0:4];
wire signed [37:0] mul_sum [0:4];

reg signed [26:0] mul_lo  [0:4];
reg signed [26:0] mul_mid [0:4];
reg signed [24:0] mul_hi  [0:4];
reg signed [31:0] mul     [0:4];

genvar g;
generate
for (g = 0; g < 5; g = g + 1) begin: g_split
    assign mb_lo[g]  = {8'd0, mul_b[g][4:0]};
    assign mb_mid[g] = {8'd0, mul_b[g][9:5]};
    assign mb_hi[g]  = {10'd0, mul_b[g][12:10]};
    assign mul_sum[g] = {{11{mul_lo[g][26]}}, mul_lo[g]} +
                        ({{11{mul_mid[g][26]}}, mul_mid[g]} <<< 5) +
                        ({{13{mul_hi[g][24]}}, mul_hi[g]} <<< 10);
end
endgenerate

integer mi;

always @ (posedge clk )
if (rst)
begin
    for (mi = 0; mi < 5; mi = mi + 1)
    begin
        mul_lo[mi]  <= 27'b0;
        mul_mid[mi] <= 27'b0;
        mul_hi[mi]  <= 25'b0;
    end
end
else
begin
    for (mi = 0; mi < 5; mi = mi + 1)
    begin
        mul_lo[mi]  <= mul_a[mi] * mb_lo[mi];
        mul_mid[mi] <= mul_a[mi] * mb_mid[mi];
        mul_hi[mi]  <= mul_a[mi] * mb_hi[mi];
    end
end

always @ (posedge clk )
if (rst)
begin
    for (mi = 0; mi < 5; mi = mi + 1)
        mul[mi] <= 32'b0;
end
else
begin
    for (mi = 0; mi < 5; mi = mi + 1)
        mul[mi] <= mul_sum[mi][31:0];
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

reg signed [31:0] o_s5;
reg signed [31:0] o_s6;
reg signed [31:0] o_s7;
reg signed [31:0] o_t0;
reg signed [31:0] o_t1;
reg signed [31:0] o_t2;
reg signed [31:0] o_t3;
reg signed [31:0] o_t4;
reg signed [31:0] o_t7;
reg signed [31:0] o_prod65;
reg signed [31:0] o_prod56;

// 1/sqrt(2) scaling: (x*181)/256 with the exact original fixed-point
// semantics - the multiply wraps at 32 bits and the signed division
// truncates toward zero (arithmetic shift + correction when p<0 and the
// remainder is non-zero). The o_t6_5/o_t5_6 butterfly sums are formed
// from pre-multiplied halves registered as soon as their operands are
// available (beats idx2 / idx5) and combined one adder deep at beat idx6,
// keeping every individual cycle shallow. Bit-identical: multiplication
// distributes over +/- modulo 2^32.
reg signed [31:0] o_prod_s76;   // (o_s7 - o_s6) * 181
reg signed [31:0] o_prod_t05;   // (o_t0 - o_s5) * 181
wire signed [31:0] prod_s76_w = (o_s7 - o_s6) * 32'sd181;
wire signed [31:0] prod_t05_w = (o_t0 - o_s5) * 32'sd181;

// truncate-toward-zero division by 256 of the wrapped product
wire signed [31:0] p65_sft  = o_prod65 >>> 8;
wire signed [31:0] p56_sft  = o_prod56 >>> 8;
wire               p65_corr = o_prod65[31] & (o_prod65[7:0] != 8'd0);
wire               p56_corr = o_prod56[31] & (o_prod56[7:0] != 8'd0);

always @ (posedge clk )
if (rst)
begin
    o_s5   <= 32'b0;
    o_s6   <= 32'b0;
    o_s7   <= 32'b0;
    o_t0   <= 32'b0;
    o_t1   <= 32'b0;
    o_t2   <= 32'b0;
    o_t3   <= 32'b0;
    o_t4   <= 32'b0;
    o_t7   <= 32'b0;
    o_prod_s76 <= 32'b0;
    o_prod_t05 <= 32'b0;
    o_prod65 <= 32'b0;
    o_prod56 <= 32'b0;
end
else
begin
    case (out_stg2_idx_q)
    3'd0:
    begin
        o_t3 <= mul[0] + mul[1]; // s3
    end
    3'd1:
    begin
        o_s7 <= mul[0] + mul[1];
        o_s6 <= mul[2] + mul[3];
        o_t0 <= mul[4];          // s0
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3; // t0
        o_t3 <= o_t0 - o_t3; // t3
        o_t7 <= o_s6 + o_s7;
        o_prod_s76 <= prod_s76_w;
    end
    3'd3:
    begin
        o_t4 <= (mul[0] - mul[1]) + (mul[2] - mul[3]);
    end
    3'd4:
    begin
        o_t0 <= mul[0] - mul[1]; // s4
        o_s5 <= mul[2] - mul[3];
    end
    3'd5:
    begin
        o_t3 <= mul[0] - mul[1]; // s2
        o_t4 <= mul[4]; // s1
        o_prod_t05 <= prod_t05_w;
    end
    3'd6:
    begin
        o_t1 <= o_t4 + o_t3;
        o_t2 <= o_t4 - o_t3;
        o_prod65 <= o_prod_s76 - o_prod_t05;
        o_prod56 <= o_prod_s76 + o_prod_t05;
    end
    default:
    begin
        o_s5 <= p65_sft + {31'b0, p65_corr}; // 1/sqrt(2)
        o_s6 <= p56_sft + {31'b0, p56_corr}; // 1/sqrt(2)
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

// Pass outputs are bounded by the l1 norm of the scaled IDCT basis
// (<= ~16k for 16-bit inputs): 16 bits are sufficient and exact.
reg signed [21:0] block_out[0:7];
reg signed [21:0] block_out_tmp;

wire signed [31:0] sft0 = (o_t0 + o_t7) >>> OUT_SHIFT;
wire signed [31:0] sft7 = (o_t0 - o_t7) >>> OUT_SHIFT;
wire signed [31:0] sft3 = (o_t3 + o_t4) >>> OUT_SHIFT;
wire signed [31:0] sft4 = (o_t3 - o_t4) >>> OUT_SHIFT;
wire signed [31:0] sft2 = (o_t2 + o_s5) >>> OUT_SHIFT;
wire signed [31:0] sft5 = (o_t2 - o_s5) >>> OUT_SHIFT;
wire signed [31:0] sft1 = (o_t1 + o_s6) >>> OUT_SHIFT;
wire signed [31:0] sft6 = (o_t1 - o_s6) >>> OUT_SHIFT;

always @ (posedge clk )
if (rst)
begin
    block_out[0] <= 22'b0;
    block_out[1] <= 22'b0;
    block_out[2] <= 22'b0;
    block_out[3] <= 22'b0;
    block_out[4] <= 22'b0;
    block_out[5] <= 22'b0;
    block_out[6] <= 22'b0;
    block_out[7] <= 22'b0;
    block_out_tmp <= 22'b0;
end
else if (out_stg3_valid_q)
begin
    if (out_stg3_idx_q == 3'd3)
    begin
        block_out[0] <= sft0[21:0];
        block_out_tmp <= sft7[21:0]; // block_out[7]
        block_out[3] <= sft3[21:0];
        block_out[4] <= sft4[21:0];
    end

    if (out_stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (out_stg3_idx_q == 3'd7)
    begin
        block_out[2] <= sft2[21:0];
        block_out[5] <= sft5[21:0];
        block_out[1] <= sft1[21:0];
        block_out[6] <= sft6[21:0];
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

wire signed [21:0] block_out_rd = block_out[ptr_q[2:0]];

assign outport_valid_o = valid_q[6];
assign outport_data_o  = {{10{block_out_rd[21]}}, block_out_rd};

generate
if (IDX_CONV != 0)
    assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
else
    assign outport_idx_o = ptr_q;
endgenerate

endmodule

//-----------------------------------------------------------------
// Corner-turn buffer entry: 32 addresses. Inter-pass values are bounded
// (<= ~16k signed), so entries are stored as 16 bits and sign-extended on
// read - bit-exact with the original 32-bit storage.
module jpeg_idct_transpose_ram
(
     input           clk0_i
    ,input  [  4:0]  addr0_i
    ,input  [ 21:0]  data0_i
    ,input           wr0_i
    ,input           clk1_i
    ,input           blk_i
    ,input  [ 15:0]  oh_i
    ,output [ 31:0]  data1_o
);

reg [21:0]   ram [31:0];
reg [21:0] ram_read1_q;

integer n;
reg [21:0] rd_mux_w;

always @ (posedge clk0_i)
begin
    if (wr0_i)
        ram[addr0_i][21:0] <= data0_i[21:0];
end

always @ *
begin
    rd_mux_w = 22'b0;
    for (n = 0; n < 16; n = n + 1)
        if (oh_i[n])
            rd_mux_w = ram[{blk_i, n[3:0]}];
end

always @ (posedge clk1_i)
begin
    ram_read1_q <= rd_mux_w;
end

assign data1_o = {{10{ram_read1_q[21]}}, ram_read1_q};

endmodule

//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
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

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
// One-hot copy of the (block-local) read address: the per-entry read select
// is then a 16:1 AND-OR mux instead of a 32:1 mux from a binary address, and
// the address sequence (+8/-8/hold/+1) maps to wiring permutations plus a
// cheap 3:1 mux on the one-hot register update.
reg [15:0]  rd_oh_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire [31:0] outport_data0_w;
wire [31:0] outport_data1_w;
wire [31:0] outport_data2_w;
wire [31:0] outport_data3_w;

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram u_ram0
(     .clk0_i(clk)
    , .addr0_i(wr_ptr_w), .data0_i(inport_data_i[21:0]), .wr0_i(wr0_w)
    , .clk1_i(clk)
    , .blk_i(block_rd_q), .oh_i(rd_oh_q)
    , .data1_o(outport_data0_w)
);

jpeg_idct_transpose_ram u_ram1
(     .clk0_i(clk)
    , .addr0_i(wr_ptr_w), .data0_i(inport_data_i[21:0]), .wr0_i(wr1_w)
    , .clk1_i(clk)
    , .blk_i(block_rd_q), .oh_i(rd_oh_q)
    , .data1_o(outport_data1_w)
);

jpeg_idct_transpose_ram u_ram2
(     .clk0_i(clk)
    , .addr0_i(wr_ptr_w), .data0_i(inport_data_i[21:0]), .wr0_i(wr2_w)
    , .clk1_i(clk)
    , .blk_i(block_rd_q), .oh_i(rd_oh_q)
    , .data1_o(outport_data2_w)
);

jpeg_idct_transpose_ram u_ram3
(     .clk0_i(clk)
    , .addr0_i(wr_ptr_w), .data0_i(inport_data_i[21:0]), .wr0_i(wr3_w)
    , .clk1_i(clk)
    , .blk_i(block_rd_q), .oh_i(rd_oh_q)
    , .data1_o(outport_data3_w)
);

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
    rd_oh_q <= 16'h0001;
else if (img_start_i || state_q == STATE_IDLE)
    rd_oh_q <= 16'h0001;           // address 0
else if (state_q == STATE_SETUP)
    rd_oh_q <= 16'h0100;           // address 8
else if (state_q == STATE_ACTIVE)
begin
    case (rd_idx_q[2:0])
    3'd6: rd_oh_q <= {rd_oh_q[14:0], rd_oh_q[15]};  // +1 (rotate)
    3'd2: rd_oh_q <= rd_oh_q;                        // hold
    default: rd_oh_q <= {rd_oh_q[7:0], rd_oh_q[15:8]}; // +8 / -8 (swap halves)
    endcase
end

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;

endmodule

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

wire          idct_x_valid_w;
wire [ 31:0]  idct_x_data_w;
wire [  5:0]  idct_x_idx_w;
wire          idct_x_accept_w;

wire [31:0] x_d0 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire [31:0] x_d1 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire [31:0] x_d2 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire [31:0] x_d3 = {{16{inport_data3_i[15]}}, inport_data3_i};

jpeg_idct_pass #(.OUT_SHIFT(11), .IDX_CONV(0), .A_W(22)) u_idct_x
(
     .clk(clk), .rst(rst), .img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(x_d0)
    ,.inport_data1_i(x_d1)
    ,.inport_data2_i(x_d2)
    ,.inport_data3_i(x_d3)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(idct_x_valid_w)
    ,.outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [  2:0]  transpose_idx_w;
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
wire          transpose_ready_w = 1'b1;

jpeg_idct_transpose u_transpose
(
     .clk(clk), .rst(rst), .img_start_i(img_start_i)
    ,.inport_valid_i(idct_x_valid_w)
    ,.inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w)
    ,.inport_accept_o(idct_x_accept_w)
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
    ,.outport_ready_i(transpose_ready_w)
);

jpeg_idct_pass #(.OUT_SHIFT(15), .IDX_CONV(1)) u_idct_y
(
     .clk(clk), .rst(rst), .img_start_i(img_start_i)
    ,.inport_valid_i(transpose_valid_w)
    ,.inport_data0_i(transpose_data0_w)
    ,.inport_data1_i(transpose_data1_w)
    ,.inport_data2_i(transpose_data2_w)
    ,.inport_data3_i(transpose_data3_w)
    ,.inport_idx_i(transpose_idx_w)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

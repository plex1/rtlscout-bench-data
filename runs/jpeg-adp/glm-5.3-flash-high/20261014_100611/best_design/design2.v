//-----------------------------------------------------------------
// jpeg_idct_2d — full 8x8 2-D IDCT (core_jpeg algorithm), optimized.
//  - narrowed datapath widths (provably sufficient, wrap semantics kept)
//  - transpose implemented as flat register array indexed by coefficient
//    index (no address decode / 32:1 read mux); read via small per-lane
//    muxes matching the fixed y-pass read schedule.
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
wire [ 20:0]  idct_x_data_w;
wire [  5:0]  idct_x_idx_w;
wire          idct_x_accept_w;

jpeg_idct_x
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(inport_data0_i)
    ,.inport_data1_i(inport_data1_i)
    ,.inport_data2_i(inport_data2_i)
    ,.inport_data3_i(inport_data3_i)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(idct_x_valid_w)
    ,.outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [ 20:0]  transpose_data0_w;
wire [ 20:0]  transpose_data1_w;
wire [ 20:0]  transpose_data2_w;
wire [ 20:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
    ,.inport_valid_i(idct_x_valid_w)
    ,.inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w)
    ,.inport_accept_o(idct_x_accept_w)
    ,.outport_ready_i(1'b1)
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
);

jpeg_idct_y
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.img_end_i(img_end_i)
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

//-----------------------------------------------------------------
// 1-D IDCT row pass (OUT_SHIFT = 11)
//-----------------------------------------------------------------
module jpeg_idct_x
#(
     parameter OUT_SHIFT        = 11
    ,parameter INPUT_WIDTH      = 16
)
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
    ,output [ 20:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

localparam [15:0] C1_16 = 4017;
localparam [15:0] C2_16 = 3784;
localparam [15:0] C3_16 = 3406;
localparam [15:0] C4_16 = 2896;
localparam [15:0] C5_16 = 2276;
localparam [15:0] C6_16 = 1567;
localparam [15:0] C7_16 = 799;

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

// signed division by 256 truncating toward zero
function signed [31:0] dv256;
    input signed [31:0] x;
    reg   signed [31:0] sh;
    begin
        sh = x >>> 8;
        if (x[31] && (|x[7:0]))
            sh = sh + 32'd1;
        dv256 = sh;
    end
endfunction

reg signed [31:0] i0;
reg signed [16:0] mul0_a;
reg signed [12:0] mul0_b;
reg signed [16:0] mul1_a;
reg signed [12:0] mul1_b;
reg signed [16:0] mul2_a;
reg signed [12:0] mul2_b;
reg signed [16:0] mul3_a;
reg signed [12:0] mul3_b;
reg signed [16:0] mul4_a;
reg signed [12:0] mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 17'b0;
    mul0_b <= 13'b0;
    mul1_a <= 17'b0;
    mul1_b <= 13'b0;
    mul2_a <= 17'b0;
    mul2_b <= 13'b0;
    mul3_a <= 17'b0;
    mul3_b <= 13'b0;
    mul4_a <= 17'b0;
    mul4_b <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3[16:0];
        mul0_b <= C2_16;
        mul1_a <= block_in_6_7[16:0];
        mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1[16:0];
        mul0_b <= C1_16;
        mul1_a <= block_in_6_7[16:0];
        mul1_b <= C7_16;
        mul2_a <= block_in_4_5[16:0];
        mul2_b <= C5_16;
        mul3_a <= block_in_2_3[16:0];
        mul3_b <= C3_16;
        mul4_a <= i0[16:0];
        mul4_b <= C4_16;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1[16:0];
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7[16:0];
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5[16:0];
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3[16:0];
        mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1[16:0];
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7[16:0];
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5[16:0];
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3[16:0];
        mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3[16:0];
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7[16:0];
        mul1_b <= C2_16;
        mul4_a <= i0[16:0];
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
end

// all true products fit in 32 bits (no wrap) -> exact
reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0;
    mul1_q <= 32'b0;
    mul2_q <= 32'b0;
    mul3_q <= 32'b0;
    mul4_q <= 32'b0;
end
else
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0;
    mul1 <= 32'b0;
    mul2 <= 32'b0;
    mul3 <= 32'b0;
    mul4 <= 32'b0;
end
else
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;

always @ (posedge clk)
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

always @ (posedge clk)
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

always @ (posedge clk)
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
reg signed [31:0] o_t5;
reg signed [31:0] o_t6;
reg signed [31:0] o_t7;
reg signed [31:0] o_t6_5;
reg signed [31:0] o_t5_6;

always @ (posedge clk)
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
    o_t5   <= 32'b0;
    o_t6   <= 32'b0;
    o_t7   <= 32'b0;
    o_t6_5 <= 32'b0;
    o_t5_6 <= 32'b0;
end
else
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
        o_t4 <= mul4; // s1
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
        o_s5 <= dv256(o_t6_5 * 181); // 1/sqrt(2)
        o_s6 <= dv256(o_t5_6 * 181); // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
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

// (32-bit)>>>11 always fits in signed 21 bits
reg signed [20:0] block_out [0:7];
reg signed [20:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 21'b0;
    block_out[1] <= 21'b0;
    block_out[2] <= 21'b0;
    block_out[3] <= 21'b0;
    block_out[4] <= 21'b0;
    block_out[5] <= 21'b0;
    block_out[6] <= 21'b0;
    block_out[7] <= 21'b0;
    block_out_tmp <= 21'b0;
end
else if (out_stg3_valid_q)
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
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = ptr_q;

endmodule

//-----------------------------------------------------------------
// Double-buffered corner turn.
// Storage: mem[buf][idx], idx = r*8+c (x-pass output coefficient index).
// Read schedule (fixed, matches the y-pass expectation):
//   at read beat (column c, beat b) lane l carries X(2l + e(b), c),
//   with e(b) = b[0] ^ b[2]  (i.e. e = 0,1,0,1,1,0,1,0 for b = 0..7).
//-----------------------------------------------------------------
module jpeg_idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 20:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i

    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 20:0]  outport_data0_o
    ,output [ 20:0]  outport_data1_o
    ,output [ 20:0]  outport_data2_o
    ,output [ 20:0]  outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg [1:0] block_ready_q;
reg       block_wr_q;
reg       block_rd_q;
reg [5:0] rd_idx_q;

reg [20:0] mem [0:1][0:63];

assign inport_accept_o = ~block_ready_q[block_wr_q];

// write port
always @ (posedge clk)
if (!rst && inport_valid_i && inport_accept_o)
    mem[block_wr_q][inport_idx_i] <= inport_data_i;

// buffer control
always @ (posedge clk)
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

    if (outport_valid_o && rd_idx_q == 6'd63)
    begin
        block_ready_q[block_rd_q] <= 1'b0;
        block_rd_q <= ~block_rd_q;
    end
end

// FSM
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

always @ (posedge clk)
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @ (posedge clk)
if (rst)
    rd_idx_q <= 6'b0;
else if (img_start_i)
    rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE)
    rd_idx_q <= rd_idx_q + 6'd1;

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o   = rd_idx_q[2:0];

// read side: prefetch beat k+1 (beat 0 during SETUP) into the output regs.
// nidx = rd_idx_q + 1 (mod 64) is the next beat to be consumed:
// lane l = X(2l + e, c), c = nidx[5:3], e = nidx[0]^nidx[2].
wire [2:0] c_w = rd_idx_q[5:3];
wire [2:0] b_w = rd_idx_q[2:0];
wire [5:0] nidx_w = rd_idx_q + 6'd1;
wire [2:0] nc_w  = state_q == STATE_SETUP ? 3'd0 : nidx_w[5:3];
wire       ne_w  = nidx_w[0] ^ nidx_w[2];

reg [20:0] rd0_q, rd1_q, rd2_q, rd3_q;

always @ (posedge clk)
if (rst)
begin
    rd0_q <= 21'b0;
    rd1_q <= 21'b0;
    rd2_q <= 21'b0;
    rd3_q <= 21'b0;
end
else if (state_q == STATE_SETUP)
begin
    rd0_q <= mem[block_rd_q][0  ];
    rd1_q <= mem[block_rd_q][16 ];
    rd2_q <= mem[block_rd_q][32 ];
    rd3_q <= mem[block_rd_q][48 ];
end
else if (state_q == STATE_ACTIVE)
begin
    if (ne_w)
    begin
        rd0_q <= mem[block_rd_q][8  + nc_w];
        rd1_q <= mem[block_rd_q][24 + nc_w];
        rd2_q <= mem[block_rd_q][40 + nc_w];
        rd3_q <= mem[block_rd_q][56 + nc_w];
    end
    else
    begin
        rd0_q <= mem[block_rd_q][0  + nc_w];
        rd1_q <= mem[block_rd_q][16 + nc_w];
        rd2_q <= mem[block_rd_q][32 + nc_w];
        rd3_q <= mem[block_rd_q][48 + nc_w];
    end
end

assign outport_data0_o = rd0_q;
assign outport_data1_o = rd1_q;
assign outport_data2_o = rd2_q;
assign outport_data3_o = rd3_q;

endmodule

//-----------------------------------------------------------------
// 1-D IDCT column pass (OUT_SHIFT = 15)
//-----------------------------------------------------------------
module jpeg_idct_y
#(
     parameter OUT_SHIFT        = 15
    ,parameter INPUT_WIDTH      = 32
)
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [ 20:0]  inport_data0_i
    ,input  [ 20:0]  inport_data1_i
    ,input  [ 20:0]  inport_data2_i
    ,input  [ 20:0]  inport_data3_i
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

// inputs are 21-bit x-pass outputs (which fit signed 21 bits exactly)
wire signed [31:0] block_in_0_1 = {{11{inport_data0_i[20]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{11{inport_data1_i[20]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{11{inport_data2_i[20]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{11{inport_data3_i[20]}}, inport_data3_i};

// signed division by 256 truncating toward zero
function signed [31:0] dv256;
    input signed [31:0] x;
    reg   signed [31:0] sh;
    begin
        sh = x >>> 8;
        if (x[31] && (|x[7:0]))
            sh = sh + 32'd1;
        dv256 = sh;
    end
endfunction

reg signed [31:0] i0;
reg signed [21:0] mul0_a;
reg signed [12:0] mul0_b;
reg signed [21:0] mul1_a;
reg signed [12:0] mul1_b;
reg signed [21:0] mul2_a;
reg signed [12:0] mul2_b;
reg signed [21:0] mul3_a;
reg signed [12:0] mul3_b;
reg signed [22:0] mul4_a;
reg signed [12:0] mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 22'b0;
    mul0_b <= 13'b0;
    mul1_a <= 22'b0;
    mul1_b <= 13'b0;
    mul2_a <= 22'b0;
    mul2_b <= 13'b0;
    mul3_a <= 22'b0;
    mul3_b <= 13'b0;
    mul4_a <= 23'b0;
    mul4_b <= 13'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3[21:0];
        mul0_b <= C2_16;
        mul1_a <= block_in_6_7[21:0];
        mul1_b <= C6_16;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1[21:0];
        mul0_b <= C1_16;
        mul1_a <= block_in_6_7[21:0];
        mul1_b <= C7_16;
        mul2_a <= block_in_4_5[21:0];
        mul2_b <= C5_16;
        mul3_a <= block_in_2_3[21:0];
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
        mul0_a <= block_in_0_1[21:0];
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7[21:0];
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5[21:0];
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3[21:0];
        mul3_b <= C5_16;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1[21:0];
        mul0_b <= C7_16;
        mul1_a <= block_in_6_7[21:0];
        mul1_b <= C1_16;
        mul2_a <= block_in_4_5[21:0];
        mul2_b <= C3_16;
        mul3_a <= block_in_2_3[21:0];
        mul3_b <= C5_16;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3[21:0];
        mul0_b <= C6_16;
        mul1_a <= block_in_6_7[21:0];
        mul1_b <= C2_16;
        mul4_a <= i0;
        mul4_b <= C4_16;
    end
    default:
        ;
    endcase
end

// 32-bit wrap semantics preserved: operands widened to 32 bits before mult
reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0;
    mul1_q <= 32'b0;
    mul2_q <= 32'b0;
    mul3_q <= 32'b0;
    mul4_q <= 32'b0;
end
else
begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0;
    mul1 <= 32'b0;
    mul2 <= 32'b0;
    mul3 <= 32'b0;
    mul4 <= 32'b0;
end
else
begin
    mul0 <= mul0_q;
    mul1 <= mul1_q;
    mul2 <= mul2_q;
    mul3 <= mul3_q;
    mul4 <= mul4_q;
end

reg        out_stg0_valid_q;
reg [2:0]  out_stg0_idx_q;

always @ (posedge clk)
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

always @ (posedge clk)
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

always @ (posedge clk)
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
reg signed [31:0] o_t5;
reg signed [31:0] o_t6;
reg signed [31:0] o_t7;
reg signed [31:0] o_t6_5;
reg signed [31:0] o_t5_6;

always @ (posedge clk)
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
    o_t5   <= 32'b0;
    o_t6   <= 32'b0;
    o_t7   <= 32'b0;
    o_t6_5 <= 32'b0;
    o_t5_6 <= 32'b0;
end
else
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
        o_t4 <= mul4; // s1
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
        o_s5 <= dv256(o_t6_5 * 181); // 1/sqrt(2)
        o_s6 <= dv256(o_t5_6 * 181); // 1/sqrt(2)
    end
    endcase
end

reg        out_stg3_valid_q;
reg [2:0]  out_stg3_idx_q;

always @ (posedge clk)
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

// (32-bit)>>>15 always fits in signed 17 bits
reg signed [16:0] block_out [0:7];
reg signed [16:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 17'b0;
    block_out[1] <= 17'b0;
    block_out[2] <= 17'b0;
    block_out[3] <= 17'b0;
    block_out[4] <= 17'b0;
    block_out[5] <= 17'b0;
    block_out[6] <= 17'b0;
    block_out[7] <= 17'b0;
    block_out_tmp <= 17'b0;
end
else if (out_stg3_valid_q)
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
assign outport_data_o  = {{15{block_out[ptr_q[2:0]][16]}}, block_out[ptr_q[2:0]]};

// ptr_conv is a bit rotation: {k[2:0], k[5:3]}
assign outport_idx_o   = {ptr_q[2:0], ptr_q[5:3]};

endmodule

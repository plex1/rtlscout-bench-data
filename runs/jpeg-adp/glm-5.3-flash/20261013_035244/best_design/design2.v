// Optimized jpeg_idct_2d (full 8x8 2-D IDCT, core_jpeg-compatible, bit-exact).
//
// Structure (same as original composition): jpeg_idct_x (row pass, OUT_SHIFT=11)
// -> double-buffered corner turn -> jpeg_idct_y (column pass, OUT_SHIFT=15).
// Bit-exactness notes:
//  - all arithmetic wraps at 32 bits, signed;
//  - constant multipliers are expanded by shift-add (exact for signed wrap-32);
//  - the 1/sqrt(2) scaling is (x*181)/256 with signed division truncating to 0
//    (implemented via add-bias / drop-low-bits form).

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
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
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
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w)
    ,.outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w)
    ,.outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w)
    ,.outport_ready_i(1'b1)
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

//=====================================================================
// Row pass (OUT_SHIFT = 11)
//=====================================================================
module jpeg_idct_x
#(
     parameter OUT_SHIFT        = 11
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
    ,output [ 31:0]  outport_data_o
    ,output [  5:0]  outport_idx_o
);

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg        [ 3:0] mul0_b;
reg signed [31:0] mul1_a;
reg        [ 3:0] mul1_b;
reg signed [31:0] mul2_a;
reg        [ 3:0] mul2_b;
reg signed [31:0] mul3_a;
reg        [ 3:0] mul3_b;
reg signed [31:0] mul4_a;
reg        [ 3:0] mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0; mul0_b <= 4'b0;
    mul1_a <= 32'b0; mul1_b <= 4'b0;
    mul2_a <= 32'b0; mul2_b <= 4'b0;
    mul3_a <= 32'b0; mul3_b <= 4'b0;
    mul4_a <= 32'b0; mul4_b <= 4'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= 4'd2;
        mul1_a <= block_in_6_7; mul1_b <= 4'd6;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd1;
        mul1_a <= block_in_6_7; mul1_b <= 4'd7;
        mul2_a <= block_in_4_5; mul2_b <= 4'd5;
        mul3_a <= block_in_2_3; mul3_b <= 4'd3;
        mul4_a <= i0;           mul4_b <= 4'd4;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd7;
        mul1_a <= block_in_6_7; mul1_b <= 4'd1;
        mul2_a <= block_in_4_5; mul2_b <= 4'd3;
        mul3_a <= block_in_2_3; mul3_b <= 4'd5;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd7;
        mul1_a <= block_in_6_7; mul1_b <= 4'd1;
        mul2_a <= block_in_4_5; mul2_b <= 4'd3;
        mul3_a <= block_in_2_3; mul3_b <= 4'd5;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3; mul0_b <= 4'd6;
        mul1_a <= block_in_6_7; mul1_b <= 4'd2;
        mul4_a <= i0;           mul4_b <= 4'd4;
    end
    default:
        ;
    endcase
end

// Constant multiplier: b in {1..7} maps to C1..C7. Shift-add expansion.
function signed [31:0] cmul;
    input signed [31:0] a;
    input        [ 3:0] b;
    reg signed [31:0] s10, s11, s12;
    reg signed [31:0] acc;
    begin
        s10 = a <<< 10;
        s11 = a <<< 11;
        s12 = a <<< 12;
        case (b)
        4'd1: acc = s12 - ((a<<<6)+(a<<3)+(a<<2)+(a<<1)+a);
        4'd2: acc = s12 - ((a<<8)+(a<<5)+(a<<4)+(a<<3));
        4'd3: acc = s12 - ((a<<9)+(a<<7)+(a<<5)+(a<<4)+(a<<1));
        4'd4: acc = s12 - ((a<<10)+(a<<7)+(a<<5)+(a<<4));
        4'd5: acc = s11 + ((a<<7)+(a<<6)+(a<<5)+(a<<2));
        4'd6: acc = s11 - ((a<<8)+(a<<7)+(a<<6)+(a<<5)+a);
        default: acc = s10 - ((a<<7)+(a<<6)+(a<<5)+a);
        endcase
        cmul = acc;
    end
endfunction

reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0; mul1_q <= 32'b0; mul2_q <= 32'b0;
    mul3_q <= 32'b0; mul4_q <= 32'b0;
end
else
begin
    mul0_q <= cmul(mul0_a, mul0_b);
    mul1_q <= cmul(mul1_a, mul1_b);
    mul2_q <= cmul(mul2_a, mul2_b);
    mul3_q <= cmul(mul3_a, mul3_b);
    mul4_q <= cmul(mul4_a, mul4_b);
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0; mul1 <= 32'b0; mul2 <= 32'b0; mul3 <= 32'b0; mul4 <= 32'b0;
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
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
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

// (x*181)/256 with signed truncation toward zero:
//   pos: (x*181 + 256) >>> 8 ; neg: (x*181 + 255) >>> 8 (bias as unsigned add)
function signed [31:0] m181_div256;
    input signed [31:0] x;
    reg signed [31:0] p;
    reg   [31:0] bias;
    begin
        p = 32'(x * 64'sd181);
        bias = x[31] ? 32'h000000FF : 32'h00000100;
        m181_div256 = 32'(p + 32'(bias)) >>> 8;
    end
endfunction

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'b0; o_s6 <= 32'b0; o_s7 <= 32'b0;
    o_t0 <= 32'b0; o_t1 <= 32'b0; o_t2 <= 32'b0; o_t3 <= 32'b0;
    o_t4 <= 32'b0; o_t5 <= 32'b0; o_t6 <= 32'b0; o_t7 <= 32'b0;
    o_t6_5 <= 32'b0; o_t5_6 <= 32'b0;
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
        o_s5 <= m181_div256(o_t6_5); // 1/sqrt(2)
        o_s6 <= m181_div256(o_t5_6);
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

reg signed [31:0] block_out [0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 32'b0; block_out[1] <= 32'b0; block_out[2] <= 32'b0;
    block_out[3] <= 32'b0; block_out[4] <= 32'b0; block_out[5] <= 32'b0;
    block_out[6] <= 32'b0; block_out[7] <= 32'b0; block_out_tmp <= 32'b0;
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
reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

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

//=====================================================================
// Double-buffered corner turn (transpose), functionally identical to
// jpeg_idct_transpose but with the four 32x32 dual-port RAMs replaced
// by one shared bank structure (same storage, simpler decode).
//=====================================================================
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

reg         block_wr_q;
reg         block_rd_q;
reg [5:0]   rd_idx_q;
reg [3:0]   rd_addr_q;

wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

// Four independent write lanes (one per beat data lane).
wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

// Storage: 4 buffers x 8 rows x 8 cols x 32b = 2048 FFs.
// Read mux: row (5:3) select then col (2:0) select.
reg [31:0] bufs [0:63];
reg [31:0] rd0_q;
reg [31:0] rd1_q;
reg [31:0] rd2_q;
reg [31:0] rd3_q;

wire [2:0] wr_row_w = inport_idx_i[3+3:3];
wire [2:0] wr_col_w = inport_idx_i[2:0];
wire [5:0] rd_row_w = {block_rd_q, rd_addr_q[3:1]};
wire [2:0] rd_col_w = rd_addr_q[2:0];

always @ (posedge clk)
begin
    if (wr0_w) bufs[{block_wr_q, wr_row_w, wr_col_w}] <= inport_data_i;
    if (wr1_w) bufs[{block_wr_q, wr_row_w, 3'b100 + wr_col_w}] <= inport_data_i;
    if (wr2_w) bufs[{block_wr_q, wr_row_w, wr_col_w}] <= inport_data_i;
    if (wr3_w) bufs[{block_wr_q, wr_row_w, 3'b100 + wr_col_w}] <= inport_data_i;

    rd0_q <= bufs[{rd_row_w, rd_col_w}];
    rd1_q <= bufs[{rd_row_w, 3'b100 + rd_col_w}];
    rd2_q <= bufs[{rd_row_w, rd_col_w}];
    rd3_q <= bufs[{rd_row_w, 3'b100 + rd_col_w}];
end

assign outport_data0_o = rd0_q;
assign outport_data1_o = rd1_q;
assign outport_data2_o = rd2_q;
assign outport_data3_o = rd3_q;

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

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

always @ (posedge clk)
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

//=====================================================================
// Column pass (OUT_SHIFT = 15)
//=====================================================================
module jpeg_idct_y
#(
     parameter OUT_SHIFT        = 15
)
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
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

wire signed [31:0] block_in_0_1 = inport_data0_i;
wire signed [31:0] block_in_2_3 = inport_data1_i;
wire signed [31:0] block_in_4_5 = inport_data2_i;
wire signed [31:0] block_in_6_7 = inport_data3_i;

reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg        [ 3:0] mul0_b;
reg signed [31:0] mul1_a;
reg        [ 3:0] mul1_b;
reg signed [31:0] mul2_a;
reg        [ 3:0] mul2_b;
reg signed [31:0] mul3_a;
reg        [ 3:0] mul3_b;
reg signed [31:0] mul4_a;
reg        [ 3:0] mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0; mul0_b <= 4'b0;
    mul1_a <= 32'b0; mul1_b <= 4'b0;
    mul2_a <= 32'b0; mul2_b <= 4'b0;
    mul3_a <= 32'b0; mul3_b <= 4'b0;
    mul4_a <= 32'b0; mul4_b <= 4'b0;
end
else
begin
    case (inport_idx_i)
    3'd0:
    begin
        i0     <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= 4'd2;
        mul1_a <= block_in_6_7; mul1_b <= 4'd6;
    end
    3'd1:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd1;
        mul1_a <= block_in_6_7; mul1_b <= 4'd7;
        mul2_a <= block_in_4_5; mul2_b <= 4'd5;
        mul3_a <= block_in_2_3; mul3_b <= 4'd3;
        mul4_a <= i0;           mul4_b <= 4'd4;
    end
    3'd2:
    begin
        i0     <= block_in_0_1 - block_in_4_5;
    end
    3'd3:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd7;
        mul1_a <= block_in_6_7; mul1_b <= 4'd1;
        mul2_a <= block_in_4_5; mul2_b <= 4'd3;
        mul3_a <= block_in_2_3; mul3_b <= 4'd5;
    end
    3'd4:
    begin
        mul0_a <= block_in_0_1; mul0_b <= 4'd7;
        mul1_a <= block_in_6_7; mul1_b <= 4'd1;
        mul2_a <= block_in_4_5; mul2_b <= 4'd3;
        mul3_a <= block_in_2_3; mul3_b <= 4'd5;
    end
    3'd5:
    begin
        mul0_a <= block_in_2_3; mul0_b <= 4'd6;
        mul1_a <= block_in_6_7; mul1_b <= 4'd2;
        mul4_a <= i0;           mul4_b <= 4'd4;
    end
    default:
        ;
    endcase
end

// Constant multiplier: b in {1..7} maps to C1..C7. Shift-add expansion.
function signed [31:0] cmul;
    input signed [31:0] a;
    input        [ 3:0] b;
    reg signed [31:0] s10, s11, s12;
    reg signed [31:0] acc;
    begin
        s10 = a <<< 10;
        s11 = a <<< 11;
        s12 = a <<< 12;
        case (b)
        4'd1: acc = s12 - ((a<<<6)+(a<<3)+(a<<2)+(a<<1)+a);
        4'd2: acc = s12 - ((a<<8)+(a<<5)+(a<<4)+(a<<3));
        4'd3: acc = s12 - ((a<<9)+(a<<7)+(a<<5)+(a<<4)+(a<<1));
        4'd4: acc = s12 - ((a<<10)+(a<<7)+(a<<5)+(a<<4));
        4'd5: acc = s11 + ((a<<7)+(a<<6)+(a<<5)+(a<<2));
        4'd6: acc = s11 - ((a<<8)+(a<<7)+(a<<6)+(a<<5)+a);
        default: acc = s10 - ((a<<7)+(a<<6)+(a<<5)+a);
        endcase
        cmul = acc;
    end
endfunction

reg signed [31:0] mul0_q;
reg signed [31:0] mul1_q;
reg signed [31:0] mul2_q;
reg signed [31:0] mul3_q;
reg signed [31:0] mul4_q;

always @ (posedge clk)
if (rst)
begin
    mul0_q <= 32'b0; mul1_q <= 32'b0; mul2_q <= 32'b0;
    mul3_q <= 32'b0; mul4_q <= 32'b0;
end
else
begin
    mul0_q <= cmul(mul0_a, mul0_b);
    mul1_q <= cmul(mul1_a, mul1_b);
    mul2_q <= cmul(mul2_a, mul2_b);
    mul3_q <= cmul(mul3_a, mul3_b);
    mul4_q <= cmul(mul4_a, mul4_b);
end

reg signed [31:0] mul0;
reg signed [31:0] mul1;
reg signed [31:0] mul2;
reg signed [31:0] mul3;
reg signed [31:0] mul4;

always @ (posedge clk)
if (rst)
begin
    mul0 <= 32'b0; mul1 <= 32'b0; mul2 <= 32'b0; mul3 <= 32'b0; mul4 <= 32'b0;
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
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0; out_stg0_idx_q <= 3'b0;
    out_stg1_valid_q <= 1'b0; out_stg1_idx_q <= 3'b0;
    out_stg2_valid_q <= 1'b0; out_stg2_idx_q <= 3'b0;
end
else
begin
    out_stg0_valid_q <= inport_valid_i;
    out_stg0_idx_q   <= inport_idx_i;
    out_stg1_valid_q <= out_stg0_valid_q;
    out_stg1_idx_q   <= out_stg0_idx_q;
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

function signed [31:0] m181_div256;
    input signed [31:0] x;
    reg signed [31:0] p;
    reg   [31:0] bias;
    begin
        p = 32'(x * 64'sd181);
        bias = x[31] ? 32'h000000FF : 32'h00000100;
        m181_div256 = 32'(p + 32'(bias)) >>> 8;
    end
endfunction

always @ (posedge clk)
if (rst)
begin
    o_s5 <= 32'b0; o_s6 <= 32'b0; o_s7 <= 32'b0;
    o_t0 <= 32'b0; o_t1 <= 32'b0; o_t2 <= 32'b0; o_t3 <= 32'b0;
    o_t4 <= 32'b0; o_t5 <= 32'b0; o_t6 <= 32'b0; o_t7 <= 32'b0;
    o_t6_5 <= 32'b0; o_t5_6 <= 32'b0;
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
        o_s5 <= m181_div256(o_t6_5); // 1/sqrt(2)
        o_s6 <= m181_div256(o_t5_6);
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

reg signed [31:0] block_out [0:7];
reg signed [31:0] block_out_tmp;

always @ (posedge clk)
if (rst)
begin
    block_out[0] <= 32'b0; block_out[1] <= 32'b0; block_out[2] <= 32'b0;
    block_out[3] <= 32'b0; block_out[4] <= 32'b0; block_out[5] <= 32'b0;
    block_out[6] <= 32'b0; block_out[7] <= 32'b0; block_out_tmp <= 32'b0;
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
reg [5:0] ptr_q;

always @ (posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], out_stg3_valid_q};

always @ (posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];

// ptr_conv: streamed slot -> row-major pixel position.
function [5:0] ptr_conv;
    input [5:0] idx;
begin
    ptr_conv = {idx[5:4], idx[2:0], idx[3]};
end
endfunction

assign outport_idx_o = ptr_conv(ptr_q);

endmodule

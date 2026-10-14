//-----------------------------------------------------------------
// jpeg_idct_2d — full 8x8 2-D IDCT (ultraembedded/core_jpeg algorithm).
// v3 = v1 (proven) + one-hot decoded RAM read (low fanout, shallow
// AND-OR instead of deep 32:1 mux) to fix the transpose read path.
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

wire         x_valid_w;
wire [ 20:0] x_data_w;
wire [  5:0] x_idx_w;
wire         x_accept_w;

jpeg_idct_1d #(
     .LANE_W    (16)
    ,.I0_W      (18)
    ,.BO_W      (21)
    ,.OUT_SHIFT (11)
    ,.SWAP_IDX  (0)
) u_idct_x (
     .clk              (clk)
    ,.rst              (rst)
    ,.img_start_i      (img_start_i)
    ,.inport_valid_i   (inport_valid_i)
    ,.inport_data0_i   (inport_data0_i)
    ,.inport_data1_i   (inport_data1_i)
    ,.inport_data2_i   (inport_data2_i)
    ,.inport_data3_i   (inport_data3_i)
    ,.inport_idx_i     (inport_idx_i)
    ,.outport_valid_o  (x_valid_w)
    ,.outport_data_o   (x_data_w)
    ,.outport_idx_o    (x_idx_w)
);

wire         t_valid_w;
wire [ 20:0] t_data0_w;
wire [ 20:0] t_data1_w;
wire [ 20:0] t_data2_w;
wire [ 20:0] t_data3_w;
wire [  2:0] t_idx_w;

jpeg_idct_transpose
u_transpose
(
     .clk               (clk)
    ,.rst               (rst)
    ,.img_start_i       (img_start_i)
    ,.inport_valid_i    (x_valid_w)
    ,.inport_data_i     (x_data_w)
    ,.inport_idx_i      (x_idx_w)
    ,.inport_accept_o   (x_accept_w)
    ,.outport_valid_o   (t_valid_w)
    ,.outport_data0_o   (t_data0_w)
    ,.outport_data1_o   (t_data1_w)
    ,.outport_data2_o   (t_data2_w)
    ,.outport_data3_o   (t_data3_w)
    ,.outport_idx_o     (t_idx_w)
    ,.outport_ready_i   (1'b1)
);

wire [16:0] y_data_w;

jpeg_idct_1d #(
     .LANE_W    (21)
    ,.I0_W      (22)
    ,.BO_W      (17)
    ,.OUT_SHIFT (15)
    ,.SWAP_IDX  (1)
) u_idct_y (
     .clk              (clk)
    ,.rst              (rst)
    ,.img_start_i      (img_start_i)
    ,.inport_valid_i   (t_valid_w)
    ,.inport_data0_i   (t_data0_w)
    ,.inport_data1_i   (t_data1_w)
    ,.inport_data2_i   (t_data2_w)
    ,.inport_data3_i   (t_data3_w)
    ,.inport_idx_i     (t_idx_w)
    ,.outport_valid_o  (outport_valid_o)
    ,.outport_data_o   (y_data_w)
    ,.outport_idx_o    (outport_idx_o)
);

assign outport_data_o = {{15{y_data_w[16]}}, y_data_w};

endmodule

//-----------------------------------------------------------------
// 1-D IDCT pass (row or column), AAN-style factorization,
// constants scaled by 4096, 32-bit wrapping arithmetic context.
//-----------------------------------------------------------------
module jpeg_idct_1d #(
     parameter LANE_W    = 16
    ,parameter I0_W      = 18
    ,parameter BO_W      = 21
    ,parameter OUT_SHIFT = 11
    ,parameter SWAP_IDX  = 0
)
(
    // Inputs
     input                        clk
    ,input                        rst
    ,input                        img_start_i
    ,input                        inport_valid_i
    ,input  [LANE_W-1:0]          inport_data0_i
    ,input  [LANE_W-1:0]          inport_data1_i
    ,input  [LANE_W-1:0]          inport_data2_i
    ,input  [LANE_W-1:0]          inport_data3_i
    ,input  [ 2:0]                inport_idx_i

    // Outputs
    ,output                       outport_valid_o
    ,output [BO_W-1:0]            outport_data_o
    ,output [ 5:0]                outport_idx_o
);

localparam signed [15:0] C1 = 4017; // cos( pi/16) x4096
localparam signed [15:0] C2 = 3784; // cos(2pi/16) x4096
localparam signed [15:0] C3 = 3406; // cos(3pi/16) x4096
localparam signed [15:0] C4 = 2896; // cos(4pi/16) x4096
localparam signed [15:0] C5 = 2276; // cos(5pi/16) x4096
localparam signed [15:0] C6 = 1567; // cos(6pi/16) x4096
localparam signed [15:0] C7 =  799; // cos(7pi/16) x4096

wire signed [LANE_W-1:0] lane0 = inport_data0_i;
wire signed [LANE_W-1:0] lane1 = inport_data1_i;
wire signed [LANE_W-1:0] lane2 = inport_data2_i;
wire signed [LANE_W-1:0] lane3 = inport_data3_i;

//-----------------------------------------------------------------
// Operand capture: one register per constant (fixed-constant
// multipliers). Even coefficient set arrives at idx {0,2,5,7},
// odd set at idx {1,3,4,6}.
//-----------------------------------------------------------------
reg signed [I0_W-1:0] i0;
reg signed [I0_W-1:0] a1, a2, a3, a4, a5, a6, a7;

always @(posedge clk)
begin
    case (inport_idx_i)
    3'd0: begin
        i0 <= lane0 + lane2;    // c0 + c4
        a2 <= lane1;            // c2 (for *C2)
        a6 <= lane3;            // c6 (for *C6)
    end
    3'd1: begin
        a1 <= lane0;            // c1 (for *C1)
        a7 <= lane3;            // c7 (for *C7)
        a5 <= lane2;            // c5 (for *C5)
        a3 <= lane1;            // c3 (for *C3)
        a4 <= i0;               // (c0+c4) (for *C4)
    end
    3'd2: begin
        i0 <= lane0 - lane2;    // c0 - c4
    end
    3'd3: begin
        a7 <= lane0;            // c1 (for *C7)
        a1 <= lane3;            // c7 (for *C1)
        a3 <= lane2;            // c5 (for *C3)
        a5 <= lane1;            // c3 (for *C5)
    end
    3'd5: begin
        a6 <= lane1;            // c2 (for *C6)
        a2 <= lane3;            // c6 (for *C2)
        a4 <= i0;               // (c0-c4) (for *C4)
    end
    default: ;
    endcase
end

//-----------------------------------------------------------------
// Constant multipliers (single register stage). Products are kept
// at 32 bits so wrapping matches the original 32-bit context.
//-----------------------------------------------------------------
reg signed [31:0] p1, p2, p3, p4, p5, p6, p7;

always @(posedge clk)
begin
    p1 <= a1 * C1;
    p2 <= a2 * C2;
    p3 <= a3 * C3;
    p4 <= a4 * C4;
    p5 <= a5 * C5;
    p6 <= a6 * C6;
    p7 <= a7 * C7;
end

//-----------------------------------------------------------------
// Pipeline tags
//-----------------------------------------------------------------
reg         stg0_valid_q, stg1_valid_q, stg2_valid_q, stg3_valid_q;
reg [ 2:0]  stg0_idx_q, stg1_idx_q, stg2_idx_q, stg3_idx_q;

always @(posedge clk)
begin
    if (rst) begin
        stg0_valid_q <= 1'b0;
        stg1_valid_q <= 1'b0;
        stg2_valid_q <= 1'b0;
        stg3_valid_q <= 1'b0;
    end
    else begin
        stg0_valid_q <= inport_valid_i;
        stg1_valid_q <= stg0_valid_q;
        stg2_valid_q <= stg1_valid_q;
        stg3_valid_q <= stg2_valid_q;
    end
    stg0_idx_q <= inport_idx_i;
    stg1_idx_q <= stg0_idx_q;
    stg2_idx_q <= stg1_idx_q;
    stg3_idx_q <= stg2_idx_q;
end

//-----------------------------------------------------------------
// Butterfly / accumulator stage
//-----------------------------------------------------------------
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

// signed division by 256 truncating toward zero
function signed [31:0] d256;
    input signed [31:0] x;
    begin
        d256 = x[31] ? ((x + 32'sd255) >>> 8) : (x >>> 8);
    end
endfunction

always @(posedge clk)
begin
    case (stg2_idx_q)
    3'd0:
    begin
        o_t3 <= p2 + p6; // s3
    end
    3'd1:
    begin
        o_s7 <= p1 + p7;
        o_s6 <= p5 + p3;
        o_t0 <= p4;      // s0
    end
    3'd2:
    begin
        o_t0 <= o_t0 + o_t3; // t0
        o_t3 <= o_t0 - o_t3; // t3
        o_t7 <= o_s6 + o_s7;
    end
    3'd3:
    begin
        o_t4 <= (p7 - p1) + (p3 - p5);
    end
    3'd4:
    begin
        o_t0 <= p7 - p1; // s4
        o_s5 <= p3 - p5;
    end
    3'd5:
    begin
        o_t3 <= p6 - p2; // s2
        o_t4 <= p4;      // s1
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
        o_s5 <= d256(o_t6_5 * 32'sd181); // 1/sqrt(2)
        o_s6 <= d256(o_t5_6 * 32'sd181); // 1/sqrt(2)
    end
    endcase
end

//-----------------------------------------------------------------
// Output stage
//-----------------------------------------------------------------
reg signed [BO_W-1:0] block_out [0:7];
reg signed [BO_W-1:0] block_out_tmp;

always @(posedge clk)
if (stg3_valid_q)
begin
    if (stg3_idx_q == 3'd3)
    begin
        block_out[0]   <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp  <= ((o_t0 - o_t7) >>> OUT_SHIFT); // block_out[7]
        block_out[3]   <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4]   <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end

    if (stg3_idx_q == 3'd6)
        block_out[7] <= block_out_tmp;

    if (stg3_idx_q == 3'd7)
    begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;

always @(posedge clk)
if (rst)
    valid_q <= 8'b0;
else if (img_start_i)
    valid_q <= 8'b0;
else
    valid_q <= {valid_q[6:0], stg3_valid_q};

reg [5:0] ptr_q;

always @(posedge clk)
if (rst)
    ptr_q <= 6'd0;
else if (img_start_i)
    ptr_q <= 6'd0;
else if (outport_valid_o)
    ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];

generate
if (SWAP_IDX) begin : idx_swap
    assign outport_idx_o = {ptr_q[2:0], ptr_q[5:3]};
end else begin : idx_direct
    assign outport_idx_o = ptr_q;
end
endgenerate

endmodule

//-----------------------------------------------------------------
// Double-buffered corner turn (four 32-entry RAMs, one-hot read)
//-----------------------------------------------------------------
module jpeg_idct_transpose #(
     parameter DATA_W = 21
)
(
    // Inputs
     input                    clk
    ,input                    rst
    ,input                    img_start_i
    ,input                    inport_valid_i
    ,input  [DATA_W-1:0]      inport_data_i
    ,input  [  5:0]           inport_idx_i
    ,input                    outport_ready_i

    // Outputs
    ,output                   inport_accept_o
    ,output                   outport_valid_o
    ,output [DATA_W-1:0]      outport_data0_o
    ,output [DATA_W-1:0]      outport_data1_o
    ,output [DATA_W-1:0]      outport_data2_o
    ,output [DATA_W-1:0]      outport_data3_o
    ,output [  2:0]           outport_idx_o
);

reg         block_wr_q;
reg         block_rd_q;
reg [ 5:0]  rd_idx_q;
reg [ 3:0]  rd_addr_q;

wire [DATA_W-1:0] outport_data0_w, outport_data1_w, outport_data2_w, outport_data3_w;

wire [ 4:0] wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

wire [ 4:0] rd_addr5_w = {block_rd_q, rd_addr_q};

// registered one-hot read enable, shared by all four RAMs: keeps the
// address decode off the RAM->IDCT data path (fanout 4 per enable)
reg [31:0] rd_en_q;

always @(posedge clk)
begin
    case (rd_addr5_w[2:0])
    3'd0: rd_en_q[ 7:0] <= 8'b00000001;
    3'd1: rd_en_q[ 7:0] <= 8'b00000010;
    3'd2: rd_en_q[ 7:0] <= 8'b00000100;
    3'd3: rd_en_q[ 7:0] <= 8'b00001000;
    3'd4: rd_en_q[ 7:0] <= 8'b00010000;
    3'd5: rd_en_q[ 7:0] <= 8'b00100000;
    3'd6: rd_en_q[ 7:0] <= 8'b01000000;
    default: rd_en_q[ 7:0] <= 8'b10000000;
    endcase
    rd_en_q[15:8]  <= 8'b0;
    rd_en_q[23:16] <= {8{rd_addr5_w[3]}};
    rd_en_q[31:24] <= {8{rd_addr5_w[4]}};
end

jpeg_idct_transpose_ram #(.DATA_W(DATA_W))
u_ram0
(
     .clk      (clk)
    ,.addr0_i  (wr_ptr_w)
    ,.data0_i  (inport_data_i)
    ,.wr0_i    (wr0_w)
    ,.addr1_i  (rd_addr5_w)
    ,.data1_o  (outport_data0_w)
);

jpeg_idct_transpose_ram #(.DATA_W(DATA_W))
u_ram1
(
     .clk      (clk)
    ,.addr0_i  (wr_ptr_w)
    ,.data0_i  (inport_data_i)
    ,.wr0_i    (wr1_w)
    ,.addr1_i  (rd_addr5_w)
    ,.data1_o  (outport_data1_w)
);

jpeg_idct_transpose_ram #(.DATA_W(DATA_W))
u_ram2
(
     .clk      (clk)
    ,.addr0_i  (wr_ptr_w)
    ,.data0_i  (inport_data_i)
    ,.wr0_i    (wr2_w)
    ,.addr1_i  (rd_addr5_w)
    ,.data1_o  (outport_data2_w)
);

jpeg_idct_transpose_ram #(.DATA_W(DATA_W))
u_ram3
(
     .clk      (clk)
    ,.addr0_i  (wr_ptr_w)
    ,.data0_i  (inport_data_i)
    ,.wr0_i    (wr3_w)
    ,.addr1_i  (rd_addr5_w)
    ,.data1_o  (outport_data3_w)
);

//-----------------------------------------------------------------
// Input Buffer
//-----------------------------------------------------------------
reg [1:0] block_ready_q;

always @(posedge clk)
if (rst || img_start_i)
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

always @(posedge clk)
if (rst)
    state_q <= STATE_IDLE;
else
    state_q <= next_state_r;

always @(posedge clk)
if (rst || img_start_i)
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

assign outport_valid_o = valid_d_q;
assign outport_idx_o   = idx_d_q;
assign outport_data0_o = data0_d_q;
assign outport_data1_o = data1_d_q;
assign outport_data2_o = data2_d_q;
assign outport_data3_o = data3_d_q;

// output pipeline stage: breaks the RAM read -> y-pass input adder
// combinational path (the stream is shifted by one cycle as a whole)
reg             valid_d_q;
reg [ 2:0]      idx_d_q;
reg [DATA_W-1:0] data0_d_q, data1_d_q, data2_d_q, data3_d_q;

always @(posedge clk)
if (rst)
    valid_d_q <= 1'b0;
else
begin
    valid_d_q   <= (state_q == STATE_ACTIVE);
    idx_d_q     <= rd_idx_q[2:0];
    data0_d_q   <= outport_data0_w;
    data1_d_q   <= outport_data1_w;
    data2_d_q   <= outport_data2_w;
    data3_d_q   <= outport_data3_w;
end

endmodule

//-----------------------------------------------------------------
// 32-entry RAM: sync write + one-hot decoded registered read
//-----------------------------------------------------------------
module jpeg_idct_transpose_ram #(
     parameter DATA_W = 21
)
(
    // Inputs
     input           clk
    ,input  [  4:0]  addr0_i
    ,input  [DATA_W-1:0]  data0_i
    ,input           wr0_i
    ,input  [  4:0]  addr1_i

    // Outputs
    ,output [DATA_W-1:0]  data1_o
);

reg [DATA_W-1:0] ram [0:31];
reg [DATA_W-1:0] rd_q;

// two-level one-hot decode of the read address (low fanout)
wire [3:0] hh = 4'b0001 << {addr1_i[4], addr1_i[3]};
wire [7:0] kk = 8'b00000001 << addr1_i[2:0];

wire [31:0] en;
genvar gi, gj;
generate
    for (gi = 0; gi < 4; gi = gi + 1)
        for (gj = 0; gj < 8; gj = gj + 1)
            assign en[gi*8+gj] = hh[gi] & kk[gj];
endgenerate

wire [DATA_W-1:0] anded [0:31];
wire [DATA_W-1:0] or4   [0:7];
genvar gq;
generate
    for (gq = 0; gq < 32; gq = gq + 1)
        assign anded[gq] = {DATA_W{en[gq]}} & ram[gq];
    for (gq = 0; gq < 8; gq = gq + 1)
        assign or4[gq] = anded[gq*4+0] | anded[gq*4+1] | anded[gq*4+2] | anded[gq*4+3];
endgenerate

wire [DATA_W-1:0] rd_c;
genvar gg;
generate
    for (gg = 0; gg < DATA_W; gg = gg + 1)
        assign rd_c[gg] = or4[0][gg] | or4[1][gg] | or4[2][gg] | or4[3][gg]
                        | or4[4][gg] | or4[5][gg] | or4[6][gg] | or4[7][gg];
endgenerate

always @(posedge clk)
begin
    if (wr0_i)
        ram[addr0_i] <= data0_i;

    rd_q <= rd_c;
end

assign data1_o = rd_q;

endmodule

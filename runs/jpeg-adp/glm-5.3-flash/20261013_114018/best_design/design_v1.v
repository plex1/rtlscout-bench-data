//=====================================================================
// jpeg_idct_2d : full 8x8 2-D IDCT (core_jpeg behavior, bit-exact)
// Optimized composition:
//   - 1-D x pass (OUT_SHIFT=11): AAN factorized algorithm, 5 time-shared
//     multipliers, 8-wide registered output (one full row per beat)
//   - corner turn: 4 single-port 64x8b RAMs (banked, double buffered)
//   - 1-D y pass (OUT_SHIFT=15)
// All arithmetic signed 32-bit wrap; >>> 11 / >>> 15; (x*181)/256 with
// signed truncating division.
//=====================================================================

//---------------------------------------------------------------------
// 1-D IDCT pass (factorized AAN algorithm, 5 time-shared multipliers)
// Emits one full row (8 values) on a single valid beat.
//---------------------------------------------------------------------
module jpeg_idct_1d #(
    parameter OUT_SHIFT = 11
)(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [15:0]   inport_data0_i
    ,input  [15:0]   inport_data1_i
    ,input  [15:0]   inport_data2_i
    ,input  [15:0]   inport_data3_i
    ,input  [ 2:0]   inport_idx_i
    ,output          outport_valid_o
    ,output [31:0]   outport_data0_o
    ,output [31:0]   outport_data1_o
    ,output [31:0]   outport_data2_o
    ,output [31:0]   outport_data3_o
    ,output [31:0]   outport_data4_o
    ,output [31:0]   outport_data5_o
    ,output [31:0]   outport_data6_o
    ,output [31:0]   outport_data7_o
);

localparam [15:0] C1_16 = 4017; // cos( pi/16) x4096
localparam [15:0] C2_16 = 3784; // cos(2pi/16) x4096
localparam [15:0] C3_16 = 3406; // cos(3pi/16) x4096
localparam [15:0] C4_16 = 2896; // cos(4pi/16) x4096
localparam [15:0] C5_16 = 2276; // cos(5pi/16) x4096
localparam [15:0] C6_16 = 1567; // cos(6pi/16) x4096
localparam [15:0] C7_16 = 799;  // cos(7pi/16) x4096

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a;
reg signed [31:0] mul0_b;
reg signed [31:0] mul1_a;
reg signed [31:0] mul1_b;
reg signed [31:0] mul2_a;
reg signed [31:0] mul2_b;
reg signed [31:0] mul3_a;
reg signed [31:0] mul3_b;
reg signed [31:0] mul4_a;
reg signed [31:0] mul4_b;

always @ (posedge clk)
if (rst)
begin
    i0     <= 32'b0;
    mul0_a <= 32'b0;
    mul0_b <= 32'b0;
    mul1_a <= 32'b0;
    mul1_b <= 32'b0;
    mul2_a <= 32'b0;
    mul2_b <= 32'b0;
    mul3_a <= 32'b0;
    mul3_b <= 32'b0;
    mul4_a <= 32'b0;
    mul4_b <= 32'b0;
end
else
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
reg        out_stg1_valid_q;
reg [2:0]  out_stg1_idx_q;
reg        out_stg2_valid_q;
reg [2:0]  out_stg2_idx_q;

always @ (posedge clk)
if (rst)
begin
    out_stg0_valid_q <= 1'b0;
    out_stg0_idx_q   <= 3'b0;
    out_stg1_valid_q <= 1'b0;
    out_stg1_idx_q   <= 3'b0;
    out_stg2_valid_q <= 1'b0;
    out_stg2_idx_q   <= 3'b0;
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
        o_s5 <= (o_t6_5 * 181) / 256; // 1/sqrt(2)
        o_s6 <= (o_t5_6 * 181) / 256; // 1/sqrt(2)
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

// Output stage: assemble the full row, emit once (idx 0 of the row)
reg signed [31:0] sh_t0_p;
reg signed [31:0] sh_t0_m;
reg signed [31:0] sh_t3_p;
reg signed [31:0] sh_t3_m;
reg signed [31:0] sh_t2_p;
reg signed [31:0] sh_t2_m;
reg signed [31:0] sh_t1_p;
reg signed [31:0] sh_t1_m;

always @ (posedge clk)
if (rst)
begin
    sh_t0_p <= 32'b0;
    sh_t0_m <= 32'b0;
    sh_t3_p <= 32'b0;
    sh_t3_m <= 32'b0;
    sh_t2_p <= 32'b0;
    sh_t2_m <= 32'b0;
    sh_t1_p <= 32'b0;
    sh_t1_m <= 32'b0;
end
else if (out_stg3_valid_q && (out_stg3_idx_q == 3'd3))
begin
    sh_t0_p <= (o_t0 + o_t7) >>> OUT_SHIFT;
    sh_t0_m <= (o_t0 - o_t7) >>> OUT_SHIFT;
    sh_t3_p <= (o_t3 + o_t4) >>> OUT_SHIFT;
    sh_t3_m <= (o_t3 - o_t4) >>> OUT_SHIFT;
end
else if (out_stg3_valid_q && (out_stg3_idx_q == 3'd7))
begin
    sh_t2_p <= (o_t2 + o_s5) >>> OUT_SHIFT;
    sh_t2_m <= (o_t2 - o_s5) >>> OUT_SHIFT;
    sh_t1_p <= (o_t1 + o_s6) >>> OUT_SHIFT;
    sh_t1_m <= (o_t1 - o_s6) >>> OUT_SHIFT;
end

assign outport_valid_o  = out_stg3_valid_q && (out_stg3_idx_q == 3'd0);
assign outport_data0_o  = sh_t0_p;
assign outport_data1_o  = sh_t1_p;
assign outport_data2_o  = sh_t2_p;
assign outport_data3_o  = sh_t3_p;
assign outport_data4_o  = sh_t3_m;
assign outport_data5_o  = sh_t2_m;
assign outport_data6_o  = sh_t1_m;
assign outport_data7_o  = sh_t0_m;

endmodule

//---------------------------------------------------------------------
// Corner turn: 4 single-port RAMs of 64 x 8b
//---------------------------------------------------------------------
module jpeg_idct_ct_ram (
     input           clk
    ,input  [ 5:0]   wr_addr_i
    ,input  [ 7:0]   wr_data_i
    ,input           wr_en_i
    ,input  [ 5:0]   rd_addr_i
    ,output reg [7:0] rd_data_o
);

reg [7:0] mem [0:63];

always @ (posedge clk)
begin
    if (wr_en_i)
        mem[wr_addr_i] <= wr_data_i;
    rd_data_o <= mem[rd_addr_i];
end

endmodule

//=====================================================================
// Top level
//=====================================================================
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

// ---- 1-D x pass -----------------------------------------------------
wire          x_valid_w;
wire [31:0]   x_d0_w, x_d1_w, x_d2_w, x_d3_w;
wire [31:0]   x_d4_w, x_d5_w, x_d6_w, x_d7_w;

jpeg_idct_1d #(.OUT_SHIFT(11)) u_x
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
    ,.outport_data0_o(x_d0_w)
    ,.outport_data1_o(x_d1_w)
    ,.outport_data2_o(x_d2_w)
    ,.outport_data3_o(x_d3_w)
    ,.outport_data4_o(x_d4_w)
    ,.outport_data5_o(x_d5_w)
    ,.outport_data6_o(x_d6_w)
    ,.outport_data7_o(x_d7_w)
);

// ---- corner turn control --------------------------------------------
reg        wr_bank_q;   // bank being written
reg        rd_bank_q;   // bank being read
reg        buf_rdy_q;   // written bank holds a complete block
reg [5:0]  ct_wr_addr_q;
reg [5:0]  ct_rd_row_q;
reg [1:0]  ct_state_q;
reg [1:0]  ct_next_state;

// row parity from the inport schedule: the row repeats every 8 beats
reg [2:0]  wr_row_q;

wire ct_wr_last = x_valid_w && (inport_idx_i == 3'd7);
wire ct_read_this_q = (rd_bank_q == wr_bank_q) & buf_rdy_q;
wire ct_rd_last = ct_read_this_q && (ct_rd_row_q == 6'd56);

wire [6:0] ct_wr_row_w = {1'b0, wr_row_q};

// split each 32-bit value into two 8-bit halves
// lane j: value j -> mem j (low half), mem j+4 (high half)
wire ct_wr0_w = x_valid_w & (ct_wr_row_w[0] == 1'b0);
wire ct_wr1_w = x_valid_w & (ct_wr_row_w[0] == 1'b1);
wire ct_wr2_w = x_valid_w & (ct_wr_row_w[1] == 1'b0);
wire ct_wr3_w = x_valid_w & (ct_wr_row_w[1] == 1'b1);

wire [7:0]  ct_d0_even, ct_d0_pear;
wire [7:0]  ct_d1_even, ct_d1_pear;
wire [7:0]  ct_d2_even, ct_d2_pear;
wire [7:0]  ct_d3_even, ct_d3_pear;

jpeg_idct_ct_ram u_ct0 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d0_w[7:0])
    ,.wr_en_i(ct_wr0_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b000})
    ,.rd_data_o(ct_d0_even)
);
jpeg_idct_ct_ram u_ct1 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d0_w[15:8])
    ,.wr_en_i(ct_wr0_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b100})
    ,.rd_data_o(ct_d0_pear)
);
jpeg_idct_ct_ram u_ct2 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d1_w[7:0])
    ,.wr_en_i(ct_wr1_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b001})
    ,.rd_data_o(ct_d1_even)
);
jpeg_idct_ct_ram u_ct3 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d1_w[15:8])
    ,.wr_en_i(ct_wr1_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b101})
    ,.rd_data_o(ct_d1_pear)
);
jpeg_idct_ct_ram u_ct4 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d2_w[7:0])
    ,.wr_en_i(ct_wr2_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b010})
    ,.rd_data_o(ct_d2_even)
);
jpeg_idct_ct_ram u_ct5 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d2_w[15:8])
    ,.wr_en_i(ct_wr2_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b110})
    ,.rd_data_o(ct_d2_pear)
);
jpeg_idct_ct_ram u_ct6 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d3_w[7:0])
    ,.wr_en_i(ct_wr3_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b011})
    ,.rd_data_o(ct_d3_even)
);
jpeg_idct_ct_ram u_ct7 (
     .clk(clk)
    ,.wr_addr_i(ct_wr_addr_q)
    ,.wr_data_i(x_d3_w[15:8])
    ,.wr_en_i(ct_wr3_w)
    ,.rd_addr_i({ct_rd_row_q[5:3], 3'b111})
    ,.rd_data_o(ct_d3_pear)
);

// write FSM: track row and block boundaries
always @ (posedge clk)
if (rst)
begin
    wr_bank_q <= 1'b0;
    rd_bank_q <= 1'b0;
    buf_rdy_q <= 1'b0;
    wr_row_q  <= 3'd0;
end
else if (img_start_i)
begin
    wr_bank_q <= 1'b0;
    rd_bank_q <= 1'b0;
    buf_rdy_q <= 1'b0;
    wr_row_q  <= 3'd0;
end
else
begin
    if (ct_wr_last)
    begin
        wr_bank_q <= ~wr_bank_q;
        buf_rdy_q <= 1'b1;
        wr_row_q  <= 3'd0;
    end
    else if (x_valid_w)
        wr_row_q <= wr_row_q + 3'd1;

    if (ct_read_this_q && ct_rd_last)
    begin
        rd_bank_q <= ~rd_bank_q;
        buf_rdy_q <= 1'b0;
    end
end

always @ (posedge clk)
if (rst)
    ct_wr_addr_q <= 6'd0;
else if (img_start_i)
    ct_wr_addr_q <= 6'd0;
else if (x_valid_w)
    ct_wr_addr_q <= ct_wr_addr_q + 6'd1;

// read FSM
always @ *
begin
    ct_next_state = ct_state_q;
    case (ct_state_q)
    2'd0: if (ct_read_this_q) ct_next_state = 2'd1;
    2'd1: ct_next_state = 2'd2;
    2'd2: if (ct_rd_last) ct_next_state = 2'd0;
    default: ;
    endcase
end

always @ (posedge clk)
if (rst)
    ct_state_q <= 2'd0;
else if (img_start_i)
    ct_state_q <= 2'd0;
else
    ct_state_q <= ct_next_state;

always @ (posedge clk)
if (rst)
    ct_rd_row_q <= 6'd0;
else if (img_start_i)
    ct_rd_row_q <= 6'd0;
else if (ct_state_q == 2'd2)
    ct_rd_row_q <= ct_rd_last ? 6'd0 : (ct_rd_row_q + 6'd8);

// corner-turn output register
reg ct_out_valid_q;
reg [2:0] ct_out_row_q;

always @ (posedge clk)
if (rst)
begin
    ct_out_valid_q <= 1'b0;
    ct_out_row_q   <= 3'd0;
end
else if (img_start_i)
begin
    ct_out_valid_q <= 1'b0;
    ct_out_row_q   <= 3'd0;
end
else
begin
    ct_out_valid_q <= (ct_state_q == 2'd2);
    ct_out_row_q   <= ct_rd_row_q[5:3];
end

// ---- y-pass ----------------------------------------------------------
reg y_valid_q;
reg [2:0] y_idx_q;
reg [31:0] y_d0_q, y_d1_q, y_d2_q, y_d3_q;
reg [31:0] y_d4_q, y_d5_q, y_d6_q, y_d7_q;

always @ (posedge clk)
if (rst)
begin
    y_valid_q <= 1'b0;
    y_idx_q   <= 3'd0;
    y_d0_q <= 32'b0; y_d1_q <= 32'b0; y_d2_q <= 32'b0; y_d3_q <= 32'b0;
    y_d4_q <= 32'b0; y_d5_q <= 32'b0; y_d6_q <= 32'b0; y_d7_q <= 32'b0;
end
else
begin
    y_valid_q <= ct_out_valid_q;
    y_idx_q   <= ct_out_row_q;
    y_d0_q <= {{24{ct_d0_even[7]}}, ct_d0_even};
    y_d1_q <= {{24{ct_d1_even[7]}}, ct_d1_even};
    y_d2_q <= {{24{ct_d2_even[7]}}, ct_d2_even};
    y_d3_q <= {{24{ct_d3_even[7]}}, ct_d3_even};
    y_d4_q <= {{24{ct_d0_pear[7]}}, ct_d0_pear};
    y_d5_q <= {{24{ct_d1_pear[7]}}, ct_d1_pear};
    y_d6_q <= {{24{ct_d2_pear[7]}}, ct_d2_pear};
    y_d7_q <= {{24{ct_d3_pear[7]}}, ct_d3_pear};
end

wire y_valid_w;
wire [31:0] y0_w, y1_w, y2_w, y3_w, y4_w, y5_w, y6_w, y7_w;

jpeg_idct_1d #(.OUT_SHIFT(15)) u_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(y_valid_q)
    ,.inport_data0_i(y_d0_q[15:0])
    ,.inport_data1_i(y_d1_q[15:0])
    ,.inport_data2_i(y_d2_q[15:0])
    ,.inport_data3_i(y_d3_q[15:0])
    ,.inport_idx_i(y_idx_q)
    ,.outport_valid_o(y_valid_w)
    ,.outport_data0_o(y0_w)
    ,.outport_data1_o(y1_w)
    ,.outport_data2_o(y2_w)
    ,.outport_data3_o(y3_w)
    ,.outport_data4_o(y4_w)
    ,.outport_data5_o(y5_w)
    ,.outport_data6_o(y6_w)
    ,.outport_data7_o(y7_w)
);

// ---- output register -------------------------------------------------
reg out_valid_q;
reg [5:0] out_idx_q;
reg [31:0] out_d0_q, out_d1_q, out_d2_q, out_d3_q;
reg [31:0] out_d4_q, out_d5_q, out_d6_q, out_d7_q;

always @ (posedge clk)
if (rst)
begin
    out_valid_q <= 1'b0;
    out_idx_q   <= 6'd0;
    out_d0_q <= 32'b0; out_d1_q <= 32'b0; out_d2_q <= 32'b0; out_d3_q <= 32'b0;
    out_d4_q <= 32'b0; out_d5_q <= 32'b0; out_d6_q <= 32'b0; out_d7_q <= 32'b0;
end
else if (img_start_i)
begin
    out_valid_q <= 1'b0;
    out_idx_q   <= 6'd0;
end
else
begin
    out_valid_q <= y_valid_w;
    out_idx_q   <= y_valid_w ? {out_idx_q[2:0] + 3'd1, out_idx_q[5:3]} : out_idx_q;
    out_d0_q <= y0_w;
    out_d1_q <= y1_w;
    out_d2_q <= y2_w;
    out_d3_q <= y3_w;
    out_d4_q <= y4_w;
    out_d5_q <= y5_w;
    out_d6_q <= y6_w;
    out_d7_q <= y7_w;
end

reg [31:0] out_mux_q;

always @ (*)
begin
    case (out_idx_q[2:0])
    3'd0: out_mux_q = out_d0_q;
    3'd1: out_mux_q = out_d1_q;
    3'd2: out_mux_q = out_d2_q;
    3'd3: out_mux_q = out_d3_q;
    3'd4: out_mux_q = out_d4_q;
    3'd5: out_mux_q = out_d5_q;
    3'd6: out_mux_q = out_d6_q;
    default: out_mux_q = out_d7_q;
    endcase
end

assign outport_data_o = out_mux_q;
assign outport_valid_o = out_valid_q;
assign outport_idx_o = out_idx_q;

endmodule

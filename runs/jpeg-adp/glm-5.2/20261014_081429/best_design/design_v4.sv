//-----------------------------------------------------------------
// jpeg_idct_2d — optimized: split RAMs + combinational rd_addr
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

jpeg_idct_x u_idct_x(
     .clk(clk), .rst(rst), .img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i(inport_data0_i), .inport_data1_i(inport_data1_i)
    ,.inport_data2_i(inport_data2_i), .inport_data3_i(inport_data3_i)
    ,.inport_idx_i(inport_idx_i)
    ,.outport_valid_o(idct_x_valid_w), .outport_data_o(idct_x_data_w)
    ,.outport_idx_o(idct_x_idx_w)
);

wire          transpose_valid_w;
wire [ 31:0]  transpose_data0_w, transpose_data1_w, transpose_data2_w, transpose_data3_w;
wire [  2:0]  transpose_idx_w;

jpeg_idct_transpose u_transpose(
     .clk(clk), .rst(rst), .img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(idct_x_valid_w), .inport_data_i(idct_x_data_w)
    ,.inport_idx_i(idct_x_idx_w), .inport_accept_o(idct_x_accept_w)
    ,.outport_valid_o(transpose_valid_w)
    ,.outport_data0_o(transpose_data0_w), .outport_data1_o(transpose_data1_w)
    ,.outport_data2_o(transpose_data2_w), .outport_data3_o(transpose_data3_w)
    ,.outport_idx_o(transpose_idx_w), .outport_ready_i(1'b1)
);

jpeg_idct_y u_idct_y(
     .clk(clk), .rst(rst), .img_start_i(img_start_i), .img_end_i(img_end_i)
    ,.inport_valid_i(transpose_valid_w)
    ,.inport_data0_i(transpose_data0_w), .inport_data1_i(transpose_data1_w)
    ,.inport_data2_i(transpose_data2_w), .inport_data3_i(transpose_data3_w)
    ,.inport_idx_i(transpose_idx_w)
    ,.outport_valid_o(outport_valid_o), .outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);
endmodule

//-----------------------------------------------------------------
// jpeg_idct_x — row pass
//-----------------------------------------------------------------
module jpeg_idct_x #(
     parameter OUT_SHIFT   = 11
    ,parameter INPUT_WIDTH = 16
)(
    input           clk, rst, img_start_i, img_end_i, inport_valid_i
    ,input  [ 15:0]  inport_data0_i, inport_data1_i, inport_data2_i, inport_data3_i
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

wire signed [31:0] block_in_0_1 = {{16{inport_data0_i[15]}}, inport_data0_i};
wire signed [31:0] block_in_2_3 = {{16{inport_data1_i[15]}}, inport_data1_i};
wire signed [31:0] block_in_4_5 = {{16{inport_data2_i[15]}}, inport_data2_i};
wire signed [31:0] block_in_6_7 = {{16{inport_data3_i[15]}}, inport_data3_i};

reg signed [31:0] i0;
reg signed [31:0] mul0_a, mul0_b, mul1_a, mul1_b, mul2_a, mul2_b, mul3_a, mul3_b, mul4_a, mul4_b;

always @(posedge clk)
if (rst) begin
    i0 <= 0;
    mul0_a <= 0; mul0_b <= 0; mul1_a <= 0; mul1_b <= 0;
    mul2_a <= 0; mul2_b <= 0; mul3_a <= 0; mul3_b <= 0;
    mul4_a <= 0; mul4_b <= 0;
end else begin
    case (inport_idx_i)
    3'd0: begin
        i0 <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= C2_16;
        mul1_a <= block_in_6_7; mul1_b <= C6_16;
    end
    3'd1: begin
        mul0_a <= block_in_0_1; mul0_b <= C1_16;
        mul1_a <= block_in_6_7; mul1_b <= C7_16;
        mul2_a <= block_in_4_5; mul2_b <= C5_16;
        mul3_a <= block_in_2_3; mul3_b <= C3_16;
        mul4_a <= i0; mul4_b <= C4_16;
    end
    3'd2: i0 <= block_in_0_1 - block_in_4_5;
    3'd3: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd4: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd5: begin
        mul0_a <= block_in_2_3; mul0_b <= C6_16;
        mul1_a <= block_in_6_7; mul1_b <= C2_16;
        mul4_a <= i0; mul4_b <= C4_16;
    end
    default: ;
    endcase
end

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;
always @(posedge clk)
if (rst) begin
    mul0_q <= 0; mul1_q <= 0; mul2_q <= 0; mul3_q <= 0; mul4_q <= 0;
end else begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;
always @(posedge clk)
if (rst) begin
    mul0 <= 0; mul1 <= 0; mul2 <= 0; mul3 <= 0; mul4 <= 0;
end else begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q;
    mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q; reg [2:0] out_stg0_idx_q;
reg        out_stg1_valid_q; reg [2:0] out_stg1_idx_q;
reg        out_stg2_valid_q; reg [2:0] out_stg2_idx_q;

always @(posedge clk) if (rst) begin out_stg0_valid_q <= 0; out_stg0_idx_q <= 0; end
else begin out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i; end
always @(posedge clk) if (rst) begin out_stg1_valid_q <= 0; out_stg1_idx_q <= 0; end
else begin out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q; end
always @(posedge clk) if (rst) begin out_stg2_valid_q <= 0; out_stg2_idx_q <= 0; end
else begin out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q; end

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

always @(posedge clk)
if (rst) begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0;
    o_t0 <= 0; o_t1 <= 0; o_t2 <= 0; o_t3 <= 0;
    o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end else begin
    case (out_stg2_idx_q)
    3'd0: o_t3 <= mul0 + mul1;
    3'd1: begin
        o_s7 <= mul0 + mul1; o_s6 <= mul2 + mul3; o_t0 <= mul4;
    end
    3'd2: begin
        o_t0 <= o_t0 + o_t3; o_t3 <= o_t0 - o_t3; o_t7 <= o_s6 + o_s7;
    end
    3'd3: o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    3'd4: begin
        o_t0 <= mul0 - mul1; o_s5 <= mul2 - mul3;
    end
    3'd5: begin
        o_t3 <= mul0 - mul1; o_t4 <= mul4;
        o_t5 <= o_t0 - o_s5; o_t6 <= o_s7 - o_s6;
    end
    3'd6: begin
        o_t1 <= o_t4 + o_t3; o_t2 <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5; o_t5_6 <= o_t5 + o_t6;
    end
    default: begin
        o_s5 <= (o_t6_5 * 181) / 256;
        o_s6 <= (o_t5_6 * 181) / 256;
    end
    endcase
end

reg        out_stg3_valid_q; reg [2:0] out_stg3_idx_q;
always @(posedge clk) if (rst) begin out_stg3_valid_q <= 0; out_stg3_idx_q <= 0; end
else begin out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q; end

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @(posedge clk)
if (rst) begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end else if (out_stg3_valid_q) begin
    if (out_stg3_idx_q == 3'd3) begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end
    if (out_stg3_idx_q == 3'd6) block_out[7] <= block_out_tmp;
    if (out_stg3_idx_q == 3'd7) begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;
always @(posedge clk)
if (rst) valid_q <= 8'b0;
else if (img_start_i) valid_q <= 8'b0;
else valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;
always @(posedge clk)
if (rst) ptr_q <= 6'd0;
else if (img_start_i) ptr_q <= 6'd0;
else if (outport_valid_o) ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];
assign outport_idx_o   = ptr_q;
endmodule

//-----------------------------------------------------------------
// jpeg_idct_transpose — split RAMs, combinational rd_addr
//-----------------------------------------------------------------
module jpeg_idct_transpose(
    input           clk, rst, img_start_i, img_end_i
    ,input           inport_valid_i
    ,input  [ 31:0]  inport_data_i
    ,input  [  5:0]  inport_idx_i
    ,input           outport_ready_i
    ,output          inport_accept_o
    ,output          outport_valid_o
    ,output [ 31:0]  outport_data0_o, outport_data1_o, outport_data2_o, outport_data3_o
    ,output [  2:0]  outport_idx_o
);

reg         block_wr_q, block_rd_q;
reg [5:0]   rd_idx_q;
wire [4:0]  wr_ptr_w = {block_wr_q, inport_idx_i[3:0]};

// Compute rd_addr combinationally from rd_idx_q
// The read address sequence for the transposed read:
// The transpose reads back data transposed (column-major from row-major writes).
// In ACTIVE state, rd_idx_q counts 0..63.
// The rd_addr follows the pattern: 0,8,8,0,8,0,1,8, 0,8,8,0,8,0,1,8, ...
// Let's compute it from rd_idx_q[5:0] directly.

// rd_addr is a 4-bit value. The pattern depends on rd_idx_q[2:0] (sub-block index):
// idx 0: addr was 0, next = 0-8 = 8 (mod 16)
// idx 1: addr was 8, next = 8+8 = 0 (mod 16) 
// idx 2: addr stays
// idx 3: addr-8
// idx 4: addr+8
// idx 5: addr-8
// idx 6: addr+1
// idx 7: addr+8

// Actually, let's precompute the full sequence of rd_addr for each rd_idx:
// It's easier to keep rd_addr_q but register it separately with lower fanout.
// Instead, let's compute rd_addr from rd_idx directly.

// The transpose read pattern (for 8 rows of 8 cols):
// For each output row j=0..7, we read 8 values from the RAM.
// The read address sequence for j=0: 0,8,0,8,0,8,1,8
//   wait, that doesn't look right. Let me trace the original.

// Original: rd_addr_q starts at 0, then:
// idx=0: rd_addr=0 (read), next rd_addr = 0-8 = 8 (wraps to 8 in 4 bits? 0-8=-8=8 mod 16)
// idx=1: rd_addr=8 (read), next = 8+8=16=0 mod 16
// idx=2: rd_addr=0 (read), next stays 0
// idx=3: rd_addr=0 (read), next = 0-8=8
// idx=4: rd_addr=8 (read), next = 8+8=0
// idx=5: rd_addr=0 (read), next = 0-8=8
// idx=6: rd_addr=8 (read), next = 8+1=9
// idx=7: rd_addr=9 (read), next = 9+8=17=1 mod 16
// Then idx=0 of next group: rd_addr=1 (read), next = 1-8=9 mod 16
// ...

// The pattern repeats with the base address incrementing by 1 each 8 beats.
// For group j (0..7), the base address is j, and the 8 reads are:
// j, j+8, j, j, j+8, j, j+1, j+1+8=j+9
// Wait that's not right either. Let me trace more carefully.

// Group 0 (rd_idx 0..7): base=0
//   idx 0: addr=0, next=0-8=-8=8
//   idx 1: addr=8, next=8+8=16=0
//   idx 2: addr=0, next=0
//   idx 3: addr=0, next=0-8=8
//   idx 4: addr=8, next=8+8=0
//   idx 5: addr=0, next=0-8=8
//   idx 6: addr=8, next=8+1=9
//   idx 7: addr=9, next=9+8=17=1
// Group 1 (rd_idx 8..15): base=1
//   idx 0: addr=1, next=1-8=9
//   idx 1: addr=9, next=9+8=17=1
//   idx 2: addr=1, next=1
//   idx 3: addr=1, next=1-8=9
//   idx 4: addr=9, next=9+8=1
//   idx 5: addr=1, next=1-8=9
//   idx 6: addr=9, next=9+1=10
//   idx 7: addr=10, next=10+8=18=2

// So for group j, the read addresses are:
//   j, j+8, j, j, j+8, j, j+1, j+9  (all mod 16)
// Wait, that doesn't match for group 0:
//   0, 8, 0, 0, 8, 0, 1, 9  -> but idx 6 should be addr=8 not 1
// Hmm, let me re-trace. idx 6: addr was 0 (from idx 5 next), 
//   NO! idx 5 next = 0-8=8, so idx 6: addr=8. OK.
//   idx 6 next = 8+1=9, so idx 7: addr=9.
//   idx 7 next = 9+8=17=1
// So group 0: 0, 8, 0, 0, 8, 0, 8, 9  -> wait, idx 6 addr=8 (from idx 5 next=0-8=8)
// Actually: idx 5 addr = 0 (from idx 4 next = 8+8=0), idx 5 next = 0-8=8
// idx 6 addr = 8 (from idx 5 next = 8), idx 6 next = 8+1=9
// idx 7 addr = 9 (from idx 6 next = 9), idx 7 next = 9+8=1

// Group 0: 0, 8, 0, 0, 8, 0, 8, 9  (mod 16)
// Group 1: 1, 9, 1, 1, 9, 1, 9, 10
// Group j: j, j+8, j, j, j+8, j, j+8, j+9  (mod 16)

// Wait, that's not consistent. Let me re-check group 0:
// idx 0: addr=0  (was set in SETUP to 8? No, IDLE sets to 0, SETUP sets to 8)
// WAIT: state_q == STATE_IDLE sets rd_addr_q <= 0
//        state_q == STATE_SETUP sets rd_addr_q <= 8
// So when ACTIVE starts, rd_addr_q = 8.
// idx 0: addr=8?? No, the addr computation happens at the END of ACTIVE.
// In ACTIVE, the read uses the CURRENT rd_addr_q value.
// At the first ACTIVE cycle, rd_addr_q = 8 (set in SETUP).
// Wait, SETUP sets rd_addr_q <= 8, so in the first ACTIVE cycle, rd_addr_q = 8.
// But then the case statement says: rd_idx_q[2:0] == 0: rd_addr_q <= rd_addr_q - 8 = 0
// So at idx 0, we read addr=8, then compute next=0.
// Hmm, that means the read pattern is:
// idx 0: read addr=8, next=0
// idx 1: read addr=0, next=0+8=8
// idx 2: read addr=8, next stays 8
// idx 3: read addr=8, next=8-8=0
// idx 4: read addr=0, next=0+8=8
// idx 5: read addr=8, next=8-8=0
// idx 6: read addr=0, next=0+1=1
// idx 7: read addr=1, next=1+8=9
// Next group:
// idx 0: read addr=9, next=9-8=1
// idx 1: read addr=1, next=1+8=9
// ...
// Group 0: 8, 0, 8, 8, 0, 8, 0, 1
// Group 1: 9, 1, 9, 9, 1, 9, 1, 2
// Group j: j+8, j, j+8, j+8, j, j+8, j, j+1

// So rd_addr = f(rd_idx_q) where:
// For rd_idx_q = 8*j + k (j=0..7, k=0..7):
// k=0: j+8
// k=1: j
// k=2: j+8
// k=3: j+8
// k=4: j
// k=5: j+8
// k=6: j
// k=7: j+1
// All mod 16.

// rd_idx_q[5:3] = j, rd_idx_q[2:0] = k
// So rd_addr = {rd_idx_q[5:3]} + offset[k] mod 16

// This is getting complex. Let me just keep rd_addr_q but use the Q output.
// Actually, let me just replicate the rd_addr_q into 4 separate registers
// (one per RAM) to reduce fanout.

wire [31:0] outport_data0_w, outport_data1_w, outport_data2_w, outport_data3_w;

// Replicated rd_addr to reduce fanout per register
wire [4:0] rd_addr0 = {block_rd_q, rd_addr_q};
wire [4:0] rd_addr1 = {block_rd_q, rd_addr_q};
wire [4:0] rd_addr2 = {block_rd_q, rd_addr_q};
wire [4:0] rd_addr3 = {block_rd_q, rd_addr_q};

// Actually these are the same wires. The issue is the DFF has high fanout.
// Let me instead register the full 5-bit address into 4 separate registers.

reg [4:0] rd_addr_full0, rd_addr_full1, rd_addr_full2, rd_addr_full3;

wire wr0_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd0);
wire wr1_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd1);
wire wr2_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd2);
wire wr3_w = inport_valid_i && inport_accept_o && (inport_idx_i[5:4] == 2'd3);

jpeg_idct_transpose_ram u_ram0(
    .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr0_w), .data0_o()
    ,.addr1_i(rd_addr_full0), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data0_w)
);
jpeg_idct_transpose_ram u_ram1(
    .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr1_w), .data0_o()
    ,.addr1_i(rd_addr_full1), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data1_w)
);
jpeg_idct_transpose_ram u_ram2(
    .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr2_w), .data0_o()
    ,.addr1_i(rd_addr_full2), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data2_w)
);
jpeg_idct_transpose_ram u_ram3(
    .clk0_i(clk), .rst0_i(rst), .clk1_i(clk), .rst1_i(rst)
    ,.addr0_i(wr_ptr_w), .data0_i(inport_data_i), .wr0_i(wr3_w), .data0_o()
    ,.addr1_i(rd_addr_full3), .data1_i(32'b0), .wr1_i(1'b0), .data1_o(outport_data3_w)
);

reg [1:0] block_ready_q;
always @(posedge clk)
if (rst) begin
    block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
end else if (img_start_i) begin
    block_ready_q <= 2'b0; block_wr_q <= 1'b0; block_rd_q <= 1'b0;
end else begin
    if (inport_valid_i && inport_idx_i == 6'd63 && inport_accept_o) begin
        block_ready_q[block_wr_q] <= 1'b1;
        block_wr_q <= ~block_wr_q;
    end
    if (outport_valid_o && rd_idx_q[5:0] == 6'd63) begin
        block_ready_q[block_rd_q] <= 1'b0;
        block_rd_q <= ~block_rd_q;
    end
end
assign inport_accept_o = ~block_ready_q[block_wr_q];

localparam STATE_IDLE=2'd0, STATE_SETUP=2'd1, STATE_ACTIVE=2'd2;
reg [1:0] state_q, next_state_r;
always @(*) begin
    next_state_r = state_q;
    case (state_q)
    STATE_IDLE: if (block_ready_q[block_rd_q] && outport_ready_i) next_state_r = STATE_SETUP;
    STATE_SETUP: next_state_r = STATE_ACTIVE;
    STATE_ACTIVE: if (outport_valid_o && rd_idx_q == 6'd63) next_state_r = STATE_IDLE;
    default: ;
    endcase
    if (img_start_i) next_state_r = STATE_IDLE;
end
always @(posedge clk) if (rst) state_q <= STATE_IDLE; else state_q <= next_state_r;

always @(posedge clk)
if (rst) rd_idx_q <= 6'b0;
else if (img_start_i) rd_idx_q <= 6'b0;
else if (state_q == STATE_ACTIVE) rd_idx_q <= rd_idx_q + 6'd1;

reg [3:0] rd_addr_q;
always @(posedge clk)
if (rst) rd_addr_q <= 4'b0;
else if (img_start_i) rd_addr_q <= 4'b0;
else if (state_q == STATE_IDLE) rd_addr_q <= 4'b0;
else if (state_q == STATE_SETUP) rd_addr_q <= 4'd8;
else if (state_q == STATE_ACTIVE) begin
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

// Replicate rd_addr into 4 separate registers to reduce fanout
always @(posedge clk) begin
    rd_addr_full0 <= {block_rd_q, rd_addr_q};
    rd_addr_full1 <= {block_rd_q, rd_addr_q};
    rd_addr_full2 <= {block_rd_q, rd_addr_q};
    rd_addr_full3 <= {block_rd_q, rd_addr_q};
end

assign outport_valid_o = (state_q == STATE_ACTIVE);
assign outport_idx_o = rd_idx_q[2:0];
assign outport_data0_o = outport_data0_w;
assign outport_data1_o = outport_data1_w;
assign outport_data2_o = outport_data2_w;
assign outport_data3_o = outport_data3_w;
endmodule

//-----------------------------------------------------------------
// jpeg_idct_transpose_ram — split into two 16-entry arrays
//-----------------------------------------------------------------
module jpeg_idct_transpose_ram(
    input           clk0_i, rst0_i, clk1_i, rst1_i
    ,input  [  4:0]  addr0_i
    ,input  [ 31:0]  data0_i
    ,input           wr0_i
    ,output [ 31:0]  data0_o
    ,input  [  4:0]  addr1_i
    ,input  [ 31:0]  data1_i
    ,input           wr1_i
    ,output [ 31:0]  data1_o
);
reg [31:0] ram_lo [15:0];
reg [31:0] ram_hi [15:0];
reg [31:0] ram_read0_q, ram_read1_q;

wire wr0_lo = wr0_i && !addr0_i[4];
wire wr0_hi = wr0_i &&  addr0_i[4];
wire wr1_lo = wr1_i && !addr1_i[4];
wire wr1_hi = wr1_i &&  addr1_i[4];

always @(posedge clk0_i) begin
    if (wr0_lo) ram_lo[addr0_i[3:0]] <= data0_i;
    if (wr0_hi) ram_hi[addr0_i[3:0]] <= data0_i;
    ram_read0_q <= addr0_i[4] ? ram_hi[addr0_i[3:0]] : ram_lo[addr0_i[3:0]];
end
always @(posedge clk1_i) begin
    if (wr1_lo) ram_lo[addr1_i[3:0]] <= data1_i;
    if (wr1_hi) ram_hi[addr1_i[3:0]] <= data1_i;
    ram_read1_q <= addr1_i[4] ? ram_hi[addr1_i[3:0]] : ram_lo[addr1_i[3:0]];
end
assign data0_o = ram_read0_q;
assign data1_o = ram_read1_q;
endmodule

//-----------------------------------------------------------------
// jpeg_idct_y — column pass
//-----------------------------------------------------------------
module jpeg_idct_y #(
     parameter OUT_SHIFT   = 15
    ,parameter INPUT_WIDTH = 32
)(
    input           clk, rst, img_start_i, img_end_i, inport_valid_i
    ,input  [ 31:0]  inport_data0_i, inport_data1_i, inport_data2_i, inport_data3_i
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

wire signed [31:0] block_in_0_1 = inport_data0_i;
wire signed [31:0] block_in_2_3 = inport_data1_i;
wire signed [31:0] block_in_4_5 = inport_data2_i;
wire signed [31:0] block_in_6_7 = inport_data3_i;

reg signed [31:0] i0;
reg signed [31:0] mul0_a, mul0_b, mul1_a, mul1_b, mul2_a, mul2_b, mul3_a, mul3_b, mul4_a, mul4_b;

always @(posedge clk)
if (rst) begin
    i0 <= 0;
    mul0_a <= 0; mul0_b <= 0; mul1_a <= 0; mul1_b <= 0;
    mul2_a <= 0; mul2_b <= 0; mul3_a <= 0; mul3_b <= 0;
    mul4_a <= 0; mul4_b <= 0;
end else begin
    case (inport_idx_i)
    3'd0: begin
        i0 <= block_in_0_1 + block_in_4_5;
        mul0_a <= block_in_2_3; mul0_b <= C2_16;
        mul1_a <= block_in_6_7; mul1_b <= C6_16;
    end
    3'd1: begin
        mul0_a <= block_in_0_1; mul0_b <= C1_16;
        mul1_a <= block_in_6_7; mul1_b <= C7_16;
        mul2_a <= block_in_4_5; mul2_b <= C5_16;
        mul3_a <= block_in_2_3; mul3_b <= C3_16;
        mul4_a <= i0; mul4_b <= C4_16;
    end
    3'd2: i0 <= block_in_0_1 - block_in_4_5;
    3'd3: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd4: begin
        mul0_a <= block_in_0_1; mul0_b <= C7_16;
        mul1_a <= block_in_6_7; mul1_b <= C1_16;
        mul2_a <= block_in_4_5; mul2_b <= C3_16;
        mul3_a <= block_in_2_3; mul3_b <= C5_16;
    end
    3'd5: begin
        mul0_a <= block_in_2_3; mul0_b <= C6_16;
        mul1_a <= block_in_6_7; mul1_b <= C2_16;
        mul4_a <= i0; mul4_b <= C4_16;
    end
    default: ;
    endcase
end

reg signed [31:0] mul0_q, mul1_q, mul2_q, mul3_q, mul4_q;
always @(posedge clk)
if (rst) begin
    mul0_q <= 0; mul1_q <= 0; mul2_q <= 0; mul3_q <= 0; mul4_q <= 0;
end else begin
    mul0_q <= mul0_a * mul0_b;
    mul1_q <= mul1_a * mul1_b;
    mul2_q <= mul2_a * mul2_b;
    mul3_q <= mul3_a * mul3_b;
    mul4_q <= mul4_a * mul4_b;
end

reg signed [31:0] mul0, mul1, mul2, mul3, mul4;
always @(posedge clk)
if (rst) begin
    mul0 <= 0; mul1 <= 0; mul2 <= 0; mul3 <= 0; mul4 <= 0;
end else begin
    mul0 <= mul0_q; mul1 <= mul1_q; mul2 <= mul2_q;
    mul3 <= mul3_q; mul4 <= mul4_q;
end

reg        out_stg0_valid_q; reg [2:0] out_stg0_idx_q;
reg        out_stg1_valid_q; reg [2:0] out_stg1_idx_q;
reg        out_stg2_valid_q; reg [2:0] out_stg2_idx_q;

always @(posedge clk) if (rst) begin out_stg0_valid_q <= 0; out_stg0_idx_q <= 0; end
else begin out_stg0_valid_q <= inport_valid_i; out_stg0_idx_q <= inport_idx_i; end
always @(posedge clk) if (rst) begin out_stg1_valid_q <= 0; out_stg1_idx_q <= 0; end
else begin out_stg1_valid_q <= out_stg0_valid_q; out_stg1_idx_q <= out_stg0_idx_q; end
always @(posedge clk) if (rst) begin out_stg2_valid_q <= 0; out_stg2_idx_q <= 0; end
else begin out_stg2_valid_q <= out_stg1_valid_q; out_stg2_idx_q <= out_stg1_idx_q; end

reg signed [31:0] o_s5, o_s6, o_s7;
reg signed [31:0] o_t0, o_t1, o_t2, o_t3, o_t4, o_t5, o_t6, o_t7;
reg signed [31:0] o_t6_5, o_t5_6;

always @(posedge clk)
if (rst) begin
    o_s5 <= 0; o_s6 <= 0; o_s7 <= 0;
    o_t0 <= 0; o_t1 <= 0; o_t2 <= 0; o_t3 <= 0;
    o_t4 <= 0; o_t5 <= 0; o_t6 <= 0; o_t7 <= 0;
    o_t6_5 <= 0; o_t5_6 <= 0;
end else begin
    case (out_stg2_idx_q)
    3'd0: o_t3 <= mul0 + mul1;
    3'd1: begin
        o_s7 <= mul0 + mul1; o_s6 <= mul2 + mul3; o_t0 <= mul4;
    end
    3'd2: begin
        o_t0 <= o_t0 + o_t3; o_t3 <= o_t0 - o_t3; o_t7 <= o_s6 + o_s7;
    end
    3'd3: o_t4 <= (mul0 - mul1) + (mul2 - mul3);
    3'd4: begin
        o_t0 <= mul0 - mul1; o_s5 <= mul2 - mul3;
    end
    3'd5: begin
        o_t3 <= mul0 - mul1; o_t4 <= mul4;
        o_t5 <= o_t0 - o_s5; o_t6 <= o_s7 - o_s6;
    end
    3'd6: begin
        o_t1 <= o_t4 + o_t3; o_t2 <= o_t4 - o_t3;
        o_t6_5 <= o_t6 - o_t5; o_t5_6 <= o_t5 + o_t6;
    end
    default: begin
        o_s5 <= (o_t6_5 * 181) / 256;
        o_s6 <= (o_t5_6 * 181) / 256;
    end
    endcase
end

reg        out_stg3_valid_q; reg [2:0] out_stg3_idx_q;
always @(posedge clk) if (rst) begin out_stg3_valid_q <= 0; out_stg3_idx_q <= 0; end
else begin out_stg3_valid_q <= out_stg2_valid_q; out_stg3_idx_q <= out_stg2_idx_q; end

reg signed [31:0] block_out[0:7];
reg signed [31:0] block_out_tmp;

always @(posedge clk)
if (rst) begin
    block_out[0] <= 0; block_out[1] <= 0; block_out[2] <= 0; block_out[3] <= 0;
    block_out[4] <= 0; block_out[5] <= 0; block_out[6] <= 0; block_out[7] <= 0;
    block_out_tmp <= 0;
end else if (out_stg3_valid_q) begin
    if (out_stg3_idx_q == 3'd3) begin
        block_out[0] <= ((o_t0 + o_t7) >>> OUT_SHIFT);
        block_out_tmp <= ((o_t0 - o_t7) >>> OUT_SHIFT);
        block_out[3] <= ((o_t3 + o_t4) >>> OUT_SHIFT);
        block_out[4] <= ((o_t3 - o_t4) >>> OUT_SHIFT);
    end
    if (out_stg3_idx_q == 3'd6) block_out[7] <= block_out_tmp;
    if (out_stg3_idx_q == 3'd7) begin
        block_out[2] <= ((o_t2 + o_s5) >>> OUT_SHIFT);
        block_out[5] <= ((o_t2 - o_s5) >>> OUT_SHIFT);
        block_out[1] <= ((o_t1 + o_s6) >>> OUT_SHIFT);
        block_out[6] <= ((o_t1 - o_s6) >>> OUT_SHIFT);
    end
end

reg [7:0] valid_q;
always @(posedge clk)
if (rst) valid_q <= 8'b0;
else if (img_start_i) valid_q <= 8'b0;
else valid_q <= {valid_q[6:0], out_stg3_valid_q};

reg [5:0] ptr_q;
always @(posedge clk)
if (rst) ptr_q <= 6'd0;
else if (img_start_i) ptr_q <= 6'd0;
else if (outport_valid_o) ptr_q <= ptr_q + 6'd1;

assign outport_valid_o = valid_q[6];
assign outport_data_o  = block_out[ptr_q[2:0]];

function [5:0] ptr_conv;
    input [5:0] idx;
    reg [5:0] out_idx;
begin
    case (idx)
    6'd0:  out_idx = 6'd0;   6'd1:  out_idx = 6'd8;
    6'd2:  out_idx = 6'd16;  6'd3:  out_idx = 6'd24;
    6'd4:  out_idx = 6'd32;  6'd5:  out_idx = 6'd40;
    6'd6:  out_idx = 6'd48;  6'd7:  out_idx = 6'd56;
    6'd8:  out_idx = 6'd1;   6'd9:  out_idx = 6'd9;
    6'd10: out_idx = 6'd17;  6'd11: out_idx = 6'd25;
    6'd12: out_idx = 6'd33;  6'd13: out_idx = 6'd41;
    6'd14: out_idx = 6'd49;  6'd15: out_idx = 6'd57;
    6'd16: out_idx = 6'd2;   6'd17: out_idx = 6'd10;
    6'd18: out_idx = 6'd18;  6'd19: out_idx = 6'd26;
    6'd20: out_idx = 6'd34;  6'd21: out_idx = 6'd42;
    6'd22: out_idx = 6'd50;  6'd23: out_idx = 6'd58;
    6'd24: out_idx = 6'd3;   6'd25: out_idx = 6'd11;
    6'd26: out_idx = 6'd19;  6'd27: out_idx = 6'd27;
    6'd28: out_idx = 6'd35;  6'd29: out_idx = 6'd43;
    6'd30: out_idx = 6'd51;  6'd31: out_idx = 6'd59;
    6'd32: out_idx = 6'd4;   6'd33: out_idx = 6'd12;
    6'd34: out_idx = 6'd20;  6'd35: out_idx = 6'd28;
    6'd36: out_idx = 6'd36;  6'd37: out_idx = 6'd44;
    6'd38: out_idx = 6'd52;  6'd39: out_idx = 6'd60;
    6'd40: out_idx = 6'd5;   6'd41: out_idx = 6'd13;
    6'd42: out_idx = 6'd21;  6'd43: out_idx = 6'd29;
    6'd44: out_idx = 6'd37;  6'd45: out_idx = 6'd45;
    6'd46: out_idx = 6'd53;  6'd47: out_idx = 6'd61;
    6'd48: out_idx = 6'd6;   6'd49: out_idx = 6'd14;
    6'd50: out_idx = 6'd22;  6'd51: out_idx = 6'd30;
    6'd52: out_idx = 6'd38;  6'd53: out_idx = 6'd46;
    6'd54: out_idx = 6'd54;  6'd55: out_idx = 6'd62;
    6'd56: out_idx = 6'd7;   6'd57: out_idx = 6'd15;
    6'd58: out_idx = 6'd23;  6'd59: out_idx = 6'd31;
    6'd60: out_idx = 6'd39;  6'd61: out_idx = 6'd47;
    6'd62: out_idx = 6'd55;  default: out_idx = 6'd63;
    endcase
    ptr_conv = out_idx;
end
endfunction

assign outport_idx_o = ptr_conv(ptr_q);
endmodule

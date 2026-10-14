//-----------------------------------------------------------------
// jpeg_idct_transpose_ram — 4-bank (8-entry each) to minimize mux depth
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
reg [31:0] ram_b0 [7:0];
reg [31:0] ram_b1 [7:0];
reg [31:0] ram_b2 [7:0];
reg [31:0] ram_b3 [7:0];
reg [31:0] ram_read0_q, ram_read1_q;

wire [1:0] bs0 = addr0_i[4:3];
wire [1:0] bs1 = addr1_i[4:3];
wire [2:0] off0 = addr0_i[2:0];
wire [2:0] off1 = addr1_i[2:0];

always @(posedge clk0_i) begin
    if (wr0_i && bs0 == 2'd0) ram_b0[off0] <= data0_i;
    if (wr0_i && bs0 == 2'd1) ram_b1[off0] <= data0_i;
    if (wr0_i && bs0 == 2'd2) ram_b2[off0] <= data0_i;
    if (wr0_i && bs0 == 2'd3) ram_b3[off0] <= data0_i;
    case (bs0)
    2'd0: ram_read0_q <= ram_b0[off0];
    2'd1: ram_read0_q <= ram_b1[off0];
    2'd2: ram_read0_q <= ram_b2[off0];
    2'd3: ram_read0_q <= ram_b3[off0];
    endcase
end
always @(posedge clk1_i) begin
    if (wr1_i && bs1 == 2'd0) ram_b0[off1] <= data1_i;
    if (wr1_i && bs1 == 2'd1) ram_b1[off1] <= data1_i;
    if (wr1_i && bs1 == 2'd2) ram_b2[off1] <= data1_i;
    if (wr1_i && bs1 == 2'd3) ram_b3[off1] <= data1_i;
    case (bs1)
    2'd0: ram_read1_q <= ram_b0[off1];
    2'd1: ram_read1_q <= ram_b1[off1];
    2'd2: ram_read1_q <= ram_b2[off1];
    2'd3: ram_read1_q <= ram_b3[off1];
    endcase
end
assign data0_o = ram_read0_q;
assign data1_o = ram_read1_q;
endmodule

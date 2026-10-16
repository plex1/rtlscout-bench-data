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
jpeg_idct_x u_idct_x(.clk(clk),.rst(rst),.img_start_i(img_start_i),.img_end_i(img_end_i),
    .inport_valid_i(inport_valid_i),.inport_data0_i(inport_data0_i),.inport_data1_i(inport_data1_i),
    .inport_data2_i(inport_data2_i),.inport_data3_i(inport_data3_i),.inport_idx_i(inport_idx_i),
    .outport_valid_o(idct_x_valid_w),.outport_data_o(idct_x_data_w),.outport_idx_o(idct_x_idx_w));
wire          transpose_valid_w;
wire [ 31:0]  transpose_data0_w;
wire [ 31:0]  transpose_data1_w;
wire [ 31:0]  transpose_data2_w;
wire [ 31:0]  transpose_data3_w;
wire [  2:0]  transpose_idx_w;
jpeg_idct_transpose u_transpose(.clk(clk),.rst(rst),.img_start_i(img_start_i),.img_end_i(img_end_i),
    .inport_valid_i(idct_x_valid_w),.inport_data_i(idct_x_data_w),.inport_idx_i(idct_x_idx_w),
    .inport_accept_o(idct_x_accept_w),
    .outport_valid_o(transpose_valid_w),.outport_data0_o(transpose_data0_w),
    .outport_data1_o(transpose_data1_w),.outport_data2_o(transpose_data2_w),
    .outport_data3_o(transpose_data3_w),.outport_idx_o(transpose_idx_w),.outport_ready_i(1'b1));
jpeg_idct_y u_idct_y(.clk(clk),.rst(rst),.img_start_i(img_start_i),.img_end_i(img_end_i),
    .inport_valid_i(transpose_valid_w),.inport_data0_i(transpose_data0_w),
    .inport_data1_i(transpose_data1_w),.inport_data2_i(transpose_data2_w),
    .inport_data3_i(transpose_data3_w),.inport_idx_i(transpose_idx_w),
    .outport_valid_o(outport_valid_o),.outport_data_o(outport_data_o),
    .outport_idx_o(outport_idx_o));
endmodule
`include "starting_point_body.v"

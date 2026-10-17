module cv32e40p_alu
(
    input logic               clk,
    input logic               rst_n,
    input logic               enable_i,
    input logic [ 6:0] operator_i,
    input logic        [31:0] operand_a_i,
    input logic        [31:0] operand_b_i,
    input logic        [31:0] operand_c_i,
    input logic [1:0] vector_mode_i,
    input logic [4:0] bmask_a_i,
    input logic [4:0] bmask_b_i,
    input logic [1:0] imm_vec_ext_i,
    input logic       is_clpx_i,
    input logic       is_subrot_i,
    input logic [1:0] clpx_shift_i,
    output logic [31:0] result_o,
    output logic        comparison_result_o,
    output logic ready_o,
    input  logic ex_ready_i
);
  assign result_o = operand_a_i + operand_b_i;
  assign comparison_result_o = 1'b0;
  assign ready_o = 1'b1;
endmodule

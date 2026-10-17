`timescale 1ns/1ps
module tb_dbg;
  logic clk=0, rst_n, enable_i, is_clpx_i, is_subrot_i;
  logic [6:0] operator_i;
  logic [31:0] operand_a_i, operand_b_i, operand_c_i;
  logic [1:0] vector_mode_i, imm_vec_ext_i, clpx_shift_i;
  logic [4:0] bmask_a_i, bmask_b_i;
  logic [31:0] result_o, ref_result;
  logic comparison_result_o, ref_cmp, ready_o, ref_ready, ex_ready_i=1;

  always #5 clk = ~clk;

  cv32e40p_alu dut(.clk(clk),.rst_n(rst_n),.enable_i(enable_i),.operator_i(operator_i),
    .operand_a_i(operand_a_i),.operand_b_i(operand_b_i),.operand_c_i(operand_c_i),
    .vector_mode_i(vector_mode_i),.bmask_a_i(bmask_a_i),.bmask_b_i(bmask_b_i),
    .imm_vec_ext_i(imm_vec_ext_i),.is_clpx_i(is_clpx_i),.is_subrot_i(is_subrot_i),
    .clpx_shift_i(clpx_shift_i),.result_o(result_o),.comparison_result_o(comparison_result_o),
    .ready_o(ready_o),.ex_ready_i(ex_ready_i));

  cv32e40p_alu ref(.clk(clk),.rst_n(rst_n),.enable_i(enable_i),.operator_i(operator_i),
    .operand_a_i(operand_a_i),.operand_b_i(operand_b_i),.operand_c_i(operand_c_i),
    .vector_mode_i(vector_mode_i),.bmask_a_i(bmask_a_i),.bmask_b_i(bmask_b_i),
    .imm_vec_ext_i(imm_vec_ext_i),.is_clpx_i(is_clpx_i),.is_subrot_i(is_subrot_i),
    .clpx_shift_i(clpx_shift_i),.result_o(ref_result),.comparison_result_o(ref_cmp),
    .ready_o(ref_ready),.ex_ready_i(ex_ready_i));

  initial begin
    rst_n = 0;
    enable_i = 1;
    repeat (3) @(negedge clk);
    rst_n = 1;

    // test 1: ADD with bmask rounding, vm=0
    operator_i = 7'b0011000; vector_mode_i = 2'b00;
    bmask_a_i = 5'h1e; bmask_b_i = 5'h0f; imm_vec_ext_i = 0;
    is_clpx_i = 0; is_subrot_i = 0; clpx_shift_i = 0;
    operand_a_i = 32'h00008000; operand_b_i = 32'h80000000; operand_c_i = 0;
    #1 $display("T1 dut=%h ref=%h", result_o, ref_result);

    // test 2: ADD clpx shift
    bmask_a_i = 0; bmask_b_i = 0; is_clpx_i = 1; clpx_shift_i = 1;
    operand_a_i = 32'hffffffff; operand_b_i = 32'h00010000;
    #1 $display("T2 dut=%h ref=%h", result_o, ref_result);

    // test 3: plain ADD
    is_clpx_i = 0; clpx_shift_i = 0;
    operand_a_i = 32'h00008000; operand_b_i = 32'h80000000;
    #1 $display("T3 dut=%h ref=%h", result_o, ref_result);

    // test 4: shift op sanity SRL
    operator_i = 7'b0100101; operand_a_i = 32'h0000ffff; operand_b_i = 32'd1;
    #1 $display("T4 dut=%h ref=%h", result_o, ref_result);

    // test 5: SLL
    operator_i = 7'b0100111;
    #1 $display("T5 dut=%h ref=%h", result_o, ref_result);

    // test 6: vm16 SRL
    operator_i = 7'b0100101; vector_mode_i = 2'b10;
    operand_a_i = 32'h0000ffff; operand_b_i = 32'h00010001;
    #1 $display("T6 dut=%h ref=%h", result_o, ref_result);

    // test 7: vm8 SRA
    operator_i = 7'b0100100; vector_mode_i = 2'b11;
    operand_a_i = 32'hff00ff00; operand_b_i = 32'h01010101;
    #1 $display("T7 dut=%h ref=%h", result_o, ref_result);

    // test 8: ROR
    operator_i = 7'b0100110; vector_mode_i = 2'b00;
    operand_a_i = 32'hdeadbeef; operand_b_i = 32'd4;
    #1 $display("T8 dut=%h ref=%h", result_o, ref_result);

    // test 9: ADDR rounding
    operator_i = 7'b0011100; vector_mode_i = 2'b00;
    bmask_a_i = 5'h03; bmask_b_i = 5'h04;
    operand_a_i = 32'h12345678; operand_b_i = 32'h9abcdef0;
    #1 $display("T9 dut=%h ref=%h", result_o, ref_result);

    // test 10: SUB vm16
    operator_i = 7'b0011001; vector_mode_i = 2'b10;
    operand_a_i = 32'h00010002; operand_b_i = 32'h00030004;
    #1 $display("T10 dut=%h ref=%h", result_o, ref_result);

    $finish;
  end
endmodule

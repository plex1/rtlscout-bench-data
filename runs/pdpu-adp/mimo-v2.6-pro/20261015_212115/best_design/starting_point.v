// PDPU dot-16, 5-stage pipelined -- readable SystemVerilog baseline.
//
// Source: github.com/qleenju/PDPU (ISCAS 2023, Apache-2.0), pdpu_top_pipelined,
// N default 4 -> 16, clock port clk -> clk to match the testbench.
//
// This is the ORIGINAL SystemVerilog, not sv2v output, so it keeps the upstream
// comments, the real generate-block names and `always_ff`. Only the changes
// yosys's front end actually requires were made:
//   * function bodies use `<name> = expr;` instead of `return expr;`
//   * `automatic` dropped from function-local variable declarations
//   * `unsigned'(x)` / `signed'(x)` -> `$unsigned(x)` / `$signed(x)`
//   * `'{a,b,c}` assignment patterns -> `{a,b,c}` (identical bit layout)
//   * simulation-only `// pragma translate_off` regions removed
//   * packed 2-D arrays lowered to flat vectors -- element e of `x` is
//     `x[e*W +: W]`, bit b of element e is `x[e*W + b]`. This is the ONE
//     construct yosys cannot handle (it cannot slice the outer dimension), and
//     the sole reason the sv2v detour existed.
//
// Module hierarchy, arithmetic and pipeline cut points are upstream verbatim.

// ---- FFNR: non-resettable flip-flop (from common_cells registers.svh) ----
`define FFNR(__q, __d, __clk)        \
  always_ff @(posedge (__clk)) begin \
    __q <= (__d);                    \
  end


package pdpu_pkg;
    // --------------
    // General Helper Function
    // ---------------
    function automatic int maximum(int a, int b);
        maximum = (a > b) ? a : b;
    endfunction

    function automatic integer clog2(input integer n);
        begin
            n = n-1;
            for(clog2=0; n>0; clog2=clog2+1)
                n = n>>1;
        end
    endfunction

endpackage

// Copyright 2016 ETH Zurich and University of Bologna.
// Copyright and related rights are licensed under the Solderpad Hardware
// License, Version 0.51 (the "License"); you may not use this file except in
// compliance with the License.  You may obtain a copy of the License at
// http://solderpad.org/licenses/SHL-0.51. Unless required by applicable law
// or agreed to in writing, software, hardware and materials distributed under
// this License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

/// cf_math_pkg: Constant Function Implementations of Mathematical Functions for HDL Elaboration
///
/// This package contains a collection of mathematical functions that are commonly used when defining
/// the value of constants in HDL code.  These functions are implemented as Verilog constants
/// functions.  Introduced in Verilog 2001 (IEEE Std 1364-2001), a constant function (§ 10.3.5) is a
/// function whose value can be evaluated at compile time or during elaboration.  A constant function
/// must be called with arguments that are constants.
package cf_math_pkg;

    /// Ceiled Division of Two Natural Numbers
    ///
    /// Returns the quotient of two natural numbers, rounded towards plus infinity.
    function automatic integer ceil_div (input longint dividend, input longint divisor);
        longint remainder;


        remainder = dividend;
        for (ceil_div = 0; remainder > 0; ceil_div++) begin
            remainder = remainder - divisor;
        end
    endfunction

    function automatic integer clog(input integer n);
        begin
            n = n-1;
            for(clog=0; n>0; clog=clog+1)
                n = n>>1;
        end
    endfunction
    /// Index width required to be able to represent up to `num_idx` indices as a binary
    /// encoded signal.
    /// Ensures that the minimum width if an index signal is `1`, regardless of parametrization.
    ///
    /// Sample usage in type definition:
    /// As parameter:
    ///   `parameter type idx_t = logic[cf_math_pkg::idx_width(NumIdx)-1:0]`
    /// As typedef:
    ///   `typedef logic [cf_math_pkg::idx_width(NumIdx)-1:0] idx_t`
    function automatic integer unsigned idx_width (input integer unsigned num_idx);
        idx_width = (num_idx > 32'd1) ? $unsigned(clog(num_idx)) : 32'd1;
    endfunction

endpackage

// full adder, the basic module that constitutes the 3:2 compressor
module fulladder(
    input logic x,y,z,
    output logic sum,carry
);
    assign sum = (x ^ y ) ^ z;
    assign carry = ((x ^ y) & z) | (x & y);
endmodule

// 5:3 counter, the basic module that constitutes the 4:2 compressor
module counter_5to3(
    input logic x1,x2,x3,x4,cin,
    output logic sum,carry,cout
);
    assign sum = x1 ^ x2 ^ x3 ^ x4 ^ cin;
    assign cout = (x1 ^ x2) & x3 | ~(x1 ^ x2) & x1;
    assign carry = (x1 ^ x2 ^ x3 ^ x4) & cin | ~(x1 ^ x2 ^ x3 ^ x4) & x4;
endmodule

// Compare to obtain the larger of two signed numbers
module comparator #(
    parameter int unsigned WIDTH = 8
)(
    input logic signed [WIDTH:0] operand_a,
    input logic signed [WIDTH:0] operand_b,
    output logic signed [WIDTH:0] result_o
);
    // unsigned max of two values (operands given in offset binary)
    logic [WIDTH:0] ua, ub;
    assign ua = operand_a;
    assign ub = operand_b;

    always_comb begin
        result_o = (ua>ub) ? operand_a : operand_b;
    end
endmodule

/*
The radix-4 booth encoding is as follows:
    code | value | operation
    000  |  +0   |     0
    001  |  +1   |    1*A
    010  |  +1   |    1*A
    011  |  +2   |    2*A
    100  |  -2   |   -2*A
    101  |  -1   |   -1*A
    110  |  -1   |   -1*A
    111  |   0   |     0
*/
module booth_encoder(
    input [2:0] code,
    output neg,
    output zero,
    output one,
    output two
);
    // Only when code = 100/101/110, neg = 1
    assign neg = code[2]&(~code[0]) | code[2]&(~code[1]);
    assign zero = ~(|code) | &code;
    assign two = ~code[2]&code[1]&code[0] | code[2]&(~code[1])&(~code[0]);
    assign one = code[1]&(~code[0]) | (~code[1])&code[0];
endmodule

// Generate partial product according to radix-4 booth encoding result
module gen_product #(
    parameter int unsigned WIDTH = 16
)(
    input logic [WIDTH-1:0] multiplicand,
    input logic [2:0] code,
    output logic [WIDTH:0] partial_prod,
    output sign
);
    logic neg,zero,one,two;
    booth_encoder u_booth_encoder(
        .code(code),
        .neg(neg),
        .zero(zero),
        .one(one),
        .two(two)
    );

    logic [WIDTH:0] temp_prod;
    always_comb begin
        if(one) begin
            temp_prod = multiplicand;
        end
        else if(two) begin
            temp_prod = multiplicand << 1;
        end
        else begin
            temp_prod = '0;
        end
    end
    assign partial_prod = neg ? (~temp_prod) : temp_prod;
    assign sign = neg;
endmodule

// Generate partial products based on modified sign extension in booth multipliers
module gen_prods #(
    parameter int unsigned WIDTH_A = 16,
    parameter int unsigned WIDTH_B = 16,
    // do not change
    parameter int unsigned COUNT = (WIDTH_B+2)/2,       // count of partial products generated in unsigned multiplication
    parameter int unsigned WIDTH_O = WIDTH_A + WIDTH_B  // bit-width of outputs
)(
    input logic [WIDTH_A-1:0] operand_a,
    input logic [WIDTH_B-1:0] operand_b,
    output [COUNT-1:0][WIDTH_O-1:0] partial_prods
);
    logic [WIDTH_B+2:0] multiplier;
    logic [COUNT-1:0][2:0] codes;
    logic [COUNT-1:0][WIDTH_A:0] temp_prods;
    logic [COUNT-1:0] signs;

    assign multiplier = {2'b00,operand_b,1'b0};
    assign codes[0] = multiplier[2:0];

    // Generate partial products
    gen_product #(
        .WIDTH(WIDTH_A)
    ) ua_gen_product(
        .multiplicand(operand_a),
        .code(codes[0]),
        .partial_prod(temp_prods[0]),
        .sign(signs[0])
    );
    assign partial_prods[0] = {~signs[0], signs[0], signs[0], temp_prods[0]};

    generate
        genvar i;
        for(i=1;i<COUNT-1;i++) begin
            assign codes[i] = multiplier[2*i+2:2*i];
            gen_product #(
                .WIDTH(WIDTH_A)
            ) ub_gen_product(
                .multiplicand(operand_a),
                .code(codes[i]),
                .partial_prod(temp_prods[i]),
                .sign(signs[i])
            );
            // modified sign extension in booth multipliers
            assign partial_prods[i] = {1'b1,~signs[i],temp_prods[i],1'b0,signs[i-1]} << (2*i-2);
        end
    endgenerate

    assign codes[COUNT-1] = multiplier[2*COUNT:2*COUNT-2];
    gen_product #(
        .WIDTH(WIDTH_A)
    ) uc_gen_product(
        .multiplicand(operand_a),
        .code(codes[COUNT-1]),
        .partial_prod(temp_prods[COUNT-1]),
        .sign(signs[COUNT-1])
    );
    assign partial_prods[COUNT-1] = {temp_prods[COUNT-1],1'b0,signs[COUNT-2]} << (2*COUNT-4);

endmodule

// 3:2 compressor, compressing 3 inputs into sum and carry
// `operands_i` is a flat vector: element k occupies bits [k*WIDTH_I +: WIDTH_I].
// (Flat rather than a 2-D packed array so yosys's Verilog front end can read it.)
module compressor_3to2 #(
    parameter int unsigned WIDTH_I = 8,
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(3)
)(
    input logic [3*WIDTH_I-1:0] operands_i,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    logic [WIDTH_I-1:0] sum;
    logic [WIDTH_I-1:0] carry;
    generate
        genvar i;
        for(i=0;i<WIDTH_I;i++) begin
            fulladder u_fulladder(
                .x(operands_i[0*WIDTH_I + i]),
                .y(operands_i[1*WIDTH_I + i]),
                .z(operands_i[2*WIDTH_I + i]),
                .sum(sum[i]),
                .carry(carry[i])
            );
        end
    endgenerate

    assign sum_o = sum;
    assign carry_o = {1'b0,carry,1'b0};
endmodule

// 4:2 compressor, compressing 4 inputs into sum and carry
// `operands_i` is a flat vector: element k occupies bits [k*WIDTH_I +: WIDTH_I].
// (Flat rather than a 2-D packed array so yosys's Verilog front end can read it.)
module compressor_4to2 #(
    parameter int unsigned WIDTH_I = 8,                             // bit-width of inputs
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(4)   // bit-width of outputs
)(
    input logic [4*WIDTH_I-1:0] operands_i,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    logic [WIDTH_I-1:0] sum;
    logic [WIDTH_I:0] cin;
    logic [WIDTH_I-1:0] cout;
    logic [WIDTH_I-1:0] carry;
    
    assign cin[0] = 1'b0;

    // Cascaded 5:3 counters according to input bit-width
    generate
        genvar i;
        for(i=0;i<WIDTH_I;i++) begin
            counter_5to3 u_counter_5to3(
                .x1(operands_i[0*WIDTH_I + i]),
                .x2(operands_i[1*WIDTH_I + i]),
                .x3(operands_i[2*WIDTH_I + i]),
                .x4(operands_i[3*WIDTH_I + i]),
                .cin(cin[i]),
                .sum(sum[i]),
                .carry(carry[i]),
                .cout(cout[i])
            );
            assign cin[i+1] = cout[i];
        end
    endgenerate

    logic [1:0] carry_temp;
    
    assign sum_o = sum;
    assign carry_temp = carry[WIDTH_I-1]+cin[WIDTH_I];
    assign carry_o = {carry_temp,carry[WIDTH_I-2:0],1'b0};
endmodule

// Recursive implementation of carry-save-adder (CSA) tree, to compress inputs with any size into sum and carry
// `operands_i` is a flat vector: element k occupies bits [k*WIDTH_I +: WIDTH_I].
// (Flat rather than a 2-D packed array so yosys's Verilog front end can read it.)
module csa_tree #(
    parameter int unsigned N = 8,           // size of inputs
    parameter int unsigned WIDTH_I = 8,     // bit-width of inputs
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(N)   // bit-width of outputs
)(
    input logic [N*WIDTH_I-1:0] operands_i,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    localparam int unsigned N_A = N/2;
    localparam int unsigned N_B = N - N_A;

    generate
        if (N==1) begin
            assign sum_o = operands_i[0 +: WIDTH_I];
            assign carry_o = '0;
        end
        else if(N==2) begin
            assign sum_o = operands_i[0 +: WIDTH_I];
            assign carry_o = operands_i[WIDTH_I +: WIDTH_I];
        end
        else if(N==3) begin
            compressor_3to2 #(
                .WIDTH_I(WIDTH_I),
                .WIDTH_O(WIDTH_O)
            ) u_compressor_3to2(
                .operands_i(operands_i),
                .sum_o(sum_o),
                .carry_o(carry_o)
            );
        end
        else if(N==4) begin
            compressor_4to2 #(
                .WIDTH_I(WIDTH_I),
                .WIDTH_O(WIDTH_O)
            ) u_compressor_4to2(
                .operands_i(operands_i),
                .sum_o(sum_o),
                .carry_o(carry_o)
            );
        end
        else begin
            logic [N_A*WIDTH_I-1:0] operands_i_A;
            logic [N_B*WIDTH_I-1:0] operands_i_B;
            logic [WIDTH_O-1:0] sum_o_A;
            logic [WIDTH_O-1:0] sum_o_B;
            logic [WIDTH_O-1:0] carry_o_A;
            logic [WIDTH_O-1:0] carry_o_B;

            // Divide the inputs into two chunks
            assign operands_i_A = operands_i[0           +: N_A*WIDTH_I];
            assign operands_i_B = operands_i[N_A*WIDTH_I +: N_B*WIDTH_I];

            csa_tree #(
                .N(N_A),
                .WIDTH_I(WIDTH_I),
                .WIDTH_O(WIDTH_O)
            ) ua_csa_tree(
                .operands_i(operands_i_A),
                .sum_o(sum_o_A),
                .carry_o(carry_o_A)
            );

            csa_tree #(
                .N(N_B),
                .WIDTH_I(WIDTH_I),
                .WIDTH_O(WIDTH_O)
            ) ub_csa_tree(
                .operands_i(operands_i_B),
                .sum_o(sum_o_B),
                .carry_o(carry_o_B)
            );

            logic [4*WIDTH_O-1:0] operands_i_C;
            assign operands_i_C = {sum_o_A, carry_o_A, sum_o_B, carry_o_B};
            
            compressor_4to2 #(
                .WIDTH_I(WIDTH_O),
                .WIDTH_O(WIDTH_O)
            ) uc_compressor_4to2(
                .operands_i(operands_i_C),
                .sum_o(sum_o),
                .carry_o(carry_o)
            );
        end
    endgenerate
endmodule

// Recursive implementation of comparator tree, to obtain the maximum value among inputs with any size
//
// `operands_i` is a flat vector of N elements, each (WIDTH+1) bits wide, element i
// occupying bits [i*(WIDTH+1) +: WIDTH+1]. (Written flat rather than as a 2-D packed
// array so the yosys Verilog front end can read this file directly.)
module comp_tree #(
    parameter N = 4,        // size of inputs
    parameter WIDTH = 8     // bit-width of inputs
)(
    input logic signed [N*(WIDTH+1)-1:0] operands_i,
    output logic signed [WIDTH:0] result_o
);
    localparam int unsigned W   = WIDTH + 1;   // bits per element
    localparam int unsigned N_A = N/2;
    localparam int unsigned N_B = N - N_A;

    generate
        if (N==1) begin : g_leaf
            assign result_o = operands_i[0 +: W];
        end

        else if (N==2) begin : g_pair
            comparator #(
                .WIDTH(WIDTH)
            ) u_comparator(
                .operand_a(operands_i[0 +: W]),
                .operand_b(operands_i[W +: W]),
                .result_o(result_o)
            );
        end

        else begin : g_split
            logic signed [N_A*W-1:0] operands_i_A;
            logic signed [N_B*W-1:0] operands_i_B;
            logic signed [WIDTH:0]   result_o_A;
            logic signed [WIDTH:0]   result_o_B;

            // Divide the inputs into two chunks
            assign operands_i_A = operands_i[0       +: N_A*W];
            assign operands_i_B = operands_i[N_A*W   +: N_B*W];

            // Module recursion for operands_i_A
            comp_tree #(
                .N(N_A),
                .WIDTH(WIDTH)
            ) ua_comp_tree(
                .operands_i(operands_i_A),
                .result_o(result_o_A)
            );

            // Module recursion for operands_i_B
            comp_tree #(
                .N(N_B),
                .WIDTH(WIDTH)
            ) ub_comp_tree(
                .operands_i(operands_i_B),
                .result_o(result_o_B)
            );

            // Compare result_o_A and result_o_B
            comparator #(
                .WIDTH(WIDTH)
            ) uc_comparator(
                .operand_a(result_o_A),
                .operand_b(result_o_B),
                .result_o(result_o)
            );
        end
    endgenerate
endmodule

// Modified radix-4 booth wallace multiplier
module radix4_booth_multiplier #(
    parameter int unsigned WIDTH_A = 16,                // bit-width of operand A
    parameter int unsigned WIDTH_B = 16,                // bit-width of operand B
    parameter int unsigned WIDTH_O = WIDTH_A + WIDTH_B  // bit-width of the multiplication result
)(
    input logic [WIDTH_A-1:0] operand_a,
    input logic [WIDTH_B-1:0] operand_b,
    // the partial products are compressed into sum_o and carry_o
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    // count of partial products generated in unsigned multiplication
    localparam int unsigned COUNT = (WIDTH_B+2)/2;
    
    // ---------------
    // Generate partial products
    // ---------------
    logic [COUNT-1:0][WIDTH_O-1:0] partial_prods;
    gen_prods #(
        .WIDTH_A(WIDTH_A),
        .WIDTH_B(WIDTH_B)
    ) u_gen_prods(
        .operand_a(operand_a),
        .operand_b(operand_b),
        .partial_prods(partial_prods)
    );

    // ---------------
    // Compress the partial products into sum and carry
    // ---------------
    csa_tree #(
        .N(COUNT),
        .WIDTH_I(WIDTH_O),
        .WIDTH_O(WIDTH_O)
    ) u_csa_tree(
        .operands_i(partial_prods),
        .sum_o(sum_o),
        .carry_o(carry_o)
    );
endmodule

// barrel shifter
module barrel_shifter #(
    parameter int unsigned WIDTH = 8,               // bit-width of input data
    parameter int unsigned SHIFT_WIDTH = 3,         // bit-width of shift amount
    // MODE=0-->shift left; MODE=1-->shift right
    parameter bit MODE = 1'b0
)(
    input logic [WIDTH-1:0] operand_i,
    input logic [SHIFT_WIDTH-1:0] shift_amount,
    output logic [WIDTH-1:0] result_o
);
    logic [SHIFT_WIDTH-1:0][WIDTH-1:0] temp_results;

    assign temp_results[SHIFT_WIDTH-1] = operand_i;

    generate
        genvar i;
        // Left Shift
        if(MODE==1'b0) begin
            for(i=SHIFT_WIDTH-1;i>0;i--) begin
                assign temp_results[i-1] = shift_amount[i] ? temp_results[i]<<(2**i) : temp_results[i];
            end
            assign result_o = shift_amount[0] ? temp_results[0]<<1 : temp_results;
        end
        // Left Right
        else begin
            for(i=SHIFT_WIDTH-1;i>0;i--) begin
                assign temp_results[i-1] = shift_amount[i] ? temp_results[i]>>(2**i) : temp_results[i];
            end
            assign result_o = shift_amount[0] ? temp_results[0]>>1 : temp_results;
        end
    endgenerate
endmodule

// Copyright (c) 2018 - 2019 ETH Zurich, University of Bologna
// All rights reserved.
//
// This code is under development and not yet released to the public.
// Until it is released, the code is under the copyright of ETH Zurich and
// the University of Bologna, and may contain confidential and/or unpublished
// work. Any reuse/redistribution is strictly forbidden without written
// permission from ETH Zurich.
//
// Bug fixes and contributions will eventually be released under the
// SolderPad open hardware license in the context of the PULP platform
// (http://www.pulp-platform.org), under the copyright of ETH Zurich and the
// University of Bologna.

/// A trailing zero counter / leading zero counter.
/// Set MODE to 0 for trailing zero counter => cnt_o is the number of trailing zeros (from the LSB)
/// Set MODE to 1 for leading zero counter  => cnt_o is the number of leading zeros  (from the MSB)
/// If the input does not contain a zero, `empty_o` is asserted. Additionally `cnt_o` contains
/// the maximum number of zeros - 1. For example:
///   in_i = 000_0000, empty_o = 1, cnt_o = 6 (mode = 0)
///   in_i = 000_0001, empty_o = 0, cnt_o = 0 (mode = 0)
///   in_i = 000_1000, empty_o = 0, cnt_o = 3 (mode = 0)
/// Furthermore, this unit contains a more efficient implementation for Verilator (simulation only).
/// This speeds up simulation significantly.
module lzc #(
  /// The width of the input vector.
  parameter int unsigned WIDTH = 2,
  /// Mode selection: 0 -> trailing zero, 1 -> leading zero
  parameter bit          MODE  = 1'b0,
  /// Dependent parameter. Do **not** change!
  ///
  /// Width of the output signal with the zero count.
  parameter int unsigned CNT_WIDTH = cf_math_pkg::idx_width(WIDTH)
) (
  /// Input vector to be counted.
  input  logic [WIDTH-1:0]     in_i,
  /// Count of the leading / trailing zeros.
  output logic [CNT_WIDTH-1:0] cnt_o,
  /// Counter is empty: Asserted if all bits in in_i are zero.
  output logic                 empty_o
);

  if (WIDTH == 1) begin : gen_degenerate_lzc

    assign cnt_o[0] = !in_i[0];
    assign empty_o = !in_i[0];

  end else begin : gen_lzc

    localparam int unsigned NumLevels = cf_math_pkg::clog(WIDTH);


    // index_lut / index_nodes are flat: entry e occupies [e*NumLevels +: NumLevels].
    logic [WIDTH*NumLevels-1:0] index_lut;
    logic [2**NumLevels-1:0] sel_nodes;
    logic [(2**NumLevels)*NumLevels-1:0] index_nodes;

    logic [WIDTH-1:0] in_tmp;

    // reverse vector if required
    always_comb begin : flip_vector
      for (int unsigned i = 0; i < WIDTH; i++) begin
        in_tmp[i] = (MODE) ? in_i[WIDTH-1-i] : in_i[i];
      end
    end

    for (genvar j = 0; $unsigned(j) < WIDTH; j++) begin : g_index_lut
      assign index_lut[(j)*NumLevels +: NumLevels] = ($unsigned(j));
    end

    for (genvar level = 0; $unsigned(level) < NumLevels; level++) begin : g_levels
      if ($unsigned(level) == NumLevels - 1) begin : g_last_level
        for (genvar k = 0; k < 2 ** level; k++) begin : g_level
          // if two successive indices are still in the vector...
          if ($unsigned(k) * 2 < WIDTH - 1) begin : g_reduce
            assign sel_nodes[2 ** level - 1 + k] = in_tmp[k * 2] | in_tmp[k * 2 + 1];
            assign index_nodes[(2 ** level - 1 + k)*NumLevels +: NumLevels] = (in_tmp[k * 2] == 1'b1)
              ? index_lut[(k * 2)*NumLevels +: NumLevels] :
                index_lut[(k * 2 + 1)*NumLevels +: NumLevels];
          end
          // if only the first index is still in the vector...
          if ($unsigned(k) * 2 == WIDTH - 1) begin : g_base
            assign sel_nodes[2 ** level - 1 + k] = in_tmp[k * 2];
            assign index_nodes[(2 ** level - 1 + k)*NumLevels +: NumLevels] = index_lut[(k * 2)*NumLevels +: NumLevels];
          end
          // if index is out of range
          if ($unsigned(k) * 2 > WIDTH - 1) begin : g_out_of_range
            assign sel_nodes[2 ** level - 1 + k] = 1'b0;
            assign index_nodes[(2 ** level - 1 + k)*NumLevels +: NumLevels] = '0;
          end
        end
      end else begin : g_not_last_level
        for (genvar l = 0; l < 2 ** level; l++) begin : g_level
          assign sel_nodes[2 ** level - 1 + l] =
              sel_nodes[2 ** (level + 1) - 1 + l * 2] | sel_nodes[2 ** (level + 1) - 1 + l * 2 + 1];
          assign index_nodes[(2 ** level - 1 + l)*NumLevels +: NumLevels] = (sel_nodes[2 ** (level + 1) - 1 + l * 2] == 1'b1)
            ? index_nodes[(2 ** (level + 1) - 1 + l * 2)*NumLevels +: NumLevels] :
              index_nodes[(2 ** (level + 1) - 1 + l * 2 + 1)*NumLevels +: NumLevels];
        end
      end
    end

    assign cnt_o = NumLevels > $unsigned(0) ? index_nodes[(0)*NumLevels +: NumLevels] : {(cf_math_pkg::clog(WIDTH)) {1'b0}};
    assign empty_o = NumLevels > $unsigned(0) ? ~sel_nodes[0] : ~(|in_i);

  end : gen_lzc

endmodule : lzc

// Mantissa normalization
module mantissa_norm #(
    parameter int unsigned WIDTH = 8,
    parameter int unsigned EXP_WIDTH = 3,
    parameter int unsigned DECIMAL_POINT = 3     // digits before decimal point
)(
    input logic [WIDTH-1:0] operand_i,
    output logic signed [EXP_WIDTH:0] exp_adjust,
    output logic [WIDTH-1:0] result_o
);
    localparam int unsigned LZC_WIDTH = pdpu_pkg::clog2(WIDTH);
    logic [LZC_WIDTH-1:0] leading_zero_count;
    logic lzc_zeroes;
    
    // leading zero count
    lzc #(
        .WIDTH(WIDTH),
        .MODE(1'b1)
    ) u_lzc(
        .in_i(operand_i),
        .cnt_o(leading_zero_count),
        .empty_o(lzc_zeroes)
    );

    // Exponent adjustment
    always_comb begin
        if(lzc_zeroes) begin
            exp_adjust = '0;
        end
        else if(leading_zero_count <= DECIMAL_POINT - 1) begin
            exp_adjust = DECIMAL_POINT - leading_zero_count - 1;
        end
        else begin
            exp_adjust = -$signed(leading_zero_count - DECIMAL_POINT + 1);
        end
    end

    // Mantissa normalization
    barrel_shifter #(
        .WIDTH(WIDTH),
        .SHIFT_WIDTH(LZC_WIDTH),
        .MODE(1'b0)
    ) u_barrel_shifter(
        .operand_i(operand_i),
        .shift_amount(leading_zero_count),
        .result_o(result_o)
    );
endmodule

// Posit decoder
module posit_decoder_generic #(
    parameter int unsigned n = 16,  // word size
    parameter int unsigned es = 1,  // exponent size
    //do not change
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd+es,       // not include sign bit
    parameter int unsigned MANT_WIDTH = n-es-3      // not include implicit bit
)(
    input logic [n-1:0] operand_i,
    output logic sign_o,
    output logic signed [EXP_WIDTH:0] rg_exp_o,  
    output logic [MANT_WIDTH:0] mant_norm_o
);

    // ---------------
    // Perform two's complement if negative
    // ---------------
    logic sign;
    logic [n-2:0] operand_value;

    assign sign = operand_i[n-1];
    logic [n-2:0] operand_neg;
    prefix_add #(
        .WIDTH(n-1)
    ) u_dec_neg(
        .a(~operand_i[n-2:0]),
        .b({(n-1){1'b0}}),
        .cin(1'b1),
        .sum(operand_neg),
        .cout()
    );
    assign operand_value = sign ? operand_neg : operand_i[n-2:0];
    
    // ---------------
    // Leading Zero Count
    // ---------------
    logic regS;
    logic [n-2:0] lzc_operand;
    logic [nd-1:0] leading_zero_count;
    logic lzc_zeroes;       // lzc_zeroes = 1 if input is all 0s

    assign regS = operand_value[n-2];
    assign lzc_operand = regS ? (~operand_value) : operand_value;

    lzc #(
        .WIDTH(n-1),
        .MODE(1'b1)    // mode=1 means "counting leading zeroes"
    ) u_lzc(
        .in_i(lzc_operand),
        .cnt_o(leading_zero_count),
        .empty_o(lzc_zeroes)
    );

    logic [nd-1:0] runlength;       // number of consecutive identical bits in regime field
    logic [nd-1:0] regime_bits;     // bitwidth of regime field
    logic signed [nd:0] regime_k;   // value represented by regime field

    assign runlength = lzc_zeroes ? (n-1) : leading_zero_count;
    assign regime_bits = lzc_zeroes ? (n-1) : (leading_zero_count+1);
    assign regime_k = regS ? ({1'b0,runlength-1}) : ({1'b1,~runlength+1});

    // ---------------
    // Left shift regime field
    // ---------------
    logic [n-2:0] op_no_rg;

    barrel_shifter #(
        .WIDTH(n-1),
        .SHIFT_WIDTH(nd),
        .MODE(1'b0)
    ) u_barrel_shifter(
        .operand_i(operand_value),
        .shift_amount(regime_bits),
        .result_o(op_no_rg)
    );

    // ---------------
    // Extract sign
    // ---------------
    assign sign_o = sign;

    // ---------------
    // Extract valid exponent
    // ---------------
    logic [es:0] exp;
    if(es==0) begin
        assign exp = 0;
    end
    else begin
        assign exp = op_no_rg[n-2:n-2-es+1];
    end
    
    assign rg_exp_o = regime_k << es | exp;

    // ---------------
    // Extract mantissa
    // ---------------
    logic implicit_bit;
    
    assign implicit_bit = |operand_i[n-2:0];
    assign mant_norm_o = {implicit_bit, op_no_rg[n-2-es:2]};
endmodule

// Posit encoder
module posit_encoder #(
    parameter int unsigned n = 16,
    parameter int unsigned es = 1,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd+es,       // not include sign bit
    parameter int unsigned MANT_WIDTH = n-es-3      // not include implicit bit
) (
    input logic sign_i,
    input logic signed [EXP_WIDTH:0] rg_exp_i,
    input logic [MANT_WIDTH:0] mant_norm_i,

    output logic [n-1:0] result_o
);
    // ---------------
    // Determine whether input is effective 0
    // ---------------
    logic input_not_zero;
    assign input_not_zero = mant_norm_i[MANT_WIDTH];

    // ---------------
    // Compute regime_k and exp
    // ---------------
    logic signed [EXP_WIDTH-es:0] regime_k;
    logic signed [es:0] exp;

    assign regime_k = rg_exp_i[EXP_WIDTH:es];
    if(es==0) begin
        assign exp = 0;
    end
    else begin
        assign exp = rg_exp_i[es-1:0];
    end

    // ---------------
    // initial regime field
    // ---------------
    logic sign_k;
    logic [n-2:0] rg_const;
    logic [n-2:0] regime;
    
    assign sign_k = rg_exp_i[EXP_WIDTH];
    assign rg_const = 1;
    assign regime = sign_k ? rg_const : ~rg_const;

    // ---------------
    // Compute regime bits
    // ---------------
    logic [EXP_WIDTH-es:0] regime_bits;
    
    assign regime_bits = (sign_k ? ~regime_k : regime_k) + 2;
    
    // ---------------
    // Combine {regime, Exp, Mantissa}
    // ---------------
    localparam int unsigned TEMP_WIDTH = (n-1)+es+MANT_WIDTH;   // regime + exponent + mantissa
    logic [TEMP_WIDTH-1:0] rg_exp_mant;

    if(es==0) begin
        assign rg_exp_mant = {regime, mant_norm_i[MANT_WIDTH-1:0]};
    end
    else begin
        assign rg_exp_mant = {regime, exp[es-1:0], mant_norm_i[MANT_WIDTH-1:0]}; 
    end

    // ---------------
    // Calculate amount of right shift
    // ---------------
    localparam int unsigned MAX_SHIFT_AMOUNT = MANT_WIDTH + es + 1;
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(MAX_SHIFT_AMOUNT+1);
    logic [SHIFT_WIDTH-1:0] shift_amount;
    
    assign shift_amount = (regime_bits>=n) ? MAX_SHIFT_AMOUNT : (regime_bits+(MANT_WIDTH+es-n+1));

    // ---------------
    // Right Shift
    // ---------------
    logic [TEMP_WIDTH+MAX_SHIFT_AMOUNT-1:0] value_before_shift, value_after_shift;
    
    assign value_before_shift = rg_exp_mant << MAX_SHIFT_AMOUNT;
    barrel_shifter #(
        .WIDTH(TEMP_WIDTH+MAX_SHIFT_AMOUNT),
        .SHIFT_WIDTH(SHIFT_WIDTH),
        .MODE(1'b1)
    ) u_barrel_shifter(
        .operand_i(value_before_shift),
        .shift_amount(shift_amount),
        .result_o(value_after_shift)
    );

    // ---------------
    // Compute abs_value and rounding_bits
    // ---------------
    logic [n-2:0] value_before_round;
    logic [MAX_SHIFT_AMOUNT-1:0] rounding_bits;

    assign {value_before_round,rounding_bits} = value_after_shift[MAX_SHIFT_AMOUNT+n-2:0];

    // ---------------
    // Perform rounding (RNE mode is applied by default)
    // ---------------
    logic round_bit;
    logic sticky_bit;
    logic round_value;
    logic [n-2:0] value_after_round;

    assign round_bit = rounding_bits[MAX_SHIFT_AMOUNT-1];
    assign sticky_bit = |rounding_bits[MAX_SHIFT_AMOUNT-2:0];
    assign round_value = round_bit & (sticky_bit | value_before_round[0]);
    // Fused rounding + final negation:
    //   out = sign ? -(v+r) : (v+r)  ==  (sign ? ~v : v) + (r ^ sign)
    logic [n-2:0] add_op;
    assign add_op = sign_i ? ~value_before_round : value_before_round;
    prefix_add #(
        .WIDTH(n-1)
    ) u_round_add(
        .a(add_op),
        .b({(n-1){1'b0}}),
        .cin(round_value ^ sign_i),
        .sum(value_after_round),
        .cout()
    );

    // ---------------
    // Output result
    // ---------------
    logic [n-1:0] normal_result;
    assign normal_result = {sign_i, value_after_round};
    assign result_o = input_not_zero ? normal_result : '0;

endmodule

/*
PDPU equipped with a fine-grained 6-stage
*/

module pdpu_top_pipelined_old #(
    parameter int unsigned N = 16,                  // dot-product size
    parameter int unsigned n_i = 8,                 // word size
    parameter int unsigned es_i = 2,                // exponent size
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14         // alignment width
)(
    input logic clk,
    input logic [N*n_i-1:0] operands_a,        // input vector (Va)
    input logic [N*n_i-1:0] operands_b,        // input vector (Vb)
    input logic [n_o-1:0] acc,                      // previous output (acc)
    output logic [n_o-1:0] result_o                 // output result (out)
);
    // ********************
    // Pipeline 0: Input Reg
    // ********************
    // group_path -name pipe0to1 -from {operands_a operands_b operands_c} -to {pipe1*}


    // ---------------
    // Posit Decoding
    // ---------------
    // the maximum exponent size (not include sign bit)
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;
    // the maximum mantissa size (not include implicit bit)
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;

    // Extract valid sign, exponent, and mantissa from input vectors
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ua_posit_decoder(
                .operand_i(operands_a[(i)*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_a[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );

            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ub_posit_decoder(
                .operand_i(operands_b[(i)*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_b[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate


    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;

    // Extract valid sign, exponent, and mantissa from previous output (acc)
    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;
    
    posit_decoder #(
        .n(n_o),
        .es(es_o)
    )uc_posit_decoder(
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    // ---------------
    // Sign Process
    // ---------------
    logic [N-1:0] signs_ab;
    logic [N:0] signs;
    
    // Perform XOR to obtain the signs of the product of Va and Vb
    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc, signs_ab};

    // ---------------
    // Exponent Addition
    // ---------------
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    
    // obtain the exponents of the product of Va and Vb through signed addition
    logic signed [N*(EXP_WIDTH+1)-1:0] rg_exp_c;
    generate
        genvar m;
        for(m=0;m<N;m++) begin
            assign rg_exp_c[(m)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(rg_exp_a[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]) + $signed(rg_exp_b[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
        end
    endgenerate

    // ********************
    // Pipeline 1
    // ********************
    // group_path -name pipe1to2 -from {pipe1*} -to {pipe2*}
    logic [N:0] pipe1_signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe1_rg_exp_c;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;

    `FFNR(pipe1_signs,signs,clk)
    `FFNR(pipe1_rg_exp_c,rg_exp_c,clk)
    `FFNR(pipe1_mants_norm_a,mants_norm_a,clk)
    `FFNR(pipe1_mants_norm_b,mants_norm_b,clk)
    `FFNR(pipe1_rg_exp_acc,rg_exp_acc,clk)
    `FFNR(pipe1_mant_norm_acc,mant_norm_acc,clk)

    // ---------------
    // Mantissa Multiplication
    // ---------------
    localparam int unsigned MUL_WIDTH = 2*(MANT_WIDTH_I+1);
    logic [N*MUL_WIDTH-1:0] mul_sum, mul_carry;

    // Mantissa multiplication is performed by a modified radix-4 booth multiplier
    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            radix4_booth_multiplier #(
                .WIDTH_A(MANT_WIDTH_I+1),
                .WIDTH_B(MANT_WIDTH_I+1)
            ) u_radix4_booth_multiplier(
                .operand_a(pipe1_mants_norm_a[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .operand_b(pipe1_mants_norm_b[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .sum_o(mul_sum[(j)*MUL_WIDTH +: MUL_WIDTH]),
                .carry_o(mul_carry[(j)*MUL_WIDTH +: MUL_WIDTH])
            );
        end
    endgenerate

    // ---------------
    // Exponent Compare
    // ---------------
    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_items;
    generate
        genvar u;
        for(u=0;u<N;u++) begin
            assign rg_exp_items[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = pipe1_rg_exp_c[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate
    assign rg_exp_items[(N)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc);

    // Obtain the maximum exponent by a recursive comparator tree
    logic signed [EXP_WIDTH:0] rg_exp_max;
    comp_tree #(
        .N(N+1),
        .WIDTH(EXP_WIDTH)
    ) u_comp_tree(
        .operands_i(rg_exp_items),
        .result_o(rg_exp_max)
    );

    // *******************
    // Pipeline 2
    // *******************
    // group_path -name pipe2to3 -from {pipe2*} -to {pipe3*}
    logic [N:0] pipe2_signs;
    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] pipe2_rg_exp_items;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_sum, pipe2_mul_carry;
    
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;


    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_items,rg_exp_items,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe2_mul_sum,mul_sum,clk)
    `FFNR(pipe2_mul_carry,mul_carry,clk)
    
    `FFNR(pipe2_mant_norm_acc,pipe1_mant_norm_acc,clk)

    // ---------------
    // The final addition
    // ---------------
    logic [N*MUL_WIDTH-1:0] mants_norm_c;
    generate
        genvar v;
        for(v=0;v<N;v++) begin
            assign mants_norm_c[(v)*MUL_WIDTH +: MUL_WIDTH] = pipe2_mul_sum[(v)*MUL_WIDTH +: MUL_WIDTH] + pipe2_mul_carry[(v)*MUL_WIDTH +: MUL_WIDTH];
        end
    endgenerate

    // ---------------
    // Mantissa Alignment
    // ---------------
    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;

    // the mantissa products of Va and Vb are kept in values of bitwidth `ALIGN_WIDTH` before alignment
    generate
        genvar k;
        if(ALIGN_WIDTH>MUL_WIDTH) begin: fixed_shift
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = mants_norm_c[(k)*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
            end
        end  
        else begin: fixed_shift
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = mants_norm_c[(k)*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH-ALIGN_WIDTH);
            end
        end
    endgenerate

    if(ALIGN_WIDTH>MANT_WIDTH_O+2) begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
    end
    else begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc >> (MANT_WIDTH_O+2-ALIGN_WIDTH);
    end

    // Align the products according to the difference between the respective exponent and the maximum exponent
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    logic [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z,s;
        // compute the respective exponent difference
        for(z=0;z<N+1;z++) begin
            assign rg_exp_diff[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $unsigned(pipe2_rg_exp_max - pipe2_rg_exp_items[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]);
        end
        // the maximum shift amount is ALIGN_WIDTH, since the bits exceeding ALIGN_WIDTH will be discarded directly.
        if(EXP_WIDTH+1>SHIFT_WIDTH) begin
            for(s=0;s<N+1;s++) begin
                assign shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH] = (|rg_exp_diff[(s)*(EXP_WIDTH+1) + SHIFT_WIDTH +: (EXP_WIDTH-SHIFT_WIDTH+1)]) ? ALIGN_WIDTH : rg_exp_diff[(s)*(EXP_WIDTH+1) +: SHIFT_WIDTH];
                barrel_shifter #(
                    .WIDTH(ALIGN_WIDTH),
                    .SHIFT_WIDTH(SHIFT_WIDTH),
                    .MODE(1'b1)
                ) u_barrel_shifter(
                    .operand_i(product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]),
                    .shift_amount(shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH]),
                    .result_o(product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH])
                );
            end
        end

        else begin
            for(s=0;s<N+1;s++) begin
                barrel_shifter #(
                    .WIDTH(ALIGN_WIDTH),
                    .SHIFT_WIDTH(EXP_WIDTH+1),
                    .MODE(1'b1)
                ) u_barrel_shifter(
                    .operand_i(product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]),
                    .shift_amount(rg_exp_diff[(s)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                    .result_o(product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH])
                );
            end
        end
    endgenerate

    // ---------------
    // Mantissa accumulation in two's complement
    // ---------------
    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;

    logic [(N+1)*(SUM_WIDTH+1)-1:0] mantissa, mantissa_comp;
    // Convert to two's complement
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin
            assign mantissa[(y)*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] = product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH];
            assign mantissa_comp[(y)*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] = pipe2_signs[y] ? (~mantissa[(y)*(SUM_WIDTH+1) +: (SUM_WIDTH+1)]+1) : mantissa[(y)*(SUM_WIDTH+1) +: (SUM_WIDTH+1)];
        end
    endgenerate

    // *******************
    // Pipeline 3
    // *******************
    // group_path -name pipe3to4 -from {pipe3*} -to {pipe4*}
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;
    logic [(N+1)*(SUM_WIDTH+1)-1:0] pipe3_mantissa_comp;

    `FFNR(pipe3_rg_exp_max,pipe2_rg_exp_max,clk)
    `FFNR(pipe3_mantissa_comp,mantissa_comp,clk)

    logic [SUM_WIDTH:0] csa_sum,csa_carry;
    csa_tree #(
        .N(N+1),
        .WIDTH_I(SUM_WIDTH+1),
        .WIDTH_O(SUM_WIDTH+1)
    ) u_csa_tree(
        .operands_i(pipe3_mantissa_comp),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );

    logic [SUM_WIDTH:0] sum_result;
    // the final addition
    assign sum_result = csa_sum + csa_carry;


    // ********************
    // Pipeline 4
    // ********************
    // group_path -name pipe4to5 -from {pipe4*} -to {pipe5*}
    logic [SUM_WIDTH:0] pipe4_sum_result;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;

    `FFNR(pipe4_sum_result,sum_result,clk)
    `FFNR(pipe4_rg_exp_max,pipe3_rg_exp_max,clk)


    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    
    assign final_sign = pipe4_sum_result[SUM_WIDTH];
    assign sum_c = final_sign ? (~pipe4_sum_result+1) : pipe4_sum_result[SUM_WIDTH-1:0];

    // ---------------
    // Mantissa normalization and exponent adjustment
    // ---------------
    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic signed [EXP_WIDTH:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;

    mantissa_norm #(
        .WIDTH(SUM_WIDTH),
        .EXP_WIDTH(EXP_WIDTH),
        .DECIMAL_POINT(CARRY_WIDTH+2)
    ) u_mantissa_norm(
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    // Note: is there any possibility of overflow?
    assign final_rg_exp = pipe4_rg_exp_max + rg_exp_adjust;

    // ---------------
    // Posit encoder and Rounding
    // ---------------
    logic [MANT_WIDTH_O+2:0] final_mant;
    
    if(SUM_WIDTH>MANT_WIDTH_O+3) begin
        logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
        logic sticky_bit;
        assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
        assign sticky_bit = |sticky_bits;
        assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2],sticky_bit};
    end
    else begin
        assign final_mant = sum_norm << (MANT_WIDTH_O+3-SUM_WIDTH);
    end


    // ********************
    // Pipeline 5
    // ********************
    // group_path -name pipe5to6 -from {pipe5*} -to {result_o}
    logic pipe5_final_sign;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;

    `FFNR(pipe5_final_sign,final_sign,clk)
    `FFNR(pipe5_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe5_final_mant,final_mant,clk)


    posit_encoder #(
        .n(n_o),
        .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );


    // ********************
    // Pipeline 6: Output Reg
    // ********************
    // No group_path
endmodule

// ---------------------------------------------------------------------------
// Parallel-prefix (Sklansky) adder: sum = a + b + cin
// Much shallower than the ripple chain yosys infers from `a + b`.
// ---------------------------------------------------------------------------
module prefix_add #(
    parameter int unsigned WIDTH = 16
)(
    input  logic [WIDTH-1:0] a,
    input  logic [WIDTH-1:0] b,
    input  logic             cin,
    output logic [WIDTH-1:0] sum,
    output logic             cout
);
    localparam int unsigned LV = pdpu_pkg::clog2(WIDTH);

    logic [WIDTH-1:0] gi, pi, pt;
    assign gi = a & b;
    assign pi = a | b;
    assign pt = a ^ b;

    logic [(LV+1)*WIDTH-1:0] gt, ptx;
    assign gt[WIDTH-1:0] = gi;
    assign ptx[WIDTH-1:0] = pi;

    generate
        genvar l, i;
        for (l = 1; l <= LV; l++) begin : levels
            for (i = 0; i < WIDTH; i++) begin : bits
                if (((i) % (2**l)) < (2**(l-1))) begin : white
                    assign gt[l*WIDTH+i] = gt[(l-1)*WIDTH+i];
                    assign ptx[l*WIDTH+i] = ptx[(l-1)*WIDTH+i];
                end
                else begin : black
                    assign gt[l*WIDTH+i] = gt[(l-1)*WIDTH+i] | (ptx[(l-1)*WIDTH+i] & gt[(l-1)*WIDTH+(i - (i % (2**l)) + (2**(l-1)) - 1)]);
                    assign ptx[l*WIDTH+i] = ptx[(l-1)*WIDTH+i] & ptx[(l-1)*WIDTH+(i - (i % (2**l)) + (2**(l-1)) - 1)];
                end
            end
        end
    endgenerate

    logic [WIDTH-1:0] go, ppo, carry_in;
    assign go = gt[LV*WIDTH +: WIDTH];
    assign ppo = ptx[LV*WIDTH +: WIDTH];

    assign carry_in[0] = cin;
    generate
        genvar j;
        for (j = 1; j < WIDTH; j++) begin : cin_chain
            assign carry_in[j] = go[j-1] | (ppo[j-1] & cin);
        end
    endgenerate

    assign sum = pt ^ carry_in;
    assign cout = go[WIDTH-1] | (ppo[WIDTH-1] & cin);
endmodule

// ===========================================================================
// Optimized 5-stage pipelined PDPU dot-16.
// Function is IDENTICAL to the upstream pdpu_top_pipelined (renamed
// pdpu_top_pipelined_old above) bit for bit:
//   * same posit decode / exponent add / ALIGN_WIDTH=14 alignment truncation /
//     accumulate / normalize / RNE encode arithmetic
// Micro-architectural changes (latency stays 5 cycles):
//   * mantissa products computed with `*`; the booth CSA tree's sum+carry that
//     the datapath used is exactly the 8-bit unsigned product
//   * sign folding: aligned terms are INVERTED instead of two's-complemented;
//     the +1 per negative term is collected as a popcount and enters the CSA
//     tree as one extra operand (removes 17 x 20-bit incrementers)
//   * exponent differences + shift clamping computed ahead of pipe2, so pipe2
//     carries 4-bit shift amounts instead of 7-bit exponents
//   * pipe2 carries products (8b) instead of booth sum+carry pairs (2x8b)
//   * the encoded posit is registered in pipe5 (result_o is a pure register)
module pdpu_top_pipelined_v1 #(
    parameter int unsigned N = 16,                  // dot-product size
    parameter int unsigned n_i = 8,                 // word size
    parameter int unsigned es_i = 2,                // exponent size
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14         // alignment width
)(
    input logic clk,
    input logic [N*n_i-1:0] operands_a,        // input vector (Va)
    input logic [N*n_i-1:0] operands_b,        // input vector (Vb)
    input logic [n_o-1:0] acc,                      // previous output (acc)
    output logic [n_o-1:0] result_o,                // output result (out)
    output logic [299:0] probe_o
);
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;
    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    localparam int unsigned MUL_WIDTH = 2*(MANT_WIDTH_I+1);
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;
    localparam int unsigned EXT_W = SUM_WIDTH + 1;   // accumulator width

    // ---------------- Cloud A: posit decoding, exp add, mantissa mul --------
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ua_posit_decoder(
                .operand_i(operands_a[(i)*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_a[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ub_posit_decoder(
                .operand_i(operands_b[(i)*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_b[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    posit_decoder #(
        .n(n_o),
        .es(es_o)
    )uc_posit_decoder(
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    logic [N-1:0] signs_ab;
    logic [N:0] signs;
    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc, signs_ab};

    logic signed [N*(EXP_WIDTH+1)-1:0] rg_exp_c;
    logic [N*MUL_WIDTH-1:0] mul_prod;
    generate
        genvar m;
        for(m=0;m<N;m++) begin: exp_add_mul
            assign rg_exp_c[(m)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(rg_exp_a[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]) + $signed(rg_exp_b[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
            assign mul_prod[(m)*MUL_WIDTH +: MUL_WIDTH] = mants_norm_a[(m)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] * mants_norm_b[(m)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)];
        end
    endgenerate

    // ********************
    // Pipeline 1
    // ********************
    logic [N:0] pipe1_signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe1_rg_exp_c;
    logic [N*MUL_WIDTH-1:0] pipe1_mul_prod;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;

    `FFNR(pipe1_signs,signs,clk)
    `FFNR(pipe1_rg_exp_c,rg_exp_c,clk)
    `FFNR(pipe1_mul_prod,mul_prod,clk)
    `FFNR(pipe1_rg_exp_acc,rg_exp_acc,clk)
    `FFNR(pipe1_mant_norm_acc,mant_norm_acc,clk)

    // ---------------- Cloud B: exponent max, diff/clamp, sign popcount ------
    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_items;
    generate
        genvar u;
        for(u=0;u<N;u++) begin: exp_items
            assign rg_exp_items[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = pipe1_rg_exp_c[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate
    assign rg_exp_items[(N)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc);

    logic signed [EXP_WIDTH:0] rg_exp_max;
    comp_tree #(
        .N(N+1),
        .WIDTH(EXP_WIDTH)
    ) u_comp_tree(
        .operands_i(rg_exp_items),
        .result_o(rg_exp_max)
    );

    logic [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    logic [N:0] kill_amount;
    generate
        genvar z;
        for(z=0;z<N+1;z++) begin: align_amount
            wire [EXP_WIDTH:0] exp_item_z = rg_exp_items[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
            prefix_add #(
                .WIDTH(EXP_WIDTH+1)
            ) u_diff_add(
                .a(rg_exp_max),
                .b(~exp_item_z),
                .cin(1'b1),
                .sum(rg_exp_diff[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .cout()
            );
            // The reference clamps the shift to ALIGN_WIDTH when the difference
            // exceeds the field; since shifting ALIGN_WIDTH bits right by
            // ALIGN_WIDTH or more already yields zero, only the >= 2**SHIFT_WIDTH
            // case must be forced to zero. The kill flag is computed in parallel
            // with the low shift bits and applied after the barrel shifter, so
            // no clamp logic sits on the shift-amount path.
            if(EXP_WIDTH+1>SHIFT_WIDTH) begin: diff_kill
                assign kill_amount[z] = |rg_exp_diff[(z)*(EXP_WIDTH+1) + SHIFT_WIDTH +: (EXP_WIDTH-SHIFT_WIDTH+1)];
            end
            else begin: no_kill
                assign kill_amount[z] = 1'b0;
            end
            assign shift_amount[(z)*SHIFT_WIDTH +: SHIFT_WIDTH] = rg_exp_diff[(z)*(EXP_WIDTH+1) +: SHIFT_WIDTH];
        end
    endgenerate

    // popcount of the N+1 term signs (balanced); its +1s replace the
    // two's-complement incrementers of the term formation
    logic [8*2-1:0] pc_l1;
    logic [4*3-1:0] pc_l2;
    logic [2*4-1:0] pc_l3;
    generate
        genvar g;
        for(g=0;g<8;g++) begin: pc_pair
            if(2*g+1 <= N) assign pc_l1[g*2 +: 2] = pipe1_signs[2*g] + pipe1_signs[2*g+1];
            else if(2*g <= N) assign pc_l1[g*2 +: 2] = pipe1_signs[2*g];
            else assign pc_l1[g*2 +: 2] = 2'b0;
        end
        for(g=0;g<4;g++) begin: pc_quad
            assign pc_l2[g*3 +: 3] = pc_l1[(2*g)*2 +: 2] + pc_l1[(2*g+1)*2 +: 2];
        end
        for(g=0;g<2;g++) begin: pc_oct
            assign pc_l3[g*4 +: 4] = pc_l2[(2*g)*3 +: 3] + pc_l2[(2*g+1)*3 +: 3];
        end
    endgenerate
    logic [CARRY_WIDTH-1:0] pipe1_pc_val;
    assign pipe1_pc_val = pc_l3[0 +: 4] + pc_l3[4 +: 4] + pipe1_signs[N];

    // ********************
    // Pipeline 2
    // ********************
    logic [N:0] pipe2_signs;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_prod;
    logic [(N+1)*SHIFT_WIDTH-1:0] pipe2_shift_amount;
    logic [N:0] pipe2_kill;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [CARRY_WIDTH-1:0] pipe2_sign_pc;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_mul_prod,pipe1_mul_prod,clk)
    `FFNR(pipe2_shift_amount,shift_amount,clk)
    `FFNR(pipe2_kill,kill_amount,clk)
    `FFNR(pipe2_mant_norm_acc,pipe1_mant_norm_acc,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe2_sign_pc,pipe1_pc_val,clk)

    // ---------------- Cloud C: align, sign fold, CSA accumulate ------------
    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;

    generate
        genvar k;
        if(ALIGN_WIDTH>MUL_WIDTH) begin: fixed_shift_hi
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
            end
        end
        else begin: fixed_shift_lo
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH-ALIGN_WIDTH);
            end
        end
    endgenerate

    generate
        if(ALIGN_WIDTH>MANT_WIDTH_O+2) begin: acc_shift_hi
            assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
        end
        else begin: acc_shift_lo
            assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc >> (MANT_WIDTH_O+2-ALIGN_WIDTH);
        end
    endgenerate

    generate
        genvar s;
        for(s=0;s<N+1;s++) begin: align
            barrel_shifter #(
                .WIDTH(ALIGN_WIDTH),
                .SHIFT_WIDTH(SHIFT_WIDTH),
                .MODE(1'b1)
            ) u_barrel_shifter(
                .operand_i(product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(pipe2_shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH])
            );
        end
    endgenerate

    // sign folding: invert instead of negate, the +1s are pipe2_sign_pc
    logic [(N+1)*EXT_W-1:0] mantissa_comp;
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin: sign_fold
            wire [EXT_W-1:0] ext_raw = {{(EXT_W-ALIGN_WIDTH){1'b0}}, product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH]};
            wire [EXT_W-1:0] ext_term = pipe2_kill[y] ? {EXT_W{1'b0}} : ext_raw;
            assign mantissa_comp[(y)*EXT_W +: EXT_W] = pipe2_signs[y] ? ~ext_term : ext_term;
        end
    endgenerate

    logic [(N+2)*EXT_W-1:0] csa_operands;
    generate
        genvar y2;
        for(y2=0;y2<N+1;y2++) begin: csa_in
            assign csa_operands[(y2)*EXT_W +: EXT_W] = mantissa_comp[(y2)*EXT_W +: EXT_W];
        end
    endgenerate
    assign csa_operands[(N+1)*EXT_W +: EXT_W] = {{(EXT_W-CARRY_WIDTH){1'b0}}, pipe2_sign_pc};

    logic [SUM_WIDTH:0] csa_sum, csa_carry;
    csa_tree #(
        .N(N+2),
        .WIDTH_I(EXT_W),
        .WIDTH_O(EXT_W)
    ) u_csa_tree(
        .operands_i(csa_operands),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );

    logic [SUM_WIDTH:0] sum_result;
    prefix_add #(
        .WIDTH(SUM_WIDTH+1)
    ) u_final_add(
        .a(csa_sum),
        .b(csa_carry),
        .cin(1'b0),
        .sum(sum_result),
        .cout()
    );

    // ********************
    // Pipeline 3
    // ********************
    logic [SUM_WIDTH:0] pipe3_sum_result;
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;

    `FFNR(pipe3_sum_result,sum_result,clk)
    `FFNR(pipe3_rg_exp_max,pipe2_rg_exp_max,clk)

    // ---------------- Cloud D: abs, normalize, exponent adjust -------------
    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;

    assign final_sign = pipe3_sum_result[SUM_WIDTH];
    logic [SUM_WIDTH:0] sum_neg;
    prefix_add #(
        .WIDTH(SUM_WIDTH+1)
    ) u_abs_add(
        .a(~pipe3_sum_result),
        .b({(SUM_WIDTH+1){1'b0}}),
        .cin(1'b1),
        .sum(sum_neg),
        .cout()
    );
    assign sum_c = final_sign ? sum_neg[SUM_WIDTH-1:0] : pipe3_sum_result[SUM_WIDTH-1:0];

    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic signed [EXP_WIDTH:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;

    mantissa_norm #(
        .WIDTH(SUM_WIDTH),
        .EXP_WIDTH(EXP_WIDTH),
        .DECIMAL_POINT(CARRY_WIDTH+2)
    ) u_mantissa_norm(
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    prefix_add #(
        .WIDTH(EXP_WIDTH+1)
    ) u_exp_add(
        .a(pipe3_rg_exp_max),
        .b(rg_exp_adjust),
        .cin(1'b0),
        .sum(final_rg_exp),
        .cout()
    );

    logic [MANT_WIDTH_O+2:0] final_mant;
    generate
        if(SUM_WIDTH>MANT_WIDTH_O+3) begin: sticky_fold
            logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
            logic sticky_bit;
            assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
            assign sticky_bit = |sticky_bits;
            assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2],sticky_bit};
        end
        else begin: final_shift
            assign final_mant = sum_norm << (MANT_WIDTH_O+3-SUM_WIDTH);
        end
    endgenerate

    // ********************
    // Pipeline 4
    // ********************
    logic pipe4_final_sign;
    logic signed [EXP_WIDTH:0] pipe4_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe4_final_mant;

    `FFNR(pipe4_final_sign,final_sign,clk)
    `FFNR(pipe4_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe4_final_mant,final_mant,clk)

    // ---------------- Cloud E: posit encoding + RNE rounding ---------------
    logic [n_o-1:0] enc_result;
    area_probe u_probe(.clk(clk), .din(sum_result), .dout(probe_o));
    posit_encoder #(
        .n(n_o),
        .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder(
        .sign_i(pipe4_final_sign),
        .rg_exp_i(pipe4_final_rg_exp),
        .mant_norm_i(pipe4_final_mant),
        .result_o(enc_result)
    );

    // ********************
    // Pipeline 5: Output Reg
    // ********************
    `FFNR(result_o,enc_result,clk)
endmodule

// ---- AREA PROBE (temporary): 300 dummy flops, used on the output ----
module area_probe(
    input logic clk,
    input logic [19:0] din,
    output logic [299:0] dout
);
    logic [299:0] dummy_q;
    wire [299:0] dummy_d = {15{din}};
    always_ff @(posedge clk) dummy_q <= dummy_d;
    assign dout = dummy_q;
endmodule

// ---------------------------------------------------------------------------
// Compact posit decoder. Bit-identical function to posit_decoder_generic
// (including the zero/NaR quirks: runlength = n-1 when the magnitude is zero,
// implicit bit = |operand_i[n-2:0]), but built from a one-hot run-length decode
// instead of lzc + barrel shifter: much smaller and shallower.
// ---------------------------------------------------------------------------
module posit_decoder #(
    parameter int unsigned n = 8,
    parameter int unsigned es = 2,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd+es,
    parameter int unsigned MANT_WIDTH = n-es-3
)(
    input  logic [n-1:0] operand_i,
    output logic sign_o,
    output logic signed [EXP_WIDTH:0] rg_exp_o,
    output logic [MANT_WIDTH:0] mant_norm_o
);
    localparam int unsigned W = n-1;
    localparam int unsigned MF_W = W-es-2;      // kept fraction bits

    assign sign_o = operand_i[n-1];

    logic [W-1:0] mag_neg;
    prefix_add #(
        .WIDTH(W)
    ) u_neg (
        .a(~operand_i[W-1:0]),
        .b({W{1'b0}}),
        .cin(1'b1),
        .sum(mag_neg),
        .cout()
    );
    logic [W-1:0] value;
    assign value = sign_o ? mag_neg : operand_i[W-1:0];
    wire regS = value[W-1];

    // one-hot run length: sel[r-1] == 1  <=>  regime run length r (r = 1..W)
    //   regS=1: run of leading ones;  regS=0: run of leading zeros (value==0 -> r=W)
    logic [W-1:0] sel;
    generate
        genvar r, i;
        for (r = 1; r < W; r++) begin : g_sel
            assign sel[r-1] = regS ? ((&value[W-1:W-r]) & ~value[W-1-r]) :
                                     ((~|value[W-1:W-r]) &  value[W-1-r]);
        end
    endgenerate
    assign sel[W-1] = regS ? &value : ~|value;

    // per-run field options: rg_exp = {k, e} with k = regS ? r-1 : -r and
    // e / mf taken from value at the run-dependent positions
    // (equivalent to op_no_rg = value << (r+1), e = op_no_rg[W-1:W-es],
    //  mf = op_no_rg[W-1-es:2])
    logic [(EXP_WIDTH+1)*W-1:0] rgexp_opts;
    logic [MF_W*W-1:0]          mf_opts;
    generate
        for (r = 1; r <= W; r++) begin : g_opts
            for (i = 0; i < es; i++) begin : g_exp
                if ((r+1) <= (W-es+i)) begin : g_e1
                    assign rgexp_opts[(r-1)*(EXP_WIDTH+1) + i] = value[W-es+i-r-1];
                end
                else begin : g_e0
                    assign rgexp_opts[(r-1)*(EXP_WIDTH+1) + i] = 1'b0;
                end
            end
            wire signed [EXP_WIDTH:0] k_part = regS ? ((r-1) << es) : ((-r) << es);
            assign rgexp_opts[(r-1)*(EXP_WIDTH+1) + es +: (nd+1)] = k_part[EXP_WIDTH:es];
            for (i = 0; i < MF_W; i++) begin : g_mf
                if (r <= (i+1)) begin : g_m1
                    assign mf_opts[(r-1)*MF_W + i] = value[i+1-r];
                end
                else begin : g_m0
                    assign mf_opts[(r-1)*MF_W + i] = 1'b0;
                end
            end
        end
    endgenerate

    // one-hot OR-mux
    logic [MF_W-1:0] mf_sel;
    always_comb begin
        rg_exp_o = '0;
        mf_sel = '0;
        for (int rr = 0; rr < W; rr++) begin
            rg_exp_o = rg_exp_o | (rgexp_opts[rr*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] & {(EXP_WIDTH+1){sel[rr]}});
            mf_sel  = mf_sel  | (mf_opts[rr*MF_W +: MF_W] & {MF_W{sel[rr]}});
        end
    end
    assign mant_norm_o = {|operand_i[W-1:0], mf_sel};
endmodule

// ---------------------------------------------------------------------------
// Dadda-style 3:2 reduction tree: reduces N W-bit operands to OUTN W-bit
// operands (all arithmetic mod 2**W). Shallow: one full-adder layer per level.
// The reduction is split across two pipeline stages (OUTN=6 then OUTN=2).
// ---------------------------------------------------------------------------
module csa_red #(
    parameter int unsigned N = 18,
    parameter int unsigned W = 20,
    parameter int unsigned OUTN = 2
)(
    input  logic [N*W-1:0] operands_i,
    output logic [OUTN*W-1:0] result_o
);
    localparam int unsigned K = (2*N+2)/3;   // ceil(2N/3)
    generate
        if (N == OUTN) begin : g_pass
            assign result_o = operands_i;
        end
        else if (N < OUTN) begin : g_pad
            assign result_o = {{(OUTN-N)*W{1'b0}}, operands_i};
        end
        else begin : g_red
            localparam int unsigned NFA = N - K;
            localparam int unsigned LEFT = 3*K - 2*N;
            logic [K*W-1:0] lv;
            genvar j;
            for (j = 0; j < NFA; j++) begin : g_fa
                compressor_3to2 #(
                    .WIDTH_I(W),
                    .WIDTH_O(W)
                ) u_c32(
                    .operands_i(operands_i[(3*j)*W +: 3*W]),
                    .sum_o(lv[(2*j)*W +: W]),
                    .carry_o(lv[(2*j+1)*W +: W])
                );
            end
            for (j = 0; j < LEFT; j++) begin : g_rest
                assign lv[(2*NFA+j)*W +: W] = operands_i[(3*NFA+j)*W +: W];
            end
            csa_red #(
                .N(K),
                .W(W),
                .OUTN(OUTN)
            ) u_next(
                .operands_i(lv),
                .result_o(result_o)
            );
        end
    endgenerate
endmodule

// ===========================================================================
// pdpu_top_pipelined (v2): same function as the upstream combinational
// pdpu_dot16, 5-cycle latency. Pipeline cuts re-balanced: flops are nearly
// free in this technology, so wide intermediate state is carried between the
// stages to keep every combinational cloud shallow:
//   cloud A: posit decode, mantissa multiply, exponent add (carry-save), popcount
//   cloud B: exponent add (final), exponent maximum tree
//   cloud C: exponent diff/kill, align shift, sign fold, CSA reduction part 1
//   cloud D: CSA reduction part 2, final add, absolute value
//   cloud E: normalize, exponent adjust, sticky / final mantissa
//   cloud F: posit encode + RNE rounding (result_o is combinational)
// ===========================================================================
module pdpu_top_pipelined #(
    parameter int unsigned N = 16,                  // dot-product size
    parameter int unsigned n_i = 8,                 // word size
    parameter int unsigned es_i = 2,                // exponent size
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14         // alignment width
)(
    input logic clk,
    input logic [N*n_i-1:0] operands_a,
    input logic [N*n_i-1:0] operands_b,
    input logic [n_o-1:0] acc,
    output logic [n_o-1:0] result_o,
    output logic [351:0] probe_o
);
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;
    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    localparam int unsigned MUL_WIDTH = 2*(MANT_WIDTH_I+1);
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;
    localparam int unsigned EXT_W = SUM_WIDTH + 1;

    // ---------------- Cloud A: decode, mul, exp add (CSA), popcount ------
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ua_posit_decoder(
                .operand_i(operands_a[(i)*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_a[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(
                .n(n_i),
                .es(es_i)
            )ub_posit_decoder(
                .operand_i(operands_b[(i)*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_b[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    posit_decoder #(
        .n(n_o),
        .es(es_o)
    )uc_posit_decoder(
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    logic [N-1:0] signs_ab;
    logic [N:0] signs;
    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc, signs_ab};

    // exponent addition in carry-save form (final add in cloud B)
    logic [N*(EXP_WIDTH+1)-1:0] exp_add_s, exp_add_c;
    logic [N*MUL_WIDTH-1:0] mul_prod;
    generate
        genvar m;
        for(m=0;m<N;m++) begin: exp_add_mul
            wire signed [EXP_WIDTH:0] ea = $signed(rg_exp_a[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
            wire signed [EXP_WIDTH:0] eb = $signed(rg_exp_b[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
            assign exp_add_s[(m)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = ea ^ eb;
            assign exp_add_c[(m)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = {ea & eb, 1'b0};
            assign mul_prod[(m)*MUL_WIDTH +: MUL_WIDTH] = mants_norm_a[(m)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] * mants_norm_b[(m)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)];
        end
    endgenerate

    // popcount of the N+1 term signs: its +1s replace the two's-complement
    // incrementers of the term formation (sign folding)
    logic [8*2-1:0] pc_l1;
    logic [4*3-1:0] pc_l2;
    logic [2*4-1:0] pc_l3;
    generate
        genvar g;
        for(g=0;g<8;g++) begin: pc_pair
            if(2*g+1 <= N) assign pc_l1[g*2 +: 2] = signs[2*g] + signs[2*g+1];
            else if(2*g <= N) assign pc_l1[g*2 +: 2] = signs[2*g];
            else assign pc_l1[g*2 +: 2] = 2'b0;
        end
        for(g=0;g<4;g++) begin: pc_quad
            assign pc_l2[g*3 +: 3] = pc_l1[(2*g)*2 +: 2] + pc_l1[(2*g+1)*2 +: 2];
        end
        for(g=0;g<2;g++) begin: pc_oct
            assign pc_l3[g*4 +: 4] = pc_l2[(2*g)*3 +: 3] + pc_l2[(2*g+1)*3 +: 3];
        end
    endgenerate
    logic [CARRY_WIDTH-1:0] sign_pc_val;
    assign sign_pc_val = pc_l3[0 +: 4] + pc_l3[4 +: 4] + signs[N];

    // ********************
    // Pipeline 1
    // ********************
    logic [N:0] pipe1_signs;
    logic [N*(EXP_WIDTH+1)-1:0] pipe1_exp_add_s, pipe1_exp_add_c;
    logic [N*MUL_WIDTH-1:0] pipe1_mul_prod;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;
    logic [CARRY_WIDTH-1:0] pipe1_sign_pc;

    `FFNR(pipe1_signs,signs,clk)
    `FFNR(pipe1_exp_add_s,exp_add_s,clk)
    `FFNR(pipe1_exp_add_c,exp_add_c,clk)
    `FFNR(pipe1_mul_prod,mul_prod,clk)
    `FFNR(pipe1_rg_exp_acc,rg_exp_acc,clk)
    `FFNR(pipe1_mant_norm_acc,mant_norm_acc,clk)
    `FFNR(pipe1_sign_pc,sign_pc_val,clk)

    // ---------------- Cloud B: exponent add (final), exponent max -------
    logic signed [N*(EXP_WIDTH+1)-1:0] rg_exp_c;
    generate
        genvar m2;
        for(m2=0;m2<N;m2++) begin: exp_add_final
            prefix_add #(
                .WIDTH(EXP_WIDTH+1)
            ) u_exp_sum(
                .a(pipe1_exp_add_s[(m2)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .b(pipe1_exp_add_c[(m2)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .cin(1'b0),
                .sum(rg_exp_c[(m2)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .cout()
            );
        end
    endgenerate

    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_items;
    generate
        genvar u;
        for(u=0;u<N;u++) begin: exp_items
            assign rg_exp_items[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = rg_exp_c[(u)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] ^ (1 << EXP_WIDTH);
        end
    endgenerate
    assign rg_exp_items[(N)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc) ^ (1 << EXP_WIDTH);

    // Partial exponent max: four groups (5,5,5,2) reduced in cloud B, the
    // final 4-input max runs in cloud C (flops are free, keeps both shallow)
    localparam int unsigned NGRP = 4;
    logic signed [NGRP*(EXP_WIDTH+1)-1:0] grp_max;
    comp_tree #(
        .N(5),
        .WIDTH(EXP_WIDTH)
    ) u_g0(
        .operands_i(rg_exp_items[0*(EXP_WIDTH+1) +: 5*(EXP_WIDTH+1)]),
        .result_o(grp_max[0*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );
    comp_tree #(
        .N(5),
        .WIDTH(EXP_WIDTH)
    ) u_g1(
        .operands_i(rg_exp_items[5*(EXP_WIDTH+1) +: 5*(EXP_WIDTH+1)]),
        .result_o(grp_max[1*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );
    comp_tree #(
        .N(5),
        .WIDTH(EXP_WIDTH)
    ) u_g2(
        .operands_i(rg_exp_items[10*(EXP_WIDTH+1) +: 5*(EXP_WIDTH+1)]),
        .result_o(grp_max[2*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );
    comp_tree #(
        .N(2),
        .WIDTH(EXP_WIDTH)
    ) u_g3(
        .operands_i(rg_exp_items[15*(EXP_WIDTH+1) +: 2*(EXP_WIDTH+1)]),
        .result_o(grp_max[3*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );
    // second max level also runs ahead of pipe2
    logic signed [2*(EXP_WIDTH+1)-1:0] pair_max;
    comp_tree #(
        .N(2),
        .WIDTH(EXP_WIDTH)
    ) u_p0(
        .operands_i(grp_max[0*(EXP_WIDTH+1) +: 2*(EXP_WIDTH+1)]),
        .result_o(pair_max[0*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );
    comp_tree #(
        .N(2),
        .WIDTH(EXP_WIDTH)
    ) u_p1(
        .operands_i(grp_max[2*(EXP_WIDTH+1) +: 2*(EXP_WIDTH+1)]),
        .result_o(pair_max[1*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
    );

    // ********************
    // Pipeline 2
    // ********************
    logic [N:0] pipe2_signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe2_rg_exp_c;
    logic signed [2*(EXP_WIDTH+1)-1:0] pipe2_pair_max;
    logic signed [EXP_WIDTH_O:0] pipe2_rg_exp_acc;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_prod;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;
    logic [CARRY_WIDTH-1:0] pipe2_sign_pc;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_c,rg_exp_c,clk)
    `FFNR(pipe2_pair_max,pair_max,clk)
    `FFNR(pipe2_rg_exp_acc,pipe1_rg_exp_acc,clk)
    `FFNR(pipe2_mul_prod,pipe1_mul_prod,clk)
    `FFNR(pipe2_mant_norm_acc,pipe1_mant_norm_acc,clk)
    `FFNR(pipe2_sign_pc,pipe1_sign_pc,clk)

    // ---------------- Cloud C: diff, kill, align, sign fold, CSA part 1 -
    logic signed [EXP_WIDTH:0] rg_exp_max_off;
    comp_tree #(
        .N(2),
        .WIDTH(EXP_WIDTH)
    ) u_comp_tree(
        .operands_i(pipe2_pair_max),
        .result_o(rg_exp_max_off)
    );
    logic signed [EXP_WIDTH:0] rg_exp_max;
    assign rg_exp_max = rg_exp_max_off ^ (1 << EXP_WIDTH);
    logic [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    logic [N:0] kill_amount;
    generate
        genvar z;
        for(z=0;z<N+1;z++) begin: align_amount
            wire signed [EXP_WIDTH:0] exp_item_z = (z<N) ? pipe2_rg_exp_c[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] : $signed(pipe2_rg_exp_acc);
            prefix_add #(
                .WIDTH(EXP_WIDTH+1)
            ) u_diff_add(
                .a(rg_exp_max),
                .b(~exp_item_z),
                .cin(1'b1),
                .sum(rg_exp_diff[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .cout()
            );
            if(EXP_WIDTH+1>SHIFT_WIDTH) begin: diff_kill
                assign kill_amount[z] = |rg_exp_diff[(z)*(EXP_WIDTH+1) + SHIFT_WIDTH +: (EXP_WIDTH-SHIFT_WIDTH+1)];
            end
            else begin: no_kill
                assign kill_amount[z] = 1'b0;
            end
            assign shift_amount[(z)*SHIFT_WIDTH +: SHIFT_WIDTH] = rg_exp_diff[(z)*(EXP_WIDTH+1) +: SHIFT_WIDTH];
        end
    endgenerate

    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;

    generate
        genvar k;
        if(ALIGN_WIDTH>MUL_WIDTH) begin: fixed_shift_hi
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
            end
        end
        else begin: fixed_shift_lo
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH-ALIGN_WIDTH);
            end
        end
    endgenerate

    generate
        if(ALIGN_WIDTH>MANT_WIDTH_O+2) begin: acc_shift_hi
            assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
        end
        else begin: acc_shift_lo
            assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc >> (MANT_WIDTH_O+2-ALIGN_WIDTH);
        end
    endgenerate

    generate
        genvar s;
        for(s=0;s<N+1;s++) begin: align
            barrel_shifter #(
                .WIDTH(ALIGN_WIDTH),
                .SHIFT_WIDTH(SHIFT_WIDTH),
                .MODE(1'b1)
            ) u_barrel_shifter(
                .operand_i(product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH])
            );
        end
    endgenerate

    logic [(N+2)*EXT_W-1:0] csa_operands;
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin: sign_fold
            wire [EXT_W-1:0] ext_raw = {{(EXT_W-ALIGN_WIDTH){1'b0}}, product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH]};
            wire [EXT_W-1:0] ext_term = kill_amount[y] ? {EXT_W{1'b0}} : ext_raw;
            assign csa_operands[(y)*EXT_W +: EXT_W] = pipe2_signs[y] ? ~ext_term : ext_term;
        end
    endgenerate
    assign csa_operands[(N+1)*EXT_W +: EXT_W] = {{(EXT_W-CARRY_WIDTH){1'b0}}, pipe2_sign_pc};

    localparam int unsigned CSA_MID = 6;
    logic [CSA_MID*EXT_W-1:0] csa_mid;
    csa_red #(
        .N(N+2),
        .W(EXT_W),
        .OUTN(CSA_MID)
    ) u_csa_top(
        .operands_i(csa_operands),
        .result_o(csa_mid)
    );

    // ********************
    // Pipeline 3
    // ********************
    logic [CSA_MID*EXT_W-1:0] pipe3_csa_mid;
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;

    `FFNR(pipe3_csa_mid,csa_mid,clk)
    `FFNR(pipe3_rg_exp_max,rg_exp_max,clk)

    // ---------------- Cloud D: CSA part 2, final add, absolute value ----
    logic [(2)*EXT_W-1:0] csa_fin;
    csa_red #(
        .N(CSA_MID),
        .W(EXT_W),
        .OUTN(2)
    ) u_csa_bot(
        .operands_i(pipe3_csa_mid),
        .result_o(csa_fin)
    );

    // sum_result = a+b; the negated form -(a+b) = ~a + ~b + 2 is computed in
    // parallel and selected by the sign (fuses final add and absolute value)
    logic [SUM_WIDTH:0] sum_result, sum_neg_t;
    prefix_add #(
        .WIDTH(SUM_WIDTH+1)
    ) u_final_add(
        .a(csa_fin[0 +: EXT_W]),
        .b(csa_fin[EXT_W +: EXT_W]),
        .cin(1'b0),
        .sum(sum_result),
        .cout()
    );
    prefix_add #(
        .WIDTH(SUM_WIDTH+1)
    ) u_neg_add(
        .a(~csa_fin[0 +: EXT_W]),
        .b(~csa_fin[EXT_W +: EXT_W]),
        .cin(1'b0),
        .sum(sum_neg_t),
        .cout()
    );

    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    assign final_sign = sum_result[SUM_WIDTH];
    wire [SUM_WIDTH:0] sum_neg = sum_neg_t + 2;
    assign sum_c = final_sign ? sum_neg[SUM_WIDTH-1:0] : sum_result[SUM_WIDTH-1:0];

    // ********************
    // Pipeline 4
    // ********************
    logic [SUM_WIDTH-1:0] pipe4_sum_c;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;
    logic pipe4_final_sign;

    `FFNR(pipe4_sum_c,sum_c,clk)
    `FFNR(pipe4_rg_exp_max,pipe3_rg_exp_max,clk)
    `FFNR(pipe4_final_sign,final_sign,clk)

    // ---------------- Cloud E: normalize, exponent adjust, sticky -------
    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic signed [EXP_WIDTH:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;

    mantissa_norm #(
        .WIDTH(SUM_WIDTH),
        .EXP_WIDTH(EXP_WIDTH),
        .DECIMAL_POINT(CARRY_WIDTH+2)
    ) u_mantissa_norm(
        .operand_i(pipe4_sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    prefix_add #(
        .WIDTH(EXP_WIDTH+1)
    ) u_exp_add(
        .a(pipe4_rg_exp_max),
        .b(rg_exp_adjust),
        .cin(1'b0),
        .sum(final_rg_exp),
        .cout()
    );

    logic [MANT_WIDTH_O+2:0] final_mant;
    generate
        if(SUM_WIDTH>MANT_WIDTH_O+3) begin: sticky_fold
            logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
            logic sticky_bit;
            assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
            assign sticky_bit = |sticky_bits;
            assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2],sticky_bit};
        end
        else begin: final_shift
            assign final_mant = sum_norm << (MANT_WIDTH_O+3-SUM_WIDTH);
        end
    endgenerate

    // ********************
    // Pipeline 5
    // ********************
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic pipe5_final_sign;

    `FFNR(pipe5_final_mant,final_mant,clk)
    `FFNR(pipe5_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe5_final_sign,pipe4_final_sign,clk)

    // ---------------- Cloud F: posit encode + RNE rounding --------------
    posit_encoder #(
        .n(n_o),
        .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );
endmodule

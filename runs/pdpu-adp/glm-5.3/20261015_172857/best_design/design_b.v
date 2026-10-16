// PDPU dot-16 pipelined, optimized: 6 comb chunks over 5 register stages.
// Key changes vs upstream: XOR+popcount sign-magnitude trick, FA-tree CSA,
// plain 4x4 multipliers, acc decoder moved into chunk 2.
`define FFNR(__q, __d, __clk)        \
  always_ff @(posedge (__clk)) begin \
    __q <= (__d);                    \
  end

package pdpu_pkg;
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

package cf_math_pkg;
    function automatic integer clog(input integer n);
        begin
            n = n-1;
            for(clog=0; n>0; clog=clog+1)
                n = n>>1;
        end
    endfunction
    function automatic integer unsigned idx_width (input integer unsigned num_idx);
        idx_width = (num_idx > 32'd1) ? $unsigned(clog(num_idx)) : 32'd1;
    endfunction
endpackage

module barrel_shifter #(
    parameter int unsigned WIDTH = 8,
    parameter int unsigned SHIFT_WIDTH = 3,
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
        if(MODE==1'b0) begin
            for(i=SHIFT_WIDTH-1;i>0;i--) begin
                assign temp_results[i-1] = shift_amount[i] ? temp_results[i]<<(2**i) : temp_results[i];
            end
            assign result_o = shift_amount[0] ? temp_results[0]<<1 : temp_results;
        end
        else begin
            for(i=SHIFT_WIDTH-1;i>0;i--) begin
                assign temp_results[i-1] = shift_amount[i] ? temp_results[i]>>(2**i) : temp_results[i];
            end
            assign result_o = shift_amount[0] ? temp_results[0]>>1 : temp_results;
        end
    endgenerate
endmodule

module lzc #(
  parameter int unsigned WIDTH = 2,
  parameter bit          MODE  = 1'b0,
  parameter int unsigned CNT_WIDTH = cf_math_pkg::idx_width(WIDTH)
) (
  input  logic [WIDTH-1:0]     in_i,
  output logic [CNT_WIDTH-1:0] cnt_o,
  output logic                 empty_o
);
  if (WIDTH == 1) begin : gen_degenerate_lzc
    assign cnt_o[0] = !in_i[0];
    assign empty_o = !in_i[0];
  end else begin : gen_lzc
    localparam int unsigned NumLevels = cf_math_pkg::clog(WIDTH);
    logic [WIDTH*NumLevels-1:0] index_lut;
    logic [2**NumLevels-1:0] sel_nodes;
    logic [(2**NumLevels)*NumLevels-1:0] index_nodes;
    logic [WIDTH-1:0] in_tmp;
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
          if ($unsigned(k) * 2 < WIDTH - 1) begin : g_reduce
            assign sel_nodes[2 ** level - 1 + k] = in_tmp[k * 2] | in_tmp[k * 2 + 1];
            assign index_nodes[(2 ** level - 1 + k)*NumLevels +: NumLevels] = (in_tmp[k * 2] == 1'b1)
              ? index_lut[(k * 2)*NumLevels +: NumLevels] :
                index_lut[(k * 2 + 1)*NumLevels +: NumLevels];
          end
          if ($unsigned(k) * 2 == WIDTH - 1) begin : g_base
            assign sel_nodes[2 ** level - 1 + k] = in_tmp[k * 2];
            assign index_nodes[(2 ** level - 1 + k)*NumLevels +: NumLevels] = index_lut[(k * 2)*NumLevels +: NumLevels];
          end
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
endmodule

module comparator #(parameter int unsigned WIDTH = 8)(
    input logic signed [WIDTH:0] operand_a,
    input logic signed [WIDTH:0] operand_b,
    output logic signed [WIDTH:0] result_o
);
    logic [1:0] sign_ab;
    logic [WIDTH-1:0] data_a;
    logic [WIDTH-1:0] data_b;
    assign sign_ab = {operand_a[WIDTH],operand_b[WIDTH]};
    assign data_a = operand_a[WIDTH-1:0];
    assign data_b = operand_b[WIDTH-1:0];
    always_comb begin
        case(sign_ab)
            2'b10: result_o = operand_b;
            2'b01: result_o = operand_a;
            default: result_o = (data_a>data_b) ? operand_a : operand_b;
        endcase
    end
endmodule

module comp_tree #(parameter N = 4, parameter WIDTH = 8)(
    input logic signed [N*(WIDTH+1)-1:0] operands_i,
    output logic signed [WIDTH:0] result_o
);
    localparam int unsigned W   = WIDTH + 1;
    localparam int unsigned N_A = N/2;
    localparam int unsigned N_B = N - N_A;
    generate
        if (N==1) begin : g_leaf
            assign result_o = operands_i[0 +: W];
        end
        else if (N==2) begin : g_pair
            comparator #(.WIDTH(WIDTH)) u_comparator(
                .operand_a(operands_i[0 +: W]),
                .operand_b(operands_i[W +: W]),
                .result_o(result_o));
        end
        else begin : g_split
            logic signed [N_A*W-1:0] operands_i_A;
            logic signed [N_B*W-1:0] operands_i_B;
            logic signed [WIDTH:0]   result_o_A;
            logic signed [WIDTH:0]   result_o_B;
            assign operands_i_A = operands_i[0       +: N_A*W];
            assign operands_i_B = operands_i[N_A*W   +: N_B*W];
            comp_tree #(.N(N_A), .WIDTH(WIDTH)) ua_comp_tree(
                .operands_i(operands_i_A), .result_o(result_o_A));
            comp_tree #(.N(N_B), .WIDTH(WIDTH)) ub_comp_tree(
                .operands_i(operands_i_B), .result_o(result_o_B));
            comparator #(.WIDTH(WIDTH)) uc_comparator(
                .operand_a(result_o_A), .operand_b(result_o_B), .result_o(result_o));
        end
    endgenerate
endmodule

module mantissa_norm #(
    parameter int unsigned WIDTH = 8,
    parameter int unsigned EXP_WIDTH = 3,
    parameter int unsigned DECIMAL_POINT = 3
)(
    input logic [WIDTH-1:0] operand_i,
    output logic signed [EXP_WIDTH:0] exp_adjust,
    output logic [WIDTH-1:0] result_o
);
    localparam int unsigned LZC_WIDTH = pdpu_pkg::clog2(WIDTH);
    logic [LZC_WIDTH-1:0] leading_zero_count;
    logic lzc_zeroes;
    lzc #(.WIDTH(WIDTH), .MODE(1'b1)) u_lzc(
        .in_i(operand_i), .cnt_o(leading_zero_count), .empty_o(lzc_zeroes));
    always_comb begin
        if(lzc_zeroes) exp_adjust = '0;
        else if(leading_zero_count <= DECIMAL_POINT - 1) exp_adjust = DECIMAL_POINT - leading_zero_count - 1;
        else exp_adjust = -$signed(leading_zero_count - DECIMAL_POINT + 1);
    end
    barrel_shifter #(.WIDTH(WIDTH), .SHIFT_WIDTH(LZC_WIDTH), .MODE(1'b0)) u_barrel_shifter(
        .operand_i(operand_i), .shift_amount(leading_zero_count), .result_o(result_o));
endmodule

module posit_decoder #(
    parameter int unsigned n = 16,
    parameter int unsigned es = 1,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd+es,
    parameter int unsigned MANT_WIDTH = n-es-3
)(
    input logic [n-1:0] operand_i,
    output logic sign_o,
    output logic signed [EXP_WIDTH:0] rg_exp_o,
    output logic [MANT_WIDTH:0] mant_norm_o
);
    logic sign;
    logic [n-2:0] operand_value;
    assign sign = operand_i[n-1];
    assign operand_value = sign ? (~operand_i[n-2:0]+1) : operand_i[n-2:0];

    logic regS;
    logic [n-2:0] lzc_operand;
    logic [nd-1:0] leading_zero_count;
    logic lzc_zeroes;
    assign regS = operand_value[n-2];
    assign lzc_operand = regS ? (~operand_value) : operand_value;
    lzc #(.WIDTH(n-1), .MODE(1'b1)) u_lzc(
        .in_i(lzc_operand), .cnt_o(leading_zero_count), .empty_o(lzc_zeroes));

    logic [nd-1:0] runlength;
    logic [nd-1:0] regime_bits;
    logic signed [nd:0] regime_k;
    assign runlength = lzc_zeroes ? (n-1) : leading_zero_count;
    assign regime_bits = lzc_zeroes ? (n-1) : (leading_zero_count+1);
    assign regime_k = regS ? ({1'b0,runlength-1}) : ({1'b1,~runlength+1});

    logic [n-2:0] op_no_rg;
    barrel_shifter #(.WIDTH(n-1), .SHIFT_WIDTH(nd), .MODE(1'b0)) u_barrel_shifter(
        .operand_i(operand_value), .shift_amount(regime_bits), .result_o(op_no_rg));

    assign sign_o = sign;
    logic [es:0] exp;
    if(es==0) begin
        assign exp = 0;
    end
    else begin
        assign exp = op_no_rg[n-2:n-2-es+1];
    end
    assign rg_exp_o = regime_k << es | exp;

    logic implicit_bit;
    assign implicit_bit = |operand_i[n-2:0];
    assign mant_norm_o = {implicit_bit, op_no_rg[n-2-es:2]};
endmodule

module posit_encoder #(
    parameter int unsigned n = 16,
    parameter int unsigned es = 1,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd+es,
    parameter int unsigned MANT_WIDTH = n-es-3
) (
    input logic sign_i,
    input logic signed [EXP_WIDTH:0] rg_exp_i,
    input logic [MANT_WIDTH:0] mant_norm_i,
    output logic [n-1:0] result_o
);
    logic input_not_zero;
    assign input_not_zero = mant_norm_i[MANT_WIDTH];

    logic signed [EXP_WIDTH-es:0] regime_k;
    logic signed [es:0] exp;
    assign regime_k = rg_exp_i[EXP_WIDTH:es];
    if(es==0) begin
        assign exp = 0;
    end
    else begin
        assign exp = rg_exp_i[es-1:0];
    end

    logic sign_k;
    logic [n-2:0] rg_const;
    logic [n-2:0] regime;
    assign sign_k = rg_exp_i[EXP_WIDTH];
    assign rg_const = 1;
    assign regime = sign_k ? rg_const : ~rg_const;

    logic [EXP_WIDTH-es:0] regime_bits;
    assign regime_bits = sign_k ? (~regime_k+2) : (regime_k+2);

    localparam int unsigned TEMP_WIDTH = (n-1)+es+MANT_WIDTH;
    logic [TEMP_WIDTH-1:0] rg_exp_mant;
    if(es==0) begin
        assign rg_exp_mant = {regime, mant_norm_i[MANT_WIDTH-1:0]};
    end
    else begin
        assign rg_exp_mant = {regime, exp[es-1:0], mant_norm_i[MANT_WIDTH-1:0]};
    end

    localparam int unsigned MAX_SHIFT_AMOUNT = MANT_WIDTH + es + 1;
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(MAX_SHIFT_AMOUNT+1);
    logic [SHIFT_WIDTH-1:0] shift_amount;
    assign shift_amount = (regime_bits>=n) ? MAX_SHIFT_AMOUNT : (regime_bits+(MANT_WIDTH+es-n+1));

    logic [TEMP_WIDTH+MAX_SHIFT_AMOUNT-1:0] value_before_shift, value_after_shift;
    assign value_before_shift = rg_exp_mant << MAX_SHIFT_AMOUNT;
    barrel_shifter #(.WIDTH(TEMP_WIDTH+MAX_SHIFT_AMOUNT), .SHIFT_WIDTH(SHIFT_WIDTH), .MODE(1'b1)) u_barrel_shifter(
        .operand_i(value_before_shift), .shift_amount(shift_amount), .result_o(value_after_shift));

    logic [n-2:0] value_before_round;
    logic [MAX_SHIFT_AMOUNT-1:0] rounding_bits;
    assign {value_before_round,rounding_bits} = value_after_shift[MAX_SHIFT_AMOUNT+n-2:0];

    logic round_bit;
    logic sticky_bit;
    logic round_value;
    logic [n-2:0] value_after_round;
    assign round_bit = rounding_bits[MAX_SHIFT_AMOUNT-1];
    assign sticky_bit = |rounding_bits[MAX_SHIFT_AMOUNT-2:0];
    assign round_value = round_bit & (sticky_bit | value_before_round[0]);
    assign value_after_round = value_before_round + round_value;

    logic [n-1:0] normal_result;
    assign normal_result = sign_i ? {1'b1,~value_after_round+1} : {1'b0,value_after_round};
    assign result_o = input_not_zero ? normal_result : '0;
endmodule

// ---- radix-4 booth multiplier (upstream structure, must match reference bit-exactly) ----
module fulladder(input logic x,y,z, output logic sum,carry);
    assign sum = (x ^ y ) ^ z;
    assign carry = ((x ^ y) & z) | (x & y);
endmodule

module booth_encoder(input [2:0] code, output neg, zero, one, two);
    assign neg = code[2]&(~code[0]) | code[2]&(~code[1]);
    assign zero = ~(|code) | &code;
    assign two = ~code[2]&code[1]&code[0] | code[2]&(~code[1])&(~code[0]);
    assign one = code[1]&(~code[0]) | (~code[1])&code[0];
endmodule

module gen_product #(parameter int unsigned WIDTH = 16)(
    input logic [WIDTH-1:0] multiplicand,
    input logic [2:0] code,
    output logic [WIDTH:0] partial_prod,
    output sign
);
    logic neg,zero,one,two;
    booth_encoder u_booth_encoder(.code(code), .neg(neg), .zero(zero), .one(one), .two(two));
    logic [WIDTH:0] temp_prod;
    always_comb begin
        if(one) temp_prod = multiplicand;
        else if(two) temp_prod = multiplicand << 1;
        else temp_prod = '0;
    end
    assign partial_prod = neg ? (~temp_prod) : temp_prod;
    assign sign = neg;
endmodule

module gen_prods #(
    parameter int unsigned WIDTH_A = 16,
    parameter int unsigned WIDTH_B = 16,
    parameter int unsigned COUNT = (WIDTH_B+2)/2,
    parameter int unsigned WIDTH_O = WIDTH_A + WIDTH_B
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
    gen_product #(.WIDTH(WIDTH_A)) ua_gen_product(
        .multiplicand(operand_a), .code(codes[0]),
        .partial_prod(temp_prods[0]), .sign(signs[0]));
    assign partial_prods[0] = {~signs[0], signs[0], signs[0], temp_prods[0]};

    generate
        genvar i;
        for(i=1;i<COUNT-1;i++) begin
            assign codes[i] = multiplier[2*i+2:2*i];
            gen_product #(.WIDTH(WIDTH_A)) ub_gen_product(
                .multiplicand(operand_a), .code(codes[i]),
                .partial_prod(temp_prods[i]), .sign(signs[i]));
            assign partial_prods[i] = {1'b1,~signs[i],temp_prods[i],1'b0,signs[i-1]} << (2*i-2);
        end
    endgenerate

    assign codes[COUNT-1] = multiplier[2*COUNT:2*COUNT-2];
    gen_product #(.WIDTH(WIDTH_A)) uc_gen_product(
        .multiplicand(operand_a), .code(codes[COUNT-1]),
        .partial_prod(temp_prods[COUNT-1]), .sign(signs[COUNT-1]));
    assign partial_prods[COUNT-1] = {temp_prods[COUNT-1],1'b0,signs[COUNT-2]} << (2*COUNT-4);
endmodule

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
                .carry(carry[i]));
        end
    endgenerate
    assign sum_o = sum;
    assign carry_o = {1'b0,carry,1'b0};
endmodule

module csa_tree #(
    parameter int unsigned N = 8,
    parameter int unsigned WIDTH_I = 8,
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(N)
)(
    input logic [N*WIDTH_I-1:0] operands_i,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
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
            compressor_3to2 #(.WIDTH_I(WIDTH_I), .WIDTH_O(WIDTH_O)) u_compressor_3to2(
                .operands_i(operands_i), .sum_o(sum_o), .carry_o(carry_o));
        end
    endgenerate
endmodule

module radix4_booth_multiplier #(
    parameter int unsigned WIDTH_A = 16,
    parameter int unsigned WIDTH_B = 16,
    parameter int unsigned WIDTH_O = WIDTH_A + WIDTH_B
)(
    input logic [WIDTH_A-1:0] operand_a,
    input logic [WIDTH_B-1:0] operand_b,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    localparam int unsigned COUNT = (WIDTH_B+2)/2;
    logic [COUNT-1:0][WIDTH_O-1:0] partial_prods;
    gen_prods #(.WIDTH_A(WIDTH_A), .WIDTH_B(WIDTH_B)) u_gen_prods(
        .operand_a(operand_a), .operand_b(operand_b), .partial_prods(partial_prods));
    csa_tree #(.N(COUNT), .WIDTH_I(WIDTH_O), .WIDTH_O(WIDTH_O)) u_csa_tree(
        .operands_i(partial_prods), .sum_o(sum_o), .carry_o(carry_o));
endmodule

// 3:2 compressor row of full adders (no ripple)
module fa_row #(parameter int unsigned W = 20)(
    input logic [W-1:0] a, b, c,
    output logic [W-1:0] s,
    output logic [W-1:0] cout
);
    assign s = a ^ b ^ c;
    assign cout = (a & b) | (c & (a ^ b));
endmodule

module pdpu_top_pipelined #(
    parameter int unsigned N = 16,
    parameter int unsigned n_i = 8,
    parameter int unsigned es_i = 2,
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14
)(
    input logic clk,
    input logic [N*n_i-1:0] operands_a,
    input logic [N*n_i-1:0] operands_b,
    input logic [n_o-1:0] acc,
    output logic [n_o-1:0] result_o
);
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;
    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    localparam int unsigned EW = EXP_WIDTH+1;
    localparam int unsigned MW = MANT_WIDTH_I+1;
    localparam int unsigned MUL_WIDTH = 2*MW;
    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;
    localparam int unsigned SW = SUM_WIDTH+1;   // 20

    // ---- chunk 1: posit decoding of the 32 inputs + product exponents ----
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*EW-1:0] rg_exp_a, rg_exp_b;
    logic [N*MW-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(.n(n_i), .es(es_i)) ua_posit_decoder(
                .operand_i(operands_a[(i)*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[(i)*EW +: EW]),
                .mant_norm_o(mants_norm_a[(i)*MW +: MW]));
            posit_decoder #(.n(n_i), .es(es_i)) ub_posit_decoder(
                .operand_i(operands_b[(i)*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[(i)*EW +: EW]),
                .mant_norm_o(mants_norm_b[(i)*MW +: MW]));
        end
    endgenerate

    logic [N:0] signs;
    logic signed [N*EW-1:0] rg_exp_c;
    generate
        genvar m;
        for(m=0;m<N;m++) begin
            assign rg_exp_c[(m)*EW +: EW] = $signed(rg_exp_a[(m)*EW +: EW]) + $signed(rg_exp_b[(m)*EW +: EW]);
        end
    endgenerate

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;
    posit_decoder #(.n(n_o), .es(es_o)) uc_posit_decoder(
        .operand_i(acc), .sign_o(sign_acc), .rg_exp_o(rg_exp_acc), .mant_norm_o(mant_norm_acc));
    assign signs = {sign_acc, signs_a ^ signs_b};

    // ---- pipe 1 ----
    logic [N:0] pipe1_signs;
    logic signed [N*EW-1:0] pipe1_rg_exp_c;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [N*MW-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;
    `FFNR(pipe1_signs,signs,clk)
    `FFNR(pipe1_rg_exp_c,rg_exp_c,clk)
    `FFNR(pipe1_rg_exp_acc,rg_exp_acc,clk)
    `FFNR(pipe1_mants_norm_a,mants_norm_a,clk)
    `FFNR(pipe1_mants_norm_b,mants_norm_b,clk)
    `FFNR(pipe1_mant_norm_acc,mant_norm_acc,clk)

    // ---- chunk 2: mantissa multiplication, acc decode, exponent max ----
    logic [N*MUL_WIDTH-1:0] mul_sum, mul_carry;
    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            radix4_booth_multiplier #(.WIDTH_A(MW), .WIDTH_B(MW)) u_radix4_booth_multiplier(
                .operand_a(pipe1_mants_norm_a[(j)*MW +: MW]),
                .operand_b(pipe1_mants_norm_b[(j)*MW +: MW]),
                .sum_o(mul_sum[(j)*MUL_WIDTH +: MUL_WIDTH]),
                .carry_o(mul_carry[(j)*MUL_WIDTH +: MUL_WIDTH]));
        end
        for(j=0;j<N;j++) begin
            assign rg_exp_items[(j)*EW +: EW] = pipe1_rg_exp_c[(j)*EW +: EW];
        end
    endgenerate
    assign rg_exp_items[(N)*EW +: EW] = $signed(pipe1_rg_exp_acc);

    logic signed [(N+1)*EW-1:0] rg_exp_items;
    logic signed [EXP_WIDTH:0] rg_exp_max;
    comp_tree #(.N(N+1), .WIDTH(EXP_WIDTH)) u_comp_tree(
        .operands_i(rg_exp_items), .result_o(rg_exp_max));

    // ---- pipe 2 ----
    logic [N:0] pipe2_signs;
    logic signed [(N+1)*EW-1:0] pipe2_rg_exp_items;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_sum, pipe2_mul_carry;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;
    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_items,rg_exp_items,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe2_mul_sum,mul_sum,clk)
    `FFNR(pipe2_mul_carry,mul_carry,clk)
    `FFNR(pipe2_mant_norm_acc,pipe1_mant_norm_acc,clk)

    // ---- chunk 3: fixed shift + alignment right shift ----
    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;
    generate
        genvar k;
        for(k=0;k<N;k++) begin
            assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = (pipe2_mul_sum[(k)*MUL_WIDTH +: MUL_WIDTH] + pipe2_mul_carry[(k)*MUL_WIDTH +: MUL_WIDTH]) << (ALIGN_WIDTH-MUL_WIDTH);
        end
    endgenerate
    assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);

    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    logic [(N+1)*EW-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z,s;
        for(z=0;z<N+1;z++) begin
            assign rg_exp_diff[(z)*EW +: EW] = $unsigned(pipe2_rg_exp_max - pipe2_rg_exp_items[(z)*EW +: EW]);
        end
        for(s=0;s<N+1;s++) begin
            assign shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH] =
                (|rg_exp_diff[(s)*EW + SHIFT_WIDTH +: (EW-SHIFT_WIDTH)]) ? ALIGN_WIDTH
                : rg_exp_diff[(s)*EW +: SHIFT_WIDTH];
            barrel_shifter #(.WIDTH(ALIGN_WIDTH), .SHIFT_WIDTH(SHIFT_WIDTH), .MODE(1'b1)) u_barrel_shifter(
                .operand_i(product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH]));
        end
    endgenerate

    // ---- pipe 3 ----
    logic [N:0] pipe3_signs;
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;
    logic [(N+1)*ALIGN_WIDTH-1:0] pipe3_product_shifted;
    `FFNR(pipe3_signs,pipe2_signs,clk)
    `FFNR(pipe3_rg_exp_max,pipe2_rg_exp_max,clk)
    `FFNR(pipe3_product_shifted,product_shifted,clk)

    // ---- chunk 4: sign-magnitude to two's complement (XOR + popcount) + CSA + final add ----
    logic [(N+1)*SW-1:0] xr;
    logic [CARRY_WIDTH-1:0] nneg;
    // DEBUG: upstream per-lane two's complement instead of XOR+popcount
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin
            assign xr[(y)*SW +: SW] = pipe3_signs[y]
                ? (~{{(SW-ALIGN_WIDTH){1'b0}}, pipe3_product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH]}+1)
                :  {{(SW-ALIGN_WIDTH){1'b0}}, pipe3_product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH]};
        end
    endgenerate
    // count of negative lanes (0..17)
    logic [4:0] pc1, pc2, pc3, pc4;
    always_comb begin
        pc1 = 5'd0;
        pc2 = 5'd0;
        pc3 = 5'd0;
        pc4 = 5'd0;
        nneg = 5'd0;
    end

    // Wallace FA tree: 17 x 20-bit operands + popcount term (18 total)
    logic [SW-1:0] a0 [0:16];
    generate
        for (genvar g=0; g<N+1; g++) begin : g_a0
            assign a0[g] = xr[(g)*SW +: SW];
        end
    endgenerate
    logic [SW-1:0] popc;
    assign popc = {{(SW-CARRY_WIDTH){1'b0}}, nneg};

    // level A: 18 -> 12
    logic [SW-1:0] A_s0,A_c0, A_s1,A_c1, A_s2,A_c2, A_s3,A_c3, A_s4,A_c4, A_s5,A_c5;
    fa_row #(.W(SW)) uA0(.a(a0[0]),  .b(a0[1]),  .c(a0[2]),  .s(A_s0), .cout(A_c0));
    fa_row #(.W(SW)) uA1(.a(a0[3]),  .b(a0[4]),  .c(a0[5]),  .s(A_s1), .cout(A_c1));
    fa_row #(.W(SW)) uA2(.a(a0[6]),  .b(a0[7]),  .c(a0[8]),  .s(A_s2), .cout(A_c2));
    fa_row #(.W(SW)) uA3(.a(a0[9]),  .b(a0[10]), .c(a0[11]), .s(A_s3), .cout(A_c3));
    fa_row #(.W(SW)) uA4(.a(a0[12]), .b(a0[13]), .c(a0[14]), .s(A_s4), .cout(A_c4));
    fa_row #(.W(SW)) uA5(.a(a0[15]), .b(a0[16]), .c(popc),   .s(A_s5), .cout(A_c5));
    // carries must be shifted left one bit position when passed to next level
    logic [SW-1:0] A_c0s,A_c1s,A_c2s,A_c3s,A_c4s,A_c5s;
    assign {A_c0s,A_c1s} = {{A_c0[SW-2:0],1'b0},{A_c1[SW-2:0],1'b0}};
    assign {A_c2s,A_c3s} = {{A_c2[SW-2:0],1'b0},{A_c3[SW-2:0],1'b0}};
    assign {A_c4s,A_c5s} = {{A_c4[SW-2:0],1'b0},{A_c5[SW-2:0],1'b0}};
    // level B: 12 -> 8
    logic [SW-1:0] B_s0,B_c0, B_s1,B_c1, B_s2,B_c2, B_s3,B_c3;
    fa_row #(.W(SW)) uB0(.a(A_s0), .b(A_c0s), .c(A_s1), .s(B_s0), .cout(B_c0));
    fa_row #(.W(SW)) uB1(.a(A_c1s), .b(A_s2), .c(A_c2s), .s(B_s1), .cout(B_c1));
    fa_row #(.W(SW)) uB2(.a(A_s3), .b(A_c3s), .c(A_s4), .s(B_s2), .cout(B_c2));
    fa_row #(.W(SW)) uB3(.a(A_c4s), .b(A_s5), .c(A_c5s), .s(B_s3), .cout(B_c3));
    // level C: 8 -> 6
    logic [SW-1:0] B_c0s,B_c1s,B_c2s,B_c3s;
    assign B_c0s = {B_c0[SW-2:0],1'b0};
    assign B_c1s = {B_c1[SW-2:0],1'b0};
    assign B_c2s = {B_c2[SW-2:0],1'b0};
    assign B_c3s = {B_c3[SW-2:0],1'b0};
    logic [SW-1:0] C_s0,C_c0, C_s1,C_c1;
    fa_row #(.W(SW)) uC0(.a(B_s0), .b(B_c0s), .c(B_s1), .s(C_s0), .cout(C_c0));
    fa_row #(.W(SW)) uC1(.a(B_c1s), .b(B_s2), .c(B_c2s), .s(C_s1), .cout(C_c1));
    // level D: 6 -> 4
    logic [SW-1:0] C_c0s,C_c1s;
    assign C_c0s = {C_c0[SW-2:0],1'b0};
    assign C_c1s = {C_c1[SW-2:0],1'b0};
    logic [SW-1:0] D_s0,D_c0, D_s1,D_c1;
    fa_row #(.W(SW)) uD0(.a(C_s0), .b(C_c0s), .c(C_s1), .s(D_s0), .cout(D_c0));
    fa_row #(.W(SW)) uD1(.a(C_c1s), .b(B_s3), .c(B_c3s), .s(D_s1), .cout(D_c1));
    // level E: 4 -> 3
    logic [SW-1:0] D_c0s,D_c1s;
    assign D_c0s = {D_c0[SW-2:0],1'b0};
    assign D_c1s = {D_c1[SW-2:0],1'b0};
    logic [SW-1:0] E_s,E_c;
    fa_row #(.W(SW)) uE0(.a(D_s0), .b(D_c0s), .c(D_s1), .s(E_s), .cout(E_c));
    // level F: 3 -> 2
    logic [SW-1:0] E_cs;
    assign E_cs = {E_c[SW-2:0],1'b0};
    logic [SW-1:0] F_s,F_c;
    fa_row #(.W(SW)) uF0(.a(E_s), .b(E_cs), .c(D_c1s), .s(F_s), .cout(F_c));

    logic [SUM_WIDTH:0] sum_result;
    assign sum_result = $signed(F_s) + $signed({F_c[SW-2:0],1'b0});

    // ---- pipe 4 ----
    logic [SUM_WIDTH:0] pipe4_sum_result;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;
    `FFNR(pipe4_sum_result,sum_result,clk)
    `FFNR(pipe4_rg_exp_max,pipe3_rg_exp_max,clk)

    // ---- chunk 5: normalize ----
    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    assign final_sign = pipe4_sum_result[SUM_WIDTH];
    assign sum_c = final_sign ? (~pipe4_sum_result+1) : pipe4_sum_result[SUM_WIDTH-1:0];

    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic signed [EXP_WIDTH:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;
    mantissa_norm #(.WIDTH(SUM_WIDTH), .EXP_WIDTH(EXP_WIDTH), .DECIMAL_POINT(CARRY_WIDTH+2)) u_mantissa_norm(
        .operand_i(sum_c), .exp_adjust(rg_exp_adjust), .result_o(sum_norm));
    assign final_rg_exp = pipe4_rg_exp_max + rg_exp_adjust;

    logic [MANT_WIDTH_O+2:0] final_mant;
    if(SUM_WIDTH>MANT_WIDTH_O+3) begin
        logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
        assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
        assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2],|sticky_bits};
    end
    else begin
        assign final_mant = sum_norm << (MANT_WIDTH_O+3-SUM_WIDTH);
    end

    // ---- pipe 5 ----
    logic pipe5_final_sign;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;
    `FFNR(pipe5_final_sign,final_sign,clk)
    `FFNR(pipe5_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe5_final_mant,final_mant,clk)

    // ---- chunk 6: posit encoder ----
    posit_encoder #(.n(n_o), .es(es_o), .EXP_WIDTH(EXP_WIDTH), .MANT_WIDTH(MANT_WIDTH_O+2)) u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o));
endmodule

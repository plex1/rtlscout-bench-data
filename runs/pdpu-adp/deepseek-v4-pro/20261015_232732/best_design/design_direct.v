// PDPU dot-16 5-stage pipelined: simplified arithmetic implementation
// Intended to be bit-exact with the original but use behavioral +/*/shift/max/sum.

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

// Behavioral barrel shifter: equivalent to the staged shifter in the original.
module barrel_shifter #(
    parameter int unsigned WIDTH = 8,
    parameter int unsigned SHIFT_WIDTH = 3,
    parameter bit MODE = 1'b0
)(
    input  logic [WIDTH-1:0]      operand_i,
    input  logic [SHIFT_WIDTH-1:0] shift_amount,
    output logic [WIDTH-1:0]      result_o
);
    // Arithmetic right shift, which is what the original staged shifter
    // synthesizes to (right shifting a signed/unsigned expression here matches).
    assign result_o = (MODE == 1'b0) ? (operand_i << shift_amount)
                                     : (operand_i >> shift_amount);
endmodule

// Original tree-priority lzc (copied from starting_point).
module lzc #(
  parameter int unsigned WIDTH = 2,
  parameter bit          MODE  = 1'b0,
  parameter int unsigned CNT_WIDTH = pdpu_pkg::clog2(WIDTH)
) (
  input  logic [WIDTH-1:0]     in_i,
  output logic [CNT_WIDTH-1:0] cnt_o,
  output logic                 empty_o
);
  if (WIDTH == 1) begin : gen_degenerate_lzc
    assign cnt_o[0] = !in_i[0];
    assign empty_o = !in_i[0];
  end else begin : gen_lzc
    localparam int unsigned NumLevels = pdpu_pkg::clog2(WIDTH);
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

    assign cnt_o = NumLevels > $unsigned(0) ? index_nodes[(0)*NumLevels +: NumLevels] : {(pdpu_pkg::clog2(WIDTH)) {1'b0}};
    assign empty_o = NumLevels > $unsigned(0) ? ~sel_nodes[0] : ~(|in_i);
  end : gen_lzc
endmodule : lzc

// Exact original comparator: sign/magnitude style used by the PDPU tree.
module comparator #(
    parameter int unsigned WIDTH = 8
)(
    input  logic signed [WIDTH:0] operand_a,
    input  logic signed [WIDTH:0] operand_b,
    output logic signed [WIDTH:0] result_o
);
    logic [1:0] sign_ab;
    logic [WIDTH-1:0] data_a;
    logic [WIDTH-1:0] data_b;

    assign sign_ab = {operand_a[WIDTH], operand_b[WIDTH]};
    assign data_a = operand_a[WIDTH-1:0];
    assign data_b = operand_b[WIDTH-1:0];

    always_comb begin
        case(sign_ab)
            2'b10: result_o = operand_b;
            2'b01: result_o = operand_a;
            default: result_o = (data_a > data_b) ? operand_a : operand_b;
        endcase
    end
endmodule

// Signed maximum tree over N elements, each (WIDTH+1) bits.
module comp_tree #(
    parameter int unsigned N = 4,
    parameter int unsigned WIDTH = 8
)(
    input  logic signed [N*(WIDTH+1)-1:0] operands_i,
    output logic signed [WIDTH:0]          result_o
);
    logic signed [WIDTH:0] tmp;
    always_comb begin
        tmp = operands_i[0 +: (WIDTH+1)];
        for (int unsigned k = 1; k < N; k++) begin
            if (operands_i[k*(WIDTH+1) +: (WIDTH+1)] > tmp)
                tmp = operands_i[k*(WIDTH+1) +: (WIDTH+1)];
        end
    end
    assign result_o = tmp;
endmodule

// Balanced modular sum tree: result is the low WIDTH bits of the sum of N inputs.
module sum_tree #(
    parameter int unsigned N = 17,
    parameter int unsigned WIDTH = 20
)(
    input  logic [N*WIDTH-1:0] operands_i,
    output logic [WIDTH-1:0]   sum_o
);
    localparam int unsigned N_A = N / 2;
    localparam int unsigned N_B = N - N_A;

    generate
        if (N == 1) begin : g_leaf
            assign sum_o = operands_i[0 +: WIDTH];
        end else if (N == 2) begin : g_pair
            assign sum_o = operands_i[0 +: WIDTH] + operands_i[WIDTH +: WIDTH];
        end else begin : g_split
            logic [N_A*WIDTH-1:0] operands_a;
            logic [N_B*WIDTH-1:0] operands_b;
            logic [WIDTH-1:0] sum_a;
            logic [WIDTH-1:0] sum_b;

            assign operands_a = operands_i[0 +: N_A*WIDTH];
            assign operands_b = operands_i[N_A*WIDTH +: N_B*WIDTH];

            sum_tree #(.N(N_A), .WIDTH(WIDTH)) ua_sum(
                .operands_i(operands_a), .sum_o(sum_a)
            );
            sum_tree #(.N(N_B), .WIDTH(WIDTH)) ub_sum(
                .operands_i(operands_b), .sum_o(sum_b)
            );
            assign sum_o = sum_a + sum_b;
        end
    endgenerate
endmodule

// Behavioral replacement for radix-4 booth multiplier.  We only use it for 4x4->8.
module radix4_booth_multiplier #(
    parameter int unsigned WIDTH_A = 16,
    parameter int unsigned WIDTH_B = 16,
    parameter int unsigned WIDTH_O = WIDTH_A + WIDTH_B
)(
    input  logic [WIDTH_A-1:0] operand_a,
    input  logic [WIDTH_B-1:0] operand_b,
    output logic [WIDTH_O-1:0] sum_o,
    output logic [WIDTH_O-1:0] carry_o
);
    logic [WIDTH_O-1:0] op_a, op_b;
    assign op_a = operand_a;
    assign op_b = operand_b;
    assign sum_o   = op_a * op_b;
    assign carry_o = '0;
endmodule

// Posit decoder (original algorithm, using behavioral lzc/barrel).
module posit_decoder #(
    parameter int unsigned n = 16,
    parameter int unsigned es = 1,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd + es,
    parameter int unsigned MANT_WIDTH = n - es - 3
)(
    input  logic [n-1:0]        operand_i,
    output logic                sign_o,
    output logic signed [EXP_WIDTH:0] rg_exp_o,
    output logic [MANT_WIDTH:0] mant_norm_o
);
    logic sign;
    logic [n-2:0] operand_value;
    assign sign = operand_i[n-1];
    assign operand_value = sign ? (~operand_i[n-2:0] + 1'b1) : operand_i[n-2:0];

    logic regS;
    logic [n-2:0] lzc_operand;
    logic [nd-1:0] leading_zero_count;
    logic lzc_zeroes;

    assign regS = operand_value[n-2];
    assign lzc_operand = regS ? (~operand_value) : operand_value;

    lzc #(
        .WIDTH(n-1),
        .MODE(1'b1)
    ) u_lzc(
        .in_i(lzc_operand),
        .cnt_o(leading_zero_count),
        .empty_o(lzc_zeroes)
    );

    logic [nd-1:0] runlength;
    logic [nd-1:0] regime_bits;
    logic signed [nd:0] regime_k;

    assign runlength   = lzc_zeroes ? (n-1) : leading_zero_count;
    assign regime_bits = lzc_zeroes ? (n-1) : (leading_zero_count + 1'd1);
    assign regime_k    = regS ? ({1'b0, runlength - 1'd1}) : ({1'b1, ~runlength + 1'd1});

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

    logic [es-1:0] exp_val;
    logic [es:0] exp;
    assign sign_o = sign;

    if (es == 0) begin : g_no_exp
        assign exp = '0;
    end else begin : g_exp
        assign exp = op_no_rg[n-2 -: es];
    end

    assign rg_exp_o = regime_k << es | exp;

    logic implicit_bit;
    assign implicit_bit = |operand_i[n-2:0];
    assign mant_norm_o = {implicit_bit, op_no_rg[n-2-es:2]};
endmodule

// Posit encoder (original algorithm, using behavioral barrel shifter).
module posit_encoder #(
    parameter int unsigned n = 16,
    parameter int unsigned es = 1,
    parameter int unsigned nd = pdpu_pkg::clog2(n-1),
    parameter int unsigned EXP_WIDTH = nd + es,
    parameter int unsigned MANT_WIDTH = n - es - 3
)(
    input  logic sign_i,
    input  logic signed [EXP_WIDTH:0] rg_exp_i,
    input  logic [MANT_WIDTH:0]       mant_norm_i,
    output logic [n-1:0]              result_o
);
    logic input_not_zero;
    assign input_not_zero = mant_norm_i[MANT_WIDTH];

    logic signed [EXP_WIDTH-es:0] regime_k;
    logic signed [es:0] exp;

    assign regime_k = rg_exp_i[EXP_WIDTH:es];
    if (es == 0) begin : g_no_exp
        assign exp = '0;
    end else begin : g_exp
        assign exp = rg_exp_i[es-1:0];
    end

    logic sign_k;
    logic [n-2:0] rg_const;
    logic [n-2:0] regime;

    assign sign_k = rg_exp_i[EXP_WIDTH];
    assign rg_const = '1;
    assign regime = sign_k ? rg_const : ~rg_const;

    logic [EXP_WIDTH-es:0] regime_bits;
    assign regime_bits = sign_k ? (~regime_k + 2'sd2) : (regime_k + 2'sd2);

    localparam int unsigned TEMP_WIDTH = (n-1) + es + MANT_WIDTH;
    logic [TEMP_WIDTH-1:0] rg_exp_mant;

    if (es == 0) begin : g_rgm0
        assign rg_exp_mant = {regime, mant_norm_i[MANT_WIDTH-1:0]};
    end else begin : g_rgm
        assign rg_exp_mant = {regime, exp[es-1:0], mant_norm_i[MANT_WIDTH-1:0]};
    end

    localparam int unsigned MAX_SHIFT_AMOUNT = MANT_WIDTH + es + 1;
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(MAX_SHIFT_AMOUNT+1);
    logic [SHIFT_WIDTH-1:0] shift_amount;

    assign shift_amount = (regime_bits >= n) ? MAX_SHIFT_AMOUNT
                                             : (regime_bits + (MANT_WIDTH + es - n + 1));

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

    logic [n-2:0] value_before_round;
    logic [MAX_SHIFT_AMOUNT-1:0] rounding_bits;
    assign {value_before_round, rounding_bits} = value_after_shift[MAX_SHIFT_AMOUNT+n-2:0];

    logic round_bit;
    logic sticky_bit;
    logic round_value;
    logic [n-2:0] value_after_round;

    assign round_bit = rounding_bits[MAX_SHIFT_AMOUNT-1];
    assign sticky_bit = |rounding_bits[MAX_SHIFT_AMOUNT-2:0];
    assign round_value = round_bit & (sticky_bit | value_before_round[0]);
    assign value_after_round = value_before_round + round_value;

    logic [n-1:0] normal_result;
    assign normal_result = sign_i ? {1'b1, ~value_after_round + 1'b1} : {1'b0, value_after_round};
    assign result_o = input_not_zero ? normal_result : '0;
endmodule

// Mantissa normalization (original algorithm, using behavioral lzc/barrel).
module mantissa_norm #(
    parameter int unsigned WIDTH = 8,
    parameter int unsigned EXP_WIDTH = 3,
    parameter int unsigned DECIMAL_POINT = 3
)(
    input  logic [WIDTH-1:0]          operand_i,
    output logic signed [EXP_WIDTH:0] exp_adjust,
    output logic [WIDTH-1:0]          result_o
);
    localparam int unsigned LZC_WIDTH = pdpu_pkg::clog2(WIDTH);
    logic [LZC_WIDTH-1:0] leading_zero_count;
    logic lzc_zeroes;

    lzc #(
        .WIDTH(WIDTH),
        .MODE(1'b1)
    ) u_lzc(
        .in_i(operand_i),
        .cnt_o(leading_zero_count),
        .empty_o(lzc_zeroes)
    );

    // Same arithmetic as the synthesizable case, but widened so -4..+3 are
    // representable in an unsigned expression (widths selected here keep the
    // behavior identical to the original).
    logic signed [EXP_WIDTH+2:0] exp_adjust_w;
    always_comb begin
        if (lzc_zeroes) begin
            exp_adjust_w = '0;
        end else if (leading_zero_count <= DECIMAL_POINT - 1) begin
            exp_adjust_w = $signed(DECIMAL_POINT) - $signed(leading_zero_count) - 1;
        end else begin
            exp_adjust_w = $signed(DECIMAL_POINT) - 1 - $signed(leading_zero_count);
        end
    end
    assign exp_adjust = exp_adjust_w;

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

// Exact CSA-tree modules from the original (kept for bit-exact accumulation).
module fulladder(
    input logic x,y,z,
    output logic sum,carry
);
    assign sum = (x ^ y ) ^ z;
    assign carry = ((x ^ y) & z) | (x & y);
endmodule

module counter_5to3(
    input logic x1,x2,x3,x4,cin,
    output logic sum,carry,cout
);
    assign sum = x1 ^ x2 ^ x3 ^ x4 ^ cin;
    assign cout = (x1 ^ x2) & x3 | ~(x1 ^ x2) & x1;
    assign carry = (x1 ^ x2 ^ x3 ^ x4) & cin | ~(x1 ^ x2 ^ x3 ^ x4) & x4;
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
                .carry(carry[i])
            );
        end
    endgenerate
    assign sum_o = sum;
    assign carry_o = {1'b0,carry,1'b0};
endmodule

module compressor_4to2 #(
    parameter int unsigned WIDTH_I = 8,
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(4)
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

module csa_tree #(
    parameter int unsigned N = 8,
    parameter int unsigned WIDTH_I = 8,
    parameter int unsigned WIDTH_O = WIDTH_I + pdpu_pkg::clog2(N)
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

// Top-level: same pipeline/arithmetic as the original, but behavioral leaf cells.
module pdpu_top_pipelined #(
    parameter int unsigned N = 16,
    parameter int unsigned n_i = 8,
    parameter int unsigned es_i = 2,
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14
)(
    input  logic clk,
    input  logic [N*n_i-1:0] operands_a,
    input  logic [N*n_i-1:0] operands_b,
    input  logic [n_o-1:0]   acc,
    output logic [n_o-1:0]   result_o
);
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;
    localparam int unsigned MANT_WIDTH_I = n_i - es_i - 3;

    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar i;
        for (i = 0; i < N; i++) begin : posit_decoding
            posit_decoder #(
                .n(n_i), .es(es_i)
            ) ua_posit_decoder(
                .operand_i(operands_a[i*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[i*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_a[i*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(
                .n(n_i), .es(es_i)
            ) ub_posit_decoder(
                .operand_i(operands_b[i*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[i*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_b[i*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate

    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o - es_o - 3;

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    posit_decoder #(
        .n(n_o), .es(es_o)
    ) uc_posit_decoder(
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    logic [N-1:0] signs_ab;
    logic [N:0] signs;
    assign signs_ab = signs_a ^ signs_b;
    assign signs    = {sign_acc, signs_ab};

    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I + 1, EXP_WIDTH_O);

    logic signed [N*(EXP_WIDTH+1)-1:0] rg_exp_c;
    generate
        genvar m;
        for (m = 0; m < N; m++) begin : exp_add
            assign rg_exp_c[m*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                $signed(rg_exp_a[m*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]) +
                $signed(rg_exp_b[m*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
        end
    endgenerate

    // pipe1
    logic [N:0] pipe1_signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe1_rg_exp_c;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;

    always_ff @(posedge clk) begin
        pipe1_signs        <= signs;
        pipe1_rg_exp_c     <= rg_exp_c;
        pipe1_mants_norm_a <= mants_norm_a;
        pipe1_mants_norm_b <= mants_norm_b;
        pipe1_rg_exp_acc   <= rg_exp_acc;
        pipe1_mant_norm_acc <= mant_norm_acc;
    end

    localparam int unsigned MUL_WIDTH = 2 * (MANT_WIDTH_I + 1);
    logic [N*MUL_WIDTH-1:0] mul_sum, mul_carry;

    generate
        genvar j;
        for (j = 0; j < N; j++) begin : multiplication
            radix4_booth_multiplier #(
                .WIDTH_A(MANT_WIDTH_I+1),
                .WIDTH_B(MANT_WIDTH_I+1)
            ) u_radix4_booth_multiplier(
                .operand_a(pipe1_mants_norm_a[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .operand_b(pipe1_mants_norm_b[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .sum_o(mul_sum[j*MUL_WIDTH +: MUL_WIDTH]),
                .carry_o(mul_carry[j*MUL_WIDTH +: MUL_WIDTH])
            );
        end
    endgenerate

    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_items;
    generate
        genvar u;
        for (u = 0; u < N; u++) begin : exp_items_prod
            assign rg_exp_items[u*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                pipe1_rg_exp_c[u*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate
    assign rg_exp_items[N*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc);

    logic signed [EXP_WIDTH:0] rg_exp_max;
    comp_tree #(
        .N(N+1), .WIDTH(EXP_WIDTH)
    ) u_comp_tree(
        .operands_i(rg_exp_items),
        .result_o(rg_exp_max)
    );

    // pipe2
    logic [N:0] pipe2_signs;
    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] pipe2_rg_exp_items;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_sum, pipe2_mul_carry;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;

    always_ff @(posedge clk) begin
        pipe2_signs         <= pipe1_signs;
        pipe2_rg_exp_items  <= rg_exp_items;
        pipe2_rg_exp_max    <= rg_exp_max;
        pipe2_mul_sum       <= mul_sum;
        pipe2_mul_carry     <= mul_carry;
        pipe2_mant_norm_acc <= pipe1_mant_norm_acc;
    end

    logic [N*MUL_WIDTH-1:0] mants_norm_c;
    generate
        genvar v;
        for (v = 0; v < N; v++) begin : final_mul_add
            assign mants_norm_c[v*MUL_WIDTH +: MUL_WIDTH] =
                pipe2_mul_sum[v*MUL_WIDTH +: MUL_WIDTH] +
                pipe2_mul_carry[v*MUL_WIDTH +: MUL_WIDTH];
        end
    endgenerate

    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;

    generate
        genvar k;
        if (ALIGN_WIDTH > MUL_WIDTH) begin : fixed_shift_prod
            for (k = 0; k < N; k++) begin : lane
                assign product[k*ALIGN_WIDTH +: ALIGN_WIDTH] =
                    mants_norm_c[k*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH - MUL_WIDTH);
            end
        end else begin : fixed_shift_prod_narrow
            for (k = 0; k < N; k++) begin : lane
                assign product[k*ALIGN_WIDTH +: ALIGN_WIDTH] =
                    mants_norm_c[k*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH - ALIGN_WIDTH);
            end
        end
    endgenerate

    if (ALIGN_WIDTH > MANT_WIDTH_O + 2) begin : acc_shift_wide
        assign product[N*ALIGN_WIDTH +: ALIGN_WIDTH] =
            pipe2_mant_norm_acc << (ALIGN_WIDTH - MANT_WIDTH_O - 2);
    end else begin : acc_shift_narrow
        assign product[N*ALIGN_WIDTH +: ALIGN_WIDTH] =
            pipe2_mant_norm_acc >> (MANT_WIDTH_O + 2 - ALIGN_WIDTH);
    end

    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    logic [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;

    generate
        genvar z, s;
        for (z = 0; z < N+1; z++) begin : exp_diff_lane
            assign rg_exp_diff[z*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                $unsigned(pipe2_rg_exp_max - pipe2_rg_exp_items[z*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]);
        end
        if (EXP_WIDTH + 1 > SHIFT_WIDTH) begin : sat_shift_wide
            for (s = 0; s < N+1; s++) begin : shifter
                assign shift_amount[s*SHIFT_WIDTH +: SHIFT_WIDTH] =
                    (|rg_exp_diff[s*(EXP_WIDTH+1) + SHIFT_WIDTH +: (EXP_WIDTH - SHIFT_WIDTH + 1)]) ?
                        ALIGN_WIDTH : rg_exp_diff[s*(EXP_WIDTH+1) +: SHIFT_WIDTH];
                barrel_shifter #(
                    .WIDTH(ALIGN_WIDTH),
                    .SHIFT_WIDTH(SHIFT_WIDTH),
                    .MODE(1'b1)
                ) u_barrel_shifter(
                    .operand_i(product[s*ALIGN_WIDTH +: ALIGN_WIDTH]),
                    .shift_amount(shift_amount[s*SHIFT_WIDTH +: SHIFT_WIDTH]),
                    .result_o(product_shifted[s*ALIGN_WIDTH +: ALIGN_WIDTH])
                );
            end
        end else begin : sat_shift_narrow
            for (s = 0; s < N+1; s++) begin : shifter
                barrel_shifter #(
                    .WIDTH(ALIGN_WIDTH),
                    .SHIFT_WIDTH(EXP_WIDTH+1),
                    .MODE(1'b1)
                ) u_barrel_shifter(
                    .operand_i(product[s*ALIGN_WIDTH +: ALIGN_WIDTH]),
                    .shift_amount(rg_exp_diff[s*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                    .result_o(product_shifted[s*ALIGN_WIDTH +: ALIGN_WIDTH])
                );
            end
        end
    endgenerate

    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;

    logic [(N+1)*(SUM_WIDTH+1)-1:0] mantissa, mantissa_comp;
    generate
        genvar y;
        for (y = 0; y < N+1; y++) begin : twos_comp
            assign mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] =
                product_shifted[y*ALIGN_WIDTH +: ALIGN_WIDTH];
            assign mantissa_comp[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] =
                pipe2_signs[y] ?
                    (~mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] + 1'b1) :
                    mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)];
        end
    endgenerate

    // pipe3
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;
    logic [(N+1)*(SUM_WIDTH+1)-1:0] pipe3_mantissa_comp;
    always_ff @(posedge clk) begin
        pipe3_rg_exp_max    <= pipe2_rg_exp_max;
        pipe3_mantissa_comp <= mantissa_comp;
    end

    logic [SUM_WIDTH:0] csa_sum, csa_carry;
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
    assign sum_result = csa_sum + csa_carry;

    // pipe4
    logic [SUM_WIDTH:0] pipe4_sum_result;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;
    always_ff @(posedge clk) begin
        pipe4_sum_result <= sum_result;
        pipe4_rg_exp_max  <= pipe3_rg_exp_max;
    end

    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    assign final_sign = pipe4_sum_result[SUM_WIDTH];
    assign sum_c = final_sign ? (~pipe4_sum_result + 1'b1) : pipe4_sum_result[SUM_WIDTH-1:0];

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

    assign final_rg_exp = pipe4_rg_exp_max + rg_exp_adjust;

    logic [MANT_WIDTH_O+2:0] final_mant;
    if (SUM_WIDTH > MANT_WIDTH_O + 3) begin : sticky_trunc
        logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
        assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
        assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2], |sticky_bits};
    end else begin : sticky_ext
        assign final_mant = sum_norm << (MANT_WIDTH_O + 3 - SUM_WIDTH);
    end

    // pipe5
    logic pipe5_final_sign;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;
    always_ff @(posedge clk) begin
        pipe5_final_sign    <= final_sign;
        pipe5_final_rg_exp  <= final_rg_exp;
        pipe5_final_mant    <= final_mant;
    end

    posit_encoder #(
        .n(n_o), .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );

    // DEBUG
    logic [127:0] opa_d [0:4];
    logic [127:0] opb_d [0:4];
    logic [15:0] acc_d [0:4];
    integer dbg_mismatch = 0;
    always_ff @(posedge clk) begin
        opa_d[4] <= opa_d[3]; opb_d[4] <= opb_d[3]; acc_d[4] <= acc_d[3];
        opa_d[3] <= opa_d[2]; opb_d[3] <= opb_d[2]; acc_d[3] <= acc_d[2];
        opa_d[2] <= opa_d[1]; opb_d[2] <= opb_d[1]; acc_d[2] <= acc_d[1];
        opa_d[1] <= opa_d[0]; opb_d[1] <= opb_d[0]; acc_d[1] <= acc_d[0];
        opa_d[0] <= operands_a; opb_d[0] <= operands_b; acc_d[0] <= acc;
    end
    always_ff @(posedge clk) begin
        if (tb.nq > 13 && result_o !== tb.q[4] && dbg_mismatch < 12) begin
            dbg_mismatch <= dbg_mismatch + 1;
            $display("MISMATCH exp=%h got=%h A=%h B=%h ACC=%h", tb.q[4], result_o, opa_d[4], opb_d[4], acc_d[4]);
            $display("  dec a0 sign=%b ex=%h m=%h", signs_a[0], rg_exp_a[0 +: 6], mants_norm_a[0 +: 4]);
            $display("  dec b0 sign=%b ex=%h m=%h", signs_b[0], rg_exp_b[0 +: 6], mants_norm_b[0 +: 4]);
            $display("  pipe1 ma0=%h mb0=%h", pipe1_mants_norm_a[0 +: 4], pipe1_mants_norm_b[0 +: 4]);
            $display("  mul_sum0=%h mul_carry0=%h mantsc0=%h", pipe2_mul_sum[0 +: 8], pipe2_mul_carry[0 +: 8], mants_norm_c[0 +: 8]);
            $display("  exp_max=%h expitems0=%h acc_ex=%h", rg_exp_max, rg_exp_items[0 +: 6], rg_exp_items[16*6 +: 6]);
            $display("  expdiff0=%h shamt0=%h prod0=%h psh0=%h", rg_exp_diff[0 +: 6], shift_amount[0 +: 4], product[0 +: 14], product_shifted[0 +: 14]);
            $display("  mant0=%h mcomp0=%h", mantissa[0 +: 20], mantissa_comp[0 +: 20]);
            $display("  csa_sum=%h csa_carry=%h sumres=%h", csa_sum, csa_carry, sum_result);
            $display("  pipe4_sum=%h final_sign=%b sum_c=%h adj=%h final_exp=%h normsum=%h final_mant=%h",
                     pipe4_sum_result, final_sign, sum_c, rg_exp_adjust, final_rg_exp, sum_norm, final_mant);
        end
    end
endmodule
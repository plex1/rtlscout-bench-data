// PDPU dot-16, 5-stage pipelined.
// Functionally identical to the upstream pdpu_top_pipelined baseline
// (starting_point.v): only the 4x4 mantissa multiply is inlined (the booth
// multiplier computed exactly a*b, and its sum/carry pair was re-added in
// stage 2 -- both replaced by a direct product).

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

// ---- basic compressors ----
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

module comparator #(
    parameter int unsigned WIDTH = 8
)(
    input logic signed [WIDTH:0] operand_a,
    input logic signed [WIDTH:0] operand_b,
    output logic signed [WIDTH:0] result_o
);
    // signed maximum (equivalent to the original sign-case + low-bit compare)
    always_comb begin
        if ($signed(operand_a) > $signed(operand_b)) result_o = operand_a;
        else                                         result_o = operand_b;
    end
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

module comp_tree #(
    parameter N = 4,
    parameter WIDTH = 8
)(
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

            assign operands_i_A = operands_i[0       +: N_A*W];
            assign operands_i_B = operands_i[N_A*W   +: N_B*W];

            comp_tree #(
                .N(N_A),
                .WIDTH(WIDTH)
            ) ua_comp_tree(
                .operands_i(operands_i_A),
                .result_o(result_o_A)
            );

            comp_tree #(
                .N(N_B),
                .WIDTH(WIDTH)
            ) ub_comp_tree(
                .operands_i(operands_i_B),
                .result_o(result_o_B)
            );

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

// barrel shifter
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

// leading zero counter
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

endmodule : lzc

// Mantissa normalization
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

    lzc #(
        .WIDTH(WIDTH),
        .MODE(1'b1)
    ) u_lzc(
        .in_i(operand_i),
        .cnt_o(leading_zero_count),
        .empty_o(lzc_zeroes)
    );

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

    assign runlength = lzc_zeroes ? (n-1) : leading_zero_count;
    assign regime_bits = lzc_zeroes ? (n-1) : (leading_zero_count+1);
    assign regime_k = regS ? ({1'b0,runlength-1}) : ({1'b1,~runlength+1});

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

// Posit encoder
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

    // Sign Process
    logic [N-1:0] signs_ab;
    logic [N:0] signs;

    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc, signs_ab};

    // Exponent Addition
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);

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
    // Mantissa Multiplication (4x4 -> 8 bit, inlined)
    // ---------------
    localparam int unsigned MUL_WIDTH = 2*(MANT_WIDTH_I+1);
    logic [N*MUL_WIDTH-1:0] mul_prod;

    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            assign mul_prod[(j)*MUL_WIDTH +: MUL_WIDTH] =
                pipe1_mants_norm_a[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] *
                pipe1_mants_norm_b[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)];
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
    logic [N:0] pipe2_signs;
    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] pipe2_rg_exp_items;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [N*MUL_WIDTH-1:0] pipe2_mul_prod;

    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;


    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_items,rg_exp_items,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe2_mul_prod,mul_prod,clk)

    `FFNR(pipe2_mant_norm_acc,pipe1_mant_norm_acc,clk)

    // ---------------
    // Mantissa Alignment
    // ---------------
    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;

    generate
        genvar k;
        if(ALIGN_WIDTH>MUL_WIDTH) begin: fixed_shift
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
            end
        end
        else begin: fixed_shift
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mul_prod[(k)*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH-ALIGN_WIDTH);
            end
        end
    endgenerate

    if(ALIGN_WIDTH>MANT_WIDTH_O+2) begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
    end
    else begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc >> (MANT_WIDTH_O+2-ALIGN_WIDTH);
    end

    // Align the products
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    logic [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_diff;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z,s;
        for(z=0;z<N+1;z++) begin
            assign rg_exp_diff[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $unsigned(pipe2_rg_exp_max - pipe2_rg_exp_items[(z)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]);
        end
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
    assign sum_result = csa_sum + csa_carry;


    // ********************
    // Pipeline 4
    // ********************
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

endmodule

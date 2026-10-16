// PDPU dot-16, 5-stage pipelined -- v7.
// v6 + stage-4 rework: balanced 3:2-compressor tree (no horizontal carry
// ripple) over 18 operands (17 lanes + popcount correction), then one 19-bit CPA.

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

    function automatic integer unsigned idx_width (input integer unsigned num_idx);
        idx_width = (num_idx > 32'd1) ? $unsigned(clog(num_idx)) : 32'd1;
    endfunction
endpackage

module fulladder(
    input logic x,y,z,
    output logic sum,carry
);
    assign sum = (x ^ y ) ^ z;
    assign carry = ((x ^ y) & z) | (x & y);
endmodule

module comparator #(
    parameter int unsigned WIDTH = 8
)(
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
            ) u_comparator (
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
            ) ua_comp_tree (
                .operands_i(operands_i_A),
                .result_o(result_o_A)
            );

            comp_tree #(
                .N(N_B),
                .WIDTH(WIDTH)
            ) ub_comp_tree (
                .operands_i(operands_i_B),
                .result_o(result_o_B)
            );

            comparator #(
                .WIDTH(WIDTH)
            ) uc_comparator (
                .operand_a(result_o_A),
                .operand_b(result_o_B),
                .result_o(result_o)
            );
        end
    endgenerate
endmodule

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

endmodule : lzc

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
    ) u_lzc (
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
    ) u_barrel_shifter (
        .operand_i(operand_i),
        .shift_amount(leading_zero_count),
        .result_o(result_o)
    );
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

    lzc #(
        .WIDTH(n-1),
        .MODE(1'b1)
    ) u_lzc (
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
    ) u_barrel_shifter (
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
    ) u_barrel_shifter (
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

// ---------------------------------------------------------------------------
// Pipelined top
// ---------------------------------------------------------------------------
module pdpu_top_pipelined #(
    parameter int unsigned N = 16,
    parameter int unsigned n_i = 8,
    parameter int unsigned es_i = 2,
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14,
    parameter int unsigned CARRY_W = 5,
    parameter int unsigned SUM_W = 19
)(
    input logic clk,
    input logic [N*n_i-1:0] operands_a,
    input logic [N*n_i-1:0] operands_b,
    input logic [n_o-1:0] acc,
    output logic [n_o-1:0] result_o
);
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;   // 5
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;                     // 3
    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;   // 6
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;                     // 11
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O); // 6
    localparam int unsigned MULW = 2*(MANT_WIDTH_I+1);                     // 8
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_W;             // 19
    localparam int unsigned SHIFT_WIDTH = pdpu_pkg::clog2(ALIGN_WIDTH+1);  // 4
    localparam int unsigned LANE_W = SUM_WIDTH;                            // 19
    localparam int unsigned L1_CNT = 18;                                   // 17 lanes + corr
    localparam int unsigned L2_CNT = 12;                                   // 6 groups of 3
    localparam int unsigned L3_CNT = 8;                                    // 4 groups of 3
    localparam int unsigned L4_CNT = 6;                                    // 2 groups of 3 + 2 pass
    localparam int unsigned L5_CNT = 4;                                    // 2 groups of 3
    localparam int unsigned L6_CNT = 3;                                    // 1 group of 3 + 1 pass

    // ================= Stage 1: decode, exp adds, pairwise max, popcount =================
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_a, mants_b;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(.n(n_i),.es(es_i)) ua_posit_decoder (
                .operand_i(operands_a[(i)*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_a[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(.n(n_i),.es(es_i)) ub_posit_decoder (
                .operand_i(operands_b[(i)*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[(i)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_b[(i)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    posit_decoder #(.n(n_o),.es(es_o)) uc_posit_decoder (
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    logic [N-1:0] signs_ab;
    logic [N:0] s1_signs;
    assign signs_ab = signs_a ^ signs_b;
    assign s1_signs = {sign_acc, signs_ab};

    // exponent sums
    logic signed [N*(EXP_WIDTH+1)-1:0] s1_rg_exp_c;
    generate
        genvar m;
        for(m=0;m<N;m++) begin
            assign s1_rg_exp_c[(m)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                $signed(rg_exp_a[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]) +
                $signed(rg_exp_b[(m)*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
        end
    endgenerate

    // pairwise maxima (tree level 1): 8 pairs of exponent sums
    logic signed [(N/2)*(EXP_WIDTH+1)-1:0] s1_exp_pairs;
    generate
        genvar h;
        for(h=0;h<N/2;h++) begin
            comparator #(.WIDTH(EXP_WIDTH)) u_pair_max (
                .operand_a(s1_rg_exp_c[(2*h)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .operand_b(s1_rg_exp_c[(2*h+1)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .result_o(s1_exp_pairs[h*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
            );
        end
    endgenerate

    // balanced popcount of signs (two's complement correction term)
    logic [15:0] pc1;
    logic [11:0] pc2;
    logic [7:0]  pc3;
    logic [4:0] s1_corr;
    generate
        genvar p;
        for(p=0;p<8;p++) begin
            assign pc1[2*p +: 2] = {1'b0, s1_signs[2*p]} + {1'b0, s1_signs[2*p+1]};
        end
        for(p=0;p<4;p++) begin
            assign pc2[3*p +: 3] = {1'b0, pc1[4*p +: 2]} + {1'b0, pc1[4*p+2 +: 2]};
        end
        for(p=0;p<2;p++) begin
            assign pc3[4*p +: 4] = {1'b0, pc2[6*p +: 3]} + {1'b0, pc2[6*p+3 +: 3]};
        end
    endgenerate
    assign s1_corr = {1'b0, pc3[3:0]} + {1'b0, pc3[7:4]} + {4'b0, s1_signs[N]};

    // ---- pipe1 ----
    logic [N:0] pipe1_signs;
    logic signed [(N/2)*(EXP_WIDTH+1)-1:0] pipe1_exp_pairs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe1_rg_exp_c;
    logic signed [EXP_WIDTH_O:0] pipe1_rg_exp_acc;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;
    logic [4:0] pipe1_corr;

    always_ff @(posedge clk) begin
        pipe1_signs         <= s1_signs;
        pipe1_exp_pairs     <= s1_exp_pairs;
        pipe1_rg_exp_c      <= s1_rg_exp_c;
        pipe1_rg_exp_acc    <= rg_exp_acc;
        pipe1_mants_norm_a  <= mants_a;
        pipe1_mants_norm_b  <= mants_b;
        pipe1_mant_norm_acc <= mant_norm_acc;
        pipe1_corr          <= s1_corr;
    end

    // ============ Stage 2: max-tree tail (9 exps) ============
    logic signed [(N/2+1)*(EXP_WIDTH+1)-1:0] exp_items2;
    assign exp_items2[(N/2)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc);
    generate
        genvar t;
        for(t=0;t<N/2;t++) begin
            assign exp_items2[t*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                pipe1_exp_pairs[t*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate

    logic signed [EXP_WIDTH:0] rg_exp_max;
    comp_tree #(
        .N(N/2+1),
        .WIDTH(EXP_WIDTH)
    ) u_comp_tree (
        .operands_i(exp_items2),
        .result_o(rg_exp_max)
    );

    // ---- pipe2 ----
    logic [N:0] pipe2_signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] pipe2_rg_exp_c;
    logic signed [EXP_WIDTH_O:0] pipe2_rg_exp_acc;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe2_mants_norm_a, pipe2_mants_norm_b;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;
    logic [4:0] pipe2_corr;

    always_ff @(posedge clk) begin
        pipe2_signs         <= pipe1_signs;
        pipe2_rg_exp_c      <= pipe1_rg_exp_c;
        pipe2_rg_exp_acc    <= pipe1_rg_exp_acc;
        pipe2_rg_exp_max    <= rg_exp_max;
        pipe2_mants_norm_a  <= pipe1_mants_norm_a;
        pipe2_mants_norm_b  <= pipe1_mants_norm_b;
        pipe2_mant_norm_acc <= pipe1_mant_norm_acc;
        pipe2_corr          <= pipe1_corr;
    end

    // ============ Stage 3: mults + diff + align + sign-XOR ============
    logic [N*MULW-1:0] s3_prod;
    generate
        genvar j;
        for(j=0;j<N;j++) begin
            assign s3_prod[(j)*MULW +: MULW] =
                pipe2_mants_norm_a[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] *
                pipe2_mants_norm_b[(j)*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)];
        end
    endgenerate

    logic [(N+1)*ALIGN_WIDTH-1:0] s3_shifted_xor;
    generate
        genvar s;
        for(s=0;s<N+1;s++) begin : align
            logic [ALIGN_WIDTH-1:0] product;
            logic [EXP_WIDTH:0] rg_exp_diff;
            logic [SHIFT_WIDTH-1:0] shift_amount;
            logic [ALIGN_WIDTH-1:0] shifted;

            if (s < N) begin
                assign product = s3_prod[(s)*MULW +: MULW] << (ALIGN_WIDTH-MULW);
            end
            else begin
                assign product = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
            end

            if (s < N) begin
                assign rg_exp_diff = $unsigned(pipe2_rg_exp_max - pipe2_rg_exp_c[(s)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]);
            end
            else begin
                assign rg_exp_diff = $unsigned(pipe2_rg_exp_max - $signed(pipe2_rg_exp_acc));
            end
            if(EXP_WIDTH+1>SHIFT_WIDTH) begin
                assign shift_amount = (|rg_exp_diff[(EXP_WIDTH+1)-1:SHIFT_WIDTH]) ? ALIGN_WIDTH[SHIFT_WIDTH-1:0]
                                                                                    : rg_exp_diff[SHIFT_WIDTH-1:0];
            end
            else begin
                assign shift_amount = rg_exp_diff[SHIFT_WIDTH-1:0];
            end
            assign shifted = product >> shift_amount;
            assign s3_shifted_xor[(s)*ALIGN_WIDTH +: ALIGN_WIDTH] = shifted ^ {ALIGN_WIDTH{pipe2_signs[s]}};
        end
    endgenerate

    // ---- pipe3 ----
    logic [N:0] pipe3_signs;
    logic [(N+1)*ALIGN_WIDTH-1:0] pipe3_shifted_xor;
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;
    logic [4:0] pipe3_corr;

    always_ff @(posedge clk) begin
        pipe3_signs         <= pipe2_signs;
        pipe3_shifted_xor   <= s3_shifted_xor;
        pipe3_rg_exp_max    <= pipe2_rg_exp_max;
        pipe3_corr          <= pipe2_corr;
    end

    // ============ Stage 4: explicit balanced 3:2 tree over 18 lanes ============
    // L1: 18 lanes (17 sign-extended aligned values + zero-extended corr)
    logic [L1_CNT*LANE_W-1:0] l1_ops;
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin
            assign l1_ops[(y)*LANE_W +: LANE_W] =
                { {(LANE_W-ALIGN_WIDTH){pipe3_signs[y]}}, pipe3_shifted_xor[(y)*ALIGN_WIDTH +: ALIGN_WIDTH] };
        end
    endgenerate
    assign l1_ops[(N+1)*LANE_W +: LANE_W] = {{(LANE_W-5){1'b0}}, pipe3_corr};

    // L1 -> L2 : 6 groups of 3
    logic [L2_CNT*LANE_W-1:0] l2_ops;
    // L2 -> L3 : 4 groups of 3
    logic [L3_CNT*LANE_W-1:0] l3_ops;
    // L3 -> L4 : 2 groups of 3 + 2 passthrough
    logic [L4_CNT*LANE_W-1:0] l4_ops;
    // L4 -> L5 : 2 groups of 3
    logic [L5_CNT*LANE_W-1:0] l5_ops;
    // L5 -> L6 : 1 group of 3 + 1 passthrough
    logic [L6_CNT*LANE_W-1:0] l6_ops;

    generate
        genvar g1, g2, g3, g4, g5, b1, b2, b3, b4, b5;
        for (g1 = 0; g1 < 6; g1++) begin : gen_l2
            logic [LANE_W-1:0] st, ct;
            for (b1 = 0; b1 < LANE_W; b1++) begin
                fulladder u_fa (
                    .x(l1_ops[(3*g1+0)*LANE_W + b1]),
                    .y(l1_ops[(3*g1+1)*LANE_W + b1]),
                    .z(l1_ops[(3*g1+2)*LANE_W + b1]),
                    .sum(st[b1]),
                    .carry(ct[b1])
                );
            end
            assign l2_ops[(2*g1+0)*LANE_W +: LANE_W] = st;
            assign l2_ops[(2*g1+1)*LANE_W +: LANE_W] = {ct[LANE_W-2:0], 1'b0};
        end
        for (g2 = 0; g2 < 4; g2++) begin : gen_l3
            logic [LANE_W-1:0] st, ct;
            for (b2 = 0; b2 < LANE_W; b2++) begin
                fulladder u_fa (
                    .x(l2_ops[(3*g2+0)*LANE_W + b2]),
                    .y(l2_ops[(3*g2+1)*LANE_W + b2]),
                    .z(l2_ops[(3*g2+2)*LANE_W + b2]),
                    .sum(st[b2]),
                    .carry(ct[b2])
                );
            end
            assign l3_ops[(2*g2+0)*LANE_W +: LANE_W] = st;
            assign l3_ops[(2*g2+1)*LANE_W +: LANE_W] = {ct[LANE_W-2:0], 1'b0};
        end
        for (g3 = 0; g3 < 2; g3++) begin : gen_l4
            logic [LANE_W-1:0] st, ct;
            for (b3 = 0; b3 < LANE_W; b3++) begin
                fulladder u_fa (
                    .x(l3_ops[(3*g3+0)*LANE_W + b3]),
                    .y(l3_ops[(3*g3+1)*LANE_W + b3]),
                    .z(l3_ops[(3*g3+2)*LANE_W + b3]),
                    .sum(st[b3]),
                    .carry(ct[b3])
                );
            end
            assign l4_ops[(2*g3+0)*LANE_W +: LANE_W] = st;
            assign l4_ops[(2*g3+1)*LANE_W +: LANE_W] = {ct[LANE_W-2:0], 1'b0};
        end
        // passthrough l3[6], l3[7]
        assign l4_ops[4*LANE_W +: LANE_W] = l3_ops[6*LANE_W +: LANE_W];
        assign l4_ops[5*LANE_W +: LANE_W] = l3_ops[7*LANE_W +: LANE_W];

        for (g4 = 0; g4 < 2; g4++) begin : gen_l5
            logic [LANE_W-1:0] st, ct;
            for (b4 = 0; b4 < LANE_W; b4++) begin
                fulladder u_fa (
                    .x(l4_ops[(3*g4+0)*LANE_W + b4]),
                    .y(l4_ops[(3*g4+1)*LANE_W + b4]),
                    .z(l4_ops[(3*g4+2)*LANE_W + b4]),
                    .sum(st[b4]),
                    .carry(ct[b4])
                );
            end
            assign l5_ops[(2*g4+0)*LANE_W +: LANE_W] = st;
            assign l5_ops[(2*g4+1)*LANE_W +: LANE_W] = {ct[LANE_W-2:0], 1'b0};
        end
        for (g5 = 0; g5 < 1; g5++) begin : gen_l6
            logic [LANE_W-1:0] st, ct;
            for (b5 = 0; b5 < LANE_W; b5++) begin
                fulladder u_fa (
                    .x(l5_ops[(3*g5+0)*LANE_W + b5]),
                    .y(l5_ops[(3*g5+1)*LANE_W + b5]),
                    .z(l5_ops[(3*g5+2)*LANE_W + b5]),
                    .sum(st[b5]),
                    .carry(ct[b5])
                );
            end
            assign l6_ops[(2*g5+0)*LANE_W +: LANE_W] = st;
            assign l6_ops[(2*g5+1)*LANE_W +: LANE_W] = {ct[LANE_W-2:0], 1'b0};
        end
        // passthrough l5[3]
        assign l6_ops[2*LANE_W +: LANE_W] = l5_ops[3*LANE_W +: LANE_W];
    endgenerate

    // final: 3 operands -> 3:2 compress then CPA
    logic [LANE_W-1:0] fs, fc, sum_result;
    generate
        genvar bf;
        for (bf = 0; bf < LANE_W; bf++) begin : gen_final_fa
            fulladder u_fa (
                .x(l6_ops[(0)*LANE_W + bf]),
                .y(l6_ops[(1)*LANE_W + bf]),
                .z(l6_ops[(2)*LANE_W + bf]),
                .sum(fs[bf]),
                .carry(fc[bf])
            );
        end
    endgenerate
    assign sum_result = fs + {fc[LANE_W-2:0], 1'b0};

    // ---- pipe4 ----
    logic [LANE_W-1:0] pipe4_sum_result;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;

    always_ff @(posedge clk) begin
        pipe4_sum_result <= sum_result;
        pipe4_rg_exp_max <= pipe3_rg_exp_max;
    end

    // ============ Stage 5: normalize ============
    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;

    assign final_sign = pipe4_sum_result[LANE_W-1];
    assign sum_c = final_sign ? (~pipe4_sum_result+1) : pipe4_sum_result;

    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic signed [EXP_WIDTH:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;

    mantissa_norm #(
        .WIDTH(SUM_WIDTH),
        .EXP_WIDTH(EXP_WIDTH),
        .DECIMAL_POINT(CARRY_W+2)
    ) u_mantissa_norm (
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );

    assign final_rg_exp = pipe4_rg_exp_max + rg_exp_adjust;

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

    // ---- pipe5 ----
    logic pipe5_final_sign;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;

    always_ff @(posedge clk) begin
        pipe5_final_sign   <= final_sign;
        pipe5_final_rg_exp <= final_rg_exp;
        pipe5_final_mant   <= final_mant;
    end

    // ============ Encoder (combinational, after pipe5) ============
    posit_encoder #(
        .n(n_o),
        .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder (
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );

endmodule

// PDPU dot-16, 5-stage pipelined -- behavioral-arithmetic rewrite.
// Bit-exact to the reference combinational pdpu_dot16.

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

    function automatic integer idx_width (input integer num_idx);
            idx_width = (num_idx > 32'd1) ? $unsigned(clog(num_idx)) : 32'd1;
    endfunction
endpackage

// ---------------------------------------------------------------------------
// Simple priority-encoder leading/trailing zero counter.
// MODE=1 -> leading zeros, MODE=0 -> trailing zeros.
// For an all-zero input the count value is irrelevant to every use in the
// datapath (empty_o gates it), so any value is fine.
// ---------------------------------------------------------------------------
module lzc #(
    parameter int unsigned WIDTH = 2,
    parameter bit          MODE  = 1'b0,
    parameter int unsigned CNT_WIDTH = cf_math_pkg::idx_width(WIDTH)
) (
    input  logic [WIDTH-1:0]     in_i,
    output logic [CNT_WIDTH-1:0] cnt_o,
    output logic                 empty_o
);
    integer i;
    always_comb begin
        cnt_o   = '0;
        empty_o = 1'b1;
        for (i = 0; i < WIDTH; i = i + 1) begin
            if (empty_o && (MODE ? in_i[WIDTH-1-i] : in_i[i])) begin
                cnt_o   = i[CNT_WIDTH-1:0];
                empty_o = 1'b0;
            end
        end
    end
endmodule

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

    assign result_o = operand_i << leading_zero_count;
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

    assign op_no_rg = operand_value << regime_bits;

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
    assign value_after_shift = value_before_shift >> shift_amount;

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

// Compare to obtain the larger of two signed numbers
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

// Recursive comparator tree: maximum of N signed (WIDTH+1)-bit values.
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

// Recursive balanced adder tree. `operands_i` is flat: element k occupies
// bits [k*W +: W] and is interpreted as a signed number. Output is W bits.
module add_tree #(
    parameter int unsigned N = 2,
    parameter int unsigned W = 8
)(
    input  logic [N*W-1:0] operands_i,
    output logic signed [W-1:0] result_o
);
    localparam int unsigned N_A = N/2;
    localparam int unsigned N_B = N - N_A;

    generate
        if (N==1) begin : g_leaf
            assign result_o = operands_i[0 +: W];
        end
        else if (N==2) begin : g_pair
            assign result_o = $signed(operands_i[0 +: W]) + $signed(operands_i[W +: W]);
        end
        else begin : g_split
            logic signed [W-1:0] sum_A;
            logic signed [W-1:0] sum_B;
            add_tree #(.N(N_A), .W(W)) ua_add_tree(
                .operands_i(operands_i[0 +: N_A*W]), .result_o(sum_A));
            add_tree #(.N(N_B), .W(W)) ub_add_tree(
                .operands_i(operands_i[N_A*W +: N_B*W]), .result_o(sum_B));
            assign result_o = sum_A + sum_B;
        end
    endgenerate
endmodule

// Wallace carry-save reduction of N two's-complement terms of W bits into a
// redundant (sum, carry) pair. All arithmetic is done modulo 2^W, which is
// exact as long as the final value fits in W bits.
module csa_reduce #(
    parameter int unsigned N = 2,
    parameter int unsigned W = 19
)(
    input  logic [N*W-1:0] ops,
    output logic [W-1:0]   sum_o,
    output logic [W-1:0]   carry_o
);
    localparam int unsigned NC = N/3;
    localparam int unsigned NR = N - 3*NC;
    localparam int unsigned NN = 2*NC + NR;

    generate
        if (N <= 2) begin : g_leaf
            assign sum_o   = ops[0 +: W];
            assign carry_o = (N==2) ? ops[W +: W] : '0;
        end
        else begin : g_reduce
            logic [NN*W-1:0] next_ops;
            for (genvar k = 0; k < NC; k++) begin : g_comp
                for (genvar b = 0; b < W; b++) begin : g_bit
                    wire a = ops[(3*k+0)*W + b];
                    wire c = ops[(3*k+1)*W + b];
                    wire d = ops[(3*k+2)*W + b];
                    assign next_ops[(2*k)*W + b] = a ^ c ^ d;
                    if (b > 0)
                        assign next_ops[(2*k+1)*W + b] = (a&c) | (a&d) | (c&d);
                    else
                        assign next_ops[(2*k+1)*W + 0] = 1'b0;
                end
            end
            for (genvar r = 0; r < NR; r++) begin : g_pass
                assign next_ops[(2*NC+r)*W +: W] = ops[(3*NC+r)*W +: W];
            end
            csa_reduce #(.N(NN), .W(W)) u_next(
                .ops(next_ops), .sum_o(sum_o), .carry_o(carry_o));
        end
    endgenerate
endmodule

/*
PDPU equipped with a fine-grained 6-stage
*/

module pdpu_top_pipelined #(
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
    // ---------------
    // Posit Decoding
    // ---------------
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

    // ---------------
    // Sign Process
    // ---------------
    logic [N-1:0] signs_ab;
    logic [N:0] signs;

    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc, signs_ab};

    // ---------------
    // Exponent Addition
    // ---------------
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
    // Mantissa Multiplication (4x4 -> 8 bits per lane)
    // ---------------
    localparam int unsigned MUL_WIDTH = 2*(MANT_WIDTH_I+1);
    logic [N*MUL_WIDTH-1:0] mants_norm_c;

    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            assign mants_norm_c[(j)*MUL_WIDTH +: MUL_WIDTH] =
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
    logic [N*MUL_WIDTH-1:0] pipe2_mants_norm_c;
    logic [MANT_WIDTH_O:0] pipe2_mant_norm_acc;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_items,rg_exp_items,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe2_mants_norm_c,mants_norm_c,clk)
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
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mants_norm_c[(k)*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
            end
        end
        else begin: fixed_shift
            for(k=0;k<N;k++) begin
                assign product[(k)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mants_norm_c[(k)*MUL_WIDTH +: MUL_WIDTH] >> (MUL_WIDTH-ALIGN_WIDTH);
            end
        end
    endgenerate

    if(ALIGN_WIDTH>MANT_WIDTH_O+2) begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);
    end
    else begin
        assign product[(N)*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mant_norm_acc >> (MANT_WIDTH_O+2-ALIGN_WIDTH);
    end

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
                assign product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH] = product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH] >> shift_amount[(s)*SHIFT_WIDTH +: SHIFT_WIDTH];
            end
        end
        else begin
            for(s=0;s<N+1;s++) begin
                assign product_shifted[(s)*ALIGN_WIDTH +: ALIGN_WIDTH] = product[(s)*ALIGN_WIDTH +: ALIGN_WIDTH] >> rg_exp_diff[(s)*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
            end
        end
    endgenerate

    // ---------------
    // Mantissa accumulation in two's complement.
    // Negation uses (-x) = (~x) + 1, so each lane contributes
    // (mag ^ sign) plus a single +1. All the +1 terms are collected into one
    // constant popcount of the signs. The XOR (cheap) is done alongside the
    // final alignment in S2 so that S3 is purely a carry-save reduction plus
    // one final add.
    // ---------------
    localparam int unsigned CARRY_WIDTH = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH = ALIGN_WIDTH + CARRY_WIDTH;
    localparam int unsigned PC_WIDTH  = CARRY_WIDTH + 1;

    logic [(N+1)*SUM_WIDTH-1:0] mantissa;
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin
            assign mantissa[(y)*SUM_WIDTH +: SUM_WIDTH] =
                product_shifted[(y)*ALIGN_WIDTH +: ALIGN_WIDTH] ^ {SUM_WIDTH{pipe2_signs[y]}};
        end
    endgenerate

    // popcount of the signs (parallel with the alignment shifter)
    logic [(N+1)*PC_WIDTH-1:0] pc_operands;
    logic [PC_WIDTH-1:0] pc_val;
    generate
        genvar p;
        for(p=0;p<N+1;p++) begin
            assign pc_operands[p*PC_WIDTH +: PC_WIDTH] = {{(PC_WIDTH-1){1'b0}}, pipe2_signs[p]};
        end
    endgenerate
    add_tree #(
        .N(N+1),
        .W(PC_WIDTH)
    ) u_pc_tree(
        .operands_i(pc_operands),
        .result_o(pc_val)
    );

    // *******************
    // Pipeline 3
    // *******************
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;
    logic [(N+1)*SUM_WIDTH-1:0] pipe3_mantissa;
    logic [PC_WIDTH-1:0] pipe3_pc;

    `FFNR(pipe3_rg_exp_max,pipe2_rg_exp_max,clk)
    `FFNR(pipe3_mantissa,mantissa,clk)
    `FFNR(pipe3_pc,pc_val,clk)

    logic [(N+2)*SUM_WIDTH-1:0] tree_operands;
    logic [SUM_WIDTH-1:0] csa_sum, csa_carry;
    logic [SUM_WIDTH-1:0] sum_result;

    assign tree_operands[(N+1)*SUM_WIDTH-1:0] = pipe3_mantissa;
    assign tree_operands[(N+1)*SUM_WIDTH +: SUM_WIDTH] = {{(SUM_WIDTH-PC_WIDTH){1'b0}}, pipe3_pc};

    csa_reduce #(
        .N(N+2),
        .W(SUM_WIDTH)
    ) u_csa_reduce(
        .ops(tree_operands),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );

    assign sum_result = chk_sum;

    logic [SUM_WIDTH-1:0] chk_sum;
    add_tree #(
        .N(N+2),
        .W(SUM_WIDTH)
    ) u_chk_tree(
        .operands_i(tree_operands),
        .result_o(chk_sum)
    );

    // ********************
    // Pipeline 4
    // ********************
    logic [SUM_WIDTH-1:0] pipe4_sum_result;
    logic signed [EXP_WIDTH:0] pipe4_rg_exp_max;

    `FFNR(pipe4_sum_result,sum_result,clk)
    `FFNR(pipe4_rg_exp_max,pipe3_rg_exp_max,clk)

    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;

    assign final_sign = pipe4_sum_result[SUM_WIDTH-1];
    assign sum_c = final_sign ? (~pipe4_sum_result+1) : pipe4_sum_result;

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

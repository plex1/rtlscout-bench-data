// pdpu_top_pipelined -- 5-stage pipelined 16-lane posit dot-product-accumulate.
// Bit-exact to the combinational reference pdpu_dot16 (ALIGN_WIDTH=14 truncation,
// NaR handling identical). Balanced pipeline + parallel-prefix arithmetic:
//   R1 in->pipe1 : decode (33x) + product-exponent add
//   R2 pipe1->2  : 16x (4x4 mul) || max tree (17x7) || popcount || prefix subtract+big flag
//   R3 pipe2->3  : align shift + sign fold + corr row + CSA compress (18 rows -> 2)
//   R4 pipe3->4  : prefix carry-propagate add
//   R5 pipe4->5  : abs + normalize + final exponent + mantissa/sticky
//   R6 pipe5->out: posit encode (rounding add merged with output negation)
// Per-lane two's-complement adders are replaced by XOR-folding plus a single
// popcount correction row folded into the compressor tree (exact in Z/(2^20)).

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

// inclusive OR-prefix scan: y[i] = |x[i:0]
module or_scan #(
    parameter int unsigned M = 7
)(
    input  logic [M-1:0] x,
    output logic [M-1:0] y
);
    localparam int unsigned L = pdpu_pkg::clog2(M);
    logic [(L+1)*M-1:0] t;
    genvar gi, gl;
    generate
        for (gi = 0; gi < M; gi++) begin : g_init
            assign t[gi] = x[gi];
        end
        for (gl = 1; gl <= L; gl++) begin : g_lvl
            localparam int unsigned SP = 1 << (gl-1);
            for (gi = 0; gi < M; gi++) begin : g_bit
                if (gi < SP) begin : g_pass
                    assign t[gl*M + gi] = t[(gl-1)*M + gi];
                end else begin : g_comb
                    assign t[gl*M + gi] = t[(gl-1)*M + gi] | t[(gl-1)*M + gi - SP];
                end
            end
        end
        assign y = t[L*M +: M];
    endgenerate
endmodule

// Kogge-Stone parallel-prefix adder: s = a + b + cin (mod 2^W)
module ks_add #(
    parameter int unsigned W = 8
)(
    input  logic [W-1:0] a,
    input  logic [W-1:0] b,
    input  logic         cin,
    output logic [W-1:0] s
);
    localparam int unsigned L = pdpu_pkg::clog2(W);
    logic [(L+1)*W-1:0] G, P;
    genvar gi, gl;
    generate
        for (gi = 0; gi < W; gi++) begin : g_init
            if (gi == 0) begin : g_cin
                assign G[0] = (a[0] & b[0]) | ((a[0] ^ b[0]) & cin);
            end else begin : g_std
                assign G[gi] = a[gi] & b[gi];
            end
            assign P[gi] = a[gi] ^ b[gi];
        end
        if (W == 1) begin : g_w1
            assign s[0] = a[0] ^ b[0] ^ cin;
        end else begin : g_main
            for (gl = 1; gl <= L; gl++) begin : g_lvl
                localparam int unsigned SP = 1 << (gl-1);
                for (gi = 0; gi < W; gi++) begin : g_bit
                    if (gi < SP) begin : g_pass
                        assign G[gl*W + gi] = G[(gl-1)*W + gi];
                        assign P[gl*W + gi] = P[(gl-1)*W + gi];
                    end else begin : g_comb
                        assign G[gl*W + gi] = G[(gl-1)*W + gi] | (P[(gl-1)*W + gi] & G[(gl-1)*W + gi - SP]);
                        assign P[gl*W + gi] = P[(gl-1)*W + gi] & P[(gl-1)*W + gi - SP];
                    end
                end
            end
            // s[i] = P0[i] ^ carry-into-bit-i; carry into bit 0 is cin
            assign s = P[W-1:0] ^ {G[L*W + W-2 -: W-1], cin};
        end
    endgenerate
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
    logic [n-2:0] z, operand_value;
    assign sign = operand_i[n-1];
    assign z = operand_i[n-2:0];
    // (-z)[i] = z[i] ^ |z[i-1:0]| (i>=1), (-z)[0] = z[0]
    logic [n-2:0] S;
    or_scan #(.M(n-1)) u_scan (.x(z), .y(S));
    assign operand_value = z ^ ({(n-1){sign}} & {S[n-3:0], 1'b0});

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
    // exp_adjust = lzc_zeroes ? 0 : (DECIMAL_POINT-1 - lz)   (both reference branches agree)
    always_comb begin
        if(lzc_zeroes) begin
            exp_adjust = '0;
        end
        else begin
            exp_adjust = $signed(DECIMAL_POINT-1) - $signed({{(EXP_WIDTH-LZC_WIDTH+1){1'b0}}, leading_zero_count});
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
    assign round_bit = rounding_bits[MAX_SHIFT_AMOUNT-1];
    assign sticky_bit = |rounding_bits[MAX_SHIFT_AMOUNT-2:0];
    assign round_value = round_bit & (sticky_bit | value_before_round[0]);

    // merged: normal_result[n-2:0] = sign_i ? -(vbr+rv) : (vbr+rv)
    //       = (vbr ^ {sign}) + (rv ^ sign)          (exact in Z/(2^(n-1)))
    logic [n-2:0] add_a, add_s;
    assign add_a = value_before_round ^ {(n-1){sign_i}};
    ks_add #(
        .W(n-1)
    ) u_round_add (
        .a(add_a),
        .b({{(n-2){1'b0}}, round_value ^ sign_i}),
        .cin(1'b0),
        .s(add_s)
    );

    logic [n-1:0] normal_result;
    assign normal_result = {sign_i, add_s};
    assign result_o = input_not_zero ? normal_result : '0;
endmodule

// balanced signed maximum tree
module maxtree #(
    parameter int unsigned N = 17,
    parameter int unsigned W = 7
)(
    input  logic [N*W-1:0] items,
    output logic signed [W-1:0] mx
);
    localparam int unsigned NA = N/2;
    localparam int unsigned NB = N - NA;
    generate
        if (N == 1) begin : leaf
            assign mx = $signed(items[0 +: W]);
        end else if (N == 2) begin : pair
            wire signed [W-1:0] a = $signed(items[0 +: W]);
            wire signed [W-1:0] b = $signed(items[W +: W]);
            assign mx = (a > b) ? a : b;
        end else begin : split
            wire signed [W-1:0] ma, mb;
            maxtree #(.N(NA), .W(W)) ua (.items(items[0 +: NA*W]), .mx(ma));
            maxtree #(.N(NB), .W(W)) ub (.items(items[NA*W +: NB*W]), .mx(mb));
            assign mx = (ma > mb) ? ma : mb;
        end
    endgenerate
endmodule

// balanced popcount tree
module popcount #(
    parameter int unsigned N = 17
)(
    input  logic [N-1:0] v,
    output logic [pdpu_pkg::clog2(N+1)-1:0] cnt
);
    localparam int unsigned NA = N/2;
    localparam int unsigned NB = N - NA;
    generate
        if (N == 1) begin : leaf
            assign cnt = v[0];
        end else if (N == 2) begin : pair
            assign cnt = v[0] + v[1];
        end else begin : split
            wire [pdpu_pkg::clog2(NA+1)-1:0] ca;
            wire [pdpu_pkg::clog2(NB+1)-1:0] cb;
            popcount #(.N(NA)) ua (.v(v[0 +: NA]), .cnt(ca));
            popcount #(.N(NB)) ub (.v(v[NA +: NB]), .cnt(cb));
            assign cnt = ca + cb;
        end
    endgenerate
endmodule

// greedy Wallace reduction of N rows (width W) down to two rows using 3:2
// compressors (full adders, no serial carry chains). Depth = 6 for N=18.
module csa_reduce #(
    parameter int unsigned N = 18,
    parameter int unsigned W = 20
)(
    input  logic [N*W-1:0] rows,
    output logic [W-1:0] s_o,
    output logic [W-1:0] c_o
);
    generate
        if (N == 1) begin : b1
            assign s_o = rows[0 +: W];
            assign c_o = '0;
        end else if (N == 2) begin : b2
            assign s_o = rows[0 +: W];
            assign c_o = rows[W +: W];
        end else begin : rec
            localparam int unsigned M = N/3;
            localparam int unsigned REM = N - 3*M;
            localparam int unsigned NN = 2*M + REM;
            wire [M*W-1:0] ts, tc;
            genvar j;
            for (j = 0; j < M; j++) begin : trip
                wire [W-1:0] ra = rows[(3*j+0)*W +: W];
                wire [W-1:0] rb = rows[(3*j+1)*W +: W];
                wire [W-1:0] rc = rows[(3*j+2)*W +: W];
                wire [W-1:0] ab = ra & rb;
                wire [W-1:0] aob = ra | rb;
                assign ts[j*W +: W] = ra ^ rb ^ rc;
                assign tc[j*W +: W] = {((ab[W-2:0]) | (rc[W-2:0] & aob[W-2:0])), 1'b0};
            end
            wire [NN*W-1:0] nxt;
            for (j = 0; j < M; j++) begin : pack
                assign nxt[(2*j)*W +: W] = ts[j*W +: W];
                assign nxt[(2*j+1)*W +: W] = tc[j*W +: W];
            end
            for (j = 0; j < REM; j++) begin : pass
                assign nxt[(2*M+j)*W +: W] = rows[(3*M+j)*W +: W];
            end
            csa_reduce #(.N(NN), .W(W)) sub (.rows(nxt), .s_o(s_o), .c_o(c_o));
        end
    endgenerate
endmodule

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
    localparam int unsigned EXP_WIDTH_I = pdpu_pkg::clog2(n_i-1) + es_i;   // 5
    localparam int unsigned MANT_WIDTH_I = n_i - es_i - 3;                  // 3
    localparam int unsigned EXP_WIDTH_O = pdpu_pkg::clog2(n_o-1) + es_o;   // 6
    localparam int unsigned MANT_WIDTH_O = n_o - es_o - 3;                  // 11
    localparam int unsigned EXP_WIDTH = pdpu_pkg::maximum(EXP_WIDTH_I+1, EXP_WIDTH_O); // 6
    localparam int unsigned MUL_W = 2*(MANT_WIDTH_I+1);                     // 8
    localparam int unsigned CARRY_W = pdpu_pkg::clog2(N+1);                 // 5
    localparam int unsigned SUM_W = ALIGN_WIDTH + CARRY_W;                  // 19
    localparam int unsigned ASH_W = pdpu_pkg::clog2(ALIGN_WIDTH+1);         // 4
    localparam int unsigned NL = N + 1;                                     // 17
    localparam int unsigned ROW_W = SUM_W + 1;                              // 20

    // ********************
    // R1: decode + product exponent add   (in -> pipe1)
    // ********************
    logic [N-1:0] signs_a, signs_b;
    logic [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_a, mants_b;

    wire sign_acc;
    wire [EXP_WIDTH_O:0] rg_exp_acc;
    wire [MANT_WIDTH_O:0] mant_acc;
    posit_decoder #(.n(n_o), .es(es_o)) u_dec_acc (
        .operand_i(acc), .sign_o(sign_acc), .rg_exp_o(rg_exp_acc), .mant_norm_o(mant_acc)
    );

    generate
        genvar gi;
        for (gi = 0; gi < N; gi++) begin : g_dec
            wire sg_a, sg_b;
            wire [EXP_WIDTH_I:0] re_a, re_b;
            wire [MANT_WIDTH_I:0] mn_a, mn_b;
            posit_decoder #(.n(n_i), .es(es_i)) ua_posit_decoder (
                .operand_i(operands_a[gi*n_i +: n_i]), .sign_o(sg_a),
                .rg_exp_o(re_a), .mant_norm_o(mn_a)
            );
            posit_decoder #(.n(n_i), .es(es_i)) ub_posit_decoder (
                .operand_i(operands_b[gi*n_i +: n_i]), .sign_o(sg_b),
                .rg_exp_o(re_b), .mant_norm_o(mn_b)
            );
            assign signs_a[gi] = sg_a;
            assign signs_b[gi] = sg_b;
            assign rg_exp_a[gi*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)] = re_a;
            assign rg_exp_b[gi*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)] = re_b;
            assign mants_a[gi*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] = mn_a;
            assign mants_b[gi*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] = mn_b;
        end
    endgenerate

    wire [N-1:0] signs_ab = signs_a ^ signs_b;
    wire [N:0] signs = {sign_acc, signs_ab};

    logic [N*(EXP_WIDTH+1)-1:0] rg_exp_c;
    generate
        genvar gj;
        for (gj = 0; gj < N; gj++) begin : g_eadd
            wire [EXP_WIDTH_I:0] ea = rg_exp_a[gj*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)];
            wire [EXP_WIDTH_I:0] eb = rg_exp_b[gj*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)];
            ks_add #(.W(EXP_WIDTH+1)) u_eadd (
                .a({ea[EXP_WIDTH_I], ea}),
                .b({eb[EXP_WIDTH_I], eb}),
                .cin(1'b0),
                .s(rg_exp_c[gj*(EXP_WIDTH+1) +: (EXP_WIDTH+1)])
            );
        end
    endgenerate

    // ---- pipe1 ----
    logic [N:0] pipe1_signs;
    logic [N*(EXP_WIDTH+1)-1:0] pipe1_exp;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_ma, pipe1_mb;
    logic [MANT_WIDTH_O:0] pipe1_macc;
    logic [EXP_WIDTH_O:0] pipe1_eacc;

    always_ff @(posedge clk) begin
        pipe1_signs <= signs;
        pipe1_exp   <= rg_exp_c;
        pipe1_ma    <= mants_a;
        pipe1_mb    <= mants_b;
        pipe1_macc  <= mant_acc;
        pipe1_eacc  <= rg_exp_acc;
    end

    // ********************
    // R2: multiply || max tree || popcount || exponent diff   (pipe1 -> pipe2)
    // ********************
    logic [N*MUL_W-1:0] prod;
    generate
        genvar gk;
        for (gk = 0; gk < N; gk++) begin : g_mul
            assign prod[gk*MUL_W +: MUL_W] =
                pipe1_ma[gk*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)] *
                pipe1_mb[gk*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)];
        end
    endgenerate

    logic [NL*(EXP_WIDTH+1)-1:0] items;
    generate
        genvar gl;
        for (gl = 0; gl < N; gl++) begin : g_it
            assign items[gl*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = pipe1_exp[gl*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate
    assign items[N*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = pipe1_eacc;

    wire signed [EXP_WIDTH:0] emax;
    maxtree #(.N(NL), .W(EXP_WIDTH+1)) u_maxtree (.items(items), .mx(emax));

    wire [CARRY_W-1:0] kpop;
    popcount #(.N(NL)) u_pop (.v(pipe1_signs), .cnt(kpop));

    // exponent differences: dv = emax - item (exact, always in [0,127])
    logic [NL*ASH_W-1:0] sav;
    logic [NL-1:0] bigv;
    generate
        genvar gn;
        for (gn = 0; gn < NL; gn++) begin : g_diff
            wire [EXP_WIDTH:0] dv;
            ks_add #(.W(EXP_WIDTH+1)) u_sub (
                .a(emax[EXP_WIDTH:0]),
                .b(~items[gn*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]),
                .cin(1'b1),
                .s(dv)
            );
            // diff >= 16: shifting by dv[3:0] would be wrong -> zero the row instead
            // (diff = 14,15 already shift all bits out naturally)
            assign bigv[gn] = dv[6] | dv[5] | dv[4];
            assign sav[gn*ASH_W +: ASH_W] = dv[ASH_W-1:0];
        end
    endgenerate

    // ---- pipe2 ----
    logic [N*MUL_W-1:0] pipe2_prod;
    logic [NL*ASH_W-1:0] pipe2_sa;
    logic [NL-1:0] pipe2_big;
    logic [N:0] pipe2_signs;
    logic [CARRY_W-1:0] pipe2_k;
    logic [MANT_WIDTH_O:0] pipe2_macc;
    logic [EXP_WIDTH:0] pipe2_max;

    always_ff @(posedge clk) begin
        pipe2_prod  <= prod;
        pipe2_sa    <= sav;
        pipe2_big   <= bigv;
        pipe2_signs <= pipe1_signs;
        pipe2_k     <= kpop;
        pipe2_macc  <= pipe1_macc;
        pipe2_max   <= emax;
    end

    // ********************
    // R3: align + fold + compress   (pipe2 -> pipe3)
    // ********************
    // pre-alignment products, left justified in ALIGN_WIDTH bits
    logic [NL*ALIGN_WIDTH-1:0] prod14;
    generate
        genvar gm;
        for (gm = 0; gm < N; gm++) begin : g_p14
            assign prod14[gm*ALIGN_WIDTH +: ALIGN_WIDTH] =
                {pipe2_prod[gm*MUL_W +: MUL_W], {(ALIGN_WIDTH-MUL_W){1'b0}}};
        end
    endgenerate
    assign prod14[N*ALIGN_WIDTH +: ALIGN_WIDTH] = {1'b0, pipe2_macc, 1'b0};

    // align shift + big-row zeroing + sign fold
    logic [NL*ALIGN_WIDTH-1:0] folded;
    generate
        genvar go;
        for (go = 0; go < NL; go++) begin : g_sh
            wire [ALIGN_WIDTH-1:0] pv = prod14[go*ALIGN_WIDTH +: ALIGN_WIDTH];
            wire [ALIGN_WIDTH-1:0] sh = pv >> pipe2_sa[go*ASH_W +: ASH_W];
            wire [ALIGN_WIDTH-1:0] shz = sh & ~{ALIGN_WIDTH{pipe2_big[go]}};
            assign folded[go*ALIGN_WIDTH +: ALIGN_WIDTH] = shz ^ {ALIGN_WIDTH{pipe2_signs[go]}};
        end
    endgenerate

    // correction row: -k*2^ALIGN_WIDTH + k (mod 2^ROW_W)
    wire [5:0] nc6 = -{1'b0, pipe2_k};
    wire [ROW_W-1:0] corr = {nc6, {(ROW_W-6-CARRY_W){1'b0}}, pipe2_k};

    logic [(NL+1)*ROW_W-1:0] csa_rows;
    generate
        genvar gp;
        for (gp = 0; gp < NL; gp++) begin : g_rows
            assign csa_rows[gp*ROW_W +: ROW_W] = {{(ROW_W-ALIGN_WIDTH){1'b0}}, folded[gp*ALIGN_WIDTH +: ALIGN_WIDTH]};
        end
    endgenerate
    assign csa_rows[NL*ROW_W +: ROW_W] = corr;

    wire [ROW_W-1:0] csa_s, csa_c;
    csa_reduce #(.N(NL+1), .W(ROW_W)) u_csa (.rows(csa_rows), .s_o(csa_s), .c_o(csa_c));

    // ---- pipe3 ----
    logic [ROW_W-1:0] pipe3_s, pipe3_c;
    logic [EXP_WIDTH:0] pipe3_max;
    always_ff @(posedge clk) begin
        pipe3_s   <= csa_s;
        pipe3_c   <= csa_c;
        pipe3_max <= pipe2_max;
    end

    // ********************
    // R4: final carry-propagate add   (pipe3 -> pipe4)
    // ********************
    wire [ROW_W-1:0] sum_full;
    ks_add #(.W(ROW_W)) u_cpa (.a(pipe3_s), .b(pipe3_c), .cin(1'b0), .s(sum_full));

    logic [ROW_W-1:0] pipe4_sum;
    logic [EXP_WIDTH:0] pipe4_max;
    always_ff @(posedge clk) begin
        pipe4_sum <= sum_full;
        pipe4_max <= pipe3_max;
    end

    // ********************
    // R5: abs + normalize   (pipe4 -> pipe5)
    // ********************
    wire p4sign = pipe4_sum[ROW_W-1];
    wire [ROW_W-1:0] negv;
    ks_add #(.W(ROW_W)) u_neg (.a(~pipe4_sum), .b('0), .cin(1'b1), .s(negv));
    wire [SUM_W-1:0] sum_c = p4sign ? negv[SUM_W-1:0] : pipe4_sum[SUM_W-1:0];

    wire signed [EXP_WIDTH:0] exp_adj;
    wire [SUM_W-1:0] sum_norm;
    mantissa_norm #(
        .WIDTH(SUM_W),
        .EXP_WIDTH(EXP_WIDTH),
        .DECIMAL_POINT(CARRY_W+2)
    ) u_mantissa_norm (
        .operand_i(sum_c),
        .exp_adjust(exp_adj),
        .result_o(sum_norm)
    );

    wire signed [EXP_WIDTH:0] fre = pipe4_max + exp_adj;
    wire sticky = |sum_norm[SUM_W-MANT_WIDTH_O-3-:6];
    wire [MANT_WIDTH_O+2:0] fman = {sum_norm[SUM_W-1:SUM_W-MANT_WIDTH_O-2], sticky};

    logic pipe5_sign;
    logic signed [EXP_WIDTH:0] pipe5_re;
    logic [MANT_WIDTH_O+2:0] pipe5_mant;
    always_ff @(posedge clk) begin
        pipe5_sign <= p4sign;
        pipe5_re   <= fre;
        pipe5_mant <= fman;
    end

    // ********************
    // R6: encode   (pipe5 -> out)
    // ********************
    wire [n_o-1:0] res;
    posit_encoder #(
        .n(n_o),
        .es(es_o),
        .EXP_WIDTH(EXP_WIDTH),
        .MANT_WIDTH(MANT_WIDTH_O+2)
    ) u_posit_encoder (
        .sign_i(pipe5_sign),
        .rg_exp_i(pipe5_re),
        .mant_norm_i(pipe5_mant),
        .result_o(res)
    );
    assign result_o = res;

endmodule

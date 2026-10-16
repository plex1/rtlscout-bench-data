// Simplified fixed-configuration pdpu_top_pipelined
// N=16, posit(8,2) inputs, posit(16,2) accumulator/result.
// Intended to be functionally identical to starting_point.v.

module posit_decoder8 (
    input  logic [7:0] operand_i,
    output logic       sign_o,
    output logic signed [5:0] rg_exp_o,
    output logic [3:0] mant_norm_o
);
    logic [6:0] operand_value;
    logic       regS;
    logic [6:0] lzc_operand;
    logic [2:0] leading_zero_count;
    logic       all_zero;

    assign sign_o = operand_i[7];
    assign operand_value = sign_o ? (~operand_i[6:0] + 1'b1) : operand_i[6:0];
    assign regS = operand_value[6];
    assign lzc_operand = regS ? ~operand_value : operand_value;
    assign all_zero = ~(|lzc_operand);

    logic found;
    always_comb begin
        leading_zero_count = 3'd0;
        found = 1'b0;
        for (int i = 6; i >= 0; i--) begin
            if (lzc_operand[i] && !found) begin
                leading_zero_count = 3'(6 - i);
                found = 1'b1;
            end
        end
    end

    logic [2:0] runlength;
    logic signed [3:0] regime_k;
    logic [6:0] op_no_rg;
    logic [2:0] exp;
    logic       implicit_bit;

    assign runlength = all_zero ? 3'd7 : leading_zero_count;
    assign regime_k = regS ? {1'b0, runlength - 3'd1} :
                             {1'b1, ~runlength + 3'd1};

    // Match the original expressions, including their 4-bit truncation of
    // (regime_k << 2) and the zero-extended 3-bit exp field.
    assign op_no_rg = operand_value << (all_zero ? 3'd7 : (leading_zero_count + 3'd1));
    assign exp = {1'b0, op_no_rg[6:5]};
    assign rg_exp_o = $signed((regime_k << 2) | exp);
    assign implicit_bit = |operand_i[6:0];
    assign mant_norm_o = {implicit_bit, op_no_rg[4:2]};
endmodule

module posit_decoder16 (
    input  logic [15:0] operand_i,
    output logic        sign_o,
    output logic signed [6:0] rg_exp_o,
    output logic [11:0] mant_norm_o
);
    logic [14:0] operand_value;
    logic        regS;
    logic [14:0] lzc_operand;
    logic [3:0]  leading_zero_count;
    logic        all_zero;

    assign sign_o = operand_i[15];
    assign operand_value = sign_o ? (~operand_i[14:0] + 1'b1) : operand_i[14:0];
    assign regS = operand_value[14];
    assign lzc_operand = regS ? ~operand_value : operand_value;
    assign all_zero = ~(|lzc_operand);

    logic found;
    always_comb begin
        leading_zero_count = 4'd0;
        found = 1'b0;
        for (int i = 14; i >= 0; i--) begin
            if (lzc_operand[i] && !found) begin
                leading_zero_count = 4'(14 - i);
                found = 1'b1;
            end
        end
    end

    logic [3:0] runlength;
    logic signed [4:0] regime_k;
    logic [14:0] op_no_rg;
    logic [2:0] exp;
    logic        implicit_bit;

    assign runlength = all_zero ? 4'd15 : leading_zero_count;
    assign regime_k = regS ? {1'b0, runlength - 4'd1} :
                             {1'b1, ~runlength + 4'd1};

    // Match the original expressions, including their 5-bit truncation of
    // (regime_k << 2) and the zero-extended 3-bit exp field.
    assign op_no_rg = operand_value << (all_zero ? 4'd15 : (leading_zero_count + 4'd1));
    assign exp = {1'b0, op_no_rg[14:13]};
    assign rg_exp_o = $signed((regime_k << 2) | exp);
    assign implicit_bit = |operand_i[14:0];
    assign mant_norm_o = {implicit_bit, op_no_rg[12:2]};
endmodule

module posit_encoder16 (
    input  logic        sign_i,
    input  logic signed [6:0] rg_exp_i,
    input  logic [13:0] mant_norm_i,
    output logic [15:0] result_o
);
    logic input_not_zero;
    logic signed [4:0] regime_k;
    logic [1:0] exp;
    logic       sign_k;
    logic [14:0] regime;
    logic signed [4:0] regime_bits;
    logic [29:0] rg_exp_mant;
    logic [45:0] value_before_shift, value_after_shift;
    logic [4:0] shift_amount;
    logic [14:0] value_before_round;
    logic [15:0] rounding_bits;
    logic       round_bit, sticky_bit, round_value;
    logic [14:0] value_after_round;
    logic [15:0] normal_result;

    assign input_not_zero = mant_norm_i[13];
    assign regime_k = rg_exp_i[6:2];
    assign exp = rg_exp_i[1:0];
    assign sign_k = rg_exp_i[6];
    assign regime = sign_k ? 15'b000000000000001 : 15'b111111111111110;

    always_comb begin
        if (sign_k) regime_bits = ~regime_k + 5'sd2;
        else        regime_bits = regime_k + 5'sd2;
    end

    assign rg_exp_mant = {regime, exp, mant_norm_i[12:0]};
    assign value_before_shift = rg_exp_mant << 16;
    assign shift_amount = (regime_bits >= 5'sd16) ? 5'd16 : regime_bits[4:0];
    assign value_after_shift = value_before_shift >> shift_amount;

    assign {value_before_round, rounding_bits} = value_after_shift[30:0];
    assign round_bit   = rounding_bits[15];
    assign sticky_bit  = |rounding_bits[14:0];
    assign round_value = round_bit & (sticky_bit | value_before_round[0]);
    assign value_after_round = value_before_round + round_value;
    assign normal_result = sign_i ? {1'b1, ~value_after_round + 1'b1} :
                                    {1'b0, value_after_round};
    assign result_o = input_not_zero ? normal_result : 16'h0000;
endmodule

module max_tree #(
    parameter int N = 17,
    parameter int W = 7
)(
    input  logic signed [N*W-1:0] x,
    output logic signed [W-1:0] y
);
    generate
        if (N == 1) begin : g1
            assign y = x;
        end else if (N == 2) begin : g2
            assign y = ($signed(x[0 +: W]) > $signed(x[W +: W])) ? x[0 +: W] : x[W +: W];
        end else begin : gs
            localparam int NA = N/2;
            localparam int NB = N - NA;
            logic signed [W-1:0] ya, yb;
            max_tree #(.N(NA), .W(W)) ua (.x(x[0 +: NA*W]), .y(ya));
            max_tree #(.N(NB), .W(W)) ub (.x(x[NA*W +: NB*W]), .y(yb));
            assign y = ($signed(ya) > $signed(yb)) ? ya : yb;
        end
    endgenerate
endmodule

module add_tree #(
    parameter int N = 17,
    parameter int W = 20
)(
    input  logic [N*W-1:0] x,
    output logic [W-1:0] y
);
    generate
        if (N == 1) begin : g1
            assign y = x;
        end else if (N == 2) begin : g2
            assign y = x[0 +: W] + x[W +: W];
        end else begin : gs
            localparam int NA = N/2;
            localparam int NB = N - NA;
            logic [W-1:0] ya, yb;
            add_tree #(.N(NA), .W(W)) ua (.x(x[0 +: NA*W]), .y(ya));
            add_tree #(.N(NB), .W(W)) ub (.x(x[NA*W +: NB*W]), .y(yb));
            assign y = ya + yb;
        end
    endgenerate
endmodule

module lzc19 (
    input  logic [18:0] operand_i,
    output logic [4:0]  cnt_o,
    output logic        empty_o
);
    logic found;
    always_comb begin
        empty_o = ~(|operand_i);
        cnt_o = 5'd0;
        found = 1'b0;
        for (int i = 18; i >= 0; i--) begin
            if (operand_i[i] && !found) begin
                cnt_o = 5'(18 - i);
                found = 1'b1;
            end
        end
    end
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
    input  logic [n_o-1:0] acc,
    output logic [n_o-1:0] result_o
);
    localparam int EW_I = 6;   // exponent width from posit8 decoder (incl. sign)
    localparam int EW_O = 7;   // exponent width from posit16 decoder / products
    localparam int MW_I = 4;   // mantissa width from posit8 decoder
    localparam int MW_O = 12;  // mantissa width from posit16 decoder
    localparam int MULW = 8;   // 4x4 product width
    localparam int EW = 7;     // common signed exponent width
    localparam int SW = 20;    // signed accumulator item width (SUM_WIDTH+1)
    localparam int SUMW = 19;  // positive sum width

    // -------------------------
    // Decode inputs (combinational)
    // -------------------------
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*EW_I-1:0] rg_exp_a, rg_exp_b;
    logic [N*MW_I-1:0] mants_norm_a, mants_norm_b;

    generate
        genvar di;
        for (di = 0; di < N; di++) begin : g_decode_in
            posit_decoder8 ua (
                .operand_i (operands_a[di*n_i +: n_i]),
                .sign_o    (signs_a[di]),
                .rg_exp_o  (rg_exp_a[di*EW_I +: EW_I]),
                .mant_norm_o(mants_norm_a[di*MW_I +: MW_I])
            );
            posit_decoder8 ub (
                .operand_i (operands_b[di*n_i +: n_i]),
                .sign_o    (signs_b[di]),
                .rg_exp_o  (rg_exp_b[di*EW_I +: EW_I]),
                .mant_norm_o(mants_norm_b[di*MW_I +: MW_I])
            );
        end
    endgenerate

    logic sign_acc_c;
    logic signed [EW_O-1:0] rg_exp_acc_c;
    logic [MW_O-1:0] mant_norm_acc_c;
    posit_decoder16 uacc (
        .operand_i (acc),
        .sign_o    (sign_acc_c),
        .rg_exp_o  (rg_exp_acc_c),
        .mant_norm_o(mant_norm_acc_c)
    );

    logic [N-1:0] signs_ab;
    logic [N:0] signs;
    logic signed [N*EW-1:0] rg_exp_c;

    assign signs_ab = signs_a ^ signs_b;
    assign signs = {sign_acc_c, signs_ab};

    generate
        genvar ex;
        for (ex = 0; ex < N; ex++) begin : g_exp_sum
            assign rg_exp_c[ex*EW +: EW] =
                $signed(rg_exp_a[ex*EW_I +: EW_I]) + $signed(rg_exp_b[ex*EW_I +: EW_I]);
        end
    endgenerate

    // -------------------------
    // Pipeline 1
    // -------------------------
    logic [N:0] p1_signs;
    logic signed [N*EW-1:0] p1_rg_exp_c;
    logic [N*MW_I-1:0] p1_mants_a, p1_mants_b;
    logic signed [EW_O-1:0] p1_rg_exp_acc;
    logic [MW_O-1:0] p1_mant_acc;

    always_ff @(posedge clk) begin
        p1_signs      <= signs;
        p1_rg_exp_c   <= rg_exp_c;
        p1_mants_a    <= mants_norm_a;
        p1_mants_b    <= mants_norm_b;
        p1_rg_exp_acc <= rg_exp_acc_c;
        p1_mant_acc   <= mant_norm_acc_c;
    end

    // -------------------------
    // Stage 2: multiply products, find max exponent
    // -------------------------
    logic [N*MULW-1:0] mants_norm_c_c;
    generate
        genvar ml;
        for (ml = 0; ml < N; ml++) begin : g_mul
            assign mants_norm_c_c[ml*MULW +: MULW] =
                p1_mants_a[ml*MW_I +: MW_I] * p1_mants_b[ml*MW_I +: MW_I];
        end
    endgenerate

    logic signed [(N+1)*EW-1:0] rg_exp_items_c;
    generate
        genvar ei;
        for (ei = 0; ei < N; ei++) begin : g_exp_items
            assign rg_exp_items_c[ei*EW +: EW] = p1_rg_exp_c[ei*EW +: EW];
        end
    endgenerate
    assign rg_exp_items_c[N*EW +: EW] = p1_rg_exp_acc;

    logic signed [EW-1:0] rg_exp_max_c;
    max_tree #(.N(N+1), .W(EW)) u_max (
        .x(rg_exp_items_c),
        .y(rg_exp_max_c)
    );

    // -------------------------
    // Pipeline 2
    // -------------------------
    logic [N:0] p2_signs;
    logic signed [(N+1)*EW-1:0] p2_exp_items;
    logic signed [EW-1:0] p2_exp_max;
    logic [N*MULW-1:0] p2_mants_c;
    logic [MW_O-1:0] p2_mant_acc;

    always_ff @(posedge clk) begin
        p2_signs    <= p1_signs;
        p2_exp_items <= rg_exp_items_c;
        p2_exp_max  <= rg_exp_max_c;
        p2_mants_c  <= mants_norm_c_c;
        p2_mant_acc <= p1_mant_acc;
    end

    // -------------------------
    // Stage 3: align products and form two's-complement terms
    // -------------------------
    logic [(N+1)*SW-1:0] mantissa_comp_c;

    generate
        genvar al;
        for (al = 0; al < N; al++) begin : g_align_lane
            logic [13:0] product;
            logic [6:0]  diff;
            logic [13:0] shifted;
            logic [19:0] ext;
            assign product = p2_mants_c[al*MULW +: MULW] << 6;
            assign diff = $unsigned(p2_exp_max - $signed(p2_exp_items[al*EW +: EW]));
            assign shifted = (diff >= 7'd14) ? 14'h0000 : (product >> diff[3:0]);
            assign ext = {6'b0, shifted};
            assign mantissa_comp_c[al*SW +: SW] =
                p2_signs[al] ? (~ext + 1'b1) : ext;
        end
    endgenerate

    logic [13:0] acc_product;
    logic [6:0]  acc_diff;
    logic [13:0] acc_shifted;
    logic [19:0] acc_ext;
    assign acc_product = p2_mant_acc << 1;
    assign acc_diff = $unsigned(p2_exp_max - $signed(p2_exp_items[N*EW +: EW]));
    assign acc_shifted = (acc_diff >= 7'd14) ? 14'h0000 : (acc_product >> acc_diff[3:0]);
    assign acc_ext = {6'b0, acc_shifted};
    assign mantissa_comp_c[N*SW +: SW] =
        p2_signs[N] ? (~acc_ext + 1'b1) : acc_ext;

    // -------------------------
    // Pipeline 3
    // -------------------------
    logic [(N+1)*SW-1:0] p3_mantissa_comp;
    logic signed [EW-1:0] p3_exp_max;

    always_ff @(posedge clk) begin
        p3_mantissa_comp <= mantissa_comp_c;
        p3_exp_max       <= p2_exp_max;
    end

    // -------------------------
    // Stage 4: adder tree -> sum magnitude
    // -------------------------
    logic [SW-1:0] sum_result_c;
    add_tree #(.N(N+1), .W(SW)) u_sumtree (
        .x(p3_mantissa_comp),
        .y(sum_result_c)
    );

    // -------------------------
    // Pipeline 4
    // -------------------------
    logic [SW-1:0] p4_sum_result;
    logic signed [EW-1:0] p4_exp_max;

    always_ff @(posedge clk) begin
        p4_sum_result <= sum_result_c;
        p4_exp_max    <= p3_exp_max;
    end

    // -------------------------
    // Stage 5: normalize and prepare encoder inputs
    // -------------------------
    logic final_sign_c;
    logic [SUMW-1:0] sum_c;
    logic [4:0] lzc_cnt;
    logic lzc_empty;
    logic signed [EW-1:0] rg_exp_adjust;
    logic signed [EW-1:0] final_rg_exp_c;
    logic [18:0] sum_norm;
    logic [13:0] final_mant_c;

    assign final_sign_c = p4_sum_result[SW-1];

    logic [SW-1:0] sum_mag;
    assign sum_mag = final_sign_c ? (~p4_sum_result + 1'b1) : p4_sum_result;
    assign sum_c = sum_mag[SUMW-1:0];

    lzc19 u_lzc (
        .operand_i(sum_c),
        .cnt_o(lzc_cnt),
        .empty_o(lzc_empty)
    );

    always_comb begin
        if (lzc_empty) rg_exp_adjust = 7'sd0;
        else if (lzc_cnt <= 5'd6) rg_exp_adjust = 7'sd7 - $signed({1'b0, lzc_cnt}) - 7'sd1;
        else rg_exp_adjust = -($signed({1'b0, lzc_cnt}) - 7'sd6);
    end

    assign sum_norm = sum_c << lzc_cnt;
    assign final_rg_exp_c = p4_exp_max + rg_exp_adjust;
    assign final_mant_c = {sum_norm[18:6], |sum_norm[5:0]};

    // -------------------------
    // Pipeline 5 + final encoder
    // -------------------------
    logic p5_final_sign;
    logic signed [EW-1:0] p5_final_rg_exp;
    logic [13:0] p5_final_mant;

    always_ff @(posedge clk) begin
        p5_final_sign   <= final_sign_c;
        p5_final_rg_exp <= final_rg_exp_c;
        p5_final_mant   <= final_mant_c;
    end

    posit_encoder16 u_enc (
        .sign_i(p5_final_sign),
        .rg_exp_i(p5_final_rg_exp),
        .mant_norm_i(p5_final_mant),
        .result_o(result_o)
    );

endmodule
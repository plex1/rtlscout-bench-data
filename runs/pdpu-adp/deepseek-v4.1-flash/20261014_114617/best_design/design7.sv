`include "starting_point.v"

// =====================================================================
// max_reduce: one level of a pairwise maximum reduction over signed
// (W+1)-bit exponents.  in_i is flat, element e is [e*(W+1) +: W+1].
// =====================================================================
module max_reduce #(
    parameter int unsigned NV = 17,
    parameter int unsigned W  = 6
)(
    input  logic signed [NV*(W+1)-1:0] in_i,
    output logic signed [((NV+1)/2)*(W+1)-1:0] out_o
);
    localparam int unsigned NO = (NV+1)/2;
    generate
        genvar p;
        for(p=0;p<NO;p++) begin: level
            if (2*p+1 < NV) begin: cmp
                comparator #(.WIDTH(W)) u_cmp(
                    .operand_a(in_i[2*p*(W+1) +: (W+1)]),
                    .operand_b(in_i[(2*p+1)*(W+1) +: (W+1)]),
                    .result_o(out_o[p*(W+1) +: (W+1)])
                );
            end else begin: pass
                assign out_o[p*(W+1) +: (W+1)] = in_i[2*p*(W+1) +: (W+1)];
            end
        end
    endgenerate
endmodule

// =====================================================================
// design7: balanced 5-stage pipeline.
//   * lane products registered (8b)
//   * exponent max tree split 3 levels (stage 2) / 2 levels (stage 3)
//   * one's complement + popcount correction instead of per-lane
//     two's-complement incrementers
//   * stage 4 does the whole CSA compression + final add
// =====================================================================
module pdpu_top_pipelined (
    input  logic              clk,
    input  logic [127:0]      operands_a,
    input  logic [127:0]      operands_b,
    input  logic [15:0]       acc,
    output logic [15:0]       result_o
);
    localparam int unsigned N  = 16;
    localparam int unsigned n_i = 8;
    localparam int unsigned es_i = 2;
    localparam int unsigned n_o = 16;
    localparam int unsigned es_o = 2;
    localparam int unsigned ALIGN_WIDTH = 14;

    localparam int unsigned EXP_WIDTH_I  = pdpu_pkg::clog2(n_i-1) + es_i;   // 5
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;                       // 3
    localparam int unsigned EXP_WIDTH_O  = pdpu_pkg::clog2(n_o-1) + es_o;    // 6
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;                       // 11
    localparam int unsigned EXP_WIDTH    = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);   // 6
    localparam int unsigned MUL_WIDTH    = 2*(MANT_WIDTH_I+1);                              // 8
    localparam int unsigned CARRY_WIDTH  = pdpu_pkg::clog2(N+1);                            // 5
    localparam int unsigned SUM_WIDTH    = ALIGN_WIDTH + CARRY_WIDTH;                       // 19
    localparam int unsigned DW           = SUM_WIDTH+1;                                     // 20
    localparam int unsigned SHIFT_WIDTH  = pdpu_pkg::clog2(ALIGN_WIDTH+1);                  // 4
    localparam int unsigned EW           = EXP_WIDTH_O+1;                                    // 7
    localparam int unsigned EW_I         = EXP_WIDTH_I+1;                                    // 6

    // ---------------- stage 1 logic: decode ----------------
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*EW_I-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    logic sign_acc;
    logic signed [EW-1:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    logic [N-1:0] signs_ab;
    logic [N:0]   signs;
    logic signed [N*EW-1:0] rg_exp_c;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(.n(n_i), .es(es_i)) ua_posit_decoder(
                .operand_i(operands_a[i*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[i*EW_I +: EW_I]),
                .mant_norm_o(mants_norm_a[i*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(.n(n_i), .es(es_i)) ub_posit_decoder(
                .operand_i(operands_b[i*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[i*EW_I +: EW_I]),
                .mant_norm_o(mants_norm_b[i*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
        end
    endgenerate

    posit_decoder #(.n(n_o), .es(es_o)) uc_posit_decoder(
        .operand_i(acc),
        .sign_o(sign_acc),
        .rg_exp_o(rg_exp_acc),
        .mant_norm_o(mant_norm_acc)
    );

    assign signs_ab = signs_a ^ signs_b;
    assign signs    = {sign_acc, signs_ab};

    generate
        genvar m;
        for(m=0;m<N;m++) begin
            assign rg_exp_c[m*EW +: EW] =
                $signed(rg_exp_a[m*EW_I +: EW_I]) + $signed(rg_exp_b[m*EW_I +: EW_I]);
        end
    endgenerate

    // ---------------- pipe1 ----------------
    logic [N:0] pipe1_signs;
    logic signed [N*EW-1:0] pipe1_exp_c;
    logic signed [EW-1:0] pipe1_rg_exp_acc;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic [MANT_WIDTH_O:0] pipe1_mant_norm_acc;

    `FFNR(pipe1_signs,signs,clk)
    `FFNR(pipe1_exp_c,rg_exp_c,clk)
    `FFNR(pipe1_rg_exp_acc,rg_exp_acc,clk)
    `FFNR(pipe1_mants_norm_a,mants_norm_a,clk)
    `FFNR(pipe1_mants_norm_b,mants_norm_b,clk)
    `FFNR(pipe1_mant_norm_acc,mant_norm_acc,clk)

    // ---------------- stage 2 logic ----------------
    logic signed [(N+1)*EW-1:0] exp_items;
    generate
        genvar u;
        for(u=0;u<N;u++) assign exp_items[u*EW +: EW] = pipe1_exp_c[u*EW +: EW];
    endgenerate
    assign exp_items[N*EW +: EW] = pipe1_rg_exp_acc;

    logic [N*MUL_WIDTH-1:0] mants_norm_c;
    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            logic [MUL_WIDTH-1:0] m_sum, m_carry;
            radix4_booth_multiplier #(.WIDTH_A(MANT_WIDTH_I+1), .WIDTH_B(MANT_WIDTH_I+1))
            u_radix4_booth_multiplier(
                .operand_a(pipe1_mants_norm_a[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .operand_b(pipe1_mants_norm_b[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .sum_o(m_sum),
                .carry_o(m_carry)
            );
            assign mants_norm_c[j*MUL_WIDTH +: MUL_WIDTH] = m_sum + m_carry;
        end
    endgenerate

    logic [ALIGN_WIDTH-1:0] acc_product;
    assign acc_product = {2'b0, pipe1_mant_norm_acc};

    // exponent max: 17 -> 9 -> 5 -> 3 (3 levels in stage 2)
    logic signed [9*EW-1:0] mx9;
    logic signed [5*EW-1:0] mx5;
    logic signed [3*EW-1:0] mx3;
    max_reduce #(.NV(17), .W(EXP_WIDTH)) r0 (.in_i(exp_items), .out_o(mx9));
    max_reduce #(.NV(9),  .W(EXP_WIDTH)) r1 (.in_i(mx9),       .out_o(mx5));
    max_reduce #(.NV(5),  .W(EXP_WIDTH)) r2 (.in_i(mx5),       .out_o(mx3));

    // ---------------- pipe2 ----------------
    logic [N:0] pipe2_signs;
    logic [N*MUL_WIDTH-1:0] pipe2_mants_norm_c;
    logic [ALIGN_WIDTH-1:0] pipe2_acc_product;
    logic signed [(N+1)*EW-1:0] pipe2_exp_items;
    logic signed [3*EW-1:0] pipe2_mx3;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_mants_norm_c,mants_norm_c,clk)
    `FFNR(pipe2_acc_product,acc_product,clk)
    `FFNR(pipe2_exp_items,exp_items,clk)
    `FFNR(pipe2_mx3,mx3,clk)

    // ---------------- stage 3 logic: max finish + align ----------------
    logic signed [2*EW-1:0] mx2;
    logic signed [EW-1:0]   exp_max;
    max_reduce #(.NV(3), .W(EXP_WIDTH)) r3 (.in_i(pipe2_mx3), .out_o(mx2));
    max_reduce #(.NV(2), .W(EXP_WIDTH)) r4 (.in_i(mx2),       .out_o(exp_max));

    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    generate
        genvar k;
        for(k=0;k<N;k++) begin
            assign product[k*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_mants_norm_c[k*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
        end
    endgenerate
    assign product[N*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe2_acc_product;

    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z;
        for(z=0;z<N+1;z++) begin: sh
            logic [EW-1:0] diff;
            assign diff = $unsigned(exp_max - pipe2_exp_items[z*EW +: EW]);
            assign shift_amount[z*SHIFT_WIDTH +: SHIFT_WIDTH] =
                (|diff[SHIFT_WIDTH +: (EW-SHIFT_WIDTH)]) ? ALIGN_WIDTH : diff[SHIFT_WIDTH-1:0];
        end
    endgenerate

    logic [(N+1)*DW-1:0] lane_op;
    logic [4:0] n_neg;
    generate
        genvar s;
        for(s=0;s<N+1;s++) begin: align
            logic [ALIGN_WIDTH-1:0] prod_sh;
            logic [ALIGN_WIDTH:0]   onec;
            barrel_shifter #(.WIDTH(ALIGN_WIDTH), .SHIFT_WIDTH(SHIFT_WIDTH), .MODE(1'b1))
            u_barrel_shifter(
                .operand_i(product[s*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(shift_amount[s*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(prod_sh)
            );
            assign onec = ~{1'b0, prod_sh};
            assign lane_op[s*DW +: DW] = pipe2_signs[s]
                ? {{(DW-ALIGN_WIDTH-1){1'b1}}, onec}
                : {{(DW-ALIGN_WIDTH){1'b0}}, prod_sh};
        end
    endgenerate

    // number of negative contributions
    logic [4:0] n0,n1,n2,n3,n4,n5,n6,n7,n8;
    assign n0 = {4'b0, pipe2_signs[0]} + {4'b0, pipe2_signs[1]};
    assign n1 = {4'b0, pipe2_signs[2]} + {4'b0, pipe2_signs[3]};
    assign n2 = {4'b0, pipe2_signs[4]} + {4'b0, pipe2_signs[5]};
    assign n3 = {4'b0, pipe2_signs[6]} + {4'b0, pipe2_signs[7]};
    assign n4 = n0 + n1;
    assign n5 = n2 + n3;
    assign n6 = n4 + n5;
    assign n7 = {4'b0, pipe2_signs[8]}  + {4'b0, pipe2_signs[9]}  + {4'b0, pipe2_signs[10]} +
                {4'b0, pipe2_signs[11]} + {4'b0, pipe2_signs[12]} + {4'b0, pipe2_signs[13]} +
                {4'b0, pipe2_signs[14]} + {4'b0, pipe2_signs[15]} + {4'b0, pipe2_signs[16]};
    assign n8 = n6 + n7;
    assign n_neg = n8;

    // ---------------- pipe3 ----------------
    logic signed [EW-1:0] pipe3_exp_max;
    logic [(N+1)*DW-1:0] pipe3_lane_op;
    logic [4:0] pipe3_n_neg;

    `FFNR(pipe3_exp_max,exp_max,clk)
    `FFNR(pipe3_lane_op,lane_op,clk)
    `FFNR(pipe3_n_neg,n_neg,clk)

    // ---------------- stage 4 logic: compression ----------------
    logic [(N+2)*DW-1:0] tree_in;
    always_comb begin
        tree_in = '0;
        for (int t = 0; t < N+1; t++) tree_in[t*DW +: DW] = pipe3_lane_op[t*DW +: DW];
        tree_in[(N+1)*DW +: DW] = {{(DW-5){1'b0}}, pipe3_n_neg};
    end

    logic [DW-1:0] csa_sum, csa_carry, sum_full;
    csa_tree #(.N(N+2), .WIDTH_I(DW), .WIDTH_O(DW)) u_csa_tree(
        .operands_i(tree_in),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );
    assign sum_full = csa_sum + csa_carry;

    // ---------------- pipe4 ----------------
    logic [DW-1:0] pipe4_sum_full;
    logic signed [EW-1:0] pipe4_exp_max;

    `FFNR(pipe4_sum_full,sum_full,clk)
    `FFNR(pipe4_exp_max,pipe3_exp_max,clk)

    // ---------------- stage 5 logic: sign/mag + normalization + round -------
    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    logic signed [EW-1:0] rg_exp_adjust;
    logic signed [EW-1:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;
    logic [MANT_WIDTH_O+2:0] final_mant;
    logic [SUM_WIDTH-MANT_WIDTH_O-3:0] sticky_bits;
    logic sticky_bit;

    assign final_sign = pipe4_sum_full[DW-1];
    assign sum_c = final_sign ? (~pipe4_sum_full[SUM_WIDTH-1:0]+1) : pipe4_sum_full[SUM_WIDTH-1:0];

    mantissa_norm #(.WIDTH(SUM_WIDTH), .EXP_WIDTH(EXP_WIDTH), .DECIMAL_POINT(CARRY_WIDTH+2))
    u_mantissa_norm(
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    assign final_rg_exp = pipe4_exp_max + rg_exp_adjust;

    assign sticky_bits = sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0];
    assign sticky_bit = |sticky_bits;
    assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2], sticky_bit};

    // ---------------- pipe5 ----------------
    logic pipe5_final_sign;
    logic signed [EW-1:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;

    `FFNR(pipe5_final_sign,final_sign,clk)
    `FFNR(pipe5_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe5_final_mant,final_mant,clk)

    // ---------------- stage 6 logic: posit encoder ----------------
    posit_encoder #(.n(n_o), .es(es_o), .EXP_WIDTH(EXP_WIDTH), .MANT_WIDTH(MANT_WIDTH_O+2))
    u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );

endmodule

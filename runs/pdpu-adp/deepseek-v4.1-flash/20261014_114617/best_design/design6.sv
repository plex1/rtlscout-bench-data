`include "starting_point.v"

// =====================================================================
// design6: aggressive register-area reduction.
//
//  * lane products registered (8b) instead of booth sum/carry (2x8b)
//  * shift amounts computed in stage 2 (registered, no exp vector needed)
//  * one's complement + popcount correction instead of 17 two's-complement
//    incrementers
//  * the whole CSA tree AND the final 20-bit addition sit in stage 3, so
//    pipeline register 3 shrinks from 17x20 bits to a single 20-bit sum
//  * normalization is split across stages 4/5
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

    localparam int unsigned EXP_WIDTH_I  = pdpu_pkg::clog2(n_i-1) + es_i;
    localparam int unsigned MANT_WIDTH_I = n_i-es_i-3;
    localparam int unsigned EXP_WIDTH_O  = pdpu_pkg::clog2(n_o-1) + es_o;
    localparam int unsigned MANT_WIDTH_O = n_o-es_o-3;
    localparam int unsigned EXP_WIDTH    = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    localparam int unsigned MUL_WIDTH    = 2*(MANT_WIDTH_I+1);
    localparam int unsigned CARRY_WIDTH  = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH    = ALIGN_WIDTH + CARRY_WIDTH;   // 19
    localparam int unsigned DW           = SUM_WIDTH+1;                 // 20
    localparam int unsigned SHIFT_WIDTH  = pdpu_pkg::clog2(ALIGN_WIDTH+1);

    // ---------------- stage 1 ----------------
    logic [N-1:0] signs_a, signs_b;
    logic signed [N*(EXP_WIDTH_I+1)-1:0] rg_exp_a, rg_exp_b;
    logic [N*(MANT_WIDTH_I+1)-1:0] mants_norm_a, mants_norm_b;

    logic sign_acc;
    logic signed [EXP_WIDTH_O:0] rg_exp_acc;
    logic [MANT_WIDTH_O:0] mant_norm_acc;

    logic [N-1:0] signs_ab;
    logic [N:0]   signs;
    logic signed [N*(EXP_WIDTH+1)-1:0] rg_exp_c;

    generate
        genvar i;
        for(i=0;i<N;i++) begin: posit_decoding
            posit_decoder #(.n(n_i), .es(es_i)) ua_posit_decoder(
                .operand_i(operands_a[i*n_i +: n_i]),
                .sign_o(signs_a[i]),
                .rg_exp_o(rg_exp_a[i*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
                .mant_norm_o(mants_norm_a[i*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)])
            );
            posit_decoder #(.n(n_i), .es(es_i)) ub_posit_decoder(
                .operand_i(operands_b[i*n_i +: n_i]),
                .sign_o(signs_b[i]),
                .rg_exp_o(rg_exp_b[i*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]),
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
            assign rg_exp_c[m*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] =
                $signed(rg_exp_a[m*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]) +
                $signed(rg_exp_b[m*(EXP_WIDTH_I+1) +: (EXP_WIDTH_I+1)]);
        end
    endgenerate

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

    // ---------------- stage 2 ----------------
    logic [N*MUL_WIDTH-1:0] mul_sum, mul_carry, mants_norm_c;
    generate
        genvar j;
        for(j=0;j<N;j++) begin: multiplication
            radix4_booth_multiplier #(.WIDTH_A(MANT_WIDTH_I+1), .WIDTH_B(MANT_WIDTH_I+1))
            u_radix4_booth_multiplier(
                .operand_a(pipe1_mants_norm_a[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .operand_b(pipe1_mants_norm_b[j*(MANT_WIDTH_I+1) +: (MANT_WIDTH_I+1)]),
                .sum_o(mul_sum[j*MUL_WIDTH +: MUL_WIDTH]),
                .carry_o(mul_carry[j*MUL_WIDTH +: MUL_WIDTH])
            );
            assign mants_norm_c[j*MUL_WIDTH +: MUL_WIDTH] =
                mul_sum[j*MUL_WIDTH +: MUL_WIDTH] + mul_carry[j*MUL_WIDTH +: MUL_WIDTH];
        end
    endgenerate

    logic signed [(N+1)*(EXP_WIDTH+1)-1:0] rg_exp_items;
    logic signed [EXP_WIDTH:0] rg_exp_max;

    generate
        genvar u;
        for(u=0;u<N;u++) begin
            assign rg_exp_items[u*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = pipe1_rg_exp_c[u*(EXP_WIDTH+1) +: (EXP_WIDTH+1)];
        end
    endgenerate
    assign rg_exp_items[N*(EXP_WIDTH+1) +: (EXP_WIDTH+1)] = $signed(pipe1_rg_exp_acc);

    comp_tree #(.N(N+1), .WIDTH(EXP_WIDTH)) u_comp_tree(
        .operands_i(rg_exp_items),
        .result_o(rg_exp_max)
    );

    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    generate
        genvar k;
        for(k=0;k<N;k++) begin
            assign product[k*ALIGN_WIDTH +: ALIGN_WIDTH] = mants_norm_c[k*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
        end
    endgenerate
    assign product[N*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe1_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);

    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z;
        for(z=0;z<N+1;z++) begin: sh
            logic [EXP_WIDTH:0] diff;
            assign diff = $unsigned(rg_exp_max - rg_exp_items[z*(EXP_WIDTH+1) +: (EXP_WIDTH+1)]);
            assign shift_amount[z*SHIFT_WIDTH +: SHIFT_WIDTH] =
                (|diff[SHIFT_WIDTH +: (EXP_WIDTH-SHIFT_WIDTH+1)]) ? ALIGN_WIDTH : diff[SHIFT_WIDTH-1:0];
        end
    endgenerate

    logic [N:0] pipe2_signs;
    logic [(N+1)*ALIGN_WIDTH-1:0] pipe2_product;
    logic [(N+1)*SHIFT_WIDTH-1:0] pipe2_shift_amount;
    logic signed [EXP_WIDTH:0] pipe2_rg_exp_max;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_product,product,clk)
    `FFNR(pipe2_shift_amount,shift_amount,clk)
    `FFNR(pipe2_rg_exp_max,rg_exp_max,clk)

    // ---------------- stage 3: align + condition + tree + add ----------------
    logic [(N+1)*DW-1:0] lane_op;
    logic [4:0] n_neg;

    generate
        genvar s;
        for(s=0;s<N+1;s++) begin: align
            logic [ALIGN_WIDTH-1:0] prod_sh;
            logic [14:0] onec;
            barrel_shifter #(.WIDTH(ALIGN_WIDTH), .SHIFT_WIDTH(SHIFT_WIDTH), .MODE(1'b1))
            u_barrel_shifter(
                .operand_i(pipe2_product[s*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(pipe2_shift_amount[s*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(prod_sh)
            );
            assign onec = ~{1'b0, prod_sh};
            assign lane_op[s*DW +: DW] = pipe2_signs[s]
                ? {{(DW-15){1'b1}}, onec}
                : {{(DW-14){1'b0}}, prod_sh};
        end
    endgenerate

    logic [4:0] n0, n1, n2, n3, n4, n5, n6, n7, n8;
    assign n0 = {4'b0, pipe2_signs[0]} + {4'b0, pipe2_signs[1]};
    assign n1 = {4'b0, pipe2_signs[2]} + {4'b0, pipe2_signs[3]};
    assign n2 = {4'b0, pipe2_signs[4]} + {4'b0, pipe2_signs[5]};
    assign n3 = {4'b0, pipe2_signs[6]} + {4'b0, pipe2_signs[7]};
    assign n4 = n0 + n1;
    assign n5 = n2 + n3;
    assign n6 = n4 + n5;
    assign n7 = {4'b0, pipe2_signs[8]} + {4'b0, pipe2_signs[9]} + {4'b0, pipe2_signs[10]} +
                {4'b0, pipe2_signs[11]} + {4'b0, pipe2_signs[12]} + {4'b0, pipe2_signs[13]} +
                {4'b0, pipe2_signs[14]} + {4'b0, pipe2_signs[15]} + {4'b0, pipe2_signs[16]};
    assign n8 = n6 + n7;
    assign n_neg = n8;

    logic [(N+2)*DW-1:0] tree_in;
    always_comb begin
        tree_in = '0;
        for (int t = 0; t < N+1; t++) tree_in[t*DW +: DW] = lane_op[t*DW +: DW];
        tree_in[(N+1)*DW +: DW] = {{(DW-5){1'b0}}, n_neg};
    end

    logic [DW-1:0] csa_sum, csa_carry, sum_full;
    csa_tree #(.N(N+2), .WIDTH_I(DW), .WIDTH_O(DW)) u_csa_tree(
        .operands_i(tree_in),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );
    assign sum_full = csa_sum + csa_carry;

    // ---------------- pipe3: just the 20-bit sum + the max exponent --------
    logic [DW-1:0] pipe3_sum_full;
    logic signed [EXP_WIDTH:0] pipe3_rg_exp_max;

    `FFNR(pipe3_sum_full,sum_full,clk)
    `FFNR(pipe3_rg_exp_max,pipe2_rg_exp_max,clk)

    // ---------------- stage 4: sign/magnitude + normalization --------------
    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    logic signed [EXP_WIDTH:0] rg_exp_adjust;
    logic [SUM_WIDTH-1:0] sum_norm;
    logic signed [EXP_WIDTH:0] final_rg_exp;

    assign final_sign = pipe3_sum_full[DW-1];
    assign sum_c = final_sign ? (~pipe3_sum_full[SUM_WIDTH-1:0]+1) : pipe3_sum_full[SUM_WIDTH-1:0];

    mantissa_norm #(.WIDTH(SUM_WIDTH), .EXP_WIDTH(EXP_WIDTH), .DECIMAL_POINT(CARRY_WIDTH+2))
    u_mantissa_norm(
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    assign final_rg_exp = pipe3_rg_exp_max + rg_exp_adjust;

    logic pipe4_final_sign;
    logic signed [EXP_WIDTH:0] pipe4_final_rg_exp;
    logic [SUM_WIDTH-1:0] pipe4_sum_norm;

    `FFNR(pipe4_final_sign,final_sign,clk)
    `FFNR(pipe4_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe4_sum_norm,sum_norm,clk)

    // ---------------- stage 5: rounding bits / sticky ----------------------
    logic [MANT_WIDTH_O+2:0] final_mant;

    assign final_mant = {pipe4_sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2],
                         |pipe4_sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0]};

    logic pipe5_final_sign;
    logic signed [EXP_WIDTH:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;

    `FFNR(pipe5_final_sign,pipe4_final_sign,clk)
    `FFNR(pipe5_final_rg_exp,pipe4_final_rg_exp,clk)
    `FFNR(pipe5_final_mant,final_mant,clk)

    // ---------------- stage 6: encoder ------------------------------------
    posit_encoder #(.n(n_o), .es(es_o), .EXP_WIDTH(EXP_WIDTH), .MANT_WIDTH(MANT_WIDTH_O+2))
    u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );

endmodule

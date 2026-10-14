`include "starting_point.v"

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

// design8: design4 + exponent max split across stages (3 levels in stage2,
// 2 levels in stage3) instead of a flat comp_tree in stage 2.
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
    localparam int unsigned EXP_WIDTH    = pdpu_pkg::maximum(EXP_WIDTH_I+1,EXP_WIDTH_O);
    localparam int unsigned MUL_WIDTH    = 2*(MANT_WIDTH_I+1);
    localparam int unsigned CARRY_WIDTH  = pdpu_pkg::clog2(N+1);
    localparam int unsigned SUM_WIDTH    = ALIGN_WIDTH + CARRY_WIDTH;
    localparam int unsigned SHIFT_WIDTH  = pdpu_pkg::clog2(ALIGN_WIDTH+1);
    localparam int unsigned EW           = EXP_WIDTH+1;   // 7
    localparam int unsigned EW_I         = EXP_WIDTH_I+1; // 6

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
    logic signed [N*EW-1:0] pipe1_rg_exp_c;
    logic [N*(MANT_WIDTH_I+1)-1:0] pipe1_mants_norm_a, pipe1_mants_norm_b;
    logic signed [EW-1:0] pipe1_rg_exp_acc;
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

    logic signed [(N+1)*EW-1:0] rg_exp_items;
    generate
        genvar u;
        for(u=0;u<N;u++) begin
            assign rg_exp_items[u*EW +: EW] = pipe1_rg_exp_c[u*EW +: EW];
        end
    endgenerate
    assign rg_exp_items[N*EW +: EW] = $signed(pipe1_rg_exp_acc);

    logic signed [9*EW-1:0] mx9;
    logic signed [5*EW-1:0] mx5;
    logic signed [3*EW-1:0] mx3;
    max_reduce #(.NV(17), .W(EXP_WIDTH)) r0 (.in_i(rg_exp_items), .out_o(mx9));
    max_reduce #(.NV(9),  .W(EXP_WIDTH)) r1 (.in_i(mx9),           .out_o(mx5));
    max_reduce #(.NV(5),  .W(EXP_WIDTH)) r2 (.in_i(mx5),           .out_o(mx3));

    logic [(N+1)*ALIGN_WIDTH-1:0] product;
    generate
        genvar k;
        for(k=0;k<N;k++) begin
            assign product[k*ALIGN_WIDTH +: ALIGN_WIDTH] = mants_norm_c[k*MUL_WIDTH +: MUL_WIDTH] << (ALIGN_WIDTH-MUL_WIDTH);
        end
    endgenerate
    assign product[N*ALIGN_WIDTH +: ALIGN_WIDTH] = pipe1_mant_norm_acc << (ALIGN_WIDTH-MANT_WIDTH_O-2);

    // ---------------- pipe2 ----------------
    logic [N:0] pipe2_signs;
    logic signed [(N+1)*EW-1:0] pipe2_rg_exp_items;
    logic signed [3*EW-1:0] pipe2_mx3;
    logic [(N+1)*ALIGN_WIDTH-1:0] pipe2_product;

    `FFNR(pipe2_signs,pipe1_signs,clk)
    `FFNR(pipe2_rg_exp_items,rg_exp_items,clk)
    `FFNR(pipe2_mx3,mx3,clk)
    `FFNR(pipe2_product,product,clk)

    // ---------------- stage 3 ----------------
    logic signed [2*EW-1:0] mx2;
    logic signed [EW-1:0]   rg_exp_max;
    max_reduce #(.NV(3), .W(EXP_WIDTH)) r3 (.in_i(pipe2_mx3), .out_o(mx2));
    max_reduce #(.NV(2), .W(EXP_WIDTH)) r4 (.in_i(mx2),       .out_o(rg_exp_max));

    logic [(N+1)*ALIGN_WIDTH-1:0] product_shifted;
    logic [(N+1)*SHIFT_WIDTH-1:0] shift_amount;
    generate
        genvar z,s;
        for(z=0;z<N+1;z++) begin: sh
            logic [EW-1:0] diff;
            assign diff = $unsigned(rg_exp_max - pipe2_rg_exp_items[z*EW +: EW]);
            assign shift_amount[z*SHIFT_WIDTH +: SHIFT_WIDTH] =
                (|diff[SHIFT_WIDTH +: (EW-SHIFT_WIDTH)]) ? ALIGN_WIDTH : diff[SHIFT_WIDTH-1:0];
        end
        for(s=0;s<N+1;s++) begin: align
            barrel_shifter #(.WIDTH(ALIGN_WIDTH), .SHIFT_WIDTH(SHIFT_WIDTH), .MODE(1'b1))
            u_barrel_shifter(
                .operand_i(pipe2_product[s*ALIGN_WIDTH +: ALIGN_WIDTH]),
                .shift_amount(shift_amount[s*SHIFT_WIDTH +: SHIFT_WIDTH]),
                .result_o(product_shifted[s*ALIGN_WIDTH +: ALIGN_WIDTH])
            );
        end
    endgenerate

    logic [(N+1)*(SUM_WIDTH+1)-1:0] mantissa, mantissa_comp;
    generate
        genvar y;
        for(y=0;y<N+1;y++) begin
            assign mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] = product_shifted[y*ALIGN_WIDTH +: ALIGN_WIDTH];
            assign mantissa_comp[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)] =
                pipe2_signs[y] ? (~mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)]+1) : mantissa[y*(SUM_WIDTH+1) +: (SUM_WIDTH+1)];
        end
    endgenerate

    // ---------------- pipe3 ----------------
    logic signed [EW-1:0] pipe3_rg_exp_max;
    logic [(N+1)*(SUM_WIDTH+1)-1:0] pipe3_mantissa_comp;

    `FFNR(pipe3_rg_exp_max,rg_exp_max,clk)
    `FFNR(pipe3_mantissa_comp,mantissa_comp,clk)

    logic [SUM_WIDTH:0] csa_sum, csa_carry;
    csa_tree #(.N(N+1), .WIDTH_I(SUM_WIDTH+1), .WIDTH_O(SUM_WIDTH+1)) u_csa_tree(
        .operands_i(pipe3_mantissa_comp),
        .sum_o(csa_sum),
        .carry_o(csa_carry)
    );

    logic [SUM_WIDTH:0] sum_result;
    assign sum_result = csa_sum + csa_carry;

    // ---------------- pipe4 ----------------
    logic [SUM_WIDTH:0] pipe4_sum_result;
    logic signed [EW-1:0] pipe4_rg_exp_max;

    `FFNR(pipe4_sum_result,sum_result,clk)
    `FFNR(pipe4_rg_exp_max,pipe3_rg_exp_max,clk)

    logic final_sign;
    logic [SUM_WIDTH-1:0] sum_c;
    logic signed [EW-1:0] rg_exp_adjust;
    logic signed [EW-1:0] final_rg_exp;
    logic [SUM_WIDTH-1:0] sum_norm;
    logic [MANT_WIDTH_O+2:0] final_mant;

    assign final_sign = pipe4_sum_result[SUM_WIDTH];
    assign sum_c = final_sign ? (~pipe4_sum_result[SUM_WIDTH-1:0]+1) : pipe4_sum_result[SUM_WIDTH-1:0];

    mantissa_norm #(.WIDTH(SUM_WIDTH), .EXP_WIDTH(EXP_WIDTH), .DECIMAL_POINT(CARRY_WIDTH+2))
    u_mantissa_norm(
        .operand_i(sum_c),
        .exp_adjust(rg_exp_adjust),
        .result_o(sum_norm)
    );
    assign final_rg_exp = pipe4_rg_exp_max + rg_exp_adjust;

    assign final_mant = {sum_norm[SUM_WIDTH-1:SUM_WIDTH-MANT_WIDTH_O-2], |sum_norm[SUM_WIDTH-MANT_WIDTH_O-3:0]};

    // ---------------- pipe5 ----------------
    logic pipe5_final_sign;
    logic signed [EW-1:0] pipe5_final_rg_exp;
    logic [MANT_WIDTH_O+2:0] pipe5_final_mant;

    `FFNR(pipe5_final_sign,final_sign,clk)
    `FFNR(pipe5_final_rg_exp,final_rg_exp,clk)
    `FFNR(pipe5_final_mant,final_mant,clk)

    posit_encoder #(.n(n_o), .es(es_o), .EXP_WIDTH(EXP_WIDTH), .MANT_WIDTH(MANT_WIDTH_O+2))
    u_posit_encoder(
        .sign_i(pipe5_final_sign),
        .rg_exp_i(pipe5_final_rg_exp),
        .mant_norm_i(pipe5_final_mant),
        .result_o(result_o)
    );

endmodule

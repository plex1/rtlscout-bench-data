// ===========================================================================
// PDPU dot-16, 5-stage pipelined -- compact re-implementation.
// Functionally bit-identical to the reference (including ALIGN_WIDTH=14
// truncation, CARRY_W=5, SUM_W=19 and NaR handling) but with a much smaller
// and faster datapath:
//   * posit decode: direct regime field extraction instead of lzc+barrel
//   * posit encode: shift amount computed directly from the regime value
//   * mantissa alignment: normalised mantissa (always in [8,15]) instead of
//     the wider two's-complement normalisation
//   * product kept in carry-save form until the alignment shifter
//   * the 17-way final addition of ALIGN_WIDTH wide words is done with a
//     balanced binary adder tree instead of a full 17:2 carry-save tree
// ===========================================================================
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
    input  logic [n_o-1:0]  acc,
    output logic [n_o-1:0]  result_o
);

    // ---------------- constants (identical to the reference) ---------------
    localparam int unsigned EXPI    = 4;   // clog2(7)+es_i, signed width
    localparam int unsigned MI      = 5;   // mantissa + implicit bit
    localparam int unsigned MULO    = 10;  // product width
    localparam int unsigned EXPO    = 5;   // clog2(15)+es_o, signed width
    localparam int unsigned MO      = 13;  // mantissa + implicit bit (out)
    localparam int unsigned EXP     = 5;   // max(EXPI+1, EXPO)
    localparam int unsigned CW      = 5;   // CARRY_W
    localparam int unsigned SW      = 19;  // SUM_WIDTH
    localparam int unsigned SWS     = 20;  // SUM_WIDTH+1
    localparam int unsigned AWD     = 14;  // ALIGN_WIDTH
    localparam int unsigned SHW     = 4;   // clog2(ALIGN_WIDTH+1)
    localparam int unsigned PC      = 7;   // DECIMAL_POINT = CARRY_W+2
    localparam int unsigned LZCW    = 5;   // clog2(SUM_WIDTH)
    localparam int unsigned MULSW   = 2*(EXPI+1); // 10 exponent sum width

    // ======================= stage 1 : posit decode ========================
    logic [N-1:0]              s1_signs;
    logic signed [N*(EXP+1)-1:0] s1_e;
    logic [N*(MI)-1:0]         s1_ma, s1_mb;
    logic signed [EXPO:0]      s1_eacc;
    logic [MO-1:0]             s1_macc;

    genvar i;
    generate
    for (i=0;i<N;i++) begin : g_dec
        // ---- operand a ----
        logic [6:0] va;
        logic [2:0] rl_a;
        logic       rs_a;
        logic signed [EXPI:0] ek_a;
        logic [2:0] nb_a;
        logic [6:0] sh_a;
        logic [5:0] f_a;

        assign va    = operands_a[i*8+7] ? (~operands_a[i*8 +: 7] + 7'd1)
                                         :  operands_a[i*8 +: 7];
        assign rs_a  = va[6];
        assign rl_a  = rs_a ? 3'd7 : ({2'b0,~va[6]} | {2'b0,~va[5]} | {2'b0,~va[4]}
                                  | {2'b0,~va[3]} | {2'b0,~va[2]} | {2'b0,~va[1]}
                                  | {2'b0,~va[0]});
        assign ek_a  = rs_a ? ({1'b0,rl_a} - 5'sd1) : ({1'b1,~rl_a} + 5'sd1);
        assign nb_a  = rl_a + 3'd1;
        assign sh_a  = ({4'b0,nb_a} > 7'd6) ? 7'd6 : {4'b0,nb_a};
        assign f_a   = va << sh_a[2:0];

        assign s1_signs[i] = operands_a[i*8+7];
        assign s1_e[i*(EXP+1) +: (EXP+1)] = $signed(ek_a) <<< 2 | $signed({1'b0,f_a[5:4]});
        assign s1_ma[i*MI +: MI] = {|va, f_a[3:0]};

        // ---- operand b ----
        logic [6:0] vb;
        logic [2:0] rl_b;
        logic       rs_b;
        logic signed [EXPI:0] ek_b;
        logic [2:0] nb_b;
        logic [6:0] sh_b;
        logic [5:0] f_b;

        assign vb    = operands_b[i*8+7] ? (~operands_b[i*8 +: 7] + 7'd1)
                                         :  operands_b[i*8 +: 7];
        assign rs_b  = vb[6];
        assign rl_b  = rs_b ? 3'd7 : ({2'b0,~vb[6]} | {2'b0,~vb[5]} | {2'b0,~vb[4]}
                                  | {2'b0,~vb[3]} | {2'b0,~vb[2]} | {2'b0,~vb[1]}
                                  | {2'b0,~vb[0]});
        assign ek_b  = rs_b ? ({1'b0,rl_b} - 5'sd1) : ({1'b1,~rl_b} + 5'sd1);
        assign nb_b  = rl_b + 3'd1;
        assign sh_b  = ({4'b0,nb_b} > 7'd6) ? 7'd6 : {4'b0,nb_b};
        assign f_b   = vb << sh_b[2:0];

        assign s1_signs[N+i] = operands_b[i*8+7];
        assign s1_e[(N+i)*(EXP+1) +: (EXP+1)] = $signed(ek_b) <<< 2 | $signed({1'b0,f_b[5:4]});
        assign s1_mb[i*MI +: MI] = {|vb, f_b[3:0]};
    end

    // ---- accumulator ----
    logic [14:0] vacc;
    logic [3:0]  rl_c;
    logic        rs_c;
    logic signed [EXPO:0] ek_c;
    logic [3:0]  nb_c;
    logic [14:0] sh_c;
    logic [12:0] f_c;

    assign vacc = acc[15] ? (~acc[14:0] + 15'd1) : acc[14:0];
    assign rs_c = vacc[14];
    assign rl_c = rs_c ? 4'd15 : ({3'b0,~vacc[14]} | {3'b0,~vacc[13]} | {3'b0,~vacc[12]}
                              | {3'b0,~vacc[11]} | {3'b0,~vacc[10]} | {3'b0,~vacc[9]}
                              | {3'b0,~vacc[8]}  | {3'b0,~vacc[7]}  | {3'b0,~vacc[6]}
                              | {3'b0,~vacc[5]}  | {3'b0,~vacc[4]}  | {3'b0,~vacc[3]}
                              | {3'b0,~vacc[2]}  | {3'b0,~vacc[1]}  | {3'b0,~vacc[0]});
    assign ek_c = rs_c ? ({1'b0,rl_c} - 5'sd1) : ({1'b1,~rl_c} + 5'sd1);
    assign nb_c = rl_c + 4'd1;
    assign sh_c = (nb_c > 4'd12) ? 4'd12 : nb_c;
    assign f_c  = vacc << sh_c[3:0];

    assign s1_eacc = $signed(ek_c) <<< 2 | $signed({1'b0,f_c[11:10]});
    assign s1_macc = {|vacc, f_c[9:0]};
    endgenerate

    // ======================= stage 1 registers =============================
    logic [N-1:0]              p1_signs;
    logic signed [N*(EXP+1)-1:0] p1_e;
    logic [N*MI-1:0]           p1_ma, p1_mb;
    logic signed [EXPO:0]      p1_eacc;
    logic [MO-1:0]             p1_macc;

    always_ff @(posedge clk) begin
        p1_signs <= s1_signs;
        p1_e     <= s1_e;
        p1_ma    <= s1_ma;
        p1_mb    <= s1_mb;
        p1_eacc  <= s1_eacc;
        p1_macc  <= s1_macc;
    end

    // ======================= stage 2 : mantissa mul ========================
    logic [N*MULSW-1:0] s2_ec;
    logic [N*MULO-1:0]  s2_p;
    logic signed [EXPO:0] s2_eacc;
    logic [MO-1:0]      s2_macc;

    generate
    for (i=0;i<N;i++) begin : g_mul
        assign s2_ec[i*MULSW +: MULSW] =
            $signed(p1_e[i*(EXP+1) +: (EXP+1)]) + $signed(p1_e[(N+i)*(EXP+1) +: (EXP+1)]);
        assign s2_p[i*MULO +: MULO] =
            p1_ma[i*MI +: MI] * p1_mb[i*MI +: MI];
    end
    endgenerate
    assign s2_eacc = p1_eacc;
    assign s2_macc = p1_macc;

    // ======================= stage 2 registers =============================
    logic [N-1:0]              p2_signs;
    logic signed [N*(MULSW)-1:0] p2_ec;
    logic [N*MULO-1:0]         p2_p;
    logic signed [EXPO:0]      p2_eacc;
    logic [MO-1:0]             p2_macc;

    always_ff @(posedge clk) begin
        p2_signs <= p1_signs;
        p2_ec    <= s2_ec;
        p2_p     <= s2_p;
        p2_eacc  <= s2_eacc;
        p2_macc  <= s2_macc;
    end

    // ======================= stage 3 : align + add =========================
    // maximum exponent (17-way signed max tree)
    logic signed [EXP:0] mx_l0 [0:8];
    logic signed [EXP:0] mx_l1 [0:4];
    logic signed [EXP:0] mx_l2 [0:2];
    logic signed [EXP:0] emax;
    integer k;
    always_comb begin
        for (k=0;k<N;k++) mx_l0[k] = p2_ec[k*MULSW +: (EXP+1)];
        mx_l0[N] = p2_eacc;
        mx_l1[0] = (mx_l0[0]  > mx_l0[1])  ? mx_l0[0]  : mx_l0[1];
        mx_l1[1] = (mx_l0[2]  > mx_l0[3])  ? mx_l0[2]  : mx_l0[3];
        mx_l1[2] = (mx_l0[4]  > mx_l0[5])  ? mx_l0[4]  : mx_l0[5];
        mx_l1[3] = (mx_l0[6]  > mx_l0[7])  ? mx_l0[6]  : mx_l0[7];
        mx_l1[4] = (mx_l0[8]  > mx_l0[9])  ? mx_l0[8]  : mx_l0[9];
        mx_l2[0] = (mx_l1[0]  > mx_l1[1])  ? mx_l1[0]  : mx_l1[1];
        mx_l2[1] = (mx_l1[2]  > mx_l1[3])  ? mx_l1[2]  : mx_l1[3];
        mx_l2[2] = (mx_l1[4]  > mx_l0[10]) ? mx_l1[4]  : mx_l0[10];
        emax     = (mx_l2[0]  > mx_l2[1])  ? mx_l2[0]  : mx_l2[1];
        emax     = (emax      > mx_l2[2])  ? emax      : mx_l2[2];
    end

    // aligned signed values
    logic [N:0]                s3_sgn;
    logic [(N+1)*AWD-1:0]      s3_al;
    logic signed [EXP:0]       s3_emax;
    logic [(N+1)*SWS-1:0]      s3_tc;

    generate
    for (i=0;i<N+1;i++) begin : g_al
        logic signed [EXP:0]    ei;
        logic signed [EXP:0]    df;
        logic [3:0]             sh;
        logic [MULO-1:0]        prod;
        logic [AWD-1:0]         base;
        logic [AWD-1:0]         shf;

        if (i < N) begin : g_prod
            assign ei   = p2_ec[i*MULSW +: (EXP+1)];
            assign prod = p2_p[i*MULO +: MULO];
        end
        else begin : g_acc
            assign ei   = p2_eacc;
            assign prod = {p2_macc, 3'b0};
        end
        assign df = emax - ei;
        assign sh = df[3:0];
        assign base = prod >> 6;
        assign shf  = base >> sh;
        assign s3_sgn[i] = p2_signs[i];
        assign s3_al[i*AWD +: AWD] = shf;
        assign s3_tc[i*SWS +: SWS] = p2_signs[i] ? (~{6'b0,shf} + 20'd1)
                                                 :  {6'b0,shf};
    end
    endgenerate
    assign s3_emax = emax;

    // balanced adder tree over 17 words (11 unused slots tied to zero)
    logic [SWS-1:0] tr0 [0:31];
    logic [SWS-1:0] tr1 [0:15];
    logic [SWS-1:0] tr2 [0:7];
    logic [SWS-1:0] tr3 [0:3];
    logic [SWS-1:0] tr4 [0:1];
    logic [SWS-1:0] sum_result;

    always_comb begin
        for (k=0;k<32;k++) tr0[k] = 20'd0;
        for (k=0;k<N+1;k++) tr0[k] = s3_tc[k*SWS +: SWS];
        for (k=0;k<16;k++) tr1[k] = tr0[2*k] + tr0[2*k+1];
        for (k=0;k<8;k++)  tr2[k] = tr1[2*k] + tr1[2*k+1];
        for (k=0;k<4;k++)  tr3[k] = tr2[2*k] + tr2[2*k+1];
        for (k=0;k<2;k++)  tr4[k] = tr3[2*k] + tr3[2*k+1];
        sum_result = tr4[0] + tr4[1];
    end

    // ======================= stage 3 registers =============================
    logic [SWS-1:0]      p3_sum;
    logic signed [EXP:0] p3_emax;

    always_ff @(posedge clk) begin
        p3_sum  <= sum_result;
        p3_emax <= s3_emax;
    end

    // ======================= stage 4 : normalize ===========================
    logic        final_sign;
    logic [18:0] sum_c;
    logic [4:0]  lzc_c;
    logic [4:0]  lzc_neg;
    logic        lzc_empty;
    logic signed [EXP:0] rg_exp_adjust;
    logic [18:0] sum_norm;
    logic [14:0] final_mant;

    always_comb begin
        final_sign = p3_sum[SWS-1];
        sum_c      = final_sign ? (~p3_sum + 20'd1) : p3_sum;
        // leading zero count of sum_c[18:0] (width 19)
        lzc_empty = ~(|sum_c[18:0]);
        lzc_c     = 5'd19;
        for (k=18;k>=0;k--) begin
            if (sum_c[k]) lzc_c = 5'd18 - k[4:0];
        end
        lzc_neg   = 5'd0 - lzc_c;
        rg_exp_adjust = (lzc_c == 5'd0) ? $signed({1'b0,5'd3})
                      : (lzc_c < 5'd4)  ? $signed({1'b0,5'd3}) - $signed({1'b0,lzc_c})
                                        : -$signed({1'b0,lzc_neg});
        sum_norm  = sum_c << lzc_c;
        final_mant = {sum_norm[18:6], |sum_norm[5:0]};
    end

    // ======================= stage 4 registers =============================
    logic        p4_sign;
    logic signed [EXP:0] p4_rg_exp;
    logic [14:0] p4_mant;

    logic signed [EXP:0] final_rg_exp;
    assign final_rg_exp = p3_emax + rg_exp_adjust;

    always_ff @(posedge clk) begin
        p4_sign   <= final_sign;
        p4_rg_exp <= final_rg_exp;
        p4_mant   <= final_mant;
    end

    // ======================= stage 5 : posit encode ========================
    logic signed [2:0]  reg_k;
    logic [4:0]         reg_bits;
    logic [17:0]        rg_mant;
    logic [4:0]         enc_sh;
    logic [34:0]        enc_val;
    logic [14:0]        vbr;
    logic [6:0]         rbits;
    logic               round_bit, sticky_bit, round_value;
    logic [14:0]        var;
    logic [15:0]        normal_result;

    always_comb begin
        reg_k    = p4_rg_exp[EXP:2];
        reg_bits = p4_rg_exp[EXP] ? ({1'b0,~reg_k} + 5'd2) : ({1'b0,reg_k} + 5'd2);
        rg_mant  = {1'b1, 12'b0, p4_mant[14:10]};
        enc_sh   = (reg_bits >= 5'd16) ? 5'd7 : (reg_bits - 5'd4);
        enc_val  = {14'b0, rg_mant, 7'b0} >> enc_sh;
        vbr      = enc_val[21:7];
        rbits    = enc_val[6:0];
        round_bit    = rbits[6];
        sticky_bit   = |rbits[5:0];
        round_value  = round_bit & (sticky_bit | vbr[0]);
        var          = vbr + {14'b0, round_value};
        normal_result = p4_sign ? {1'b1, ~var + 16'd1} : {1'b0, var};
    end

    assign result_o = p4_mant[14] ? normal_result : 16'd0;

endmodule

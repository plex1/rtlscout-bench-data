// 5-stage pipelined 16-lane posit(8,2) dot-product-accumulate unit.
// Bit-exact rewrite of the reference pdpu_top_pipelined datapath.
//
//   L1: decode 16 a/b lanes + acc, exponent add, mantissa multiply, pack regs
//   L2: max-exponent tree, final product add, exponent differences, pack regs
//   L3: barrel-align, two's complement, pack regs
//   L4: CSA tree + final add, sign/magnitude, leading-zero count, pack regs
//   L5: posit encode -> result_o

// ---------------------------------------------------------------------------
// Leading-zero counter: cnt = number of leading zeros of in, capped at W.
//   in != 0 -> index (from the MSB) of the highest set bit
//   in == 0 -> W
module lzc_lead #(
    parameter int W  = 19,
    parameter int CW = 5
)(
    input  logic [W-1:0]  in,
    output logic [CW-1:0] cnt
);
    localparam int LV = $clog2(W);        // tree levels
    localparam int NP = 1 << LV;          // padded to a power of two
    localparam int IW = LV + 1;           // index width inside the tree

    // leaf k corresponds to in[W-1-k]
    logic [(LV+1)*NP-1:0]    sel_f;       // sel_f[l*NP + k]
    logic [(LV+1)*NP*IW-1:0] idx_f;       // idx_f[(l*NP+k)*IW +: IW]

    genvar l, k;
    generate
        for (k = 0; k < NP; k++) begin : g_leaf
            if (k < W) begin : g_real
                assign sel_f[0*NP + k] = in[W-1-k];
                assign idx_f[(0*NP + k)*IW +: IW] = ~in[W-1-k];
            end else begin : g_pad
                assign sel_f[0*NP + k] = 1'b0;
                assign idx_f[(0*NP + k)*IW +: IW] = {IW{1'b0}};
            end
        end
        for (l = 1; l <= LV; l++) begin : g_lvl
            for (k = 0; k < (NP >> l); k++) begin : g_node
                assign sel_f[l*NP + k] =
                    sel_f[(l-1)*NP + 2*k] | sel_f[(l-1)*NP + 2*k + 1];
                assign idx_f[(l*NP + k)*IW +: IW] =
                    sel_f[(l-1)*NP + 2*k]
                    ? {1'b0, idx_f[((l-1)*NP + 2*k)*IW +: IW-1]}
                    : {1'b1, idx_f[((l-1)*NP + 2*k + 1)*IW +: IW-1]};
            end
        end
    endgenerate

    assign cnt = sel_f[LV*NP] ? idx_f[(LV*NP)*IW +: CW] : W[CW-1:0];
endmodule

// ---------------------------------------------------------------------------
module pdpu_top_pipelined #(
    parameter int unsigned N = 16,
    parameter int unsigned n_i = 8,
    parameter int unsigned es_i = 2,
    parameter int unsigned n_o = 16,
    parameter int unsigned es_o = 2,
    parameter int unsigned ALIGN_WIDTH = 14
)(
    input  logic              clk,
    input  logic [N*n_i-1:0]  operands_a,
    input  logic [N*n_i-1:0]  operands_b,
    input  logic [n_o-1:0]    acc,
    output logic [n_o-1:0]    result_o
);
    // ---------------------------------------------------------------- widths
    localparam int unsigned EW_I = $clog2(n_i-1) + es_i;   // 3
    localparam int unsigned EW_O = $clog2(n_o-1) + es_o;   // 5
    localparam int unsigned EW   = (EW_I+1 > EW_O) ? EW_I+1 : EW_O;  // 5
    localparam int unsigned AW   = EW + 1;                 // 6 : signed rg_exp
    localparam int unsigned MW   = n_i - es_i - 2;         // 4 : product mantissa
    localparam int unsigned PW   = 2*MW;                   // 8
    localparam int unsigned PA   = 8;                      // encoded lane width
    localparam int unsigned TB   = ALIGN_WIDTH + PA + 1;   // 23 : packed lane reg
    localparam int unsigned CW   = $clog2(N+1);            // 5
    localparam int unsigned SW   = ALIGN_WIDTH + CW;       // 19 : SUM_WIDTH
    localparam int unsigned AN   = 2*N + 1;                // 33
    localparam int unsigned TW   = SW + 1;                 // 20
    localparam int unsigned OW   = ALIGN_WIDTH + 3;        // 17 : acc window
    localparam int unsigned MW5  = 12;                     // 11 mantissa + sticky
    localparam int unsigned REG5 = AW + 1 + MW5;           // 18 : L5 operand reg

    // ================================================================ L1
    // posit(8,2) decode of the a/b lanes
    //   v     = (s ? ~x+1 : x)[6:0]
    //   c     = v[6] ? ~v : v               7-bit lzc input (0 only when v==0)
    //   k     = {!v[6], v[6] ? ~lz : lz}    leading-zero count of c (0..7)
    //   e     = v[6] ? 4*k : -4*k+4         real exponent (6-bit signed)
    //   D[7]  = 1                           implicit (non-zero) bit
    //   D[4:0]= v[6] ? v[4:0]<<1 : v[5:1]   normalised significand
    //   D[7] |= v[6] & v[1]                 rounding bit
    //   D[8]  = v[6] & v[0]                 sticky bit
    logic [N*AW-1:0] dec_e;
    logic [N*PA-1:0] dec_d;
    logic [N-1:0]    dec_nz, dec_nar;

    genvar g;
    generate
        for (g = 0; g < N; g++) begin : g_dec
            logic [7:0]  xa, xb;
            logic        sa, sb;
            logic [6:0]  va, vb;
            logic [6:0]  ca, cb;
            logic [2:0]  la, lb;
            logic [3:0]  ka, kb;
            logic [3:0]  ea, eb;

            assign xa = operands_a[g*n_i +: n_i];
            assign xb = operands_b[g*n_i +: n_i];
            assign sa = xa[7];
            assign sb = xb[7];
            assign va = sa ? (~xa[6:0] + 7'd1) : xa[6:0];
            assign vb = sb ? (~xb[6:0] + 7'd1) : xb[6:0];
            assign ca = va[6] ? ~va : va;
            assign cb = vb[6] ? ~vb : vb;

            lzc_lead #(.W(7), .CW(3)) u_lza (.in(ca), .cnt(la));
            lzc_lead #(.W(7), .CW(3)) u_lzb (.in(cb), .cnt(lb));

            assign ka = {~va[6], la};
            assign kb = {~vb[6], lb};
            assign ea = va[6] ? ({2'b00, ka} << 2) : (4'd4 - {2'b00, ka});
            assign eb = vb[6] ? ({2'b00, kb} << 2) : (4'd4 - {2'b00, kb});

            assign dec_e[g*AW +: AW] = {{(AW-6){ea[3]}}, ea, 2'b00};
            assign dec_d[g*PA +: PA] = {5'b0,
                                        va[6] & va[0],
                                        va[6] & va[1],
                                        1'b1,
                                        va[6] ? va[4:0] : va[5:1]};
            assign dec_nz[g]  = |xa;
            assign dec_nar[g] = ~((|xa) & (|xb));
        end
    endgenerate

    // posit(16,2) decode of the accumulator
    //   v = (s ? ~x+1 : x)[14:0]
    //   l = v[14] ? 14-lzc(~v) : lzc(v)      regime run length
    //   E = (v[14] ? -l : l) << 2
    //   M = 1<<16 | (v[14] ? (~v<<(l+1)) : (v<<(l+1))) & 0x1FFFF
    //   acc_exp  = (E << 1) | (v[14] ? 0 : 2)          (Q1 exponent)
    //   acc_mant = 1<<16 | (M << (15-l))               (Q1 mantissa, 1<<16 implicit)
    logic [AW-1:0]     acc_e;
    logic [MW5-1:0]    acc_m;
    logic              acc_nar;
    generate
        if (1) begin : g_acc
            logic [15:0]       ax;
            logic              as;
            logic [14:0]       av, ac;
            logic [3:0]        al, rl, lpc;
            logic [14:0]       t0, t1, t2;
            logic [18:0]       M, sh;
            logic signed [5:0] E;
            logic [18:0]       mant_l1;

            assign ax = acc;
            assign as = ax[15];
            assign av = as ? (~ax[14:0] + 15'd1) : ax[14:0];
            assign ac = av[14] ? ~av : av;

            lzc_lead #(.W(15), .CW(4)) u_lzacc (.in(ac), .cnt(al));

            assign rl  = av[14] ? (4'd14 - al) : al;
            assign lpc = (rl >= 4'd14) ? 4'd15 : (rl + 4'd1);
            assign t0  = {15'd0, (av[14] & ~as)};
            assign t1  = av[14] ? ~av : av;
            assign t2  = t1 << lpc;
            assign M   = {1'b1, 2'b00, t2};
            assign sh  = M << 15;
            assign mant_l1 = sh >> (4'd16 - rl);
            assign E   = av[14] ? {2'b00, rl} : (6'd0 - {2'b00, rl});
            assign acc_e = {{(AW-8){E[5]}}, E, t0[0], 1'b0};
            assign acc_m = mant_l1[16 +: MW5];
            assign acc_nar = ~(|acc);
        end
    endgenerate

    // ---- L1 register pack ------------------------------------------------
    logic [N*TB-1:0] l1_pk;
    logic [N*TB-1:0] l2_pk;
    logic [N*TB-1:0] l3_pk;

    generate
        for (g = 0; g < N; g++) begin : g_pack
            logic [AW-1:0]    eprod;
            logic [PA-1:0]    da, db;
            logic [PA+1:0]    mul;
            logic [PA+AW+1:0] item;
            assign eprod = dec_e[g*AW +: AW] + dec_e[g*AW +: AW];
            assign da = dec_d[g*PA +: PA];
            assign db = dec_d[g*PA +: PA];
            assign mul = da[7:0] * db[7:0];
            assign item = {mul[9], mul[7:4], da[8] | db[8], eprod, 1'b0};
            assign l1_pk[g*TB +: TB] = {{(TB-PA-AW-1){1'b0}}, dec_nz[g], item};
        end
    endgenerate

    always_ff @(posedge clk) begin
        l1_pa   <= dec_d;
        l1_pb   <= dec_d;
        l1_ea   <= dec_e;
        l1_eb   <= dec_e;
        l1_nz   <= dec_nz;
        l1_nar  <= dec_nar;
        l1_macc <= acc_m;
    end

    // ================================================================ L2
    logic [AW-1:0] l2_eacc;
    logic [N*AW-1:0] l2_ea, l2_eb;
    logic [MW5-1:0] l2_macc;
    logic [N-1:0]   l2_nar;

    always_ff @(posedge clk) begin
        l2_ea   <= l1_ea;
        l2_eb   <= l1_eb;
        l2_eacc <= acc_e;
        l2_macc <= acc_m;
        l2_nar  <= l1_nar;
    end

    // ---- max exponent tree (33 values) -----------------------------------
    localparam int unsigned TREE_N = AN;
    logic [TREE_N*AW-1:0]    tree_val;
    logic [(TREE_N+1)*AW-1:0] tree_rd;

    generate
        for (g = 0; g < N; g++) begin : g_tree_leaf
            assign tree_val[g*AW +: AW] = l2_ea[g*AW +: AW] + l2_eb[g*AW +: AW];
        end
    endgenerate
    assign tree_val[N*AW +: AW] = l2_eacc;
    assign tree_val[(N+1)*AW +: AW] = {AW{1'b1}};   // -inf sentinel

    generate
        for (g = 0; g < TREE_N; g++) begin : g_tree_l0
            assign tree_rd[g*AW +: AW] = tree_val[g*AW +: AW];
        end
    endgenerate

    // radix-2 compare tree, ceil(N/2) -> ... -> 1
    localparam int unsigned LVL0 = (TREE_N + 1) / 2;    // 17
    logic [(LVL0+2)*AW-1:0] lvl_a, lvl_b;
    generate
        for (g = 0; g < LVL0; g++) begin : g_lvl_a
            if (2*g + 1 < TREE_N) begin : g_two
                assign lvl_a[g*AW +: AW] =
                    ($signed(tree_rd[(2*g)*AW +: AW]) > $signed(tree_rd[(2*g+1)*AW +: AW]))
                    ? tree_rd[(2*g)*AW +: AW] : tree_rd[(2*g+1)*AW +: AW];
            end else begin : g_one
                assign lvl_a[g*AW +: AW] = tree_rd[(2*g)*AW +: AW];
            end
        end
    endgenerate
    assign lvl_a[LVL0*AW +: AW] = {AW{1'b1}};
    assign lvl_a[(LVL0+1)*AW +: AW] = {AW{1'b1}};

    localparam int unsigned LVL1 = (LVL0 + 1) / 2;      // 9
    generate
        for (g = 0; g < LVL1; g++) begin : g_lvl_b
            assign lvl_b[g*AW +: AW] =
                ($signed(lvl_a[(2*g)*AW +: AW]) > $signed(lvl_a[(2*g+1)*AW +: AW]))
                ? lvl_a[(2*g)*AW +: AW] : lvl_a[(2*g+1)*AW +: AW];
        end
    endgenerate
    assign lvl_b[LVL1*AW +: AW] = {AW{1'b1}};
    assign lvl_b[(LVL1+1)*AW +: AW] = {AW{1'b1}};

    localparam int unsigned LVL2 = (LVL1 + 1) / 2;      // 5
    logic [(LVL2+2)*AW-1:0] lvl_c;
    generate
        for (g = 0; g < LVL2; g++) begin : g_lvl_c
            assign lvl_c[g*AW +: AW] =
                ($signed(lvl_b[(2*g)*AW +: AW]) > $signed(lvl_b[(2*g+1)*AW +: AW]))
                ? lvl_b[(2*g)*AW +: AW] : lvl_b[(2*g+1)*AW +: AW];
        end
    endgenerate
    assign lvl_c[LVL2*AW +: AW] = {AW{1'b1}};
    assign lvl_c[(LVL2+1)*AW +: AW] = {AW{1'b1}};

    localparam int unsigned LVL3 = (LVL2 + 1) / 2;      // 3
    logic [(LVL3+2)*AW-1:0] lvl_d;
    generate
        for (g = 0; g < LVL3; g++) begin : g_lvl_d
            assign lvl_d[g*AW +: AW] =
                ($signed(lvl_c[(2*g)*AW +: AW]) > $signed(lvl_c[(2*g+1)*AW +: AW]))
                ? lvl_c[(2*g)*AW +: AW] : lvl_c[(2*g+1)*AW +: AW];
        end
    endgenerate
    assign lvl_d[LVL3*AW +: AW] = {AW{1'b1}};
    assign lvl_d[(LVL3+1)*AW +: AW] = {AW{1'b1}};

    logic [AW-1:0] l2_emax;
    assign l2_emax =
        ($signed(lvl_d[0*AW +: AW]) > $signed(lvl_d[1*AW +: AW]))
        ? lvl_d[0*AW +: AW] : lvl_d[1*AW +: AW];

    // ================================================================ L3
    logic [AW-1:0]   l3_emax;
    logic [N*TB-1:0] l3_s;
    logic [AN*SW-1:0] l3_sh;

    always_ff @(posedge clk) begin
        l3_emax <= l2_emax;
        l3_pk   <= l2_pk;
    end

    // ---- align -----------------------------------------------------------
    logic [AN*SW-1:0] aligned;
    generate
        for (g = 0; g < AN; g++) begin : g_align
            if (g < N) begin : g_lane
                logic [AW-1:0] ed;
                logic [3:0]    sh4;
                logic [15:0]   x;
                logic [17:0]   prod;
                logic [17:0]   prod_shifted;
                assign ed  = l2_emax - l2_ea[g*AW +: AW] - l2_eb[g*AW +: AW];
                assign sh4 = ed[3:0];
                assign prod = l2_pk[g*TB + 1 + AW +: 18];
                assign prod_shifted = (ed[AW-1:4] != 0) ? 18'd0 : (prod >> sh4);
                assign aligned[g*SW +: SW] = {prod_shifted, 1'b0};
            end else begin : g_acc_lane
                logic [AW-1:0] ed;
                logic [3:0]    sh4;
                logic [17:0]   prod;
                logic [17:0]   prod_shifted;
                assign ed  = l2_emax - l2_eacc;
                assign sh4 = ed[3:0];
                assign prod = {1'b0, l2_macc, 5'b0};
                assign prod_shifted = (ed[AW-1:4] != 0) ? 18'd0 : (prod >> sh4);
                assign aligned[g*SW +: SW] = {prod_shifted, 1'b0};
            end
        end
    endgenerate

    logic [AN*TW-1:0] l4_op;
    generate
        for (g = 0; g < AN; g++) begin : g_tc
            assign l4_op[g*TW +: TW] = {2'b00, aligned[g*SW +: SW]};
        end
    endgenerate

    // ================================================================ L4
    logic [AN*TW-1:0] l4_reg;
    logic [AW-1:0]    l4_emax;
    logic [N-1:0]     l4_nar;

    always_ff @(posedge clk) begin
        l4_reg  <= l4_op;
        l4_emax <= l3_emax;
        l4_nar  <= l2_nar;
    end

    // ---- CSA tree (33 x 20) ---------------------------------------------
    // level 0: 33 -> 22
    localparam int unsigned L1N = 22;
    logic [L1N*TW-1:0] t1;
    generate
        for (g = 0; g < 11; g++) begin : g_csa1
            logic [TW-1:0] a, b, c, s, cy;
            assign a = l4_reg[(3*g+0)*TW +: TW];
            assign b = l4_reg[(3*g+1)*TW +: TW];
            assign c = l4_reg[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t1[(2*g+0)*TW +: TW] = s;
            assign t1[(2*g+1)*TW +: TW] = cy;
        end
    endgenerate

    // level 1: 22 -> 15
    localparam int unsigned L2N = 15;
    logic [L2N*TW-1:0] t2;
    generate
        for (g = 0; g < 7; g++) begin : g_csa2
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t1[(3*g+0)*TW +: TW];
            assign b = t1[(3*g+1)*TW +: TW];
            assign c = t1[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t2[(2*g+0)*TW +: TW] = s;
            assign t2[(2*g+1)*TW +: TW] = cy;
        end
        assign t2[14*TW +: TW] = t1[21*TW +: TW];
    endgenerate

    // level 2: 15 -> 10
    localparam int unsigned L3N = 10;
    logic [L3N*TW-1:0] t3;
    generate
        for (g = 0; g < 5; g++) begin : g_csa3
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t2[(3*g+0)*TW +: TW];
            assign b = t2[(3*g+1)*TW +: TW];
            assign c = t2[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t3[(2*g+0)*TW +: TW] = s;
            assign t3[(2*g+1)*TW +: TW] = cy;
        end
    endgenerate

    // level 3: 10 -> 7
    localparam int unsigned L4N = 7;
    logic [L4N*TW-1:0] t4;
    generate
        for (g = 0; g < 3; g++) begin : g_csa4
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t3[(3*g+0)*TW +: TW];
            assign b = t3[(3*g+1)*TW +: TW];
            assign c = t3[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t4[(2*g+0)*TW +: TW] = s;
            assign t4[(2*g+1)*TW +: TW] = cy;
        end
        assign t4[6*TW +: TW] = t3[9*TW +: TW];
    endgenerate

    // level 4: 7 -> 5
    localparam int unsigned L5N = 5;
    logic [L5N*TW-1:0] t5;
    generate
        for (g = 0; g < 2; g++) begin : g_csa5
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t4[(3*g+0)*TW +: TW];
            assign b = t4[(3*g+1)*TW +: TW];
            assign c = t4[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t5[(2*g+0)*TW +: TW] = s;
            assign t5[(2*g+1)*TW +: TW] = cy;
        end
        assign t5[4*TW +: TW] = t4[6*TW +: TW];
    endgenerate

    // level 5: 5 -> 4
    localparam int unsigned L6N = 4;
    logic [L6N*TW-1:0] t6;
    generate
        for (g = 0; g < 1; g++) begin : g_csa6
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t5[(3*g+0)*TW +: TW];
            assign b = t5[(3*g+1)*TW +: TW];
            assign c = t5[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t6[(2*g+0)*TW +: TW] = s;
            assign t6[(2*g+1)*TW +: TW] = cy;
        end
        assign t6[2*TW +: TW] = t5[3*TW +: TW];
        assign t6[3*TW +: TW] = t5[4*TW +: TW];
    endgenerate

    // level 6: 4 -> 3
    localparam int unsigned L7N = 3;
    logic [L7N*TW-1:0] t7;
    generate
        for (g = 0; g < 1; g++) begin : g_csa7
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t6[(3*g+0)*TW +: TW];
            assign b = t6[(3*g+1)*TW +: TW];
            assign c = t6[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t7[(2*g+0)*TW +: TW] = s;
            assign t7[(2*g+1)*TW +: TW] = cy;
        end
        assign t7[2*TW +: TW] = t6[3*TW +: TW];
    endgenerate

    // level 7: 3 -> 2
    logic [2*TW-1:0] t8;
    generate
        for (g = 0; g < 1; g++) begin : g_csa8
            logic [TW-1:0] a, b, c, s, cy;
            assign a = t7[(3*g+0)*TW +: TW];
            assign b = t7[(3*g+1)*TW +: TW];
            assign c = t7[(3*g+2)*TW +: TW];
            assign s  = a ^ b ^ c;
            assign cy = ((a ^ b) & c) | (a & b);
            assign t8[(2*g+0)*TW +: TW] = s;
            assign t8[(2*g+1)*TW +: TW] = cy;
        end
    endgenerate

    // ---- final add ------------------------------------------------------
    logic [TW-1:0] sum_result;
    assign sum_result = t8[0 +: TW] + t8[TW +: TW];

    logic        final_sign;
    logic [SW-1:0] sum_c;
    assign final_sign = sum_result[SW];
    assign sum_c = final_sign ? ((~sum_result[SW-1:0]) + 1'b1) : sum_result[SW-1:0];

    // ---- normalisation --------------------------------------------------
    logic [3:0]  lz;
    logic [3:0]  lzc_zero;
    logic [4:0]  shift;
    logic [AW-1:0] rg_exp_adjust;
    logic [AW-1:0] final_rg_exp;
    logic [16:0] sum_norm;
    logic [MW5-1:0] final_mant;

    lzc_lead #(.W(SW), .CW(5)) u_lzsum (.in(sum_c), .cnt({lzc_zero, lz}));

    assign shift = lzc_zero ? 5'd14 : {1'b0, lz};
    assign sum_norm = {1'b0, sum_c} << shift;
    assign rg_exp_adjust = lzc_zero ? {AW{1'b0}}
                                    : {{(AW-6){1'b0}}, 3'b100} - {1'b0, lz, 1'b0};
    assign final_rg_exp = l4_emax + rg_exp_adjust;

    logic sticky_bit;
    assign sticky_bit = |sum_norm[4:0];
    assign final_mant = {sum_norm[15:5], sticky_bit};

    // ================================================================ L5
    logic [REG5-1:0] l5_reg;
    always_ff @(posedge clk) begin
        l5_reg <= {l4_emax, final_sign, final_mant};
    end

    // ---- posit(16,2) encode ---------------------------------------------
    generate
        if (1) begin : g_enc
            logic [AW-1:0]       rg_exp;
            logic                sign_i;
            logic [MW5-1:0]      mant_norm;
            logic signed [AW-es_o-1:0] regime_k;
            logic signed [es_o:0]      exp;
            logic                sign_k;
            logic [14:0]         regime;
            logic [AW-es_o-1:0]  regime_bits;
            logic [21:0]         rg_exp_mant;
            logic [21:0]         value_before_shift;
            logic [21:0]         value_after_shift;
            logic [3:0]          shift_amount;
            logic [15:0]         value_before_round;
            logic [2:0]          rounding_bits;
            logic                round_bit, sticky_bit_e, round_value;
            logic [15:0]         value_after_round;
            logic [15:0]         normal_result;
            logic                input_not_zero;

            assign rg_exp    = l5_reg[MW5+1 +: AW];
            assign sign_i    = l5_reg[MW5];
            assign mant_norm = l5_reg[MW5-1:0];

            assign input_not_zero = mant_norm[MW5-1];
            assign regime_k = rg_exp[AW-1:es_o];
            assign exp      = rg_exp[es_o-1:0];
            assign sign_k   = rg_exp[AW-1];
            assign regime   = sign_k ? 15'd1 : ~15'd1;
            assign regime_bits = sign_k ? ((~regime_k) + 2) : (regime_k + 2);

            assign rg_exp_mant = {regime, exp[es_o-1:0], mant_norm[MW5-2:0]};
            assign value_before_shift = rg_exp_mant << 3;
            assign shift_amount = (regime_bits >= 16) ? 4'd11
                                 : (regime_bits[3:0] + 4'd4);
            assign value_after_shift = value_before_shift >> shift_amount;

            assign {value_before_round, rounding_bits} = value_after_shift[18:0];
            assign round_bit = rounding_bits[2];
            assign sticky_bit_e = |rounding_bits[1:0];
            assign round_value = round_bit & (sticky_bit_e | value_before_round[0]);
            assign value_after_round = value_before_round + {15'd0, round_value};
            assign normal_result = sign_i ? {1'b1, (~value_after_round[14:0]) + 15'd1}
                                          : {1'b0, value_after_round[14:0]};
            assign result_o = input_not_zero ? normal_result : 16'd0;
        end
    endgenerate

endmodule

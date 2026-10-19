// ---------------------------------------------------------------------------
// 5-stage pipelined 16-lane posit(8,2) dot-product-accumulate unit.
// Bit-exact re-implementation of the reference pdpu_top_pipelined datapath.
//
//   A : 32x posit(8,2) decode + acc decode + exp add + max-tree level 1
//   B : 16x 4x4 mantissa multiply + max-tree levels 2..4
//   C : max selection + exponent difference + barrel alignment
//   D : sign fold + popcount + Wallace reduction
//   E : final add, |.|, leading-zero count, normalisation
//   F : posit(16,2) encode
// ---------------------------------------------------------------------------

module smax #(
    parameter int W = 6
)(
    input  logic [W-1:0] a,
    input  logic [W-1:0] b,
    output logic [W-1:0] o
);
    assign o = ($signed(a) > $signed(b)) ? a : b;
endmodule

// leading-zero count, W when the input is zero
module lzc_p #(
    parameter int W  = 19,
    parameter int CW = 5
)(
    input  logic [W-1:0]  in,
    output logic [CW-1:0] cnt
);
    integer i;
    always_comb begin
        cnt = W;
        for (i = 0; i < W; i++)
            if (in[i]) cnt = W - 1 - i;
    end
endmodule

// One 3:2 Wallace reduction level: NP inputs -> NP-NT outputs, the first 2*NT
// being (sum, carry<<1) pairs, the rest passed through.
module wallace_lvl #(
    parameter int NP = 18,
    parameter int NT = 6,
    parameter int W  = 20
)(
    input  logic [NP*W-1:0]      i,
    output logic [(NP-NT)*W-1:0] o
);
    generate
        genvar g;
        for (g = 0; g < NT; g++) begin : g_fa
            logic [W-1:0] a, b, c;
            assign a = i[(3*g+0)*W +: W];
            assign b = i[(3*g+1)*W +: W];
            assign c = i[(3*g+2)*W +: W];
            assign o[(2*g+0)*W +: W] = a ^ b ^ c;
            assign o[(2*g+1)*W +: W] = ((a & b) | (b & c) | (a & c)) << 1;
        end
        for (g = 0; g < NP-3*NT; g++) begin : g_pass
            assign o[(2*NT+g)*W +: W] = i[(3*NT+g)*W +: W];
        end
    endgenerate
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
    localparam int unsigned PW = ALIGN_WIDTH;   // 14 : aligned product width
    localparam int unsigned TW = 20;            // signed accumulate width
    localparam int unsigned NW = N + 1;         // 17 accumulate operands
    localparam int unsigned EW = 6;             // exponent field width (signed)

    // ============================================================== block A
    // ---- posit(8,2) decode ------------------------------------------------
    //  v   = (s ? ~x+1 : x)[6:0]                  magnitude
    //  c   = v[6] ? ~v : v                        lzc input (0 only if v==0)
    //  lz  = leading-zero count of c
    //  sh  = {1'b0,v} << lz                       leading 1 lands on bit 6
    //  k   = v[6] ? lz-1 : -lz                    regime value
    //  e   = 4*k + sh[5:4]                        signed exponent
    //  mant= {|v, sh[3:1]}                        implicit bit + 3 fraction
    logic [N*EW-1:0] t_val;
    logic [N*4-1:0]  m_val;
    logic [N-1:0]    s_val;
    logic [EW-1:0]   e_acc;
    logic [11:0]     m_acc;
    logic            s_acc;

    genvar g;
    generate
        for (g = 0; g < N; g++) begin : g_dec
            logic [7:0] xa, xb, sha, shb;
            logic [6:0] va, vb, ca, cb;
            logic [2:0] la, lb;
            logic [4:0] ka, kb;
            logic [5:0] ea, eb;

            assign xa = operands_a[g*n_i +: n_i];
            assign xb = operands_b[g*n_i +: n_i];
            assign va = xa[7] ? (~xa[6:0] + 7'd1) : xa[6:0];
            assign vb = xb[7] ? (~xb[6:0] + 7'd1) : xb[6:0];
            assign ca = va[6] ? ~va : va;
            assign cb = vb[6] ? ~vb : vb;

            lzc_p #(.W(7), .CW(3)) u_la (.in(ca), .cnt(la));
            lzc_p #(.W(7), .CW(3)) u_lb (.in(cb), .cnt(lb));

            assign sha = {1'b0, va} << la;
            assign shb = {1'b0, vb} << lb;

            assign ka = va[6] ? ({1'b0, la} - 5'd1) : (5'd0 - {2'b00, la});
            assign kb = vb[6] ? ({1'b0, lb} - 5'd1) : (5'd0 - {2'b00, lb});

            assign ea = {ka, 2'b00} + {4'd0, sha[5:4]};
            assign eb = {kb, 2'b00} + {4'd0, shb[5:4]};

            assign t_val[g*EW +: EW] = ea;
            assign m_val[g*4 +: 4]   = {|xa[6:0], sha[3:1]};
            assign s_val[g]          = xa[7] ^ xb[7];
        end
    endgenerate

    // ---- posit(16,2) decode of the accumulator ----------------------------
    generate
        if (1) begin : g_accdec
            logic [14:0] av, ac, avs;
            logic [3:0]  lz;
            logic [4:0]  kk;
            logic [5:0]  ee;

            assign av = acc[15] ? (~acc[14:0] + 15'd1) : acc[14:0];
            assign ac = av[14] ? ~av : av;

            lzc_p #(.W(15), .CW(4)) u_lacc (.in(ac), .cnt(lz));

            assign avs = {1'b0, av} << lz;
            assign kk  = av[14] ? ({1'b0, lz} - 5'd1) : (5'd0 - {2'b00, lz});
            assign ee  = {kk, 2'b00} + {4'd0, avs[13:12]};

            assign e_acc = ee;
            assign m_acc = {|acc[14:0], avs[11:1]};
            assign s_acc = acc[15];
        end
    endgenerate

    // ---- max tree level 1 -------------------------------------------------
    logic [N*EW-1:0] tsum;
    logic [8*EW-1:0] tmx1;

    generate
        for (g = 0; g < N; g++) begin : g_tsum
            assign tsum[g*EW +: EW] = $signed(t_val[g*EW +: EW])
                                    + $signed(t_val[g*EW +: EW]);
        end
        for (g = 0; g < 8; g++) begin : g_mx1
            smax #(.W(EW)) u_m (.a(tsum[(2*g)*EW +: EW]), .b(tsum[(2*g+1)*EW +: EW]),
                                .o(tmx1[g*EW +: EW]));
        end
    endgenerate

    // ============================================================== pipe1
    logic [N*EW-1:0] p1_t;
    logic [N*4-1:0]  p1_m;
    logic [N-1:0]    p1_s;
    logic [8*EW-1:0] p1_tmx;
    logic [EW-1:0]   p1_eacc;
    logic [11:0]     p1_macc;
    logic            p1_sacc;

    always_ff @(posedge clk) begin
        p1_t    <= tsum;
        p1_m    <= m_val;
        p1_s    <= s_val;
        p1_tmx  <= tmx1;
        p1_eacc <= e_acc;
        p1_macc <= m_acc;
        p1_sacc <= s_acc;
    end

    // ============================================================== block B
    logic [NW*PW-1:0] prod;

    generate
        for (g = 0; g < N; g++) begin : g_mul
            assign prod[g*PW +: PW] = (p1_m[g*4 +: 4] * p1_m[g*4 +: 4]) << 6;
        end
    endgenerate
    assign prod[N*PW +: PW] = {1'b0, p1_macc, 1'b0};

    // ---- max tree levels 2..4 --------------------------------------------
    logic [4*EW-1:0] tmx2;
    logic [2*EW-1:0] tmx3;
    logic [EW-1:0]   tmx4;

    generate
        for (g = 0; g < 4; g++) begin : g_mx2
            smax #(.W(EW)) u_m (.a(p1_tmx[(2*g)*EW +: EW]), .b(p1_tmx[(2*g+1)*EW +: EW]),
                                .o(tmx2[g*EW +: EW]));
        end
        for (g = 0; g < 2; g++) begin : g_mx3
            smax #(.W(EW)) u_m (.a(tmx2[(2*g)*EW +: EW]), .b(tmx2[(2*g+1)*EW +: EW]),
                                .o(tmx3[g*EW +: EW]));
        end
        smax #(.W(EW)) u_mx4 (.a(tmx3[0 +: EW]), .b(tmx3[EW +: EW]), .o(tmx4));
    endgenerate

    // ============================================================== pipe2
    logic [N*EW-1:0]  p2_t;
    logic [N-1:0]     p2_s;
    logic [NW*PW-1:0] p2_prod;
    logic [EW-1:0]    p2_eacc;
    logic [EW-1:0]    p2_tmx;
    logic             p2_sacc;

    always_ff @(posedge clk) begin
        p2_t    <= p1_t;
        p2_s    <= p1_s;
        p2_prod <= prod;
        p2_eacc <= p1_eacc;
        p2_tmx  <= tmx4;
        p2_sacc <= p1_sacc;
    end

    // ============================================================== block C
    logic [EW-1:0] e_max;

    smax #(.W(EW)) u_emax (.a(p2_tmx), .b(p2_eacc), .o(e_max));

    logic [NW*PW-1:0] aligned;

    generate
        for (g = 0; g < NW; g++) begin : g_aln
            logic [EW-1:0] ed;
            logic [13:0]   sh1, sh2, sh3, sh4, pv;
            assign ed = (g == N) ? (e_max - p2_eacc)
                                 : (e_max - p2_t[g*EW +: EW]);
            assign pv = p2_prod[g*PW +: PW];
            assign sh1 = ed[0] ? {1'b0, pv[13:1]} : pv;
            assign sh2 = ed[1] ? {2'b0, sh1[13:2]} : sh1;
            assign sh3 = ed[2] ? {4'b0, sh2[13:4]} : sh2;
            assign sh4 = ed[3] ? {8'b0, sh3[13:8]} : sh3;
            assign aligned[g*PW +: PW] = (|ed[5:4]) ? 14'd0 : sh4;
        end
    endgenerate

    // ============================================================== pipe3
    logic [NW*PW-1:0] p3_al;
    logic [NW-1:0]    p3_s;
    logic [EW-1:0]    p3_emax;

    always_ff @(posedge clk) begin
        p3_al   <= aligned;
        p3_s    <= {p2_sacc, p2_s};
        p3_emax <= e_max;
    end

    // ============================================================== block D
    // signed term = s ? -a : a.  In a full-width two's complement this is
    // (a ^ s) + s, so the column sum is sum(a ^ s) + popcount(s).
    logic [NW*TW-1:0] wt0;
    logic [4:0]       pop;

    generate
        for (g = 0; g < NW; g++) begin : g_fold
            assign wt0[g*TW +: TW] = {{(TW-PW){1'b0}}, p3_al[g*PW +: PW]}
                                   ^ {TW{p3_s[g]}};
        end
    endgenerate

    logic [1:0] pc0,pc1,pc2,pc3,pc4,pc5,pc6,pc7;
    logic [2:0] pd0,pd1,pd2,pd3;
    logic [3:0] pe0,pe1;
    assign pc0 = {1'b0,p3_s[0]}  + {1'b0,p3_s[1]};
    assign pc1 = {1'b0,p3_s[2]}  + {1'b0,p3_s[3]};
    assign pc2 = {1'b0,p3_s[4]}  + {1'b0,p3_s[5]};
    assign pc3 = {1'b0,p3_s[6]}  + {1'b0,p3_s[7]};
    assign pc4 = {1'b0,p3_s[8]}  + {1'b0,p3_s[9]};
    assign pc5 = {1'b0,p3_s[10]} + {1'b0,p3_s[11]};
    assign pc6 = {1'b0,p3_s[12]} + {1'b0,p3_s[13]};
    assign pc7 = {1'b0,p3_s[14]} + {1'b0,p3_s[15]};
    assign pd0 = {1'b0,pc0} + {1'b0,pc1};
    assign pd1 = {1'b0,pc2} + {1'b0,pc3};
    assign pd2 = {1'b0,pc4} + {1'b0,pc5};
    assign pd3 = {1'b0,pc6} + {1'b0,pc7};
    assign pe0 = {1'b0,pd0} + {1'b0,pd1};
    assign pe1 = {1'b0,pd2} + {1'b0,pd3};
    assign pop = {1'b0,pe0} + {1'b0,pe1} + {4'b0,p3_s[16]};

    // Wallace reduction: 18 operands -> 2
    localparam int unsigned W1N = 12;
    localparam int unsigned W2N = 8;
    localparam int unsigned W3N = 6;
    localparam int unsigned W4N = 4;
    localparam int unsigned W5N = 3;

    logic [18*TW-1:0]  wv0;
    logic [W1N*TW-1:0] wv1;
    logic [W2N*TW-1:0] wv2;
    logic [W3N*TW-1:0] wv3;
    logic [W4N*TW-1:0] wv4;
    logic [W5N*TW-1:0] wv5;
    logic [2*TW-1:0]   wv6;

    assign wv0 = {wt0, {{(TW-5){1'b0}}, pop}};

    wallace_lvl #(.NP(18),  .NT(6), .W(TW)) u_w1 (.i(wv0), .o(wv1));
    wallace_lvl #(.NP(W1N), .NT(4), .W(TW)) u_w2 (.i(wv1), .o(wv2));
    wallace_lvl #(.NP(W2N), .NT(2), .W(TW)) u_w3 (.i(wv2), .o(wv3));
    wallace_lvl #(.NP(W3N), .NT(2), .W(TW)) u_w4 (.i(wv3), .o(wv4));
    wallace_lvl #(.NP(W4N), .NT(1), .W(TW)) u_w5 (.i(wv4), .o(wv5));
    wallace_lvl #(.NP(W5N), .NT(1), .W(TW)) u_w6 (.i(wv5), .o(wv6));

    // ============================================================== pipe4
    logic [2*TW-1:0] p4_w;
    logic [EW-1:0]   p4_emax;

    always_ff @(posedge clk) begin
        p4_w    <= wv6;
        p4_emax <= p3_emax;
    end

    // ============================================================== block E
    logic [TW-1:0] a_sum;
    logic          f_sign;
    logic [18:0]   sum_c;

    assign a_sum  = p4_w[0 +: TW] + p4_w[TW +: TW];
    assign f_sign = a_sum[TW-1];
    assign sum_c  = f_sign ? ((~a_sum + 1'b1) & 19'h7FFFF) : a_sum[18:0];

    // ---- leading-zero count and normalisation -----------------------------
    logic [4:0]    lz;
    logic [4:0]    shl;
    logic [18:0]   sum_norm;
    logic [13:0]   f_mant;
    logic [EW-1:0] f_exp;

    lzc_p #(.W(19), .CW(5)) u_lzs (.in(sum_c), .cnt(lz));

    assign shl      = (lz > 5'd18) ? 5'd0 : lz;
    assign sum_norm = sum_c << shl;
    assign f_mant   = {sum_norm[18:6], |sum_norm[5:0]};
    assign f_exp    = p4_emax + (6'sd6 - {1'b0, shl});

    // ============================================================== pipe5
    logic [EW-1:0] p5_exp;
    logic          p5_sign;
    logic [13:0]   p5_mant;

    always_ff @(posedge clk) begin
        p5_exp  <= f_exp;
        p5_sign <= f_sign;
        p5_mant <= f_mant;
    end

    // ============================================================== block F
    generate
        if (1) begin : g_enc
            logic signed [3:0] rk;
            logic [1:0]        ex;
            logic              sk;
            logic [14:0]       regime;
            logic [3:0]        rbits;
            logic [29:0]       rg_exp_mant;
            logic [45:0]       vbs, vas;
            logic [14:0]       vbr;
            logic [15:0]       rlow;
            logic              rbit, stky, rval;
            logic [14:0]       rounded;
            logic [15:0]       res16;

            assign rk = p5_exp[5:2];
            assign ex = p5_exp[1:0];
            assign sk = p5_exp[5];

            assign regime = sk ? 15'd1 : ~15'd1;
            assign rbits  = sk ? ((~rk) + 4'd2) : (rk + 4'd2);
            assign rg_exp_mant = {regime, ex, p5_mant[12:0]};

            assign vbs = {16'd0, rg_exp_mant} << 16;
            assign vas = vbs >> rbits;

            assign vbr   = vas[30:16];
            assign rlow  = vas[15:0];
            assign rbit  = rlow[15];
            assign stky  = |rlow[14:0];
            assign rval  = rbit & (stky | vbr[0]);
            assign rounded = vbr + {14'd0, rval};

            assign res16 = p5_sign ? {1'b1, (~rounded[14:0]) + 15'd1}
                                   : {1'b0, rounded[14:0]};

            assign result_o = p5_mant[13] ? res16 : 16'd0;
        end
    endgenerate

endmodule

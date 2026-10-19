// ===========================================================================
//  PDPU dot-16 : 5-stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model, including the
//  ALIGN_WIDTH=14 truncation, the CARRY_W=5 / SUM_WIDTH=19 wrap-around, the
//  reference decoder/encoder conventions and the NaR handling.
//
//  stage A : posit decode + exponent add + mantissa multiply   (-> reg p1)
//  stage B : 17-way signed maximum exponent tree               (-> reg p2)
//  stage C : alignment barrel shift + sign folding             (-> reg p3)
//  stage D : balanced adder tree + two's complement correction (-> reg p4)
//  stage E : normalisation and posit encoding                  (-> result_o)
//
//  Cost reductions vs. the reference:
//   * posit decode without a 15-bit leading-zero counter and 15-bit shifter
//   * two's complement of the 17 aligned words replaced by a bit inversion
//     plus one constant (K - K<<14) added once at the end
//   * balanced binary adder tree instead of a 17:2 carry-save tree plus a
//     20-bit carry-propagate adder
//   * posit encoder shifts a 31-bit window by s-1 instead of a 46-bit one
// ===========================================================================
module pdpu_top_pipelined (
    input  logic         clk,
    input  logic [127:0] operands_a,
    input  logic [127:0] operands_b,
    input  logic [15:0]  acc,
    output logic [15:0]  result_o
);
    integer i;

    // ------------------------------------------------------------------
    //  posit(8,2) decode  ->  {sign, rg_exp[5:0], mant[3:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] hsb;
        logic [3:0] rl, ek;
        logic       rs, emp;
        va  = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        rs  = va[6];
        x   = rs ? ~va : va;
        emp = (x == 7'd0);
        hsb = 3'd0;
        for (int t = 0; t < 7; t++) if (x[t]) hsb = t[2:0];
        rl  = emp ? 4'd7 : (4'd6 - {1'b0, hsb});
        ek  = rs ? (rl - 4'd1) : (4'd0 - rl);
        al  = va << (emp ? 3'd7 : (3'd7 - hsb));
        dec8 = {p[7], ek, al[6:5], |p[6:0], al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode  ->  {sign, rg_exp[6:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [19:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  hsb;
        logic [4:0]  rl, ek;
        logic        rs, emp;
        va  = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        rs  = va[14];
        x   = rs ? ~va : va;
        emp = (x == 15'd0);
        hsb = 4'd0;
        for (int t = 0; t < 15; t++) if (x[t]) hsb = t[3:0];
        rl  = emp ? 5'd15 : (5'd14 - {1'b0, hsb});
        ek  = rs ? (rl - 5'd1) : (5'd0 - rl);
        al  = va << (emp ? 4'd15 : (4'd15 - hsb));
        dec16 = {p[15], ek, al[14:13], |p[14:0], al[12:2]};
    endfunction

    // ------------------------------------------------------------------
    //  Exact replica of the reference 4x4 radix-4 Booth multiplier with
    //  3:2 compression (including its bit truncations).
    // ------------------------------------------------------------------
    function automatic [7:0] mul44(input [3:0] a, input [3:0] b);
        logic [2:0] c0, c1, c2;
        logic       n0, o0, w0, n1, o1, w1, n2, o2, w2;
        logic [4:0] t0, t1, t2, p0, p1, p2;
        logic       s0, s1;
        logic [7:0] pp0, pp1, pp2, sm, cy;
        logic [6:0] tmp2;
        c0 = {b[1], b[0], 1'b0};
        c1 = {b[3], b[2], b[1]};
        c2 = {1'b0, 1'b0, b[3]};

        n0 = c0[2] & (~c0[0] | ~c0[1]);
        o0 = (c0[1] & ~c0[0]) | (~c0[1] & c0[0]);
        w0 = (~c0[2] & c0[1] & c0[0]) | (c0[2] & ~c0[1] & ~c0[0]);
        t0 = o0 ? {1'b0, a} : (w0 ? {a, 1'b0} : 5'd0);
        p0 = n0 ? ~t0 : t0;
        s0 = n0;
        pp0 = {~s0, s0, s0, p0};

        n1 = c1[2] & (~c1[0] | ~c1[1]);
        o1 = (c1[1] & ~c1[0]) | (~c1[1] & c1[0]);
        w1 = (~c1[2] & c1[1] & c1[0]) | (c1[2] & ~c1[1] & ~c1[0]);
        t1 = o1 ? {1'b0, a} : (w1 ? {a, 1'b0} : 5'd0);
        p1 = n1 ? ~t1 : t1;
        s1 = n1;
        pp1 = {~s1, p1, 1'b0, s0};

        n2 = c2[2] & (~c2[0] | ~c2[1]);
        o2 = (c2[1] & ~c2[0]) | (~c2[1] & c2[0]);
        w2 = (~c2[2] & c2[1] & c2[0]) | (c2[2] & ~c2[1] & ~c2[0]);
        t2 = o2 ? {1'b0, a} : (w2 ? {a, 1'b0} : 5'd0);
        p2 = n2 ? ~t2 : t2;
        tmp2 = {p2, 1'b0, s1};
        pp2 = {tmp2[5:0], 2'b00};

        sm = pp0 ^ pp1 ^ pp2;
        cy = (pp0 & pp1) | (pp0 & pp2) | (pp1 & pp2);
        mul44 = sm + {cy[6:0], 1'b0};
    endfunction

    function automatic [6:0] smax(input [6:0] a, input [6:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // ==================================================================
    //  stage A : decode, exponent add, mantissa multiply
    // ==================================================================
    logic [16:0] a_sg;
    logic [6:0]  a_ep [0:16];
    logic [13:0] a_wd [0:16];

    genvar g;
    generate
    for (g = 0; g < 16; g++) begin : gdec
        logic [10:0] da, db;
        assign da = dec8(operands_a[g*8 +: 8]);
        assign db = dec8(operands_b[g*8 +: 8]);
        assign a_sg[g] = da[10] ^ db[10];
        assign a_ep[g] = {{da[9]}, da[9:4]} + {{db[9]}, db[9:4]};
        assign a_wd[g] = {mul44(da[3:0], db[3:0]), 6'b0};
    end
    endgenerate

    logic [19:0] dacc;
    assign dacc = dec16(acc);
    assign a_sg[16] = dacc[19];
    assign a_ep[16] = dacc[18:12];
    assign a_wd[16] = {1'b0, dacc[11:0], 1'b0};

    logic [16:0] p1_sg;
    logic [6:0]  p1_ep [0:16];
    logic [13:0] p1_wd [0:16];

    always_ff @(posedge clk) begin
        p1_sg <= a_sg;
        for (i = 0; i < 17; i++) begin
            p1_ep[i] <= a_ep[i];
            p1_wd[i] <= a_wd[i];
        end
    end

    // ==================================================================
    //  stage B : 17-way signed maximum exponent
    // ==================================================================
    logic [6:0] b_emax;
    logic [6:0] lv1 [0:7];
    logic [6:0] lv2 [0:3];
    logic [6:0] lv3 [0:1];

    always_comb begin
        for (i = 0; i < 8; i++) lv1[i] = smax(p1_ep[2*i], p1_ep[2*i+1]);
        for (i = 0; i < 4; i++) lv2[i] = smax(lv1[2*i], lv1[2*i+1]);
        lv3[0] = smax(lv2[0], lv2[1]);
        lv3[1] = smax(lv2[2], lv2[3]);
        b_emax = smax(smax(lv3[0], lv3[1]), p1_ep[16]);
    end

    logic [16:0] p2_sg;
    logic [6:0]  p2_ep [0:16];
    logic [13:0] p2_wd [0:16];
    logic [6:0]  p2_emax;

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_emax <= b_emax;
        for (i = 0; i < 17; i++) begin
            p2_ep[i] <= p1_ep[i];
            p2_wd[i] <= p1_wd[i];
        end
    end

    // ==================================================================
    //  stage C : alignment shift + sign folding (bit inversion)
    // ==================================================================
    logic [13:0] c_cw [0:16];

    generate
    for (g = 0; g < 17; g++) begin : gal
        logic [6:0] df;
        logic [3:0] sh;
        assign df = p2_emax - p2_ep[g];
        assign sh = (|df[6:4]) ? 4'd15 : df[3:0];
        assign c_cw[g] = (p2_wd[g] >> sh) ^ {14{p2_sg[g]}};
    end
    endgenerate

    logic [4:0] c_k;
    logic [1:0] k1 [0:7];
    logic [2:0] k2 [0:3];
    logic [3:0] k3 [0:1];
    logic [4:0] k4;

    always_comb begin
        for (i = 0; i < 8; i++) k1[i] = {1'b0, p2_sg[2*i]} + {1'b0, p2_sg[2*i+1]};
        for (i = 0; i < 4; i++) k2[i] = {1'b0, k1[2*i]} + {1'b0, k1[2*i+1]};
        for (i = 0; i < 2; i++) k3[i] = {1'b0, k2[2*i]} + {1'b0, k2[2*i+1]};
        k4   = {1'b0, k3[0]} + {1'b0, k3[1]};
        c_k  = k4 + {4'b0, p2_sg[16]};
    end

    logic [13:0] p3_cw [0:16];
    logic [4:0]  p3_k;
    logic [6:0]  p3_emax;

    always_ff @(posedge clk) begin
        p3_k    <= c_k;
        p3_emax <= p2_emax;
        for (i = 0; i < 17; i++) p3_cw[i] <= c_cw[i];
    end

    // ==================================================================
    //  stage D : balanced adder tree + two's complement correction
    // ==================================================================
    logic [14:0] d1 [0:7];
    logic [15:0] d2 [0:3];
    logic [16:0] d3 [0:1];
    logic [17:0] d4;
    logic [19:0] d_corr, d_side, d_tot;

    always_comb begin
        for (i = 0; i < 8; i++) d1[i] = {1'b0, p3_cw[2*i]} + {1'b0, p3_cw[2*i+1]};
        for (i = 0; i < 4; i++) d2[i] = {1'b0, d1[2*i]} + {1'b0, d1[2*i+1]};
        for (i = 0; i < 2; i++) d3[i] = {1'b0, d2[2*i]} + {1'b0, d2[2*i+1]};
        d4     = {1'b0, d3[0]} + {1'b0, d3[1]};
        d_corr = {15'b0, p3_k} - {1'b0, p3_k, 14'b0};
        d_side = {6'b0, p3_cw[16]} + d_corr;
        d_tot  = {2'b0, d4} + d_side;
    end

    logic [19:0] p4_sum;
    logic [6:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= d_tot;
        p4_emax <= p3_emax;
    end

    // ==================================================================
    //  stage E : normalisation and posit encoding
    // ==================================================================
    logic        e_fsign;
    logic [18:0] e_mag;
    logic [4:0]  e_lzc;
    logic [6:0]  e_rgexp;
    logic [13:0] e_fmant;

    logic [4:0]  e_sreg, e_s1;
    logic [30:0] e_win, e_sh;
    logic [14:0] e_vbr;
    logic        e_rbit, e_stk, e_rval;
    logic [14:0] e_var, e_enc;

    always_comb begin
        logic [18:0] sn;
        logic [6:0]  adj;
        logic        sk;
        logic [4:0]  k5;
        logic [30:0] msk;

        // ---- magnitude and leading zero count ----
        e_fsign = p4_sum[19];
        e_mag   = e_fsign ? (~p4_sum[18:0] + 19'd1) : p4_sum[18:0];
        e_lzc   = 5'd19;
        for (i = 0; i < 19; i++) if (e_mag[i]) e_lzc = 5'd18 - i[4:0];
        adj     = (|e_mag) ? (7'sd6 - $signed({2'b0, e_lzc})) : 7'sd0;
        e_rgexp = $signed({1'b0, p4_emax}) + adj;
        sn      = e_mag << e_lzc;
        e_fmant = {sn[18:6], |sn[5:0]};

        // ---- posit encoding ----
        sk     = e_rgexp[6];
        k5     = e_rgexp[6:2];
        e_sreg = sk ? ((5'd0 - k5) + 5'd2) : (k5 + 5'd2);
        e_s1   = (e_sreg > 5'd16) ? 5'd16 : e_sreg;
        e_s1   = e_s1 - 5'd1;
        e_win  = {1'b0, (sk ? 15'h0001 : 15'h7FFE), e_rgexp[1:0], e_fmant};
        e_sh   = e_win >> e_s1;
        e_vbr  = e_sh[15:1];
        e_rbit = e_sh[0];
        msk    = (31'd1 << e_s1) - 31'd1;
        e_stk  = |(e_win & msk);
        e_rval = e_rbit & (e_stk | e_vbr[0]);
        e_var  = e_vbr + {14'b0, e_rval};
        e_enc  = e_fsign ? (~e_var + 15'd1) : e_var;
    end

    always_ff @(posedge clk) begin
        result_o <= e_fmant[13] ? {e_fsign, e_enc} : 16'd0;
    end

endmodule

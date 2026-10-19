// ===========================================================================
//  PDPU dot-16 : 5-stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model, including the
//  ALIGN_WIDTH=14 truncation, CARRY_W=5 / SUM_WIDTH=19 wrap-around, the
//  reference decoder's exponent convention and the NaR handling.
//
//  stage 1 : posit decode + exponent add + 4x4 mantissa multiply
//  stage 2 : 17-way signed maximum exponent tree
//  stage 3 : alignment barrel shifters + sign folding (bit inversion)
//  stage 4 : balanced adder tree + two's complement correction
//  stage 5 : normalisation (lzc + shift) and posit encoding
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
    //  posit(8,2) decode -> {sign, rg_exp[5:0], mant[3:0]}
    //    rg_exp = (4*regime_k + exp) mod 64   (signed)
    //    mant   = {implicit, frac[2:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] hsb, rl, rb;
        logic [3:0] ek;
        logic       regS, emp;
        va   = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        regS = va[6];
        x    = regS ? ~va : va;
        emp  = (x == 7'd0);
        hsb  = 3'd0;
        for (int t = 0; t < 7; t++) if (x[t]) hsb = t[2:0];
        rl   = emp ? 3'd7 : (3'd6 - hsb);
        rb   = emp ? 3'd7 : (3'd7 - hsb);
        ek   = regS ? ({1'b0, rl} - 4'd1) : (4'd0 - {1'b0, rl});
        al   = va << rb;
        dec8 = {p[7], ek, al[6:5], |p[6:0], al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode -> {sign, rg_exp[6:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [19:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  hsb, rl, rb;
        logic [4:0]  ek;
        logic        regS, emp;
        va    = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        regS  = va[14];
        x     = regS ? ~va : va;
        emp   = (x == 15'd0);
        hsb   = 4'd0;
        for (int t = 0; t < 15; t++) if (x[t]) hsb = t[3:0];
        rl    = emp ? 4'd15 : (4'd14 - hsb);
        rb    = emp ? 4'd15 : (4'd15 - hsb);
        ek    = regS ? ({1'b0, rl} - 5'd1) : (5'd0 - {1'b0, rl});
        al    = va << rb;
        dec16 = {p[15], ek, al[14:13], |p[14:0], al[12:2]};
    endfunction

    // ------------------------------------------------------------------
    //  4x4 mantissa product, exactly as the reference's narrow radix-4
    //  Booth multiplier with 3:2 compression (including its truncation).
    // ------------------------------------------------------------------
    function automatic [7:0] mul44(input [3:0] a, input [3:0] b);
        logic [4:0]  a5, t0, t1, t2, p0, p1, p2;
        logic        s0, s1, o0, o1, w1;
        logic [7:0]  q0, q1, q2, sm, cy;
        a5 = {1'b0, a};
        // code0 = {b1, b0, 0}
        s0 = b[1];
        o0 = b[0];
        t0 = o0 ? a5 : ((b[1] & ~b[0]) ? (a5 << 1) : 5'd0);
        p0 = s0 ? ~t0 : t0;
        q0 = {~s0, s0, s0, p0};
        // code1 = {b3, b2, b1}
        s1 = b[3] & ~(b[2] & b[1]);
        o1 = b[2] ^ b[1];
        w1 = (~b[3] & b[2] & b[1]) | (b[3] & ~b[2] & ~b[1]);
        t1 = o1 ? a5 : (w1 ? (a5 << 1) : 5'd0);
        p1 = s1 ? ~t1 : t1;
        q1 = {~s1, p1, 1'b0, s0};
        // code2 = {0, 0, b3}
        t2 = b[3] ? a5 : 5'd0;
        q2 = {1'b0, t2, 1'b0, s1} << 2;
        // 3:2 compression, carry truncated to 8 bits
        sm = q0 ^ q1 ^ q2;
        cy = (q0 & q1) | (q0 & q2) | (q1 & q2);
        mul44 = sm + {cy[5:0], 1'b0};
    endfunction

    function automatic [6:0] smax(input [6:0] a, input [6:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // ==================================================================
    //  stage 1 : decode + multiply + exponent add
    // ==================================================================
    logic [16:0] s1_sg;
    logic [6:0]  s1_ep [0:16];
    logic [13:0] s1_wd [0:16];

    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [10:0] da, db;
        assign da = dec8(operands_a[i*8 +: 8]);
        assign db = dec8(operands_b[i*8 +: 8]);
        assign s1_sg[i] = da[10] ^ db[10];
        assign s1_ep[i] = {{da[9]}, da[9:4]} + {{db[9]}, db[9:4]};
        assign s1_wd[i] = {mul44(da[3:0], db[3:0]), 6'b0};
    end
    endgenerate

    logic [19:0] dacc;
    assign dacc = dec16(acc);
    assign s1_sg[16] = dacc[19];
    assign s1_ep[16] = dacc[18:12];
    assign s1_wd[16] = {1'b0, dacc[11:0], 1'b0};

    logic [16:0] p1_sg;
    logic [6:0]  p1_ep [0:16];
    logic [13:0] p1_wd [0:16];

    always_ff @(posedge clk) begin
        p1_sg <= s1_sg;
        for (i = 0; i < 17; i++) begin
            p1_ep[i] <= s1_ep[i];
            p1_wd[i] <= s1_wd[i];
        end
    end

    // ==================================================================
    //  stage 2 : 17-way signed maximum exponent
    // ==================================================================
    logic [6:0] s2_emax;
    logic [6:0] t1 [0:7];
    logic [6:0] t2 [0:3];
    logic [6:0] t3 [0:1];

    always_comb begin
        for (i = 0; i < 8; i++) t1[i] = smax(p1_ep[2*i], p1_ep[2*i+1]);
        for (i = 0; i < 4; i++) t2[i] = smax(t1[2*i], t1[2*i+1]);
        for (i = 0; i < 2; i++) t3[i] = smax(t2[2*i], t2[2*i+1]);
        s2_emax = smax(smax(t3[0], t3[1]), p1_ep[16]);
    end

    logic [16:0] p2_sg;
    logic [6:0]  p2_ep [0:16];
    logic [13:0] p2_wd [0:16];
    logic [6:0]  p2_emax;

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_emax <= s2_emax;
        for (i = 0; i < 17; i++) begin
            p2_ep[i] <= p1_ep[i];
            p2_wd[i] <= p1_wd[i];
        end
    end

    // ==================================================================
    //  stage 3 : alignment shift + sign folding
    // ==================================================================
    logic [13:0] s3_cw [0:16];
    logic [4:0]  s3_k;

    generate
    for (i = 0; i < 17; i++) begin : g_al
        logic [6:0] df;
        logic [3:0] sh;
        assign df = p2_emax - p2_ep[i];
        assign sh = (|df[6:4]) ? 4'd14 : df[3:0];
        assign s3_cw[i] = (p2_wd[i] >> sh) ^ {14{p2_sg[i]}};
    end
    endgenerate

    logic [1:0] c1 [0:7];
    logic [2:0] c2 [0:3];
    logic [3:0] c3 [0:1];
    logic [4:0] c4;

    always_comb begin
        for (i = 0; i < 8; i++) c1[i] = {1'b0, p2_sg[2*i]} + {1'b0, p2_sg[2*i+1]};
        for (i = 0; i < 4; i++) c2[i] = {1'b0, c1[2*i]} + {1'b0, c1[2*i+1]};
        for (i = 0; i < 2; i++) c3[i] = {1'b0, c2[2*i]} + {1'b0, c2[2*i+1]};
        c4   = {1'b0, c3[0]} + {1'b0, c3[1]};
        s3_k = c4 + {4'b0, p2_sg[16]};
    end

    logic [13:0] p3_cw [0:16];
    logic [4:0]  p3_k;
    logic [6:0]  p3_emax;

    always_ff @(posedge clk) begin
        p3_k    <= s3_k;
        p3_emax <= p2_emax;
        for (i = 0; i < 17; i++) p3_cw[i] <= s3_cw[i];
    end

    // ==================================================================
    //  stage 4 : balanced adder tree + two's complement correction
    // ==================================================================
    logic [14:0] s4_a [0:7];
    logic [15:0] s4_b [0:3];
    logic [16:0] s4_c [0:1];
    logic [17:0] s4_d;
    logic [18:0] s4_e;
    logic [19:0] s4_sum;

    always_comb begin
        for (i = 0; i < 8; i++) s4_a[i] = {1'b0, p3_cw[2*i]} + {1'b0, p3_cw[2*i+1]};
        for (i = 0; i < 4; i++) s4_b[i] = {1'b0, s4_a[2*i]} + {1'b0, s4_a[2*i+1]};
        s4_c[0] = {1'b0, s4_b[0]} + {1'b0, s4_b[1]};
        s4_c[1] = {1'b0, s4_b[2]} + {1'b0, s4_b[3]};
        s4_d    = {1'b0, s4_c[0]} + {1'b0, s4_c[1]};
        s4_e    = {1'b0, s4_d} + {5'b0, p3_cw[16]};
        s4_sum  = {1'b0, s4_e} + {14'b0, p3_k} - {1'b0, p3_k, 14'b0};
    end

    logic [19:0] p4_sum;
    logic [6:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= s4_sum;
        p4_emax <= p3_emax;
    end

    // ==================================================================
    //  stage 5 : normalise + posit encode
    // ==================================================================
    logic        s5_sign;
    logic [18:0] s5_mag;
    logic [4:0]  s5_lzc;
    logic [6:0]  s5_rgexp;
    logic [13:0] s5_fmant;
    logic        s5_nz;

    logic [4:0]  e_rb, e_sh;
    logic [30:0] e_m31, e_y;
    logic [14:0] e_mag, e_fin;
    logic        e_rbit, e_stky, e_rval;
    logic [15:0] e_out;

    always_comb begin
        logic [18:0] sn;
        logic [6:0]  adj;
        logic [4:0]  k5;
        logic        sk;
        s5_sign = p4_sum[19];
        s5_mag  = s5_sign ? (~p4_sum[18:0] + 19'd1) : p4_sum[18:0];
        s5_nz   = |s5_mag;
        s5_lzc  = 5'd19;
        for (i = 0; i < 19; i++) if (s5_mag[i]) s5_lzc = 5'd18 - i[4:0];
        adj     = s5_nz ? (7'sd6 - $signed({2'b0, s5_lzc})) : 7'sd0;
        s5_rgexp = $signed({1'b0, p4_emax}) + adj;
        sn      = s5_mag << s5_lzc;
        s5_fmant = {sn[18:6], |sn[5:0]};

        // ---- posit encode ----
        sk    = s5_rgexp[6];
        k5    = s5_rgexp[6:2];
        e_rb  = sk ? ((5'd0 - k5) + 5'd2) : (k5 + 5'd2);
        e_sh  = (e_rb >= 5'd16) ? 5'd16 : e_rb;
        e_m31 = {1'b0, sk ? 15'h0001 : 15'h7FFE, s5_rgexp[1:0], s5_fmant[12:0]};
        e_y   = e_m31 >> e_sh;
        e_rbit = e_y[15];
        e_stky = |e_y[14:0];
        e_rval = e_rbit & e_stky;
        e_mag  = e_y[14:0] + {14'b0, e_rval};
        e_fin  = sk ? (~e_mag + 15'd1) : e_mag;
        e_out  = {sk, e_fin};
    end

    always_ff @(posedge clk) begin
        result_o <= s5_fmant[13] ? e_out : 16'd0;
    end

endmodule

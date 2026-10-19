// ===========================================================================
//  PDPU dot-16 : 5-stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model, including the
//  ALIGN_WIDTH=14 truncation, the CARRY_W=5 / SUM_WIDTH=19 wrap-around, the
//  reference decoder/encoder conventions and the NaR handling.
//
//  stage 1 : posit decode + exponent add + mantissa multiply
//  stage 2 : 17-way signed maximum exponent tree
//  stage 3 : alignment barrel shift + sign folding (bit inversion)
//  stage 4 : balanced adder tree + two's complement correction
//  stage 5 : normalisation and posit encoding
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
    //    rg_exp = (4*regime_k + exp)          (6-bit signed)
    //    mant   = {implicit, frac[2:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] hsb, rb;
        logic [3:0] rl, ek;
        logic       rs, emp;
        va  = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        rs  = va[6];
        x   = rs ? ~va : va;
        emp = (x == 7'd0);
        hsb = 3'd0;
        for (int t = 0; t < 7; t++) if (x[t]) hsb = t[2:0];
        rl  = emp ? 4'd7  : (4'd6 - {1'b0, hsb});
        rb  = emp ? 3'd7  : (3'd7 - hsb);
        ek  = rs  ? (rl - 4'd1) : (4'd0 - rl);
        al  = va << rb;
        dec8 = {p[7], ek, al[6:5], |va, al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode -> {sign, rg_exp[6:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [19:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  hsb, rb;
        logic [4:0]  rl, ek;
        logic        rs, emp;
        va  = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        rs  = va[14];
        x   = rs ? ~va : va;
        emp = (x == 15'd0);
        hsb = 4'd0;
        for (int t = 0; t < 15; t++) if (x[t]) hsb = t[3:0];
        rl  = emp ? 5'd15 : (5'd14 - {1'b0, hsb});
        rb  = emp ? 4'd15 : (4'd15 - hsb);
        ek  = rs  ? (rl - 5'd1) : (5'd0 - rl);
        al  = va << rb;
        dec16 = {p[15], ek, al[14:13], |va, al[12:2]};
    endfunction

    function automatic [6:0] smax(input [6:0] a, input [6:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // ==================================================================
    //  stage 1 : decode, exponent add, mantissa multiply
    // ==================================================================
    logic [16:0] s1_sg;
    logic [6:0]  s1_ep [0:16];
    logic [13:0] s1_wd [0:16];

    genvar g;
    generate
    for (g = 0; g < 16; g++) begin : gdec
        logic [10:0] da, db;
        logic [7:0]  pr;
        assign da = dec8(operands_a[g*8 +: 8]);
        assign db = dec8(operands_b[g*8 +: 8]);
        assign s1_sg[g] = da[10] ^ db[10];
        assign s1_ep[g] = {{da[9]}, da[9:4]} + {{db[9]}, db[9:4]};
        assign pr       = da[3:0] * db[3:0];
        assign s1_wd[g] = {pr, 6'b0};
    end
    endgenerate

    logic [19:0] dacc;
    assign dacc = dec16(acc);
    assign s1_sg[16] = dacc[19];
    assign s1_ep[16] = {{dacc[18]}, dacc[18:12]};
    assign s1_wd[16] = {2'b0, dacc[11:0]};

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
    //  stage 3 : alignment shift + sign folding (bit inversion)
    // ==================================================================
    logic [13:0] c_cw [0:16];

    generate
    for (g = 0; g < 17; g++) begin : gal
        logic [6:0] df;
        logic [3:0] sh;
        assign df = p2_emax - p2_ep[g];
        assign sh = (|df[6:4]) ? 4'd14 : df[3:0];
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
    //  stage 4 : balanced adder tree + two's complement correction
    // ==================================================================
    logic [14:0] d1 [0:7];
    logic [15:0] d2 [0:3];
    logic [16:0] d3 [0:1];
    logic [17:0] d4;
    logic [18:0] d5;
    logic [19:0] d_tot;

    always_comb begin
        for (i = 0; i < 8; i++) d1[i] = {1'b0, p3_cw[2*i]} + {1'b0, p3_cw[2*i+1]};
        for (i = 0; i < 4; i++) d2[i] = {1'b0, d1[2*i]} + {1'b0, d1[2*i+1]};
        d3[0] = {1'b0, d2[0]} + {1'b0, d2[1]};
        d3[1] = {1'b0, d2[2]} + {1'b0, d2[3]};
        d4    = {1'b0, d3[0]} + {1'b0, d3[1]};
        d5    = {1'b0, d4} + {5'b0, p3_cw[16]};
        d_tot = {1'b0, d5} + {15'b0, p3_k} - {1'b0, p3_k, 14'b0};
    end

    logic [19:0] p4_sum;
    logic [6:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= d_tot;
        p4_emax <= p3_emax;
    end

    // ==================================================================
    //  stage 5 : normalisation and posit encoding
    // ==================================================================
    logic        e_fsign;
    logic [18:0] e_mag;
    logic [4:0]  e_lzc;
    logic [6:0]  e_rgexp;
    logic [14:0] e_fmant;

    logic [4:0]  e_rb, e_sh;
    logic [29:0] e_win, e_y;
    logic [14:0] e_vbr;
    logic        e_rbit, e_stk, e_rval;
    logic [14:0] e_var, e_enc;

    always_comb begin
        logic [18:0] sn;
        logic [6:0]  adj;
        logic        sk;
        logic [4:0]  k6;
        logic [28:0] msk;

        // ---- magnitude, leading zero count, normalisation ----
        e_fsign = p4_sum[19];
        e_mag   = e_fsign ? (~p4_sum[18:0] + 19'd1) : p4_sum[18:0];
        e_lzc   = 5'd19;
        for (i = 0; i < 19; i++) if (e_mag[i]) e_lzc = 5'd18 - i[4:0];
        adj     = (|e_mag) ? (7'sd6 - $signed({2'b0, e_lzc})) : 7'sd0;
        e_rgexp = $signed({1'b0, p4_emax}) + adj;
        sn      = e_mag << e_lzc;
        e_fmant = {sn[18:5], |sn[4:0]};

        // ---- posit encoding ----
        sk     = e_rgexp[6];
        k6     = e_rgexp[6:2];
        e_rb   = sk ? (5'd1 - k6) : (k6 + 5'd2);
        e_sh   = (e_rb > 5'd16) ? 5'd16 : e_rb;
        e_win  = {1'b0, (sk ? 15'h0001 : 15'h7FFE), e_rgexp[1:0], e_fmant[12:0]};
        e_y    = e_win >> e_sh;
        e_vbr  = e_y[14:0];
        e_rbit = e_win[e_sh - 5'd1];
        msk    = (29'd1 << (e_sh - 5'd1)) - 29'd1;
        e_stk  = |(e_win[28:0] & msk);
        e_rval = e_rbit & (e_stk | e_vbr[0]);
        e_var  = e_vbr + {14'b0, e_rval};
        e_enc  = e_fsign ? (~e_var + 15'd1) : e_var;
    end

    always_ff @(posedge clk) begin
        result_o <= e_fmant[13] ? {e_fsign, e_enc} : 16'd0;
    end

endmodule

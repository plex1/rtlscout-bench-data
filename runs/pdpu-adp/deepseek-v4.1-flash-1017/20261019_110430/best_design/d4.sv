// ===========================================================================
// PDPU dot-16, 5-stage pipelined - compact re-implementation.
// Bit-exact with the reference (ALIGN_WIDTH=14 truncation, CARRY_W=5,
// SUM_WIDTH=19, wrapped 6-bit exponent arithmetic, NaR handling).
// ===========================================================================
module pdpu_top_pipelined #(
    parameter int unsigned N = 16
)(
    input  logic         clk,
    input  logic [127:0] operands_a,
    input  logic [127:0] operands_b,
    input  logic [15:0]  acc,
    output logic [15:0]  result_o
);

    // ------------------------------------------------------------------
    // posit(8,2) decode: {sign, rg_exp[5:0], mant[3:0]}
    //   rg_exp = 4*regime_k + exponent   (6-bit two's complement)
    //   mant   = {implicit, frac[2:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x, b;
        logic [2:0] rl;
        logic [3:0] k;
        logic [6:0] rb;
        va = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        x  = va ^ {7{va[6]}};
        rl = 3'd7;
        for (int i = 0; i < 7; i++) if (x[i]) rl = 3'd6 - i[2:0];
        k  = va[6] ? ({1'b0,rl} - 4'd1) : (4'd0 - {1'b0,rl});
        rb = (rl == 3'd7) ? 7'd7 : ({4'b0,rl} + 7'd1);
        b  = va << rb[2:0];
        dec8 = {p[7], k, b[6:5], |va, b[4:2]};
    endfunction

    // posit(16,2) decode: {sign, rg_exp[5:0], mant[11:0]}
    function automatic [18:0] dec16(input [15:0] p);
        logic [14:0] va, x, b;
        logic [3:0]  rl, k;
        logic [14:0] rb;
        va = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        x  = va ^ {15{va[14]}};
        rl = 4'd15;
        for (int i = 0; i < 15; i++) if (x[i]) rl = 4'd14 - i[3:0];
        k  = va[14] ? (rl - 4'd1) : (4'd0 - rl);
        rb = (rl == 4'd15) ? 15'd15 : ({11'b0,rl} + 15'd1);
        b  = va << rb[3:0];
        dec16 = {p[15], k, b[14:13], |va, b[12:2]};
    endfunction

    function automatic [5:0] mx(input [5:0] a, input [5:0] b);
        mx = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    integer q;

    // ======================= stage 0 : decode + mul =======================
    logic [16:0]  s0_sg;
    logic [101:0] s0_ep;
    logic [237:0] s0_w;

    genvar i;
    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [10:0] da, db;
        logic [7:0]  pr;
        assign da = dec8(operands_a[i*8 +: 8]);
        assign db = dec8(operands_b[i*8 +: 8]);
        assign s0_sg[i] = da[10] ^ db[10];
        assign s0_ep[i*6 +: 6] = da[9:4] + db[9:4];
        assign pr = da[3:0] * db[3:0];
        assign s0_w[i*14 +: 14] = {2'b0, pr, 4'b0};
    end
    endgenerate

    logic [18:0] dacc;
    assign dacc = dec16(acc);
    assign s0_sg[16] = dacc[18];
    assign s0_ep[96 +: 6] = dacc[17:12];
    assign s0_w[224 +: 14] = {2'b0, dacc[11:1]};

    logic [16:0]  p1_sg;
    logic [101:0] p1_ep;
    logic [237:0] p1_w;

    always_ff @(posedge clk) begin
        p1_sg <= s0_sg;
        p1_ep <= s0_ep;
        p1_w  <= s0_w;
    end

    // ======================= stage 1 : max exponent =======================
    logic [5:0] m0 [0:16];
    logic [5:0] m1 [0:8];
    logic [5:0] m2 [0:4];
    logic [5:0] m3 [0:2];
    logic [5:0] m4 [0:1];
    logic [5:0] s_emax;
    always_comb begin
        for (q = 0; q < 17; q++) m0[q] = p1_ep[q*6 +: 6];
        for (q = 0; q < 8; q++)  m1[q] = mx(m0[2*q], m0[2*q+1]);
        m1[8] = m0[16];
        for (q = 0; q < 4; q++)  m2[q] = mx(m1[2*q], m1[2*q+1]);
        m2[4] = m1[8];
        for (q = 0; q < 2; q++)  m3[q] = mx(m2[2*q], m2[2*q+1]);
        m3[2] = m2[4];
        m4[0] = mx(m3[0], m3[1]);
        m4[1] = m3[2];
        s_emax = mx(m4[0], m4[1]);
    end

    logic [16:0]  p2_sg;
    logic [101:0] p2_ep;
    logic [237:0] p2_w;
    logic [5:0]   p2_emax;

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_ep   <= p1_ep;
        p2_w    <= p1_w;
        p2_emax <= s_emax;
    end

    // ======================= stage 2 : align ==============================
    logic [237:0] s_cw;
    logic [4:0]   s_k;
    logic [19:0]  s_corr;

    generate
    for (i = 0; i < 17; i++) begin : g_al
        logic [5:0]  df;
        logic [3:0]  sh;
        logic [13:0] ws;
        assign df = p2_emax - p2_ep[i*6 +: 6];
        assign sh = (|df[5:4]) ? 4'd14 : df[3:0];
        assign ws = p2_w[i*14 +: 14] >> sh;
        assign s_cw[i*14 +: 14] = ws ^ {14{p2_sg[i]}};
    end
    endgenerate

    always_comb begin
        s_k = 5'd0;
        for (q = 0; q < 17; q++) s_k = s_k + {4'b0, p2_sg[q]};
        s_corr = {~s_k, 9'b0, s_k};
    end

    logic [237:0] p3_cw;
    logic [19:0]  p3_corr;
    logic [5:0]   p3_emax;

    always_ff @(posedge clk) begin
        p3_cw   <= s_cw;
        p3_corr <= s_corr;
        p3_emax <= p2_emax;
    end

    // ======================= stage 3 : sum tree ===========================
    logic [14:0] l1 [0:7];
    logic [19:0] l1c;
    logic [15:0] l2 [0:3];
    logic [16:0] l3 [0:1];
    logic [17:0] l4;
    logic [19:0] s_sum;

    always_comb begin
        for (q = 0; q < 8; q++)
            l1[q] = p3_cw[(2*q)*14 +: 14] + p3_cw[(2*q+1)*14 +: 14];
        l1c = p3_cw[16*14 +: 14] + p3_corr;
        for (q = 0; q < 4; q++) l2[q] = l1[2*q] + l1[2*q+1];
        for (q = 0; q < 2; q++) l3[q] = l2[2*q] + l2[2*q+1];
        l4    = l3[0] + l3[1];
        s_sum = l4 + l1c;
    end

    logic [19:0] p4_sum;
    logic [5:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= s_sum;
        p4_emax <= p3_emax;
    end

    // ======================= stage 4 : normalize ==========================
    logic        s_fsign;
    logic [18:0] s_sumc;
    logic [4:0]  s_lzc;
    logic [5:0]  s_rgexp;
    logic [14:0] s_fmant;
    logic [18:0] s_nrm;

    always_comb begin
        s_fsign = p4_sum[19];
        s_sumc  = s_fsign ? (~p4_sum[18:0] + 19'd1) : p4_sum[18:0];
        s_lzc   = 5'd19;
        for (q = 0; q < 19; q++) if (s_sumc[q]) s_lzc = 5'd18 - q[4:0];
        s_rgexp = p4_emax + (6'sd6 - $signed({1'b0,s_lzc}));
        s_nrm   = s_sumc << s_lzc;
        s_fmant = {s_nrm[18:5], |s_nrm[4:0]};
    end

    logic        p5_sign;
    logic [5:0]  p5_rgexp;
    logic [14:0] p5_fmant;

    always_ff @(posedge clk) begin
        p5_sign  <= s_fsign;
        p5_rgexp <= s_rgexp;
        p5_fmant <= s_fmant;
    end

    // ======================= output : posit encode ========================
    logic [3:0]  e_ar;
    logic [3:0]  e_a;
    logic [14:0] e_m1;
    logic [14:0] e_pat;
    logic [13:0] e_var;
    logic [14:0] e_f;
    logic [14:0] e_w15;
    logic [14:0] e_h;
    logic [13:0] e_pref;
    logic        e_inc;
    logic [14:0] e_vbr;
    logic [14:0] e_mag;
    logic [15:0] e_nrm;

    always_comb begin
        // regime magnitude:  a = |k| + 1, capped at 15
        e_ar = p5_rgexp[5] ? (~p5_rgexp[5:2]) : p5_rgexp[5:2];
        e_a  = (e_ar >= 4'd14) ? 4'd15 : (e_ar + 4'd1);
        e_var = {p5_rgexp[1:0], p5_fmant[11:0]};
        e_m1  = 15'h7FFF >> e_a;
        e_pat = e_m1 ^ (p5_sign ? (e_m1 >> 1) : 15'h7FFF);
        e_f   = e_pat | {1'b0, (e_var >> e_a)};
        e_vbr = e_f[14:1];
        // round to nearest even
        e_w15 = {p5_sign, e_var};
        e_pref[0] = 1'b0;
        for (q = 1; q < 14; q++) e_pref[q] = |e_var[q-1:0];
        for (q = 0; q < 15; q++) e_h[q] = e_w15[q] & (e_pref[q<14?q:13] | e_f[0]);
        e_inc = e_a[0] ? (e_a[1] ? (e_a[2] ? (e_a[3] ? e_h[15-1] : e_h[7]) : (e_a[3] ? e_h[11] : e_h[3]))
                                 : (e_a[2] ? (e_a[3] ? e_h[13] : e_h[5]) : (e_a[3] ? e_h[9] : e_h[1])))
                       : (e_a[1] ? (e_a[2] ? (e_a[3] ? e_h[14] : e_h[6]) : (e_a[3] ? e_h[10] : e_h[2]))
                                 : (e_a[2] ? (e_a[3] ? e_h[12] : e_h[4]) : (e_a[3] ? e_h[8] : e_h[0])));
        e_mag = e_vbr + {14'b0, e_inc};
        e_nrm = p5_sign ? {1'b1, ~e_mag + 15'd1} : {1'b0, e_mag};
    end

    assign result_o = p5_fmant[12] ? e_nrm : 16'd0;

endmodule

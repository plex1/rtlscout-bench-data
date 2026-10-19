// ===========================================================================
// PDPU dot-16, 5-stage pipelined posit DPA unit.
// Bit-exact with the reference model (ALIGN_WIDTH=14 truncation, CARRY_W=5,
// SUM_WIDTH=19, wrapped exponent arithmetic and NaR handling).
//
//  1 : posit decode (direct regime extraction, no lzc / no 15-bit shifter)
//  2 : mantissa multiply (4x4) and exponent addition
//  3 : maximum-exponent tree
//  4 : alignment barrel shifters + sign folding (bit inversion + popcount)
//  5 : balanced adder tree (17 aligned words + sign correction)
//  6 : normalisation and posit encoding
// ===========================================================================
module pdpu_top_pipelined #(
    parameter int N = 16
)(
    input  logic         clk,
    input  logic [127:0] operands_a,
    input  logic [127:0] operands_b,
    input  logic [15:0]  acc,
    output logic [15:0]  result_o
);

    integer i, k;

    // ------------------------------------------------------------------
    // posit(8,2) decode -> {rg_exp[5:0], mant[3:0]}
    //   rg_exp = 4*regime_k + exponent      (6 bit two's complement)
    //   mant   = {implicit, frac[2:0]}
    // ------------------------------------------------------------------
    function automatic [9:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] h, rl, sh;
        logic [3:0] kk;
        logic       rs;
        va = p[7] ? (~p[6:0] + 7'd1) : p[6:0];
        rs = va[6];
        x  = rs ? ~va : va;
        h  = 3'd7;
        for (int t = 0; t < 7; t++) if (x[t]) h = t[2:0];
        rl = (h == 3'd7) ? 3'd6 : (3'd6 - h);
        sh = (h == 3'd7) ? 3'd6 : (3'd7 - h);
        kk = rs ? ({1'b0, rl} - 4'd1) : (4'd0 - {1'b0, rl});
        al = va << sh;
        dec8 = {kk, al[6:5], |p[6:0], al[4:2]};
    endfunction

    // posit(16,2) decode -> {rg_exp[6:0], mant[11:0]}
    function automatic [18:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  h, rl, sh, kk;
        logic        rs;
        va = p[15] ? (~p[14:0] + 15'd1) : p[14:0];
        rs = va[14];
        x  = rs ? ~va : va;
        h  = 4'd15;
        for (int t = 0; t < 15; t++) if (x[t]) h = t[3:0];
        rl = (h == 4'd15) ? 4'd14 : (4'd14 - h);
        sh = (h == 4'd15) ? 4'd14 : (4'd15 - h);
        kk = rs ? (rl - 4'd1) : (4'd0 - rl);
        al = va << sh;
        dec16 = {kk, al[14:13], |p[14:0], al[12:2]};
    endfunction

    function automatic [6:0] smax(input [6:0] a, input [6:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // ======================= 1 : decode =================================
    logic [16:0] s1_sg;
    logic [6:0]  s1_ep  [0:16];
    logic [13:0] s1_w   [0:16];

    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [9:0] da, db;
        logic [7:0] pr;
        assign da = dec8(operands_a[i*8 +: 8]);
        assign db = dec8(operands_b[i*8 +: 8]);
        assign s1_sg[i]   = operands_a[i*8+7] ^ operands_b[i*8+7];
        assign s1_ep[i]   = {{(6){da[9]}}, da[9:4]} + {{(6){db[9]}}, db[9:4]};
        assign pr = da[3:0] * db[3:0];
        assign s1_w[i] = {2'b0, pr, 4'b0};
    end
    endgenerate

    logic [18:0] dacc;
    assign dacc = dec16(acc);
    assign s1_sg[16] = acc[15];
    assign s1_ep[16] = {{(6){dacc[18]}}, dacc[18:12]};
    assign s1_w[16]  = {1'b0, dacc[11:0], 1'b0};

    logic [16:0] p1_sg;
    logic [6:0]  p1_ep [0:16];
    logic [13:0] p1_w  [0:16];

    always_ff @(posedge clk) begin
        p1_sg <= s1_sg;
        for (i = 0; i < 17; i++) begin
            p1_ep[i] <= s1_ep[i];
            p1_w[i]  <= s1_w[i];
        end
    end

    // ======================= 2 : exponent add / mul =====================
    logic [6:0]  s2_emax;
    logic [16:0] p2_sg;
    logic [6:0]  p2_ep [0:16];
    logic [13:0] p2_w  [0:16];

    always_comb begin
        s2_emax = smax(smax(smax(smax(smax(smax(smax(smax(p1_ep[0], p1_ep[1]),
                                                          p1_ep[2]), p1_ep[3]),
                                                    smax(p1_ep[4], p1_ep[5])),
                                               smax(smax(p1_ep[6], p1_ep[7]),
                                                    smax(p1_ep[8], p1_ep[9]))),
                                          smax(smax(smax(p1_ep[10], p1_ep[11]),
                                                    smax(p1_ep[12], p1_ep[13])),
                                               smax(smax(p1_ep[14], p1_ep[15]),
                                                    p1_ep[16])));
    end

    always_ff @(posedge clk) begin
        p2_sg <= p1_sg;
        for (i = 0; i < 17; i++) begin
            p2_ep[i] <= p1_ep[i];
            p2_w[i]  <= p1_w[i];
        end
    end

    // ======================= 3 : max exponent ==========================
    logic [6:0]  s3_emax;
    logic [16:0] p3_sg;
    logic [6:0]  p3_ep [0:16];
    logic [13:0] p3_w  [0:16];

    assign s3_emax = s2_emax;

    always_ff @(posedge clk) begin
        p3_sg <= p2_sg;
        for (i = 0; i < 17; i++) begin
            p3_ep[i] <= p2_ep[i];
            p3_w[i]  <= p2_w[i];
        end
    end

    // ======================= 4 : align =================================
    logic [13:0] s4_cw [0:16];
    logic [4:0]  s4_k;
    logic [6:0]  s4_emax;

    generate
    for (i = 0; i < 17; i++) begin : g_al
        logic [6:0]  df;
        logic [3:0]  sh;
        logic [13:0] ws;
        assign df = s3_emax - p3_ep[i];
        assign sh = (|df[6:4]) ? 4'd14 : df[3:0];
        assign ws = p3_w[i] >> sh;
        assign s4_cw[i] = ws ^ {14{p3_sg[i]}};
    end
    endgenerate

    assign s4_emax = s3_emax;

    always_comb begin
        logic [4:0] c;
        c = 5'd0;
        for (k = 0; k < 17; k++) c = c + {4'b0, p3_sg[k]};
        s4_k = c;
    end

    logic [13:0] p4_cw [0:16];
    logic [4:0]  p4_k;
    logic [6:0]  p4_emax;

    always_ff @(posedge clk) begin
        for (i = 0; i < 17; i++) p4_cw[i] <= s4_cw[i];
        p4_k    <= s4_k;
        p4_emax <= s4_emax;
    end

    // ======================= 5 : adder tree ============================
    logic [14:0] s5_a [0:7];
    logic [15:0] s5_b [0:3];
    logic [16:0] s5_c [0:1];
    logic [17:0] s5_d;
    logic [18:0] s5_s;
    logic [19:0] s5_tot;

    always_comb begin
        for (k = 0; k < 8; k++) s5_a[k] = {1'b0, p4_cw[2*k]} + {1'b0, p4_cw[2*k+1]};
        s5_b[0] = {1'b0, s5_a[0]} + {1'b0, s5_a[1]};
        s5_b[1] = {1'b0, s5_a[2]} + {1'b0, s5_a[3]};
        s5_b[2] = {1'b0, s5_a[4]} + {1'b0, s5_a[5]};
        s5_b[3] = {1'b0, s5_a[6]} + {1'b0, s5_a[7]};
        s5_c[0] = {1'b0, s5_b[0]} + {1'b0, s5_b[1]};
        s5_c[1] = {1'b0, s5_b[2]} + {1'b0, s5_b[3]};
        s5_d    = {1'b0, s5_c[0]} + {1'b0, s5_c[1]};
        s5_s    = {1'b0, s5_d} + {1'b0, p4_cw[16]};
        // total = s5_s + K - K*2^14  (mod 2^20)
        s5_tot  = {1'b0, s5_s} + {15'b0, p4_k} - {p4_k, 14'b0};
    end

    logic [19:0] p5_tot;
    logic [6:0]  p5_emax;

    always_ff @(posedge clk) begin
        p5_tot  <= s5_tot;
        p5_emax <= p4_emax;
    end

    // ======================= 6 : normalise + encode =====================
    logic        s6_sign;
    logic [18:0] s6_mag;
    logic [4:0]  s6_lzc;
    logic [6:0]  s6_rgexp;
    logic [13:0] s6_fmant;
    logic        s6_nz;

    always_comb begin
        logic [18:0] sm;
        logic [18:0] sn;
        logic [6:0]  adj;
        s6_sign = p5_tot[19];
        sm      = s6_sign ? (~p5_tot[18:0] + 19'd1) : p5_tot[18:0];
        s6_lzc  = 5'd19;
        for (k = 0; k < 19; k++) if (sm[k]) s6_lzc = 5'd18 - k[4:0];
        adj     = (s6_lzc == 5'd19) ? 7'sd0 : (7'sd6 - $signed({2'b0, s6_lzc}));
        s6_rgexp = p5_emax + adj;
        sn      = sm << s6_lzc;
        s6_nz   = sn[18];
        s6_fmant = {sn[17:6], |sn[5:0]};
    end

    // ---- posit encoder ----
    logic        s6_kneg;
    logic [4:0]  s6_rb;
    logic [4:0]  s6_sm1;
    logic [29:0] s6_ext;
    logic [29:0] s6_sh;
    logic [14:0] s6_var;
    logic [15:0] s6_msk;
    logic        s6_stk;
    logic        s6_rbit;
    logic        s6_rval;
    logic [14:0] s6_mag2;
    logic [15:0] s6_out;

    always_comb begin
        logic [4:0] kk;
        logic [14:0] rg;
        kk       = s6_rgexp[6:2];
        s6_kneg  = s6_rgexp[6];
        s6_rb    = s6_kneg ? (~kk + 5'd2) : (kk + 5'd2);
        rg       = (s6_rb >= 5'd16) ? 5'd16 : s6_rb;
        s6_ext   = {s6_kneg ? 15'h0001 : 15'h7FFE, s6_rgexp[1:0], s6_fmant};
        s6_sm1   = rg - 5'd1;
        s6_sh    = s6_ext >> s6_sm1;
        s6_var   = s6_sh[15:1];
        s6_rbit  = s6_sh[0];
        s6_msk   = (16'd1 << s6_sm1) - 16'd1;
        s6_stk   = |(s6_ext[15:0] & s6_msk);
        s6_rval  = s6_rbit & (s6_stk | s6_var[0]);
        s6_mag2  = s6_var + {14'b0, s6_rval};
        s6_out   = s6_sign ? {1'b1, ~s6_mag2 + 15'd1} : {1'b0, s6_mag2};
    end

    assign result_o = s6_nz ? s6_out : 16'd0;

endmodule

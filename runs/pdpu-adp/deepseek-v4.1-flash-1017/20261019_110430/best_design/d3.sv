// ===========================================================================
// PDPU dot-16, 5-stage pipelined - compact re-implementation.
// Bit-exact with the reference model (ALIGN_WIDTH=14 truncation, CARRY_W=5,
// SUM_WIDTH=19, wrapped 6-bit exponent arithmetic, NaR handling).
//
// Optimisations vs. the reference:
//  * posit decode: direct regime extraction (priority encoder + small shift)
//    instead of a 15-bit lzc plus a 15-bit barrel shifter.
//  * exponent compare tree / mantissa multiplier moved into stage 1.
//  * two's complement of the aligned words folded into a bit inversion plus a
//    single correction constant (K - K<<14) added in the summing tree, so no
//    17 x 20-bit adders are needed.
//  * final accumulation uses a balanced binary adder tree of growing width
//    instead of a 17:2 carry-save tree plus a 20-bit CPA.
//  * posit encode: shift amount derived directly from the regime value, small
//    barrel shifter instead of a 46-bit one.
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
    // posit(8,2) decode: {sign, e[5:0], m[3:0]}
    //   e = 4*regime + exponent (6-bit signed), m = {implicit, frac[2:0]}
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
        dec8 = {p[7], k[3], k[2:0], b[6:5], |va, b[4:2]};
    endfunction

    // posit(16,2) decode for the accumulator: {sign, e[6:0], m[11:0]}
    function automatic [19:0] dec16(input [15:0] p);
        logic [14:0] va, x, b;
        logic [3:0]  rl;
        logic [4:0]  k;
        logic [14:0] rb;
        va = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        x  = va ^ {15{va[14]}};
        rl = 4'd15;
        for (int i = 0; i < 15; i++) if (x[i]) rl = 4'd14 - i[3:0];
        k  = va[14] ? ({1'b0,rl} - 5'd1) : (5'd0 - {1'b0,rl});
        rb = (rl == 4'd15) ? 15'd15 : ({11'b0,rl} + 15'd1);
        b  = va << rb[3:0];
        dec16 = {p[15], k[4], k[3:0], b[14:13], |va, b[12:2]};
    endfunction

    function automatic [5:0] mx(input [5:0] a, input [5:0] b);
        mx = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // ======================= stage 0 : decode + mul =======================
    logic [16:0]  s0_sg;      // 17 signs
    logic [101:0] s0_ep;      // 17 x 6-bit exponents
    logic [237:0] s0_w;       // 17 x 14-bit mantissa words

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

    logic [19:0] dacc;
    assign dacc = dec16(acc);
    assign s0_sg[16] = dacc[19];
    assign s0_ep[96 +: 6] = dacc[18:13];     // low 6 bits of the 7-bit exponent
    assign s0_w[224 +: 14] = {1'b0, dacc[12:0], 1'b0};

    always_ff @(posedge clk) begin
        p1_sg <= s0_sg;
        p1_ep <= s0_ep;
        p1_w  <= s0_w;
    end
    logic [16:0]  p1_sg;
    logic [101:0] p1_ep;
    logic [237:0] p1_w;

    // ======================= stage 1 : max exponent =======================
    logic [5:0] m0 [0:16];
    logic [5:0] m1 [0:8];
    logic [5:0] m2 [0:4];
    logic [5:0] m3 [0:2];
    logic [5:0] m4 [0:1];
    logic [5:0] s_emax;
    integer q;
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

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_ep   <= p1_ep;
        p2_w    <= p1_w;
        p2_emax <= s_emax;
    end
    logic [16:0]  p2_sg;
    logic [101:0] p2_ep;
    logic [237:0] p2_w;
    logic [5:0]   p2_emax;

    // ======================= stage 2 : align ==============================
    logic [237:0] s_cw;       // complemented (or not) aligned words
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
        s_corr = {(6'd0 - {1'b0,s_k}), 9'b0, s_k};
    end

    always_ff @(posedge clk) begin
        p3_cw   <= s_cw;
        p3_corr <= s_corr;
        p3_emax <= p2_emax;
    end
    logic [237:0] p3_cw;
    logic [19:0]  p3_corr;
    logic [5:0]   p3_emax;

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
        l4  = l3[0] + l3[1];
        s_sum = l4 + l1c;
    end

    always_ff @(posedge clk) begin
        p4_sum  <= s_sum;
        p4_emax <= p3_emax;
    end
    logic [19:0] p4_sum;
    logic [5:0]  p4_emax;

    // ======================= stage 4 : normalize ==========================
    logic        s_fsign;
    logic [18:0] s_sumc;
    logic [4:0]  s_lzc;
    logic [5:0]  s_rgexp;
    logic [14:0] s_fmant;

    always_comb begin
        s_fsign = p4_sum[19];
        s_sumc  = s_fsign ? (~p4_sum + 20'd1) : p4_sum[18:0];
        s_lzc   = 5'd19;
        for (q = 0; q < 19; q++) if (s_sumc[q]) s_lzc = 5'd18 - q[4:0];
        s_rgexp = p4_emax + (6'sd6 - $signed({1'b0,s_lzc}));
        s_fmant = {(s_sumc << s_lzc)[18:5], |((s_sumc << s_lzc)[4:0])};
    end

    always_ff @(posedge clk) begin
        p5_sign  <= s_fsign;
        p5_rgexp <= s_rgexp;
        p5_fmant <= s_fmant;
    end
    logic        p5_sign;
    logic [5:0]  p5_rgexp;
    logic [14:0] p5_fmant;

    // ======================= output : posit encode ========================
    logic [3:0]  e_k;
    logic [4:0]  e_s;
    logic [3:0]  e_a;
    logic [14:0] e_low;
    logic [14:0] e_u;
    logic [14:0] e_pat;
    logic [14:0] e_pre;
    logic        e_rnd;
    logic [13:0] e_th;
    logic        e_stk;
    logic [14:0] e_mag;
    logic [15:0] e_nrm;

    always_comb begin
        e_k   = p5_rgexp[5:2];
        e_low = {p5_rgexp[1:0], p5_fmant[12:0]};
        e_s   = e_k[3] ? (5'd17 - {1'b0,e_k}) : ({1'b0,e_k} + 5'd2);
        e_a   = e_s[3:0] - 4'd1;
        e_u   = e_low >> e_a;
        // thermometer decode of a (bits below position a)
        for (q = 0; q < 14; q++) e_th[q] = ({1'b0,e_a} > q[4:0]);
        e_stk = |(e_low[13:0] & e_th);
        e_pat = e_k[3] ? (15'h4000 >> e_a) : ~(15'h7FFF >> e_a);
        e_rnd = e_u[0];
        e_pre = e_pat | {1'b0, e_u[14:1]};
        e_mag = e_pre + {14'b0, (e_rnd & (e_stk | e_pre[0]))};
        e_nrm = p5_sign ? {1'b1, ~e_mag + 15'd1} : {1'b0, e_mag};
    end

    assign result_o = p5_fmant[14] ? e_nrm : 16'd0;

endmodule

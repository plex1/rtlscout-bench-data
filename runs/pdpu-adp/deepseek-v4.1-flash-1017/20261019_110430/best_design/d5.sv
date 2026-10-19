// ===========================================================================
// PDPU dot-16, 5-stage pipelined -- compact re-implementation.
// Bit-exact with the reference model.
//
//  stage 1 : posit decode (direct regime extraction) + mantissa multiply
//  stage 2 : maximum exponent tree
//  stage 3 : alignment shifts, bit complement, sign count
//  stage 4 : balanced adder tree (17 words + sign correction)
//  stage 5 : normalisation (lzc + shift)
//  output  : posit encoding (regime pattern + exponent + rounded fraction)
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

    localparam int EW = 7;      // exponent width

    // ------------------------------------------------------------------
    // posit(8,2) decode -> {sign, rg_exp[5:0], mant[3:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x, b;
        logic [2:0] rl, rb;
        logic [3:0] k;
        va = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        x  = va ^ {7{va[6]}};
        rl = 3'd7;
        for (int i = 0; i < 7; i++) if (x[i]) rl = 3'd6 - i[2:0];
        k  = va[6] ? ({1'b0, rl} - 4'd1) : (4'd0 - {1'b0, rl});
        rb = (rl == 3'd7) ? 3'd7 : (rl + 3'd1);
        b  = va << rb;
        dec8 = {p[7], k[3:0], b[6:5], |va, b[4:2]};
    endfunction

    // ------------------------------------------------------------------
    // posit(16,2) decode -> {sign, rg_exp[5:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [18:0] dec16(input [15:0] p);
        logic [14:0] va, x, b;
        logic [3:0]  rl, k;
        logic [14:0] rb;
        va = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        x  = va ^ {15{va[14]}};
        rl = 4'd15;
        for (int i = 0; i < 15; i++) if (x[i]) rl = 4'd14 - i[3:0];
        k  = va[14] ? (rl - 4'd1) : (4'd0 - rl);
        rb = (rl == 4'd15) ? 4'd15 : ({11'b0, rl} + 15'd1);
        b  = va << rb[3:0];
        dec16 = {p[15], k[3:0], b[14:13], |va, b[12:2]};
    endfunction

    function automatic [EW-1:0] mx(input [EW-1:0] a, input [EW-1:0] b);
        mx = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    integer q;

    // ======================= stage 1 : decode + multiply =================
    logic [16:0]      s1_sg;
    logic [EW*17-1:0] s1_ep;
    logic [13:0]      s1_w [0:16];

    genvar i;
    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [10:0] da, db;
        logic [7:0]  pr;
        assign da = dec8(operands_a[i*8 +: 8]);
        assign db = dec8(operands_b[i*8 +: 8]);
        assign s1_sg[i]         = da[10] ^ db[10];
        assign s1_ep[i*EW +: EW] = $signed(da[9:4]) + $signed(db[9:4]);
        assign pr = da[3:0] * db[3:0];
        assign s1_w[i] = {2'b0, pr, 4'b0};
    end
    endgenerate

    logic [18:0] dacc;
    assign dacc = dec16(acc);
    assign s1_sg[16]           = dacc[18];
    assign s1_ep[16*EW +: EW]  = {{(EW-6){dacc[17]}}, dacc[17:12]};
    assign s1_w[16]            = {1'b0, dacc[11:0], 1'b0};

    logic [16:0]      p1_sg;
    logic [EW*17-1:0] p1_ep;
    logic [13:0]      p1_w [0:16];

    always_ff @(posedge clk) begin
        p1_sg <= s1_sg;
        p1_ep <= s1_ep;
        for (int j = 0; j < 17; j++) p1_w[j] <= s1_w[j];
    end

    // ======================= stage 2 : max exponent ======================
    logic [EW-1:0] s2_e [0:16];
    logic [EW-1:0] s2_a [0:7];
    logic [EW-1:0] s2_b [0:3];
    logic [EW-1:0] s2_c [0:1];
    logic [EW-1:0] s2_emax;

    always_comb begin
        for (q = 0; q < 17; q++) s2_e[q] = p1_ep[q*EW +: EW];
        for (q = 0; q < 8; q++)  s2_a[q] = mx(s2_e[2*q], s2_e[2*q+1]);
        s2_b[0] = mx(s2_a[0], s2_a[1]);
        s2_b[1] = mx(s2_a[2], s2_a[3]);
        s2_b[2] = mx(s2_a[4], s2_a[5]);
        s2_b[3] = mx(s2_a[6], s2_a[7]);
        s2_c[0] = mx(s2_b[0], s2_b[1]);
        s2_c[1] = mx(s2_b[2], s2_b[3]);
        s2_emax = mx(mx(s2_c[0], s2_c[1]), s2_e[16]);
    end

    logic [EW-1:0] p2_emax;
    logic [16:0]   p2_sg;
    logic [EW*17-1:0] p2_ep;
    logic [13:0]   p2_w [0:16];

    always_ff @(posedge clk) begin
        p2_emax <= s2_emax;
        p2_sg   <= p1_sg;
        p2_ep   <= p1_ep;
        for (int j = 0; j < 17; j++) p2_w[j] <= p1_w[j];
    end

    // ======================= stage 3 : align =============================
    logic [13:0] s3_cw [0:16];
    logic [4:0]  s3_k;

    generate
    for (i = 0; i < 17; i++) begin : g_al
        logic [EW-1:0] df;
        logic [3:0]    sh;
        logic [13:0]   ws;
        assign df = p2_emax - p2_ep[i*EW +: EW];
        assign sh = (|df[EW-1:4]) ? 4'd14 : df[3:0];
        assign ws = p2_w[i] >> sh;
        assign s3_cw[i] = ws ^ {14{p2_sg[i]}};
    end
    endgenerate

    always_comb begin
        logic [1:0] c1 [0:7];
        logic [2:0] c2 [0:3];
        logic [3:0] c3 [0:1];
        logic [4:0] c4, c5;
        for (q = 0; q < 8; q++) c1[q] = {1'b0,p2_sg[2*q]} + {1'b0,p2_sg[2*q+1]};
        c2[0] = {1'b0,c1[0]} + {1'b0,c1[1]};
        c2[1] = {1'b0,c1[2]} + {1'b0,c1[3]};
        c2[2] = {1'b0,c1[4]} + {1'b0,c1[5]};
        c2[3] = {1'b0,c1[6]} + {1'b0,c1[7]};
        c3[0] = {1'b0,c2[0]} + {1'b0,c2[1]};
        c3[1] = {1'b0,c2[2]} + {1'b0,c2[3]};
        c4 = {1'b0,c3[0]} + {1'b0,c3[1]};
        c5 = c4 + {4'b0, p2_sg[16]};
        s3_k = c5;
    end

    logic [13:0] p3_cw [0:16];
    logic [4:0]  p3_k;
    logic [EW-1:0] p3_emax;

    always_ff @(posedge clk) begin
        for (int j = 0; j < 17; j++) p3_cw[j] <= s3_cw[j];
        p3_k    <= s3_k;
        p3_emax <= p2_emax;
    end

    // ======================= stage 4 : adder tree ========================
    logic [14:0] s4_a [0:7];
    logic [19:0] s4_ka;
    logic [15:0] s4_b [0:3];
    logic [16:0] s4_c [0:1];
    logic [17:0] s4_d;
    logic [19:0] s4_sum;

    always_comb begin
        for (q = 0; q < 8; q++) s4_a[q] = {1'b0,p3_cw[2*q]} + {1'b0,p3_cw[2*q+1]};
        s4_ka = {1'b0, p3_cw[16]} + {{(20-5-14){1'b0}}, {~p3_k, 9'b0, p3_k}};
        for (q = 0; q < 4; q++) s4_b[q] = {1'b0,s4_a[2*q]} + {1'b0,s4_a[2*q+1]};
        s4_c[0] = {1'b0,s4_b[0]} + {1'b0,s4_b[1]};
        s4_c[1] = {1'b0,s4_b[2]} + {1'b0,s4_b[3]};
        s4_d = {1'b0,s4_c[0]} + {1'b0,s4_c[1]};
        s4_sum = s4_ka + {2'b0, s4_d};
    end

    logic [19:0] p4_sum;
    logic [EW-1:0] p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= s4_sum;
        p4_emax <= p3_emax;
    end

    // ======================= stage 5 : normalise =========================
    logic        s5_sign;
    logic [18:0] s5_sumc;
    logic [4:0]  s5_lzc;
    logic [5:0]  s5_sh;
    logic [18:0] s5_nrm;
    logic [13:0] s5_fmant;
    logic [EW-1:0] s5_rgexp;

    always_comb begin
        logic [4:0] lz;
        logic [18:0] shv;
        s5_sign = p4_sum[19];
        s5_sumc = s5_sign ? (~p4_sum[18:0] + 19'd1) : p4_sum[18:0];
        lz = 5'd19;
        for (q = 0; q < 19; q++) if (s5_sumc[q]) lz = 5'd18 - q[4:0];
        s5_lzc = lz;
        s5_sh  = {1'b0, lz} - 6'd6;             // lzc - 6
        shv    = (lz >= 5'd6) ? (s5_sumc >> (lz - 5'd6)) : (s5_sumc << (5'd6 - lz));
        s5_nrm = shv;
        s5_fmant = {shv[17:5], |shv[4:0]};
        s5_rgexp = p4_emax + (6'sd6 - $signed({1'b0, lz}));
    end

    logic        p5_sign;
    logic [EW-1:0] p5_rgexp;
    logic [13:0] p5_fmant;

    always_ff @(posedge clk) begin
        p5_sign  <= s5_sign;
        p5_rgexp <= s5_rgexp;
        p5_fmant <= s5_fmant;
    end

    // ======================= output : posit encode =======================
    logic [3:0]  e_kmag;
    logic [3:0]  e_a;
    logic [13:0] e_low;
    logic [14:0] e_pat;
    logic [14:0] e_fld;
    logic [14:0] e_mag;
    logic        e_stk, e_rnd, e_inc;
    logic [15:0] e_nrm;

    always_comb begin
        logic [14:0] th;
        logic [14:0] below;
        logic [14:0] shf;
        e_kmag = p5_rgexp[5] ? (~p5_rgexp[5:2] + 4'd1) : p5_rgexp[5:2];
        // a = regime field length:  k>=0 -> k+2 ; k<0 -> |k|+1
        e_a = p5_rgexp[5] ? ({1'b0, e_kmag} - 4'd1) : ({1'b0, e_kmag} + 4'd2);
        // low 15 bits: exponent then fraction
        e_low = {p5_rgexp[1:0], p5_fmant[12:0]};
        // pattern
        if (p5_rgexp[5])
            e_pat = 15'h4000 >> (e_a - 4'd1);
        else
            e_pat = 15'h7FFF << (4'd16 - e_a);
        shf = e_low >> e_a;
        e_fld = e_pat | shf;
        // rounding: round bit = e_low[a-1], sticky = |e_low[a-2:0]
        for (q = 0; q < 15; q++) th[q] = ({1'b0, e_a} > q[4:0] + 5'd1);
        below = e_low & th;
        e_stk = |below;
        e_rnd = e_low >> (e_a - 4'd1);
        e_inc = e_rnd & (e_stk | e_fld[0]);
        e_mag = e_fld + {14'b0, e_inc};
        e_nrm = p5_sign ? {1'b1, ~e_mag + 15'd1} : {1'b0, e_mag};
    end

    assign result_o = p5_fmant[13] ? e_nrm : 16'd0;

endmodule

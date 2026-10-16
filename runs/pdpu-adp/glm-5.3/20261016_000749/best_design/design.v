// Pipelined posit(8,2) x16 dot-product-accumulate -> posit(16,2)
// Bit-exact re-implementation of pdpu_top_pipelined (5 flop stages).
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

    // ---------------- helpers ----------------
    function automatic [2:0] msb7(input [6:0] v);
        msb7 = v[6] ? 3'd0 : v[5] ? 3'd1 : v[4] ? 3'd2 : v[3] ? 3'd3 :
               v[2] ? 3'd4 : v[1] ? 3'd5 : v[0] ? 3'd6 : 3'd7;
    endfunction

    function automatic [4:0] msb15(input [14:0] v);
        msb15 = v[14] ? 5'd0 : v[13] ? 5'd1 : v[12] ? 5'd2 : v[11] ? 5'd3 :
                v[10] ? 5'd4 : v[9]  ? 5'd5 : v[8]  ? 5'd6 : v[7]  ? 5'd7 :
                v[6]  ? 5'd8 : v[5]  ? 5'd9 : v[4]  ? 5'd10 : v[3] ? 5'd11 :
                v[2]  ? 5'd12 : v[1] ? 5'd13 : v[0] ? 5'd14 : 5'd15;
    endfunction

    function automatic [4:0] msb19(input [18:0] v);
        msb19 = v[18] ? 5'd0 : v[17] ? 5'd1 : v[16] ? 5'd2 : v[15] ? 5'd3 :
                v[14] ? 5'd4 : v[13] ? 5'd5 : v[12] ? 5'd6 : v[11] ? 5'd7 :
                v[10] ? 5'd8 : v[9]  ? 5'd9 : v[8]  ? 5'd10 : v[7] ? 5'd11 :
                v[6]  ? 5'd12 : v[5] ? 5'd13 : v[4] ? 5'd14 : v[3] ? 5'd15 :
                v[2]  ? 5'd16 : v[1] ? 5'd17 : v[0] ? 5'd18 : 5'd19;
    endfunction

    // posit(8,2) decode -> {sign, rg_exp[5:0] signed, mant[3:0]}
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] mag, t, sh;
        logic [2:0] mi, l;
        logic [5:0] k, re;
        mag = p[7] ? (~p[6:0] + 1'b1) : p[6:0];
        t   = mag[6] ? ~mag : mag;
        l   = msb7(t);
        sh  = mag << ({1'b0, l} + 3'd1);
        k   = mag[6] ? ({3'b000, l} - 6'd1) : (6'd1 - {3'b000, l} - 6'd1);
        re  = {k[3:0], sh[6:5]};
        dec8 = {p[7], re, |mag, sh[4:2]};
    endfunction

    // posit(16,2) decode -> {sign, rg_exp[7:0] signed, mant[11:0]}
    function automatic [20:0] dec16(input [15:0] p);
        logic [14:0] mag, t, sh;
        logic [4:0]  mi, l;
        logic [6:0]  k;
        logic [6:0]  tmp;
        logic [7:0]  re;
        mag = p[15] ? (~p[14:0] + 1'b1) : p[14:0];
        t   = mag[14] ? ~mag : mag;
        l   = msb15(t);
        sh  = mag << ({1'b0, l} + 5'd1);
        k   = mag[14] ? ({2'b00, l} - 7'd1) : (7'd1 - {2'b00, l} - 7'd1);
        tmp = {k[4:0], sh[14:13]};
        re  = {tmp[6], tmp};
        dec16 = {p[15], re, |mag, sh[12:2]};
    endfunction

    function automatic [7:0] smax8(input [7:0] a, input [7:0] b);
        smax8 = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    function automatic [13:0] rs14(input [13:0] v, input [3:0] s);
        logic [13:0] a, b, c;
        a = s[3] ? {8'b0, v[13:8]} : v;
        b = s[2] ? {4'b0, a[13:4]} : a;
        c = s[1] ? {2'b0, b[13:2]} : b;
        rs14 = s[0] ? {1'b0, c[13:1]} : c;
    endfunction

    function automatic [18:0] ls19(input [18:0] v, input [4:0] s);
        logic [18:0] a, b, c, d;
        a = s[4] ? {v[2:0], 16'b0} : v;
        b = s[3] ? {a[10:0], 8'b0} : a;
        c = s[2] ? {b[14:0], 4'b0} : b;
        d = s[1] ? {c[16:0], 2'b0} : c;
        ls19 = s[0] ? {d[17:0], 1'b0} : d;
    endfunction

    function automatic [30:0] ls31(input [30:0] v, input [4:0] s);
        logic [30:0] a, b, c, d;
        a = s[4] ? {v[14:0], 16'b0} : v;
        b = s[3] ? {a[22:0], 8'b0} : a;
        c = s[2] ? {b[26:0], 4'b0} : b;
        d = s[1] ? {c[28:0], 2'b0} : c;
        ls31 = s[0] ? {d[29:0], 1'b0} : d;
    endfunction

    function automatic [4:0] popc17(input [16:0] v);
        logic [1:0] q0, q1, q2, q3, q4, q5, q6, q7;
        logic [2:0] r0, r1, r2, r3;
        logic [3:0] s0, s1;
        q0 = {1'b0, v[0]}  + {1'b0, v[1]};
        q1 = {1'b0, v[2]}  + {1'b0, v[3]};
        q2 = {1'b0, v[4]}  + {1'b0, v[5]};
        q3 = {1'b0, v[6]}  + {1'b0, v[7]};
        q4 = {1'b0, v[8]}  + {1'b0, v[9]};
        q5 = {1'b0, v[10]} + {1'b0, v[11]};
        q6 = {1'b0, v[12]} + {1'b0, v[13]};
        q7 = {1'b0, v[14]} + {1'b0, v[15]};
        r0 = {1'b0, q0} + {1'b0, q1};
        r1 = {1'b0, q2} + {1'b0, q3};
        r2 = {1'b0, q4} + {1'b0, q5};
        r3 = {1'b0, q6} + {1'b0, q7};
        s0 = {1'b0, r0} + {1'b0, r1};
        s1 = {1'b0, r2} + {1'b0, r3};
        popc17 = {1'b0, s0} + {1'b0, s1} + {4'b0, v[16]};
    endfunction

    // ================= A: decode =================
    logic [16:0]     signs_c;
    logic [17*8-1:0] ec;
    logic [16*4-1:0] ma_v, mb_v;
    logic [11:0]     macc_v;

    generate
    for (genvar gi = 0; gi < N; gi++) begin : dec
        wire [10:0] da = dec8(operands_a[gi*n_i +: n_i]);
        wire [10:0] db = dec8(operands_b[gi*n_i +: n_i]);
        assign ma_v[gi*4 +: 4] = da[3:0];
        assign mb_v[gi*4 +: 4] = db[3:0];
        assign signs_c[gi]     = da[10] ^ db[10];
        assign ec[gi*8 +: 8]   = {{2{da[9]}}, da[9:4]} + {{2{db[9]}}, db[9:4]};
    end
    endgenerate

    wire [20:0] dac = dec16(acc);
    assign signs_c[16]   = dac[20];
    assign macc_v        = dac[11:0];
    assign ec[16*8 +: 8] = dac[19:12];

    // -------- pipe1 --------
    logic [16:0]     p1_signs;
    logic [17*8-1:0] p1_e;
    logic [16*4-1:0] p1_ma, p1_mb;
    logic [11:0]     p1_macc;

    always_ff @(posedge clk) begin
        p1_signs <= signs_c;
        p1_e     <= ec;
        p1_ma    <= ma_v;
        p1_mb    <= mb_v;
        p1_macc  <= macc_v;
    end

    // ================= B: multiply, exp max, shift amounts =================
    logic [16*8-1:0] prod_v;
    logic [17*4-1:0] shf_v;
    logic [7:0]      emax;
    logic [4:0]      pc;

    generate
    for (genvar gi = 0; gi < N; gi++) begin : mul
        assign prod_v[gi*8 +: 8] = p1_ma[gi*4 +: 4] * p1_mb[gi*4 +: 4];
    end
    endgenerate

    logic [8*8-1:0] mx1;
    logic [4*8-1:0] mx2;
    logic [2*8-1:0] mx3;
    logic [7:0]     mx4;
    generate
    for (genvar g = 0; g < 8; g++) assign mx1[g*8 +: 8] = smax8(p1_e[g*16 +: 8], p1_e[g*16+8 +: 8]);
    for (genvar g = 0; g < 4; g++) assign mx2[g*8 +: 8] = smax8(mx1[g*16 +: 8], mx1[g*16+8 +: 8]);
    for (genvar g = 0; g < 2; g++) assign mx3[g*8 +: 8] = smax8(mx2[g*16 +: 8], mx2[g*16+8 +: 8]);
    endgenerate
    assign mx4  = smax8(mx3[0 +: 8], mx3[8 +: 8]);
    assign emax = smax8(mx4, p1_e[16*8 +: 8]);

    generate
    for (genvar g = 0; g < 17; g++) begin : dif
        wire [7:0] d = emax - p1_e[g*8 +: 8];
        assign shf_v[g*4 +: 4] = (|d[7:4]) ? 4'd14 : d[3:0];
    end
    endgenerate

    assign pc = popc17(p1_signs);

    // -------- pipe2: multiply, exp max, shift amounts --------
    logic [16:0]     p2_signs;
    logic [16*8-1:0] p2_prod;
    logic [17*4-1:0] p2_shf;
    logic [7:0]      p2_emax;
    logic [4:0]      p2_pc;
    logic [13:0]     p2_accal;

    always_ff @(posedge clk) begin
        p2_signs <= p1_signs;
        p2_prod  <= prod_v;
        p2_shf   <= shf_v;
        p2_emax  <= emax;
        p2_pc    <= pc;
        p2_accal <= {p1_macc, 1'b0};
    end

    // ================= C: align, sign fold, CSA levels 1-2 =================
    logic [17*20-1:0] csa_in;
    generate
    for (genvar g = 0; g < 16; g++) begin : alnp
        wire [13:0] sh = rs14({p2_prod[g*8 +: 8], 6'b0}, p2_shf[g*4 +: 4]);
        assign csa_in[g*20 +: 20] = {{6{p2_signs[g]}}, sh ^ {14{p2_signs[g]}}};
    end
    endgenerate
    wire [13:0] sh16 = rs14(p2_accal, p2_shf[16*4 +: 4]);
    assign csa_in[16*20 +: 20] = {{6{p2_signs[16]}}, sh16 ^ {14{p2_signs[16]}}};

    // reduction: 17 -> 12 -> 8
    logic [12*20-1:0] L1;
    logic [8*20-1:0]  L2;
    generate
        for (genvar g = 0; g < 5; g++) begin : r1
            wire [19:0] x = csa_in[g*60 +: 20];
            wire [19:0] y = csa_in[g*60+20 +: 20];
            wire [19:0] z = csa_in[g*60+40 +: 20];
            assign L1[g*40 +: 20]     = x ^ y ^ z;
            assign L1[g*40+20 +: 20]  = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
        assign L1[10*20 +: 20] = csa_in[15*20 +: 20];
        assign L1[11*20 +: 20] = csa_in[16*20 +: 20];

        for (genvar g = 0; g < 4; g++) begin : r2
            wire [19:0] x = L1[g*60 +: 20];
            wire [19:0] y = L1[g*60+20 +: 20];
            wire [19:0] z = L1[g*60+40 +: 20];
            assign L2[g*40 +: 20]    = x ^ y ^ z;
            assign L2[g*40+20 +: 20] = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
    endgenerate

    // -------- pipe3 --------
    logic [8*20-1:0] p3_L;
    logic [7:0]      p3_emax;
    logic [4:0]      p3_pc;
    always_ff @(posedge clk) begin
        p3_L    <= L2;
        p3_emax <= p2_emax;
        p3_pc   <= p2_pc;
    end

    // ================= D: CSA levels 3-6 + final add =================
    logic [6*20-1:0] L3;
    logic [4*20-1:0] L4;
    logic [3*20-1:0] L5;
    logic [19:0]     csa_s, csa_c;
    generate
        for (genvar g = 0; g < 2; g++) begin : r3
            wire [19:0] x = p3_L[g*60 +: 20];
            wire [19:0] y = p3_L[g*60+20 +: 20];
            wire [19:0] z = p3_L[g*60+40 +: 20];
            assign L3[g*40 +: 20]    = x ^ y ^ z;
            assign L3[g*40+20 +: 20] = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
        assign L3[4*20 +: 20] = p3_L[6*20 +: 20];
        assign L3[5*20 +: 20] = p3_L[7*20 +: 20];

        for (genvar g = 0; g < 2; g++) begin : r4
            wire [19:0] x = L3[g*60 +: 20];
            wire [19:0] y = L3[g*60+20 +: 20];
            wire [19:0] z = L3[g*60+40 +: 20];
            assign L4[g*40 +: 20]    = x ^ y ^ z;
            assign L4[g*40+20 +: 20] = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
        begin : r5
            wire [19:0] x = L4[0 +: 20];
            wire [19:0] y = L4[20 +: 20];
            wire [19:0] z = L4[40 +: 20];
            assign L5[0 +: 20]    = x ^ y ^ z;
            assign L5[20 +: 20]   = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
        assign L5[2*20 +: 20] = L4[3*20 +: 20];
        begin : r6
            wire [19:0] x = L5[0 +: 20];
            wire [19:0] y = L5[20 +: 20];
            wire [19:0] z = L5[40 +: 20];
            assign csa_s = x ^ y ^ z;
            assign csa_c = {x[18:0] & y[18:0] | x[18:0] & z[18:0] | y[18:0] & z[18:0], 1'b0};
        end
    endgenerate

    wire [19:0] sum_result = csa_s + csa_c + {15'd0, p3_pc};

    // -------- pipe4 --------
    logic [19:0] p4_sum;
    logic [7:0]  p4_emax;
    always_ff @(posedge clk) begin
        p4_sum  <= sum_result;
        p4_emax <= p3_emax;
    end

    // ================= E: normalize =================
    wire        final_sign = p4_sum[19];
    wire [19:0] sum_neg    = ~p4_sum + 20'd1;
    wire [18:0] sum_c      = final_sign ? sum_neg[18:0] : p4_sum[18:0];
    wire [4:0]  lzc        = msb19(sum_c);
    wire [18:0] sum_norm   = ls19(sum_c, lzc);
    wire [7:0]  final_exp  = p4_emax + 8'sd6 - {3'b000, lzc};
    wire        sticky     = |sum_norm[5:0];
    wire [13:0] final_mant = {sum_norm[18:6], sticky};

    // -------- pipe5 --------
    logic        p5_sign;
    logic [7:0]  p5_exp;
    logic [13:0] p5_mant;
    always_ff @(posedge clk) begin
        p5_sign <= final_sign;
        p5_exp  <= final_exp;
        p5_mant <= final_mant;
    end

    // ================= F: encode (combinational) =================
    wire        not_zero  = p5_mant[13];
    wire [5:0]  regime_k  = p5_exp[7:2];
    wire [1:0]  eexp      = p5_exp[1:0];
    wire        sign_k    = p5_exp[7];
    wire [14:0] regime    = sign_k ? 15'h0001 : 15'h7FFE;
    wire [5:0]  rbits     = sign_k ? (~regime_k + 6'd2) : (regime_k + 6'd2);
    wire [4:0]  sft       = (rbits >= 6'd16) ? 5'd16 : rbits[4:0];
    wire [29:0] rgm       = {regime, eexp, p5_mant[12:0]};
    wire [30:0] Y         = ls31({1'b0, rgm}, 5'd16 - sft);
    wire [14:0] vbr       = Y[30:16];
    wire        rnd       = Y[15] & ((|Y[14:0]) | vbr[0]);
    wire [14:0] valrnd    = vbr + {14'b0, rnd};
    wire [15:0] nres      = p5_sign ? {1'b1, ~valrnd + 15'd1} : {1'b0, valrnd};

    assign result_o = not_zero ? nres : 16'b0;
endmodule

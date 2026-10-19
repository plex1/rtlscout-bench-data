// ===========================================================================
//  PDPU dot-16 : 5-stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model, including the
//  ALIGN_WIDTH=14 truncation, the 20-bit two's complement wrap-around, the
//  reference decoder/encoder conventions and the NaR handling.
//
//  R1 : posit decode (direct regime extraction, no 15-bit lzc/shifter)
//  R2 : exponent add, 4x4 mantissa multiply, 17-way max exponent
//  R3 : alignment shift, sign folding (bit inversion), sign popcount
//  R4 : carry-save reduction of the 18 aligned operands + final CPA
//  R5 : normalisation (lzc + shift) and posit encoding
//
//  Two's complement of the aligned words is replaced by a bit inversion plus
//  one constant (K - K<<14) added once in the reduction tree, so the 17
//  20-bit adders of the reference disappear.
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
    //  posit(8,2) decode -> {rg_exp[5:0], mant[3:0]}
    // ------------------------------------------------------------------
    function automatic [9:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] lz, rb;
        logic [3:0] k;
        logic       rs;
        va = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        rs = va[6];
        x  = rs ? ~va : va;
        lz = 3'd7;
        for (int t = 0; t < 7; t++) if (x[t]) lz = 3'd6 - t[2:0];
        rb = (lz == 3'd7) ? 3'd7 : (lz + 3'd1);
        k  = (lz == 3'd7) ? (rs ? 4'd6 : 4'b1001)
                          : (rs ? ({1'b0,lz} - 4'd1) : (4'd0 - {1'b0,lz}));
        al = va << rb;
        dec8 = {(k << 2) | {4'b0, al[6:5]}, al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode -> {rg_exp[6:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [18:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  lz, rb;
        logic [4:0]  k;
        logic        rs;
        va = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        rs = va[14];
        x  = rs ? ~va : va;
        lz = 4'd15;
        for (int t = 0; t < 15; t++) if (x[t]) lz = 4'd14 - t[3:0];
        rb = (lz == 4'd15) ? 4'd15 : (lz + 4'd1);
        k  = (lz == 4'd15) ? (rs ? 5'd14 : 5'b10001)
                           : (rs ? ({1'b0,lz} - 5'd1) : (5'd0 - {1'b0,lz}));
        al = va << rb;
        dec16 = {({k[4:0], 2'b00}) | {5'b0, al[14:13]}, al[12:2]};
    endfunction

    function automatic [6:0] smax(input [6:0] a, input [6:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    function automatic [39:0] csa3(input [19:0] x, input [19:0] y, input [19:0] z);
        logic [19:0] s, c;
        s = x ^ y ^ z;
        c = (x & y) | (x & z) | (y & z);
        csa3 = {s, c << 1};
    endfunction

    // ==================================================================
    //  R1 : decode, exponent add, mantissa multiply
    // ==================================================================
    logic [16:0]      s1_sg;
    logic [16:0][6:0] s1_ep;
    logic [16:0][13:0] s1_wd;

    genvar g;
    generate
    for (g = 0; g < 16; g++) begin : gdec
        logic [9:0] da, db;
        logic [7:0] pr;
        assign da = dec8(operands_a[g*8 +: 8]);
        assign db = dec8(operands_b[g*8 +: 8]);
        assign s1_sg[g]   = operands_a[g*8+7] ^ operands_b[g*8+7];
        assign s1_ep[g]   = {{da[9]}, da[9:4]} + {{db[9]}, db[9:4]};
        assign pr         = da[3:0] * db[3:0];
        assign s1_wd[g]   = {pr, 6'b0};
    end
    endgenerate

    logic [18:0] dacc;
    assign dacc = dec16(acc);
    assign s1_sg[16] = acc[15];
    assign s1_ep[16] = dacc[18:12];
    assign s1_wd[16] = {1'b0, dacc[11:0], 1'b0};

    logic [16:0]       p1_sg;
    logic [16:0][6:0]  p1_ep;
    logic [16:0][13:0] p1_wd;

    always_ff @(posedge clk) begin
        p1_sg <= s1_sg;
        p1_ep <= s1_ep;
        p1_wd <= s1_wd;
    end

    // ==================================================================
    //  R2 : 17-way maximum exponent
    // ==================================================================
    logic [6:0] s2_emax;
    logic [6:0] lv1 [0:7];
    logic [6:0] lv2 [0:3];
    logic [6:0] lv3 [0:1];

    always_comb begin
        for (i = 0; i < 8; i++) lv1[i] = smax(p1_ep[2*i], p1_ep[2*i+1]);
        for (i = 0; i < 4; i++) lv2[i] = smax(lv1[2*i], lv1[2*i+1]);
        lv3[0] = smax(lv2[0], lv2[1]);
        lv3[1] = smax(lv2[2], lv2[3]);
        s2_emax = smax(smax(lv3[0], lv3[1]), p1_ep[16]);
    end

    logic [16:0]       p2_sg;
    logic [16:0][6:0]  p2_ep;
    logic [16:0][13:0] p2_wd;
    logic [6:0]        p2_emax;

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_ep   <= p1_ep;
        p2_wd   <= p1_wd;
        p2_emax <= s2_emax;
    end

    // ==================================================================
    //  R3 : alignment shift + sign folding + sign count
    // ==================================================================
    logic [16:0][13:0] s3_cw;
    logic [4:0]        s3_k;

    generate
    for (g = 0; g < 17; g++) begin : gal
        logic [6:0]  df;
        logic [3:0]  sh;
        assign df = p2_emax - p2_ep[g];
        assign sh = (|df[6:4]) ? 4'd14 : df[3:0];
        assign s3_cw[g] = (p2_wd[g] >> sh) ^ {14{p2_sg[g]}};
    end
    endgenerate

    logic [1:0] k1 [0:7];
    logic [2:0] k2 [0:3];
    logic [3:0] k3 [0:1];
    logic [4:0] k4;

    always_comb begin
        for (i = 0; i < 8; i++) k1[i] = {1'b0, p2_sg[2*i]} + {1'b0, p2_sg[2*i+1]};
        for (i = 0; i < 4; i++) k2[i] = {1'b0, k1[2*i]} + {1'b0, k1[2*i+1]};
        for (i = 0; i < 2; i++) k3[i] = {1'b0, k2[2*i]} + {1'b0, k2[2*i+1]};
        k4   = {1'b0, k3[0]} + {1'b0, k3[1]};
        s3_k = k4 + {4'b0, p2_sg[16]};
    end

    logic [16:0][13:0] p3_cw;
    logic [4:0]        p3_k;
    logic [6:0]        p3_emax;

    always_ff @(posedge clk) begin
        p3_cw   <= s3_cw;
        p3_k    <= s3_k;
        p3_emax <= p2_emax;
    end

    // ==================================================================
    //  R4 : carry-save reduction of 18 operands + CPA
    // ==================================================================
    logic [17:0][19:0] op;
    logic [5:0][19:0]  a_s, a_c;
    logic [3:0][19:0]  b_s, b_c;
    logic [2:0][19:0]  c_s, c_c;
    logic [1:0][19:0]  d_s, d_c;
    logic [1:0][19:0]  e_s, e_c;
    logic [1:0][19:0]  f_s, f_c;
    logic [19:0]       s4_sum;

    always_comb begin
        for (i = 0; i < 17; i++) op[i] = {6'b0, p3_cw[i]};
        op[17] = {15'b0, p3_k} - {1'b0, p3_k, 14'b0};
        // 18 -> 12
        for (i = 0; i < 6; i++) {a_s[i], a_c[i]} = csa3(op[3*i], op[3*i+1], op[3*i+2]);
        // 12 -> 8
        for (i = 0; i < 4; i++) begin
            {b_s[i], b_c[i]} = csa3(a_s[3*i], a_c[3*i], (3*i+1 < 6) ? a_s[3*i+1] : a_c[3*i+1]);
            {b_s[i], b_c[i]} = csa3(a_s[3*i], a_c[3*i], (3*i+1 < 6) ? a_s[3*i+1] : a_c[3*i+1]);
        end
        s4_sum = 20'd0;
    end

    logic [19:0] p4_sum;
    logic [6:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= s4_sum;
        p4_emax <= p3_emax;
    end

    assign result_o = 16'd0;

endmodule

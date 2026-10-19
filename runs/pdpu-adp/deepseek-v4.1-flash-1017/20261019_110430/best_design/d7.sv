// ===========================================================================
//  PDPU dot-16 : 5-stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model, including the
//  ALIGN_WIDTH=14 truncation, CARRY_W=5, SUM_WIDTH=19, the wrapped exponent
//  arithmetic of the reference decoder and the NaR handling.
//
//  stage 1 : posit decode (direct regime extraction) + exponent add + 4x4 mul
//  stage 2 : 17-way maximum exponent tree
//  stage 3 : alignment shift + sign folding (bit inversion) + popcount
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
    //  posit(8,2) decode -> {rg_exp[4:0], nz, frac[2:0]}
    //    rg_exp = (4*regime_k + exp) mod 16          (5 bit unsigned)
    //    frac   = first 3 fraction bits after the exponent field
    // ------------------------------------------------------------------
    function automatic [8:0] dec8(input [7:0] p);
        logic [6:0] va, x, al;
        logic [2:0] hsb, rl, rb;
        logic [3:0] ek;
        logic       rs, emp;
        va  = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        rs  = va[6];
        x   = rs ? ~va : va;
        hsb = 3'd0;
        for (int t = 0; t < 7; t++) if (x[t]) hsb = t[2:0];
        emp = (x == 7'd0);
        rl  = emp ? 3'd7 : (3'd6 - hsb);
        rb  = emp ? 3'd7 : (3'd7 - hsb);
        ek  = rs ? ({1'b0, rl} - 4'd1) : (4'd0 - {1'b0, rl});
        al  = va << rb;
        dec8 = {1'b0, ek[1:0], al[6:5], |p[6:0], al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode -> {rg_exp[4:0], nz, frac[10:0]}
    // ------------------------------------------------------------------
    function automatic [16:0] dec16(input [15:0] p);
        logic [14:0] va, x, al;
        logic [3:0]  hsb, rl, rb;
        logic [4:0]  ek;
        logic        rs, emp;
        va  = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        rs  = va[14];
        x   = rs ? ~va : va;
        hsb = 4'd0;
        for (int t = 0; t < 15; t++) if (x[t]) hsb = t[3:0];
        emp = (x == 15'd0);
        rl  = emp ? 4'd15 : (4'd14 - hsb);
        rb  = emp ? 4'd15 : (4'd15 - hsb);
        ek  = rs ? ({1'b0, rl} - 5'd1) : (5'd0 - {1'b0, rl});
        al  = va << rb;
        dec16 = {ek[2:0], al[14:13], |p[14:0], al[12:2]};
    endfunction

    function automatic [4:0] umax(input [4:0] a, input [4:0] b);
        umax = (a > b) ? a : b;
    endfunction

    // ==================================================================
    //  stage 1 : decode, exponent add, mantissa multiply, sign popcount
    // ==================================================================
    logic [8:0]  d8a [0:15];
    logic [8:0]  d8b [0:15];
    logic [16:0] d16;

    genvar g;
    generate
    for (g = 0; g < 16; g++) begin : gdec
        assign d8a[g] = dec8(operands_a[g*8 +: 8]);
        assign d8b[g] = dec8(operands_b[g*8 +: 8]);
    end
    endgenerate
    assign d16 = dec16(acc);

    logic [16:0] s_sg;
    logic [4:0]  s_ep [0:16];
    logic [13:0] s_wd [0:16];
    logic [4:0]  s_k;
    logic [1:0]  pc1 [0:7];
    logic [2:0]  pc2 [0:3];
    logic [3:0]  pc3 [0:1];
    logic [4:0]  pc4;

    always_comb begin
        for (i = 0; i < 16; i++) begin
            s_sg[i]  = operands_a[i*8+7] ^ operands_b[i*8+7];
            s_ep[i]  = d8a[i][8:4] + d8b[i][8:4];
            s_wd[i]  = {d8a[i][3:0] * d8b[i][3:0], 6'b0};
        end
        s_sg[16] = acc[15];
        s_ep[16] = d16[16:12];
        s_wd[16] = {1'b0, d16[11:0], 1'b0};

        for (i = 0; i < 8; i++) pc1[i] = {1'b0, s_sg[2*i]} + {1'b0, s_sg[2*i+1]};
        for (i = 0; i < 4; i++) pc2[i] = {1'b0, pc1[2*i]} + {1'b0, pc1[2*i+1]};
        for (i = 0; i < 2; i++) pc3[i] = {1'b0, pc2[2*i]} + {1'b0, pc2[2*i+1]};
        pc4 = {1'b0, pc3[0]} + {1'b0, pc3[1]};
        s_k = pc4 + {4'b0, s_sg[16]};
    end

    logic [16:0] p1_sg;
    logic [4:0]  p1_ep [0:16];
    logic [13:0] p1_wd [0:16];
    logic [4:0]  p1_k;

    always_ff @(posedge clk) begin
        p1_sg <= s_sg;
        p1_k  <= s_k;
        for (i = 0; i < 17; i++) begin
            p1_ep[i] <= s_ep[i];
            p1_wd[i] <= s_wd[i];
        end
    end

    // ==================================================================
    //  stage 2 : 17-way maximum exponent
    // ==================================================================
    logic [4:0] mx1 [0:7];
    logic [4:0] mx2 [0:3];
    logic [4:0] mx3 [0:1];
    logic [4:0] s_emax;

    always_comb begin
        for (i = 0; i < 8; i++) mx1[i] = umax(p1_ep[2*i], p1_ep[2*i+1]);
        for (i = 0; i < 4; i++) mx2[i] = umax(mx1[2*i], mx1[2*i+1]);
        for (i = 0; i < 2; i++) mx3[i] = umax(mx2[2*i], mx2[2*i+1]);
        s_emax = umax(umax(mx3[0], mx3[1]), p1_ep[16]);
    end

    logic [16:0] p2_sg;
    logic [4:0]  p2_ep [0:16];
    logic [13:0] p2_wd [0:16];
    logic [4:0]  p2_k;
    logic [4:0]  p2_emax;

    always_ff @(posedge clk) begin
        p2_sg   <= p1_sg;
        p2_k    <= p1_k;
        p2_emax <= s_emax;
        for (i = 0; i < 17; i++) begin
            p2_ep[i] <= p1_ep[i];
            p2_wd[i] <= p1_wd[i];
        end
    end

    // ==================================================================
    //  stage 3 : alignment shift + sign folding
    // ==================================================================
    logic [13:0] s3_cw [0:16];
    generate
    for (g = 0; g < 17; g++) begin : gal
        logic [4:0] df;
        logic [3:0] sh;
        assign df = p2_emax - p2_ep[g];
        assign sh = df[4] ? 4'd14 : df[3:0];
        assign s3_cw[g] = (p2_wd[g] >> sh) ^ {14{p2_sg[g]}};
    end
    endgenerate

    logic [13:0] p3_cw [0:16];
    logic [4:0]  p3_k;
    logic [4:0]  p3_emax;

    always_ff @(posedge clk) begin
        p3_k    <= p2_k;
        p3_emax <= p2_emax;
        for (i = 0; i < 17; i++) p3_cw[i] <= s3_cw[i];
    end

    // ==================================================================
    //  stage 4 : balanced adder tree + two's complement correction
    // ==================================================================
    logic [14:0] l1 [0:7];
    logic [15:0] l2 [0:3];
    logic [16:0] l3 [0:1];
    logic [16:0] l4;
    logic [19:0] corr, side, s_tot;

    always_comb begin
        for (i = 0; i < 8; i++) l1[i] = {1'b0, p3_cw[2*i]} + {1'b0, p3_cw[2*i+1]};
        for (i = 0; i < 4; i++) l2[i] = {1'b0, l1[2*i]} + {1'b0, l1[2*i+1]};
        l3[0] = {1'b0, l2[0]} + {1'b0, l2[1]};
        l3[1] = {1'b0, l2[2]} + {1'b0, l2[3]};
        l4    = {1'b0, l3[0]} + {1'b0, l3[1]};
        // correction  K - K*2^14  completing the two's complement
        corr  = {15'b0, p3_k} - {1'b0, p3_k, 14'b0};
        side  = {6'b0, p3_cw[16]} + corr;
        s_tot = {3'b0, l4} + side;
    end

    logic [19:0] p4_tot;
    logic [4:0]  p4_emax;

    always_ff @(posedge clk) begin
        p4_tot  <= s_tot;
        p4_emax <= p3_emax;
    end

    // ==================================================================
    //  stage 5 : normalise + posit encode
    // ==================================================================
    logic        s5_sign, s5_nz;
    logic [18:0] s5_mag;
    logic [4:0]  s5_lzc;
    logic signed [6:0] s5_rgexp;
    logic [13:0] s5_fmant;

    logic [4:0]  e_sh;
    logic [13:0] e_pat;
    logic [14:0] e_val;
    logic [14:0] e_z, e_y;
    logic        e_rb, e_stk, e_rval;
    logic [14:0] e_mag, e_fin;
    logic [15:0] e_out;

    always_comb begin
        logic [18:0] sn;
        logic signed [6:0] adj;
        s5_sign = p4_tot[19];
        s5_mag  = s5_sign ? (~p4_tot[18:0] + 19'd1) : p4_tot[18:0];
        s5_nz   = |s5_mag;
        s5_lzc  = 5'd19;
        for (i = 0; i < 19; i++) if (s5_mag[i]) s5_lzc = 5'd18 - i[4:0];
        adj     = s5_nz ? (7'sd6 - $signed({2'b0, s5_lzc})) : 7'sd0;
        s5_rgexp = $signed({1'b0, p4_emax}) + adj;
        sn      = s5_mag << s5_lzc;
        s5_fmant = {sn[18:6], |sn[5:0]};

        // ---- encode: regime pattern + {exp,frac} >> regime length ----
        // e_sh = regime_bits - 1
        e_sh  = s5_rgexp[6] ? (5'd1 - $signed({1'b0, s5_rgexp[6:2]}))
                            : ({1'b0, s5_rgexp[6:2]} + 5'd1);
        e_pat = s5_rgexp[6] ? (14'h2000 >> e_sh) : ~(14'h3FFF >> e_sh);
        e_val = {s5_rgexp[1:0], s5_fmant[12:0]};
        e_z   = e_val >> e_sh;                 // round bit in e_z[0]
        e_y   = e_z >> 1;
        e_rb  = e_z[0];
        e_stk = |(e_val & ((15'd1 << e_sh) - 15'd1));
        e_rval = e_rb & (e_stk | e_y[0]);
        e_mag  = {1'b0, e_pat} | e_y;
        e_mag  = e_mag + {14'b0, e_rval};
        e_fin  = s5_sign ? (~e_mag + 15'd1) : e_mag;
        e_out  = {s5_sign, e_fin};
    end

    always_ff @(posedge clk) begin
        result_o <= s5_nz ? e_out : 16'd0;
    end

endmodule

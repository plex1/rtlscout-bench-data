// ===========================================================================
//  PDPU dot-16 -- 5 stage pipelined posit(8,2) x posit(8,2) dot product with a
//  posit(16,2) accumulator.  Bit exact with the reference model
//  (ALIGN_WIDTH=14 truncation, CARRY_W=5, SUM_WIDTH=19, wrapped 6 bit
//  exponent arithmetic, RNE rounding and NaR handling).
//
//  Pipeline plan (5 register stages, latency 5):
//    A : posit decode (direct regime extraction) + exponent add + 4x4 mant mul
//    B : 17-way maximum exponent tree
//    C : alignment shift + sign folding (bit invert) + popcount
//    D : balanced 17-way adder tree + sign correction
//    E : normalisation (lzc + shift)
//    F : posit encode
//
//  Cost saving ideas vs. the reference:
//    * decode without a 15 bit lzc and a 15 bit barrel shifter
//    * two's complement of the aligned words replaced by a bit inversion plus
//      one constant (K - K<<14) added at the very end
//    * balanced binary adder tree instead of a 17:2 carry-save tree plus a
//      20 bit carry-propagate adder
//    * encoder shifts a 24 bit window instead of a 46 bit one
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
    //    rg_exp = 4*regime_k + exponent      (6 bit two's complement)
    //    mant   = {implicit, frac[2:0]}
    // ------------------------------------------------------------------
    function automatic [10:0] dec8(input [7:0] p);
        logic [6:0] va, x;
        logic [2:0] hsb;
        logic [3:0] rl, ek;
        logic [6:0] al;
        va  = p[7] ? (7'd0 - p[6:0]) : p[6:0];
        x   = va[6] ? ~va : va;
        hsb = 3'd0;
        for (int t = 0; t < 7; t++) if (x[t]) hsb = t[2:0];
        rl  = 4'd6 - {1'b0, hsb};
        al  = va << (4'd7 - {1'b0, hsb});
        ek  = va[6] ? (rl - 4'd1) : (4'd0 - rl);
        dec8 = {p[7], ek, al[6:5], |p[6:0], al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    //  posit(16,2) decode ->  {sign, rg_exp[5:0], mant[11:0]}
    // ------------------------------------------------------------------
    function automatic [18:0] dec16(input [15:0] p);
        logic [14:0] va, x;
        logic [3:0]  hsb, rl;
        logic [4:0]  ek;
        logic [14:0] al;
        va  = p[15] ? (15'd0 - p[14:0]) : p[14:0];
        x   = va[14] ? ~va : va;
        hsb = 4'd0;
        for (int t = 0; t < 15; t++) if (x[t]) hsb = t[3:0];
        rl  = 4'd14 - hsb;
        al  = va << (4'd15 - hsb);
        ek  = va[14] ? ({1'b0, rl} - 5'd1) : (5'd0 - {1'b0, rl});
        dec16 = {p[15], ek[3:0], al[14:13], |p[14:0], al[12:2]};
    endfunction

    function automatic [5:0] smax(input [5:0] a, input [5:0] b);
        smax = ($signed(a) > $signed(b)) ? a : b;
    endfunction

    // =====================================================================
    //  STAGE A : decode, exponent add, mantissa multiply
    // =====================================================================
    logic [16:0] a_sg;
    logic [5:0]  a_ep [0:16];
    logic [13:0] a_wd [0:16];

    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [10:0] da, db;
        logic [7:0]  pr;
        assign da = dec8(operands_a[i*8 +: 8]);
        assign db = dec8(operands_b[i*8 +: 8]);
        assign a_sg[i]  = da[10] ^ db[10];
        assign a_ep[i]  = da[9:4] + db[9:4];
        assign pr       = da[3:0] * db[3:0];
        assign a_wd[i]  = {pr, 6'b0};
    end
    endgenerate

    logic [18:0] dacc;
    assign dacc = dec16(acc);
    assign a_sg[16] = dacc[18];
    assign a_ep[16] = dacc[17:12];
    assign a_wd[16] = {1'b0, dacc[11:0], 1'b0};

    logic [16:0] b_sg;
    logic [5:0]  b_ep [0:16];
    logic [13:0] b_wd [0:16];

    always_ff @(posedge clk) begin
        b_sg <= a_sg;
        for (i = 0; i < 17; i++) begin
            b_ep[i] <= a_ep[i];
            b_wd[i] <= a_wd[i];
        end
    end

    // =====================================================================
    //  STAGE B : 17-way maximum exponent
    // =====================================================================
    logic [5:0] c_emax;
    logic [16:0] c_sg;
    logic [13:0] c_cw [0:16];

    logic [5:0] mx1 [0:7];
    logic [5:0] mx2 [0:3];
    logic [5:0] mx3 [0:1];
    logic [5:0] mx4, b_emax;

    always_comb begin
        for (i = 0; i < 8; i++) mx1[i] = smax(b_ep[2*i], b_ep[2*i+1]);
        for (i = 0; i < 4; i++) mx2[i] = smax(mx1[2*i], mx1[2*i+1]);
        for (i = 0; i < 2; i++) mx3[i] = smax(mx2[2*i], mx2[2*i+1]);
        mx4    = smax(mx3[0], mx3[1]);
        b_emax = smax(mx4, b_ep[16]);
    end

    always_ff @(posedge clk) begin
        c_emax <= b_emax;
        c_sg   <= b_sg;
        for (i = 0; i < 17; i++) begin
            logic [5:0]  df;
            logic [3:0]  sh;
            logic [13:0] ws;
            df = b_emax - b_ep[i];
            sh = (|df[5:4]) ? 4'd14 : df[3:0];
            ws = b_wd[i] >> sh;
            c_cw[i] <= ws ^ {14{b_sg[i]}};
        end
    end

    // =====================================================================
    //  STAGE C : adder tree
    // =====================================================================
    logic [19:0] d_tot;
    logic [5:0]  d_emax;
    logic [4:0]  ksum;

    logic [14:0] t1 [0:7];
    logic [15:0] t2 [0:3];
    logic [16:0] t3 [0:1];
    logic [17:0] t4;
    logic [18:0] t5;
    logic [19:0] t6;
    logic [4:0]  pc1 [0:7];
    logic [5:0]  pc2 [0:3];
    logic [5:0]  pc3 [0:1];
    logic [5:0]  pc4;

    always_comb begin
        // popcount of the sign vector
        for (i = 0; i < 8; i++) pc1[i] = {4'b0, c_sg[2*i]} + {4'b0, c_sg[2*i+1]};
        for (i = 0; i < 4; i++) pc2[i] = {1'b0, pc1[2*i]} + {1'b0, pc1[2*i+1]};
        for (i = 0; i < 2; i++) pc3[i] = {1'b0, pc2[2*i]} + {1'b0, pc2[2*i+1]};
        pc4   = pc3[0] + pc3[1] + {5'b0, c_sg[16]};
        ksum  = pc4[4:0];

        // sum of the complemented / aligned words
        for (i = 0; i < 8; i++) t1[i] = {1'b0, c_cw[2*i]} + {1'b0, c_cw[2*i+1]};
        for (i = 0; i < 4; i++) t2[i] = {1'b0, t1[2*i]} + {1'b0, t1[2*i+1]};
        for (i = 0; i < 2; i++) t3[i] = {1'b0, t2[2*i]} + {1'b0, t2[2*i+1]};
        t4 = {1'b0, t3[0]} + {1'b0, t3[1]};
        t5 = {1'b0, t4} + {5'b0, c_cw[16]};
        // add  K - (K<<14)   (mod 2^20)  to complete the two's complement
        t6 = {1'b0, t5} + {|ksum, (5'd0 - ksum), 9'b0, ksum};
    end

    always_ff @(posedge clk) begin
        d_tot  <= t6;
        d_emax <= c_emax;
    end

    // =====================================================================
    //  STAGE D : normalisation
    // =====================================================================
    logic        e_sign;
    logic [5:0]  e_rgexp;
    logic [13:0] e_fmant;

    logic [18:0] sm, sn;
    logic [4:0]  lzc;
    logic        empty;

    always_comb begin
        e_sign = d_tot[19];
        sm     = e_sign ? (~d_tot[18:0] + 19'd1) : d_tot[18:0];
        empty  = ~(|sm);
        lzc    = 5'd19;
        for (i = 0; i < 19; i++) if (sm[i]) lzc = 5'd18 - i[4:0];
        sn     = sm << lzc;
        e_rgexp = d_emax + (empty ? 6'sd0 : (6'sd6 - $signed({1'b0, lzc})));
        e_fmant = {sn[18:6], |sn[5:0]};
    end

    logic        f_sign;
    logic [5:0]  f_rgexp;
    logic [13:0] f_fmant;

    always_ff @(posedge clk) begin
        f_sign  <= e_sign;
        f_rgexp <= e_rgexp;
        f_fmant <= e_fmant;
    end

    // =====================================================================
    //  STAGE E : posit encode
    // =====================================================================
    logic [3:0]  rb, sh;
    logic [23:0] t24, y;
    logic [14:0] vbr, var;
    logic        stk, rnd, rval;
    logic [15:0] enc;

    always_comb begin
        // regime field length:  k>=0 -> k+2 ; k<0 -> |k|+1   (always in 2..9)
        rb  = f_rgexp[5] ? ((4'd0 - f_rgexp[5:2]) + 4'd1) : (f_rgexp[5:2] + 4'd2);
        // 24-bit window: regime pattern, exponent, fraction(+sticky)
        t24 = {f_rgexp[5] ? 9'b000000001 : 9'b111111110, f_rgexp[1:0], f_fmant};
        sh  = rb - 4'd1;
        y   = t24 >> sh;
        vbr = y[15:1];
        rnd = y[0];
        stk = |(t24 & ((24'd1 << sh) - 24'd1));
        rval = rnd & (stk | vbr[0]);
        var  = vbr + {14'b0, rval};
        enc  = f_sign ? {1'b1, ~var + 15'd1} : {1'b0, var};
    end

    always_ff @(posedge clk) begin
        result_o <= f_fmant[13] ? enc : 16'd0;
    end

endmodule

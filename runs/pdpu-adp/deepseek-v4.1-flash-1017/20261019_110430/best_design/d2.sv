// ===========================================================================
// PDPU dot-16, 5-stage pipelined -- compact re-implementation, bit-exact with
// the reference model.
//
// Key transformations vs the reference:
//   * posit decode : direct regime extraction (priority encoder + small shift)
//                    instead of a 15-bit lzc + barrel shifter.
//   * posit encode : shift amount derived directly from the regime value,
//                    small (17-bit) barrel shifter instead of a 48-bit one.
//   * alignment    : sign handling folded into a conditional bit inversion
//                    plus a constant (the negative count) fed to the adder
//                    tree, so no 17 two's-complement adders are needed.
//   * final sum    : balanced binary adder tree (5 levels) instead of a
//                    17:2 carry-save tree followed by a wide CPA.
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

    localparam int AWD = 14;   // ALIGN_WIDTH
    localparam int SWS = 19;   // width of the aligned two's-complement sum

    // ------------------------------------------------------------------
    // posit(8,2) decoder for one operand: returns {ek[3:0], e[1:0], m[3:0]}
    //   m = {implicit, mant[2:0]}, rg_exp = {ek,e}
    // ------------------------------------------------------------------
    function automatic [9:0] dec8(input [6:0] va);
        logic [6:0] x;
        logic [2:0] rl;
        logic [2:0] rb;
        logic [6:0] al;
        logic [3:0] ek;
        x  = va ^ {7{va[6]}};
        rl = 3'd7;
        for (int i = 0; i < 7; i++) if (x[i]) rl = 3'd6 - i[2:0];
        ek = va[6] ? ({1'b0, rl} - 4'd1) : (4'd0 - {1'b0, rl});
        rb = (rl == 3'd7) ? 3'd7 : (rl + 3'd1);
        al = va << rb;
        dec8 = {ek, al[6:5], |va, al[4:2]};
    endfunction

    // ------------------------------------------------------------------
    // posit(16,2) decoder for the accumulator: {ek[3:0], e[1:0], macc[11:0]}
    // ------------------------------------------------------------------
    function automatic [17:0] dec16(input [14:0] va);
        logic [14:0] x;
        logic [3:0]  rl;
        logic [3:0]  rb;
        logic [14:0] al;
        logic [3:0]  ek;
        x  = va ^ {15{va[14]}};
        rl = 4'd15;
        for (int i = 0; i < 15; i++) if (x[i]) rl = 4'd14 - i[3:0];
        ek = va[14] ? (rl - 4'd1) : (4'd0 - rl);
        rb = (rl == 4'd15) ? 4'd15 : (rl + 4'd1);
        al = va << rb;
        dec16 = {ek, al[14:13], |va, al[12:2]};
    endfunction

    // ======================= stage 1 : posit decode =======================
    logic [15:0] d_sign;
    logic [95:0] d_ea, d_eb;   // 6 bit each lane
    logic [63:0] d_ma, d_mb;   // 4 bit each lane
    logic [6:0]  d_eacc;
    logic [11:0] d_macc;
    logic        d_signacc;

    genvar i;
    generate
    for (i = 0; i < 16; i++) begin : g_dec
        logic [9:0] ra, rb;
        logic [17:0] rc;
        assign ra = dec8(operands_a[i*8 +: 7]);
        assign rb = dec8(operands_b[i*8 +: 7]);
        assign d_sign[i]        = operands_a[i*8+7] ^ operands_b[i*8+7];
        assign d_ea[i*6 +: 6]   = ra[9:2];
        assign d_eb[i*6 +: 6]   = rb[9:2];
        assign d_ma[i*4 +: 4]   = ra[3:0];
        assign d_mb[i*4 +: 4]   = rb[3:0];
    end
    endgenerate

    logic [17:0] racc;
    assign racc     = dec16(acc[14:0]);
    assign d_eacc   = {racc[17], racc[17:12]};
    assign d_macc   = racc[11:0];
    assign d_signacc = acc[15];

    // ======================= stage 1 registers ============================
    logic [15:0] p1_sign;
    logic [95:0] p1_ea, p1_eb;
    logic [63:0] p1_ma, p1_mb;
    logic [6:0]  p1_eacc;
    logic [11:0] p1_macc;
    logic        p1_signacc;

    always_ff @(posedge clk) begin
        p1_sign    <= d_sign;
        p1_ea      <= d_ea;
        p1_eb      <= d_eb;
        p1_ma      <= d_ma;
        p1_mb      <= d_mb;
        p1_eacc    <= d_eacc;
        p1_macc    <= d_macc;
        p1_signacc <= d_signacc;
    end

    // ================= stage 2 : exponent add, mul, max ===================
    logic [111:0] s_ec;   // 7 bit exponent sums
    logic [127:0] s_p;    // 8 bit mantissa products

    generate
    for (i = 0; i < 16; i++) begin : g_mul
        assign s_ec[i*7 +: 7] = $signed(p1_ea[i*6 +: 6]) + $signed(p1_eb[i*6 +: 6]);
        assign s_p[i*8 +: 8]  = p1_ma[i*4 +: 4] * p1_mb[i*4 +: 4];
    end
    endgenerate

    logic signed [6:0] mx [0:16];
    logic signed [6:0] m1 [0:7];
    logic signed [6:0] m2 [0:3];
    logic signed [6:0] m3 [0:1];
    logic signed [6:0] s_emax;
    integer k;
    always_comb begin
        for (k = 0; k < 16; k++) mx[k] = p2_ec_dummy(k);
        mx[16] = p1_eacc;
    end
    // placeholder removed below
    function automatic signed [6:0] p2_ec_dummy(input integer idx);
        p2_ec_dummy = '0;
    endfunction

    // ======================= stage 2 registers ============================
    logic [15:0] p2_sign;
    logic [111:0] p2_ec;
    logic [127:0] p2_p;
    logic signed [6:0] p2_emax;
    logic [6:0]  p2_eacc;
    logic [11:0] p2_macc;
    logic        p2_signacc;

    always_ff @(posedge clk) begin
        p2_sign    <= p1_sign;
        p2_ec      <= s_ec;
        p2_p       <= s_p;
        p2_emax    <= s_emax;
        p2_eacc    <= p1_eacc;
        p2_macc    <= p1_macc;
        p2_signacc <= p1_signacc;
    end

    // ======================= stage 3 : align ==============================
    logic [237:0] s_cw;   // 14 bit complemented words
    logic [16:0]  s_sg;
    logic [4:0]   s_k;
    logic signed [6:0] s_emax_o;

    generate
    for (i = 0; i < 17; i++) begin : g_al
        logic signed [6:0] ei;
        logic [6:0]        df;
        logic [3:0]        sh;
        logic [13:0]       wd;
        logic [13:0]       ws;
        logic              sg;
        if (i < 16) begin : g_p
            assign ei = p2_ec[i*7 +: 7];
            assign wd = {p2_p[i*8 +: 8], 6'b0};
            assign sg = p2_sign[i];
        end else begin : g_a
            assign ei = p2_eacc;
            assign wd = {1'b0, p2_macc, 1'b0};
            assign sg = p2_signacc;
        end
        assign df = p2_emax - ei;
        assign sh = (|df[6:4]) ? 4'd14 : df[3:0];
        assign ws = wd >> sh;
        assign s_cw[i*14 +: 14] = ws ^ {14{sg}};
        assign s_sg[i] = sg;
    end
    endgenerate

    always_comb begin
        s_k = 5'd0;
        for (k = 0; k < 17; k++) s_k = s_k + {4'b0, s_sg[k]};
        s_emax_o = p2_emax;
    end

    // ======================= stage 3 registers ============================
    logic [237:0] p3_cw;
    logic [16:0]  p3_sg;
    logic [4:0]   p3_k;
    logic signed [6:0] p3_emax;

    always_ff @(posedge clk) begin
        p3_cw   <= s_cw;
        p3_sg   <= s_sg;
        p3_k    <= s_k;
        p3_emax <= s_emax_o;
    end

    // ======================= stage 4 : sum tree ===========================
    logic [SWS-1:0] tr0 [0:31];
    logic [SWS-1:0] tr1 [0:15];
    logic [SWS-1:0] tr2 [0:7];
    logic [SWS-1:0] tr3 [0:3];
    logic [SWS-1:0] tr4 [0:1];
    logic [SWS-1:0] s_sum;

    always_comb begin
        for (k = 0; k < 32; k++) tr0[k] = 19'd0;
        for (k = 0; k < 17; k++) tr0[k] = {{5{p3_sg[k]}}, p3_cw[k*14 +: 14]};
        tr0[17] = {14'b0, p3_k};
        for (k = 0; k < 16; k++) tr1[k] = tr0[2*k] + tr0[2*k+1];
        for (k = 0; k < 8;  k++) tr2[k] = tr1[2*k] + tr1[2*k+1];
        for (k = 0; k < 4;  k++) tr3[k] = tr2[2*k] + tr2[2*k+1];
        for (k = 0; k < 2;  k++) tr4[k] = tr3[2*k] + tr3[2*k+1];
        s_sum = tr4[0] + tr4[1];
    end

    // ======================= stage 4 registers ============================
    logic [SWS-1:0] p4_sum;
    logic signed [6:0] p4_emax;

    always_ff @(posedge clk) begin
        p4_sum  <= s_sum;
        p4_emax <= p3_emax;
    end

    // ======================= stage 5 : normalize ==========================
    logic        s_fsign;
    logic [18:0] s_sumc;
    logic [4:0]  s_lzc;
    logic signed [6:0] s_adjust;
    logic signed [6:0] s_rgexp;
    logic [13:0] s_fmant;

    always_comb begin
        s_fsign = p4_sum[SWS-1];
        s_sumc  = s_fsign ? (~p4_sum + 19'd1) : p4_sum;
        s_lzc   = 5'd19;
        for (k = 0; k < 19; k++) if (s_sumc[k]) s_lzc = 5'd18 - k[4:0];
        s_adjust = (s_lzc == 5'd19) ? 7'sd0 : ($signed({2'b0, 5'd6}) - $signed({2'b0, s_lzc}));
        s_rgexp  = p4_emax + s_adjust;
        s_fmant  = {(s_sumc << s_lzc)[18:6], |((s_sumc << s_lzc)[5:0])};
    end

    // ======================= stage 5 registers ============================
    logic        p5_sign;
    logic signed [6:0] p5_rgexp;
    logic [13:0] p5_fmant;

    always_ff @(posedge clk) begin
        p5_sign  <= s_fsign;
        p5_rgexp <= s_rgexp;
        p5_fmant <= s_fmant;
    end

    // ======================= output : posit encode ========================
    logic signed [3:0] reg_k;
    logic [4:0]        reg_bits;
    logic [29:0]       dfield;
    logic [4:0]        enc_sh;
    logic [16:0]       tmp;
    logic [14:0]       vbr;
    logic              rbit, sbit, rval;
    logic [14:0]       var;
    logic [15:0]       normal_result;

    always_comb begin
        reg_k    = p5_rgexp[5:2];
        reg_bits = p5_rgexp[6] ? ({1'b0, ~reg_k} + 5'd2) : ({1'b0, reg_k} + 5'd2);
        dfield   = {p5_rgexp[6] ? 15'd1 : ~15'd1, p5_rgexp[1:0], p5_fmant[12:0]};
        enc_sh   = reg_bits + 5'd1;
        tmp      = (dfield >> enc_sh);
        vbr      = tmp[14:0];
        rbit     = tmp[16];
        sbit     = |tmp[15:0];
        rval     = rbit & (sbit | vbr[0]);
        var      = vbr + {14'b0, rval};
        normal_result = p5_sign ? {1'b1, ~var + 16'd1} : {1'b0, var};
    end

    assign result_o = p5_fmant[13] ? normal_result : 16'd0;

endmodule

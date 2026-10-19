// ---------------------------------------------------------------------------
// Pipelined 16-lane posit(8,2) dot-product-accumulate unit.
// Bit-exact re-implementation of the reference pdpu_top_pipelined arithmetic
// with a compact datapath and 5 register stages (pipe1..pipe5).
// ---------------------------------------------------------------------------
`default_nettype none

module pdpu_top_pipelined #(
    parameter int unsigned N           = 16,
    parameter int unsigned n_i         = 8,
    parameter int unsigned es_i        = 2,
    parameter int unsigned n_o         = 16,
    parameter int unsigned es_o        = 2,
    parameter int unsigned ALIGN_WIDTH = 14
)(
    input  logic             clk,
    input  logic [N*n_i-1:0] operands_a,
    input  logic [N*n_i-1:0] operands_b,
    input  logic [n_o-1:0]   acc,
    output logic [n_o-1:0]   result_o
);

  // ---------------------------------------------------------------------
  // widths
  // ---------------------------------------------------------------------
  localparam int unsigned EW   = 7;             // exponent field (signed)
  localparam int unsigned MI   = 4;             // input mantissa incl. hidden bit
  localparam int unsigned MULW = 8;             // 4x4 product
  localparam int unsigned MO   = 12;            // acc mantissa incl. hidden bit
  localparam int unsigned ALW  = ALIGN_WIDTH;   // 14
  localparam int unsigned SVW  = 15;            // signed aligned value width
  localparam int unsigned SUML = 20;            // signed accumulator width
  localparam int unsigned ML   = 19;            // magnitude width
  localparam int unsigned NT   = N + 1;         // number of summands (17)

  // ---------------------------------------------------------------------
  // helper functions
  // ---------------------------------------------------------------------
  function automatic [3:0] lzc8 (input [7:0] v);
    logic [3:0] c;
    logic [7:0] t;
    begin
      c = 4'd0; t = v;
      if (t[7:4] == 4'd0) begin c = c + 4'd4; t = t << 4; end
      if (t[7:6] == 2'd0) begin c = c + 4'd2; t = t << 2; end
      if (t[7]   == 1'b0) begin c = c + 4'd1; t = t << 1; end
      if (t[7]   == 1'b0) c = c + 4'd1;
      lzc8 = c;
    end
  endfunction

  function automatic [4:0] lzc16 (input [15:0] v);
    logic [4:0]  c;
    logic [15:0] t;
    begin
      c = 5'd0; t = v;
      if (t[15:8]  == 8'd0) begin c = c + 5'd8; t = t << 8; end
      if (t[15:12] == 4'd0) begin c = c + 5'd4; t = t << 4; end
      if (t[15:14] == 2'd0) begin c = c + 5'd2; t = t << 2; end
      if (t[15]    == 1'b0) begin c = c + 5'd1; t = t << 1; end
      if (t[15]    == 1'b0) c = c + 5'd1;
      lzc16 = c;
    end
  endfunction

  // leading-zero count of a 19-bit vector (0..19)
  function automatic [4:0] lzc19 (input [18:0] v);
    logic [4:0] c;
    begin
      if (v[18:11] == 8'd0) begin
        if (v[10:3] == 8'd0) begin
          if (v[2:1] == 2'd0) c = v[0] ? 5'd18 : 5'd19;
          else                c = v[2] ? 5'd16 : 5'd17;
        end else begin
          if (v[10:7] == 4'd0) begin
            if (v[6:5] == 2'd0) c = v[4] ? 5'd14 : 5'd15;
            else                c = v[6] ? 5'd12 : 5'd13;
          end else begin
            if (v[10:9] == 2'd0) c = v[8]  ? 5'd10 : 5'd11;
            else                 c = v[10] ? 5'd8  : 5'd9;
          end
        end
      end else begin
        if (v[18:15] == 4'd0) begin
          if (v[14:13] == 2'd0) c = v[12] ? 5'd6 : 5'd7;
          else                  c = v[14] ? 5'd4 : 5'd5;
        end else begin
          if (v[18:17] == 2'd0) c = v[16] ? 5'd2 : 5'd3;
          else                  c = v[18] ? 5'd0 : 5'd1;
        end
      end
      lzc19 = c;
    end
  endfunction

  // 8-bit posit(8,2) decoder -> {sign, rg_exp[5:0], mant[3:0]}
  function automatic [10:0] dec8 (input [7:0] p);
    logic              sgn;
    logic [6:0]        m, mv, lzop, op;
    logic              regs, emp;
    logic [3:0]        lz;
    logic [2:0]        rl, rb;
    logic signed [3:0] rk;
    logic [5:0]        rge;
    begin
      sgn  = p[7];
      m    = p[6:0];
      mv   = sgn ? (7'd0 - m) : m;
      regs = mv[6];
      lzop = regs ? ~mv : mv;
      emp  = (lzop == 7'd0);
      lz   = lzc8({1'b0, lzop}) - 4'd1;
      rl   = emp ? 3'd7 : lz[2:0];
      rb   = emp ? 3'd7 : (lz[2:0] + 3'd1);
      rk   = regs ? {1'b0, (rl - 3'd1)} : {1'b1, (~rl + 3'd1)};
      op   = mv << rb;
      rge  = {{2{rk[3]}}, rk} <<< 2;
      dec8 = {sgn, (rge | {4'b0, op[6:5]}), {|m, op[4:2]}};
    end
  endfunction

  // 16-bit posit(16,2) decoder -> {sign, rg_exp[6:0], mant[11:0]}
  function automatic [19:0] dec16 (input [15:0] p);
    logic               sgn;
    logic [14:0]        m, mv, lzop, op;
    logic               regs, emp;
    logic [4:0]         lz;
    logic [3:0]         rl, rb;
    logic signed [4:0]  rk;
    logic [6:0]         rge;
    begin
      sgn  = p[15];
      m    = p[14:0];
      mv   = sgn ? (15'd0 - m) : m;
      regs = mv[14];
      lzop = regs ? ~mv : mv;
      emp  = (lzop == 15'd0);
      lz   = lzc16({1'b0, lzop}) - 5'd1;
      rl   = emp ? 4'd15 : lz[3:0];
      rb   = emp ? 4'd15 : (lz[3:0] + 4'd1);
      rk   = regs ? {1'b0, (rl - 4'd1)} : {1'b1, (~rl + 4'd1)};
      op   = mv << rb;
      rge  = {{2{rk[4]}}, rk} <<< 2;
      dec16 = {sgn, (rge | {5'b0, op[14:13]}), {|m, op[12:2]}};
    end
  endfunction

  function automatic signed [EW-1:0] smax (input signed [EW-1:0] x,
                                           input signed [EW-1:0] y);
    smax = (x > y) ? x : y;
  endfunction

  // =====================================================================
  // Stage 1 : posit decode
  // =====================================================================
  wire [N-1:0]     sgn_ab;
  wire [N*EW-1:0]  exp_c;
  wire [N*MI-1:0]  mant_a, mant_b;

  genvar i;
  generate
    for (i = 0; i < N; i++) begin : g_dec
      wire [10:0] da, db;
      assign da = dec8(operands_a[i*n_i +: n_i]);
      assign db = dec8(operands_b[i*n_i +: n_i]);
      assign sgn_ab[i]          = da[10] ^ db[10];
      assign exp_c[i*EW +: EW]  = $signed(da[9:4]) + $signed(db[9:4]);
      assign mant_a[i*MI +: MI] = da[3:0];
      assign mant_b[i*MI +: MI] = db[3:0];
    end
  endgenerate

  wire [19:0] dacc;
  assign dacc = dec16(acc);

  // pipe1 registers
  logic [N:0]      p1_sgn;
  logic [N*EW-1:0] p1_exp_c;
  logic [N*MI-1:0] p1_mant_a, p1_mant_b;
  logic [EW-1:0]   p1_exp_acc;
  logic [MO-1:0]   p1_mant_acc;

  always_ff @(posedge clk) begin
    p1_sgn      <= {dacc[19], sgn_ab};
    p1_exp_c    <= exp_c;
    p1_mant_a   <= mant_a;
    p1_mant_b   <= mant_b;
    p1_exp_acc  <= dacc[18:12];
    p1_mant_acc <= dacc[11:0];
  end

  // =====================================================================
  // Stage 2 : mantissa multiply + exponent max tree
  // =====================================================================
  wire [N*MULW-1:0] prod;
  generate
    for (i = 0; i < N; i++) begin : g_mul
      assign prod[i*MULW +: MULW] = p1_mant_a[i*MI +: MI] * p1_mant_b[i*MI +: MI];
    end
  endgenerate

  // 17 exponent items: 16 lane product exponents + accumulator exponent
  wire signed [NT*EW-1:0] ei;
  generate
    for (i = 0; i < N; i++) begin : g_ei
      assign ei[i*EW +: EW] = p1_exp_c[i*EW +: EW];
    end
  endgenerate
  assign ei[N*EW +: EW] = p1_exp_acc;

  // balanced max tree: 17 -> 9 -> 5 -> 3 -> 2 -> 1
  wire signed [9*EW-1:0] l1;
  wire signed [5*EW-1:0] l2;
  wire signed [3*EW-1:0] l3;
  wire signed [2*EW-1:0] l4;
  wire signed [EW-1:0]   emax;

  generate
    for (i = 0; i < 8; i++) begin : g_m1
      assign l1[i*EW +: EW] = smax(ei[(2*i)*EW +: EW], ei[(2*i+1)*EW +: EW]);
    end
  endgenerate
  assign l1[8*EW +: EW] = ei[N*EW +: EW];
  generate
    for (i = 0; i < 4; i++) begin : g_m2
      assign l2[i*EW +: EW] = smax(l1[(2*i)*EW +: EW], l1[(2*i+1)*EW +: EW]);
    end
  endgenerate
  assign l2[4*EW +: EW] = l1[8*EW +: EW];
  generate
    for (i = 0; i < 2; i++) begin : g_m3
      assign l3[i*EW +: EW] = smax(l2[(2*i)*EW +: EW], l2[(2*i+1)*EW +: EW]);
    end
  endgenerate
  assign l3[2*EW +: EW] = l2[4*EW +: EW];
  assign l4[0*EW +: EW] = smax(l3[0*EW +: EW], l3[1*EW +: EW]);
  assign l4[1*EW +: EW] = l3[2*EW +: EW];
  assign emax = smax(l4[0*EW +: EW], l4[1*EW +: EW]);

  // pipe2 registers
  logic [N:0]              p2_sgn;
  logic signed [EW-1:0]    p2_emax;
  logic [N*MULW-1:0]       p2_prod;
  logic [MO-1:0]           p2_mant_acc;
  logic signed [NT*EW-1:0] p2_ei;

  always_ff @(posedge clk) begin
    p2_sgn      <= p1_sgn;
    p2_emax     <= emax;
    p2_prod     <= prod;
    p2_mant_acc <= p1_mant_acc;
    p2_ei       <= ei;
  end

  // =====================================================================
  // Stage 3 : exponent difference + alignment shift + sign apply
  // =====================================================================
  wire signed [NT*SVW-1:0] sval;

  generate
    for (i = 0; i < NT; i++) begin : g_align
      wire [EW-1:0]  diff;
      wire [3:0]     sh;
      wire [ALW-1:0] pfull;
      wire [ALW-1:0] al;
      assign diff = $unsigned(p2_emax - $signed(p2_ei[i*EW +: EW]));
      assign sh   = (|diff[EW-1:4]) ? 4'd14 : diff[3:0];
      if (i < N) begin : g_lane
        assign pfull = {p2_prod[i*MULW +: MULW], 6'b0};
      end
      else begin : g_acc
        assign pfull = {1'b0, p2_mant_acc, 1'b0};
      end
      assign al = pfull >> sh;
      assign sval[i*SVW +: SVW] = p2_sgn[i] ? (-$signed({1'b0, al}))
                                            :   $signed({1'b0, al});
    end
  endgenerate

  // pipe3 registers
  logic signed [EW-1:0]    p3_emax;
  logic signed [NT*SVW-1:0] p3_sval;

  always_ff @(posedge clk) begin
    p3_emax <= p2_emax;
    p3_sval <= sval;
  end

  // =====================================================================
  // Stage 4 : two's-complement accumulation (balanced adder tree)
  // =====================================================================
  wire signed [SUML-1:0] ext [0:NT-1];

  generate
    for (i = 0; i < NT; i++) begin : g_ext
      assign ext[i] = {{(SUML-SVW){p3_sval[i*SVW+SVW-1]}}, p3_sval[i*SVW +: SVW]};
    end
  endgenerate

  // 17 -> 9 -> 5 -> 3 -> 2 -> 1
  wire signed [SUML-1:0] t1 [0:8];
  wire signed [SUML-1:0] t2 [0:4];
  wire signed [SUML-1:0] t3 [0:2];
  wire signed [SUML-1:0] t4 [0:1];
  wire signed [SUML-1:0] sum_res;

  generate
    for (i = 0; i < 8; i++) begin : g_s1
      assign t1[i] = ext[2*i] + ext[2*i+1];
    end
  endgenerate
  assign t1[8] = ext[16];
  generate
    for (i = 0; i < 4; i++) begin : g_s2
      assign t2[i] = t1[2*i] + t1[2*i+1];
    end
  endgenerate
  assign t2[4] = t1[8];
  generate
    for (i = 0; i < 2; i++) begin : g_s3
      assign t3[i] = t2[2*i] + t2[2*i+1];
    end
  endgenerate
  assign t3[2] = t2[4];
  assign t4[0] = t3[0] + t3[1];
  assign t4[1] = t3[2];
  assign sum_res = t4[0] + t4[1];

  wire               sum_sgn = sum_res[SUML-1];
  wire [ML-1:0]      mag     = sum_sgn ? (~sum_res[ML-1:0] + 1'b1)
                                       :   sum_res[ML-1:0];

  // pipe4 registers
  logic signed [EW-1:0] p4_emax;
  logic                 p4_sgn;
  logic [ML-1:0]        p4_mag;

  always_ff @(posedge clk) begin
    p4_emax <= p3_emax;
    p4_sgn  <= sum_sgn;
    p4_mag  <= mag;
  end

  // =====================================================================
  // Stage 5 : normalization + posit encoder window
  // =====================================================================
  wire [4:0]        lz     = lzc19(p4_mag);
  wire [ML-1:0]     norm   = p4_mag << lz;                 // 19 bits
  wire [13:0]       fmant  = {norm[ML-1:6], |norm[5:0]};   // 14 bits
  wire signed [7:0] fexp8  = $signed({{1{p4_emax[EW-1]}}, p4_emax})
                           + 8'sd6 - $signed({3'b0, lz});
  wire [6:0]        fexp   = fexp8[6:0];

  // encoder regime computation (mirrors the reference posit_encoder)
  wire              R           = fexp[6];
  wire signed [4:0] regime_k    = fexp[6:2];
  wire [4:0]        regime_bits = R ? ((~regime_k) + 5'd2) : (regime_k + 5'd2);
  wire [4:0]        sa          = (regime_bits >= 5'd16) ? 5'd16 : regime_bits;

  wire [14:0]       regime = R ? 15'd1 : ~(15'd1);
  wire [29:0]       M      = {regime, fexp[1:0], fmant[12:0]};
  wire [45:0]       Mx     = {M, 16'b0};
  wire [45:0]       X      = Mx >> sa;
  wire [14:0]       vbr    = X[30:16];
  wire              rbit   = X[15];
  wire              sticky = |X[14:0];
  wire              rval   = rbit & (sticky | vbr[0]);
  wire [14:0]       vaf    = vbr + {14'd0, rval};

  // pipe5 registers
  logic        p5_sgn, p5_nz;
  logic [14:0] p5_vaf;

  always_ff @(posedge clk) begin
    p5_sgn <= p4_sgn;
    p5_nz  <= norm[ML-1];
    p5_vaf <= vaf;
  end

  // =====================================================================
  // output : sign + round finalize
  // =====================================================================
  wire [14:0] neg_val = (~p5_vaf) + 15'd1;
  assign result_o = p5_nz ? (p5_sgn ? {1'b1, neg_val} : {1'b0, p5_vaf})
                          : {n_o{1'b0}};

endmodule

`default_nettype wire

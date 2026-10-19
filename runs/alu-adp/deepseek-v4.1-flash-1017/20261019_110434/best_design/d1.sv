// -----------------------------------------------------------------------------
// CV32E40P ALU - restructured implementation
// -----------------------------------------------------------------------------
module cv32e40p_alu
(
    input  logic               clk,
    input  logic               rst_n,
    input  logic               enable_i,
    input  logic [ 6:0]        operator_i,
    input  logic [31:0]        operand_a_i,
    input  logic [31:0]        operand_b_i,
    input  logic [31:0]        operand_c_i,
    input  logic [1:0]         vector_mode_i,
    input  logic [4:0]         bmask_a_i,
    input  logic [4:0]         bmask_b_i,
    input  logic [1:0]         imm_vec_ext_i,
    input  logic               is_clpx_i,
    input  logic               is_subrot_i,
    input  logic [1:0]         clpx_shift_i,
    output logic [31:0]        result_o,
    output logic               comparison_result_o,
    output logic               ready_o,
    input  logic               ex_ready_i
);

  // ---- opcodes ----
  localparam logic [6:0] ALU_ADD   = 7'h18;
  localparam logic [6:0] ALU_SUB   = 7'h19;
  localparam logic [6:0] ALU_ADDU  = 7'h1a;
  localparam logic [6:0] ALU_SUBU  = 7'h1b;
  localparam logic [6:0] ALU_ADDR  = 7'h1c;
  localparam logic [6:0] ALU_SUBR  = 7'h1d;
  localparam logic [6:0] ALU_ADDUR = 7'h1e;
  localparam logic [6:0] ALU_SUBUR = 7'h1f;
  localparam logic [6:0] ALU_AND   = 7'h15;
  localparam logic [6:0] ALU_OR    = 7'h2e;
  localparam logic [6:0] ALU_XOR   = 7'h2f;
  localparam logic [6:0] ALU_SRA   = 7'h24;
  localparam logic [6:0] ALU_SRL   = 7'h25;
  localparam logic [6:0] ALU_ROR   = 7'h26;
  localparam logic [6:0] ALU_SLL   = 7'h27;
  localparam logic [6:0] ALU_BEXT  = 7'h28;
  localparam logic [6:0] ALU_BEXTU = 7'h29;
  localparam logic [6:0] ALU_BINS  = 7'h2a;
  localparam logic [6:0] ALU_BCLR  = 7'h2b;
  localparam logic [6:0] ALU_BSET  = 7'h2c;
  localparam logic [6:0] ALU_BREV  = 7'h49;
  localparam logic [6:0] ALU_FF1   = 7'h36;
  localparam logic [6:0] ALU_FL1   = 7'h37;
  localparam logic [6:0] ALU_CNT   = 7'h34;
  localparam logic [6:0] ALU_CLB   = 7'h35;
  localparam logic [6:0] ALU_EXTS  = 7'h3e;
  localparam logic [6:0] ALU_EXT   = 7'h3f;
  localparam logic [6:0] ALU_INS   = 7'h2d;
  localparam logic [6:0] ALU_LTS   = 7'h00;
  localparam logic [6:0] ALU_LTU   = 7'h01;
  localparam logic [6:0] ALU_SLTS  = 7'h02;
  localparam logic [6:0] ALU_SLTU  = 7'h03;
  localparam logic [6:0] ALU_LES   = 7'h04;
  localparam logic [6:0] ALU_LEU   = 7'h05;
  localparam logic [6:0] ALU_SLETS = 7'h06;
  localparam logic [6:0] ALU_SLETU = 7'h07;
  localparam logic [6:0] ALU_GTS   = 7'h08;
  localparam logic [6:0] ALU_GTU   = 7'h09;
  localparam logic [6:0] ALU_GES   = 7'h0a;
  localparam logic [6:0] ALU_GEU   = 7'h0b;
  localparam logic [6:0] ALU_EQ    = 7'h0c;
  localparam logic [6:0] ALU_NE    = 7'h0d;
  localparam logic [6:0] ALU_MIN   = 7'h10;
  localparam logic [6:0] ALU_MINU  = 7'h11;
  localparam logic [6:0] ALU_MAX   = 7'h12;
  localparam logic [6:0] ALU_MAXU  = 7'h13;
  localparam logic [6:0] ALU_ABS   = 7'h14;
  localparam logic [6:0] ALU_CLIP  = 7'h16;
  localparam logic [6:0] ALU_CLIPU = 7'h17;
  localparam logic [6:0] ALU_SHUF  = 7'h3a;
  localparam logic [6:0] ALU_SHUF2 = 7'h3b;
  localparam logic [6:0] ALU_PCKLO = 7'h38;
  localparam logic [6:0] ALU_PCKHI = 7'h39;
  localparam logic [6:0] ALU_DIVU  = 7'h30;
  localparam logic [6:0] ALU_DIV   = 7'h31;
  localparam logic [6:0] ALU_REMU  = 7'h32;
  localparam logic [6:0] ALU_REM   = 7'h33;

  localparam logic [1:0] VEC_MODE32 = 2'b00;
  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8  = 2'b11;

  // ---------------------------------------------------------------------------
  // Kogge-Stone parallel-prefix adders (log-depth carry chains)
  // ---------------------------------------------------------------------------
  function automatic logic [32:0] ks_add32(input logic [31:0] x, input logic [31:0] y);
    logic [31:0] p0, p, g, pp, gg, cin;
    int          i, d;
    p0 = x ^ y;
    p  = p0;
    g  = x & y;
    for (d = 0; d < 5; d++) begin
      for (i = 0; i < 32; i++) begin
        if (i >= (1 << d)) begin
          pp[i] = p[i] & p[i-(1<<d)];
          gg[i] = g[i] | (p[i] & g[i-(1<<d)]);
        end else begin
          pp[i] = p[i];
          gg[i] = g[i];
        end
      end
      p = pp;
      g = gg;
    end
    cin = {g[30:0], 1'b0};
    ks_add32 = p0 ^ cin;
  endfunction

  function automatic logic [36:0] ks_add36(input logic [35:0] x, input logic [35:0] y);
    logic [35:0] p0, p, g, pp, gg, cin;
    int          i, d;
    p0 = x ^ y;
    p  = p0;
    g  = x & y;
    for (d = 0; d < 6; d++) begin
      for (i = 0; i < 36; i++) begin
        if (i >= (1 << d)) begin
          pp[i] = p[i] & p[i-(1<<d)];
          gg[i] = g[i] | (p[i] & g[i-(1<<d)]);
        end else begin
          pp[i] = p[i];
          gg[i] = g[i];
        end
      end
      p = pp;
      g = gg;
    end
    cin = {g[34:0], 1'b0};
    ks_add36 = {x[35] ^ y[35] ^ g[35], p0 ^ cin};
  endfunction

  logic [31:0] a, b, c;
  assign a = operand_a_i;
  assign b = operand_b_i;
  assign c = operand_c_i;

  logic [31:0] a_rev, a_neg_rev;
  genvar gi;
  generate
    for (gi = 0; gi < 32; gi++) begin : gen_a_rev
      assign a_rev[gi]     = a[31-gi];
      assign a_neg_rev[gi] = ~a[31-gi];
    end
  endgenerate

  // ---------------------------------------------------------------------------
  // operator predecode
  // ---------------------------------------------------------------------------  logic is_abs, is_clip, is_clipu, b_negate, carry_b;
  logic is_addfam, is_round, is_div;
  logic shift_left, shift_arith, shift_use_round, is_ror;

  assign is_abs   = (operator_i == ALU_ABS);
  assign is_clip  = (operator_i == ALU_CLIP);
  assign is_clipu = (operator_i == ALU_CLIPU);
  assign is_div   = (operator_i[6:2] == 5'b01100);

  assign is_addfam = (operator_i[6:3] == 4'b0011);
  assign is_round  = is_addfam & operator_i[2];

  assign b_negate = (operator_i == ALU_SUB)  | (operator_i == ALU_SUBR) |
                    (operator_i == ALU_SUBU) | (operator_i == ALU_SUBUR) | is_subrot_i;

  assign carry_b  = b_negate | is_abs | is_clip;

  assign shift_left  = (operator_i == ALU_SLL) | (operator_i == ALU_BINS) |
                       (operator_i == ALU_FL1) | (operator_i == ALU_CLB)  |
                       (operator_i == ALU_BREV) | is_div;
  assign shift_arith = (operator_i == ALU_SRA) | (operator_i == ALU_BEXT) |
                       (operator_i == ALU_ADD) | (operator_i == ALU_SUB) |
                       (operator_i == ALU_ADDR)| (operator_i == ALU_SUBR);
  assign shift_use_round = is_addfam;
  assign is_ror = (operator_i == ALU_ROR);

  // ---------------------------------------------------------------------------
  // bit mask
  // ---------------------------------------------------------------------------
  // bmask = run of (bmask_a_i+1) ones starting at bit bmask_b_i:
  //   bit j is set iff  bmask_b_i <= j <= bmask_a_i + bmask_b_i
  // implemented with two thermometer decodes instead of two barrel shifters.
  logic [31:0] bmask, bmask_inv;
  logic [ 5:0] bmask_sum;
  assign bmask_sum = {1'b0, bmask_a_i} + {1'b0, bmask_b_i};
  generate
    for (gi = 0; gi < 32; gi++) begin : gen_bmask
      assign bmask[gi] = (bmask_b_i <= gi) & (bmask_sum >= gi);
    end
  endgenerate
  assign bmask_inv = ~bmask;

  // ---------------------------------------------------------------------------
  // adder
  // ---------------------------------------------------------------------------
  logic [31:0] op_a_eff, op_b_eff;
  logic [35:0] add_a, add_b;
  logic signed [36:0] add_full;
  logic [31:0] add_out, round_val, sum;

  assign op_a_eff = is_abs ? ~a : (is_subrot_i ? {b[15:0], a[31:16]} : a);
  assign op_b_eff = b_negate ? (is_subrot_i ? ~{a[15:0], b[31:16]} : ~b) : b;

  always_comb begin
    add_a = {op_a_eff[31:24], 1'b1, op_a_eff[23:16], 1'b1,
             op_a_eff[15:8],  1'b1, op_a_eff[7:0],   1'b1};
    add_b = {op_b_eff[31:24], 1'b0, op_b_eff[23:16], 1'b0,
             op_b_eff[15:8],  1'b0, op_b_eff[7:0],   1'b0};

    if (carry_b) begin
      add_b[0] = 1'b1;
      if (vector_mode_i == VEC_MODE16) begin
        add_b[18] = 1'b1;
      end else if (vector_mode_i == VEC_MODE8) begin
        add_b[9]  = 1'b1;
        add_b[18] = 1'b1;
        add_b[27] = 1'b1;
      end
    end else begin
      if (vector_mode_i == VEC_MODE16) begin
        add_a[18] = 1'b0;
      end else if (vector_mode_i == VEC_MODE8) begin
        add_a[9]  = 1'b0;
        add_a[18] = 1'b0;
        add_a[27] = 1'b0;
      end
    end
  end

  assign add_full = ks_add36 (add_a, add_b);
  assign add_out  = {add_full[35:28], add_full[26:19], add_full[17:10], add_full[8:1]};

  assign round_val = is_round ? {1'b0, bmask[31:1]} : 32'b0;
  assign sum       = ks_add32 (add_out, round_val);

  // ---------------------------------------------------------------------------
  // shifter
  // ---------------------------------------------------------------------------
  logic [31:0] shift_amt_in, shift_amt_left, shift_amt_eff, shift_amt_int;
  logic [31:0] shift_amt_norm;
  logic [31:0] shift_op_a;
  logic [63:0] shift_op_a_64;
  logic [31:0] shift_right_result, shift_left_result, shift_result;
  logic [ 5:0] div_shift;
  logic        div_valid;

  assign shift_amt_in = div_valid ? {26'b0, div_shift} : b;

  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_amt_left[15:0]  = shift_amt_in[31:16];
        shift_amt_left[31:16] = shift_amt_in[15:0];
      end
      VEC_MODE8: begin
        shift_amt_left[7:0]   = shift_amt_in[31:24];
        shift_amt_left[15:8]  = shift_amt_in[23:16];
        shift_amt_left[23:16] = shift_amt_in[15:8];
        shift_amt_left[31:24] = shift_amt_in[7:0];
      end
      default: shift_amt_left = shift_amt_in;
    endcase
  end

  assign shift_amt_norm = is_clpx_i ? {14'b0, clpx_shift_i, 14'b0, clpx_shift_i}
                                    : {4{3'b000, bmask_b_i}};
  assign shift_amt_eff  = shift_left ? shift_amt_left : shift_amt_in;
  assign shift_amt_int  = shift_use_round ? shift_amt_norm : shift_amt_eff;

  assign shift_op_a = shift_left ? a_rev : (shift_use_round ? sum : a);

  assign shift_op_a_64 = is_ror ? {shift_op_a, shift_op_a}
                                : {{32{shift_arith & shift_op_a[31]}}, shift_op_a};

  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_right_result[31:16] = $signed({shift_arith & shift_op_a[31], shift_op_a[31:16]})
                                    >>> shift_amt_int[19:16];
        shift_right_result[15:0]  = $signed({shift_arith & shift_op_a[15], shift_op_a[15:0]})
                                    >>> shift_amt_int[3:0];
      end
      VEC_MODE8: begin
        shift_right_result[31:24] = $signed({shift_arith & shift_op_a[31], shift_op_a[31:24]})
                                    >>> shift_amt_int[26:24];
        shift_right_result[23:16] = $signed({shift_arith & shift_op_a[23], shift_op_a[23:16]})
                                    >>> shift_amt_int[18:16];
        shift_right_result[15:8]  = $signed({shift_arith & shift_op_a[15], shift_op_a[15:8]})
                                    >>> shift_amt_int[10:8];
        shift_right_result[7:0]   = $signed({shift_arith & shift_op_a[7], shift_op_a[7:0]})
                                    >>> shift_amt_int[2:0];
      end
      default: shift_right_result = shift_op_a_64 >> shift_amt_int[4:0];
    endcase
  end

  generate
    for (gi = 0; gi < 32; gi++) begin : gen_shift_left_result
      assign shift_left_result[gi] = shift_right_result[31-gi];
    end
  endgenerate

  assign shift_result = shift_left ? shift_left_result : shift_right_result;

  // ---------------------------------------------------------------------------
  // comparisons
  // ---------------------------------------------------------------------------
  logic        cmp_sgn;
  logic [ 3:0] sgn_mask;
  logic [ 3:0] byte_eq, byte_gt, lane_eq, lane_gt;
  logic        cmp_use_eq, cmp_eq_only, cmp_inv;
  logic [ 3:0] lane_flag;

  assign cmp_sgn = (operator_i == ALU_GTS)  | (operator_i == ALU_GES)  |
                   (operator_i == ALU_LTS)  | (operator_i == ALU_LES)  |
                   (operator_i == ALU_SLTS) | (operator_i == ALU_SLETS)|
                   (operator_i == ALU_MIN)  | (operator_i == ALU_MAX)  |
                   is_abs | is_clip | is_clipu;

  always_comb begin
    if (cmp_sgn) begin
      case (vector_mode_i)
        VEC_MODE8:  sgn_mask = 4'b1111;
        VEC_MODE16: sgn_mask = 4'b1010;
        default:    sgn_mask = 4'b1000;
      endcase
    end else begin
      sgn_mask = 4'b0000;
    end
  end

  generate
    for (gi = 0; gi < 4; gi++) begin : gen_cmp_byte
      assign byte_eq[gi] = (a[8*gi +: 8] == b[8*gi +: 8]);
      assign byte_gt[gi] = $signed({a[8*gi+7] & sgn_mask[gi], a[8*gi +: 8]}) >
                           $signed({b[8*gi+7] & sgn_mask[gi], b[8*gi +: 8]});
    end
  endgenerate

  always_comb begin
    lane_eq = {4{byte_eq[3] & byte_eq[2] & byte_eq[1] & byte_eq[0]}};
    lane_gt = {4{byte_gt[3] | (byte_eq[3] & (byte_gt[2] | (byte_eq[2] &
                    (byte_gt[1] | (byte_eq[1] & byte_gt[0])))))}};
    case (vector_mode_i)
      VEC_MODE16: begin
        lane_eq[1:0] = {2{byte_eq[0] & byte_eq[1]}};
        lane_eq[3:2] = {2{byte_eq[2] & byte_eq[3]}};
        lane_gt[1:0] = {2{byte_gt[1] | (byte_eq[1] & byte_gt[0])}};
        lane_gt[3:2] = {2{byte_gt[3] | (byte_eq[3] & byte_gt[2])}};
      end
      VEC_MODE8: begin
        lane_eq = byte_eq;
        lane_gt = byte_gt;
      end
      default: ;
    endcase
  end

  assign cmp_eq_only = (operator_i == ALU_EQ) | (operator_i == ALU_NE);
  assign cmp_use_eq  = (operator_i == ALU_GES)  | (operator_i == ALU_GEU)  |
                       (operator_i == ALU_LTS)  | (operator_i == ALU_LTU)  |
                       (operator_i == ALU_SLTS) | (operator_i == ALU_SLTU);
  assign cmp_inv     = (operator_i == ALU_LES)  | (operator_i == ALU_LEU)  |
                       (operator_i == ALU_LTS)  | (operator_i == ALU_LTU)  |
                       (operator_i == ALU_SLTS) | (operator_i == ALU_SLTU) |
                       (operator_i == ALU_SLETS)| (operator_i == ALU_SLETU)|
                       (operator_i == ALU_NE);

  generate
    for (gi = 0; gi < 4; gi++) begin : gen_lane_flag
      assign lane_flag[gi] = (cmp_eq_only ? lane_eq[gi]
                                          : (lane_gt[gi] | (lane_eq[gi] & cmp_use_eq))) ^ cmp_inv;
    end
  endgenerate

  assign comparison_result_o = lane_flag[3];

  // ---------------------------------------------------------------------------
  // min / max / abs / clip
  // ---------------------------------------------------------------------------
  logic        do_min;
  logic [ 3:0] sel_mm;
  logic [31:0] mm_b, mm_res, clip_res;
  logic        clip_lo;

  assign do_min = (operator_i == ALU_MIN) | (operator_i == ALU_MINU) | is_clip | is_clipu;
  assign sel_mm = lane_gt ^ {4{do_min}};
  assign mm_b   = is_abs ? add_out : b;

  generate
    for (gi = 0; gi < 4; gi++) begin : gen_mm
      assign mm_res[8*gi +: 8] = sel_mm[gi] ? a[8*gi +: 8] : mm_b[8*gi +: 8];
    end
  endgenerate

  assign clip_lo  = is_clipu ? (a[31] | (a == 32'b0)) : (add_full[36] | (a == ~b));
  assign clip_res = clip_lo ? (is_clipu ? 32'b0 : ~b) : mm_res;

  // ---------------------------------------------------------------------------
  // bit manipulation
  // ---------------------------------------------------------------------------
  logic [31:0] bextins_res, bclr_res, bset_res;
  logic        ext_sign;

  assign ext_sign    = (operator_i == ALU_BEXT) & shift_result[bmask_a_i];
  assign bextins_res = (bmask & shift_result) |
                       (((operator_i == ALU_BINS) ? c : {32{ext_sign}}) & bmask_inv);
  assign bclr_res    = a & bmask_inv;
  assign bset_res    = a | bmask;

  // ---------------------------------------------------------------------------
  // bit reverse
  // ---------------------------------------------------------------------------
  logic [31:0] radix2_rev, radix4_rev, radix8_rev, reverse_result;

  generate
    for (gi = 0; gi < 32; gi++) begin : gen_radix2
      assign radix2_rev[gi] = shift_result[31-gi];
    end
    for (gi = 0; gi < 16; gi++) begin : gen_radix4
      assign radix4_rev[2*gi+1:2*gi] = shift_result[31-2*gi -: 2];
    end
    for (gi = 0; gi < 10; gi++) begin : gen_radix8
      assign radix8_rev[3*gi+2:3*gi] = shift_result[31-3*gi -: 3];
    end
    assign radix8_rev[31:30] = 2'b0;
  endgenerate

  always_comb begin
    case (bmask_a_i[1:0])
      2'b01:   reverse_result = radix4_rev;
      2'b10:   reverse_result = radix8_rev;
      default: reverse_result = radix2_rev;
    endcase
  end

  // ---------------------------------------------------------------------------
  // shuffle / pack / extract / insert
  // ---------------------------------------------------------------------------
  logic [3:0][1:0] shuf_byte_sel;
  logic [3:0]      shuf_reg_sel;
  logic [1:0]      shuf_reg1_sel, shuf_reg0_sel;
  logic [3:0]      shuf_through;
  logic [31:0]     shuf_r0_in, shuf_r1_in, shuf_r0, shuf_r1, shuf_res, pack_result;

  always_comb begin
    shuf_reg_sel  = 4'b0;
    shuf_reg1_sel = 2'b01;
    shuf_reg0_sel = 2'b10;
    shuf_through  = 4'b1111;

    case (operator_i)
      ALU_EXT, ALU_EXTS: begin
        if (operator_i == ALU_EXTS) shuf_reg1_sel = 2'b11;
        if (vector_mode_i == VEC_MODE8) begin
          shuf_reg_sel[3:1] = 3'b111;
          shuf_reg_sel[0]   = 1'b0;
        end else begin
          shuf_reg_sel[3:2] = 2'b11;
          shuf_reg_sel[1:0] = 2'b00;
        end
      end

      ALU_PCKLO: begin
        shuf_reg1_sel = 2'b00;
        if (vector_mode_i == VEC_MODE8) begin
          shuf_through = 4'b0011;
          shuf_reg_sel = 4'b0001;
        end else begin
          shuf_reg_sel = 4'b0011;
        end
      end

      ALU_PCKHI: begin
        shuf_reg1_sel = 2'b00;
        if (vector_mode_i == VEC_MODE8) begin
          shuf_through = 4'b1100;
          shuf_reg_sel = 4'b0100;
        end else begin
          shuf_reg_sel = 4'b0011;
        end
      end

      ALU_SHUF2: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_reg_sel[3] = ~b[26];
          shuf_reg_sel[2] = ~b[18];
          shuf_reg_sel[1] = ~b[10];
          shuf_reg_sel[0] = ~b[2];
        end else begin
          shuf_reg_sel[3] = ~b[17];
          shuf_reg_sel[2] = ~b[17];
          shuf_reg_sel[1] = ~b[1];
          shuf_reg_sel[0] = ~b[1];
        end
      end

      ALU_INS: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_reg0_sel = 2'b00;
          case (imm_vec_ext_i)
            2'b00:   shuf_reg_sel = 4'b1110;
            2'b01:   shuf_reg_sel = 4'b1101;
            2'b10:   shuf_reg_sel = 4'b1011;
            default: shuf_reg_sel = 4'b0111;
          endcase
        end else begin
          shuf_reg0_sel   = 2'b01;
          shuf_reg_sel[3] = ~imm_vec_ext_i[0];
          shuf_reg_sel[2] = ~imm_vec_ext_i[0];
          shuf_reg_sel[1] = imm_vec_ext_i[0];
          shuf_reg_sel[0] = imm_vec_ext_i[0];
        end
      end

      default: ;
    endcase
  end

  always_comb begin
    shuf_byte_sel = '0;
    case (operator_i)
      ALU_EXTS, ALU_EXT: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_byte_sel[3] = imm_vec_ext_i;
          shuf_byte_sel[2] = imm_vec_ext_i;
          shuf_byte_sel[1] = imm_vec_ext_i;
          shuf_byte_sel[0] = imm_vec_ext_i;
        end else begin
          shuf_byte_sel[3] = {imm_vec_ext_i[0], 1'b1};
          shuf_byte_sel[2] = {imm_vec_ext_i[0], 1'b1};
          shuf_byte_sel[1] = {imm_vec_ext_i[0], 1'b1};
          shuf_byte_sel[0] = {imm_vec_ext_i[0], 1'b0};
        end
      end

      ALU_PCKLO: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_byte_sel = 4'b0;
        end else begin
          shuf_byte_sel[3] = 2'b01;
          shuf_byte_sel[2] = 2'b00;
          shuf_byte_sel[1] = 2'b01;
          shuf_byte_sel[0] = 2'b00;
        end
      end

      ALU_PCKHI: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_byte_sel = 4'b0;
        end else begin
          shuf_byte_sel[3] = 2'b11;
          shuf_byte_sel[2] = 2'b10;
          shuf_byte_sel[1] = 2'b11;
          shuf_byte_sel[0] = 2'b10;
        end
      end

      ALU_SHUF2, ALU_SHUF: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_byte_sel[3] = b[25:24];
          shuf_byte_sel[2] = b[17:16];
          shuf_byte_sel[1] = b[9:8];
          shuf_byte_sel[0] = b[1:0];
        end else begin
          shuf_byte_sel[3] = {b[16], 1'b1};
          shuf_byte_sel[2] = {b[16], 1'b0};
          shuf_byte_sel[1] = {b[0], 1'b1};
          shuf_byte_sel[0] = {b[0], 1'b0};
        end
      end

      ALU_INS: begin
        shuf_byte_sel[3] = 2'b11;
        shuf_byte_sel[2] = 2'b10;
        shuf_byte_sel[1] = 2'b01;
        shuf_byte_sel[0] = 2'b00;
      end

      default: ;
    endcase
  end

  assign shuf_r0_in = shuf_reg0_sel[1] ? a :
                      (shuf_reg0_sel[0] ? {2{a[15:0]}} : {4{a[7:0]}});
  assign shuf_r1_in = shuf_reg1_sel[1] ?
                      {{8{a[31]}}, {8{a[23]}}, {8{a[15]}}, {8{a[7]}}} :
                      (shuf_reg1_sel[0] ? c : b);

  generate
    for (gi = 0; gi < 4; gi++) begin : gen_shuf
      assign shuf_r0[8*gi +: 8] = shuf_byte_sel[gi][1] ?
                                  (shuf_byte_sel[gi][0] ? shuf_r0_in[31:24] : shuf_r0_in[23:16]) :
                                  (shuf_byte_sel[gi][0] ? shuf_r0_in[15:8]  : shuf_r0_in[7:0]);
      assign shuf_r1[8*gi +: 8] = shuf_byte_sel[gi][1] ?
                                  (shuf_byte_sel[gi][0] ? shuf_r1_in[31:24] : shuf_r1_in[23:16]) :
                                  (shuf_byte_sel[gi][0] ? shuf_r1_in[15:8]  : shuf_r1_in[7:0]);
      assign shuf_res[8*gi +: 8] = shuf_reg_sel[gi] ? shuf_r1[8*gi +: 8] : shuf_r0[8*gi +: 8];
      assign pack_result[8*gi +: 8] = shuf_through[gi] ? shuf_res[8*gi +: 8] : c[8*gi +: 8];
    end
  endgenerate

  // ---------------------------------------------------------------------------
  // bit counting
  // ---------------------------------------------------------------------------
  logic [31:0] ff_input;
  logic [ 4:0] ff1_result;
  logic        ff_no_one;
  logic [ 4:0] clb5_result;
  logic [ 5:0] clb6_result;
  logic [ 4:0] fl1_result;
  logic [ 5:0] bitop_result;
  logic [ 5:0] cnt_result;

  assign ff_input = (operator_i == ALU_FF1) ? a :
                    ((operator_i == ALU_FL1) | (operator_i == ALU_DIVU) |
                     (operator_i == ALU_REMU)) ? a_rev :
                    (a[31] ? a_neg_rev : a_rev);

  cv32e40p_ff_one ff_one_i (
      .in_i        (ff_input),
      .first_one_o (ff1_result),
      .no_ones_o   (ff_no_one)
  );

  cv32e40p_popcnt popcnt_i (
      .in_i     (a),
      .result_o (cnt_result)
  );

  assign clb5_result = ff1_result - 5'd1;
  assign clb6_result = {1'b0, ff1_result} - 6'd1;
  assign fl1_result  = 5'd31 - ff1_result;

  always_comb begin
    bitop_result = 6'b0;
    case (operator_i)
      ALU_FF1: bitop_result = ff_no_one ? 6'd32 : {1'b0, ff1_result};
      ALU_FL1: bitop_result = ff_no_one ? 6'd32 : {1'b0, fl1_result};
      ALU_CNT: bitop_result = cnt_result;
      ALU_CLB: bitop_result = ff_no_one ? (a[31] ? 6'd31 : 6'd0) : clb6_result;
      default: ;
    endcase
  end

  // ---------------------------------------------------------------------------
  // serial divider
  // ---------------------------------------------------------------------------
  logic [31:0] result_div;
  logic        div_ready, div_signed, div_op_a_signed;
  logic [ 5:0] div_shift_int;

  assign div_signed      = operator_i[0];
  assign div_op_a_signed = a[31] & div_signed;
  assign div_shift_int   = ff_no_one ? 6'd31 : clb6_result;
  assign div_shift       = div_shift_int + (div_op_a_signed ? 6'd0 : 6'd1);
  assign div_valid       = enable_i & is_div;

  cv32e40p_alu_div alu_div_i (
      .Clk_CI      (clk),
      .Rst_RBI     (rst_n),
      .OpA_DI      (b),
      .OpB_DI      (shift_left_result),
      .OpBShift_DI (div_shift),
      .OpBIsZero_SI(cnt_result == 6'd0),
      .OpBSign_SI  (div_op_a_signed),
      .OpCode_SI   (operator_i[1:0]),
      .Res_DO      (result_div),
      .InVld_SI    (div_valid),
      .OutRdy_SI   (ex_ready_i),
      .OutVld_SO   (div_ready)
  );

  assign ready_o = div_ready;

  // ---------------------------------------------------------------------------
  // result mux
  // ---------------------------------------------------------------------------
  always_comb begin
    case (operator_i)
      ALU_AND: result_o = a & b;
      ALU_OR:  result_o = a | b;
      ALU_XOR: result_o = a ^ b;

      ALU_ADD, ALU_SUB, ALU_ADDU, ALU_SUBU,
      ALU_ADDR, ALU_SUBR, ALU_ADDUR, ALU_SUBUR,
      ALU_SLL, ALU_SRL, ALU_SRA, ALU_ROR: result_o = shift_result;

      ALU_BINS, ALU_BEXT, ALU_BEXTU: result_o = bextins_res;
      ALU_BCLR:                      result_o = bclr_res;
      ALU_BSET:                      result_o = bset_res;
      ALU_BREV:                      result_o = reverse_result;

      ALU_SHUF, ALU_SHUF2, ALU_PCKLO, ALU_PCKHI,
      ALU_EXT, ALU_EXTS, ALU_INS: result_o = pack_result;

      ALU_MIN, ALU_MINU, ALU_MAX, ALU_MAXU: result_o = mm_res;

      ALU_ABS: result_o = is_clpx_i ? {add_out[31:16], a[15:0]} : mm_res;

      ALU_CLIP, ALU_CLIPU: result_o = clip_res;

      ALU_EQ, ALU_NE, ALU_GTS, ALU_GTU, ALU_GES, ALU_GEU,
      ALU_LTS, ALU_LTU, ALU_LES, ALU_LEU:
        result_o = {{8{lane_flag[3]}}, {8{lane_flag[2]}}, {8{lane_flag[1]}}, {8{lane_flag[0]}}};

      ALU_SLTS, ALU_SLTU, ALU_SLETS, ALU_SLETU: result_o = {31'b0, lane_flag[3]};

      ALU_FF1, ALU_FL1, ALU_CLB, ALU_CNT: result_o = {26'b0, bitop_result};

      ALU_DIV, ALU_DIVU, ALU_REM, ALU_REMU: result_o = result_div;

      default: result_o = 32'b0;
    endcase
  end

endmodule


// -----------------------------------------------------------------------------
// popcount
// -----------------------------------------------------------------------------
module cv32e40p_popcnt (
    input  logic [31:0] in_i,
    output logic [ 5:0] result_o
);
  logic [15:0][1:0] cnt_l1;
  logic [ 7:0][2:0] cnt_l2;
  logic [ 3:0][3:0] cnt_l3;
  logic [ 1:0][4:0] cnt_l4;
  genvar l;
  generate
    for (l = 0; l < 16; l++) begin : gen_l1
      assign cnt_l1[l] = {1'b0, in_i[2*l]} + {1'b0, in_i[2*l+1]};
    end
    for (l = 0; l < 8; l++) begin : gen_l2
      assign cnt_l2[l] = {1'b0, cnt_l1[2*l]} + {1'b0, cnt_l1[2*l+1]};
    end
    for (l = 0; l < 4; l++) begin : gen_l3
      assign cnt_l3[l] = {1'b0, cnt_l2[2*l]} + {1'b0, cnt_l2[2*l+1]};
    end
    for (l = 0; l < 2; l++) begin : gen_l4
      assign cnt_l4[l] = {1'b0, cnt_l3[2*l]} + {1'b0, cnt_l3[2*l+1]};
    end
  endgenerate
  assign result_o = {1'b0, cnt_l4[0]} + {1'b0, cnt_l4[1]};
endmodule


// -----------------------------------------------------------------------------
// find first one
// -----------------------------------------------------------------------------
module cv32e40p_ff_one (
    input  logic [31:0] in_i,
    output logic [ 4:0] first_one_o,
    output logic        no_ones_o
);
  logic [15:0] v1, b0;
  logic [ 7:0] v2, b1;
  logic [ 3:0] v3, b2;
  logic [ 1:0] v4, b3;
  logic        v5, b4;
  genvar k;
  generate
    for (k = 0; k < 16; k++) begin : gen_l1
      assign v1[k] = in_i[2*k] | in_i[2*k+1];
      assign b0[k] = ~in_i[2*k];
    end
    for (k = 0; k < 8; k++) begin : gen_l2
      assign v2[k] = v1[2*k] | v1[2*k+1];
      assign b1[k] = ~v1[2*k];
    end
    for (k = 0; k < 4; k++) begin : gen_l3
      assign v3[k] = v2[2*k] | v2[2*k+1];
      assign b2[k] = ~v2[2*k];
    end
    for (k = 0; k < 2; k++) begin : gen_l4
      assign v4[k] = v3[2*k] | v3[2*k+1];
      assign b3[k] = ~v3[2*k];
    end
    assign v5 = v4[0] | v4[1];
    assign b4 = ~v4[0];
  endgenerate

  assign first_one_o = {b4,
                        v4[0] ? b3[0] : b3[1],
                        v3[0] ? b2[0] : (v3[1] ? b2[1] : (v3[2] ? b2[2] : b2[3])),
                        v2[0] ? b1[0] : (v2[1] ? b1[1] : (v2[2] ? b1[2] : (v2[3] ? b1[3] :
                                (v2[4] ? b1[4] : (v2[5] ? b1[5] : (v2[6] ? b1[6] : b1[7])))))),
                        v1[0] ? b0[0] : (v1[1] ? b0[1] : (v1[2] ? b0[2] : (v1[3] ? b0[3] :
                                (v1[4] ? b0[4] : (v1[5] ? b0[5] : (v1[6] ? b0[6] : (v1[7] ? b0[7] :
                                (v1[8] ? b0[8] : (v1[9] ? b0[9] : (v1[10] ? b0[10] : (v1[11] ? b0[11] :
                                (v1[12] ? b0[12] : (v1[13] ? b0[13] : (v1[14] ? b0[14] : b0[15]))))))))))))))};
  assign no_ones_o   = ~v5;
endmodule


// -----------------------------------------------------------------------------
// simple serial divider (cycle behaviour identical to upstream)
// -----------------------------------------------------------------------------
module cv32e40p_alu_div (
    input  logic        Clk_CI,
    input  logic        Rst_RBI,
    input  logic [31:0] OpA_DI,
    input  logic [31:0] OpB_DI,
    input  logic [ 5:0] OpBShift_DI,
    input  logic        OpBIsZero_SI,
    input  logic        OpBSign_SI,
    input  logic [ 1:0] OpCode_SI,
    input  logic        InVld_SI,
    input  logic        OutRdy_SI,
    output logic        OutVld_SO,
    output logic [31:0] Res_DO
);
  logic [31:0] ResReg_DP, ResReg_DN, ResReg_DP_rev;
  logic [31:0] AReg_DP, AReg_DN;
  logic [31:0] BReg_DP, BReg_DN;

  logic RemSel_SP, CompInv_SP, ResInv_SP;
  logic RemSel_SN, CompInv_SN, ResInv_SN;

  logic [31:0] AddMux_D, AddOut_D, AddTmp_D, BMux_D, OutMux_D;
  logic [ 5:0] Cnt_DP, Cnt_DN;
  logic        CntZero_S;

  logic ARegEn_S, BRegEn_S, ResRegEn_S, ABComp_S, PmSel_S, LoadEn_S;

  typedef enum logic [1:0] {IDLE, DIVIDE, FINISH} state_t;
  state_t State_SN, State_SP;

  assign PmSel_S  = LoadEn_S & ~(OpCode_SI[0] & (OpA_DI[31] ^ OpBSign_SI));
  assign AddMux_D = (LoadEn_S) ? OpA_DI : BReg_DP;
  assign BMux_D   = (LoadEn_S) ? OpB_DI : {CompInv_SP, BReg_DP[31:1]};

  genvar index;
  generate
    for (index = 0; index < 32; index++) begin : gen_bit_swapping
      assign ResReg_DP_rev[index] = ResReg_DP[31-index];
    end
  endgenerate

  assign OutMux_D = (RemSel_SP) ? AReg_DP : ResReg_DP_rev;
  assign Res_DO   = (ResInv_SP) ? -$signed(OutMux_D) : OutMux_D;

  assign ABComp_S = ((AReg_DP == BReg_DP) | ((AReg_DP > BReg_DP) ^ CompInv_SP)) &
                    ((|AReg_DP) | OpBIsZero_SI);

  assign AddTmp_D = (LoadEn_S) ? 32'b0 : AReg_DP;
  assign AddOut_D = (PmSel_S) ? AddTmp_D + AddMux_D : AddTmp_D - $signed(AddMux_D);

  assign Cnt_DN    = (LoadEn_S) ? OpBShift_DI : ((~CntZero_S) ? Cnt_DP - 6'd1 : Cnt_DP);
  assign CntZero_S = ~(|Cnt_DP);

  always_comb begin
    State_SN   = State_SP;
    OutVld_SO  = 1'b0;
    LoadEn_S   = 1'b0;
    ARegEn_S   = 1'b0;
    BRegEn_S   = 1'b0;
    ResRegEn_S = 1'b0;

    case (State_SP)
      IDLE: begin
        OutVld_SO = 1'b1;
        if (InVld_SI) begin
          OutVld_SO = 1'b0;
          ARegEn_S  = 1'b1;
          BRegEn_S  = 1'b1;
          LoadEn_S  = 1'b1;
          State_SN  = DIVIDE;
        end
      end
      DIVIDE: begin
        ARegEn_S   = ABComp_S;
        BRegEn_S   = 1'b1;
        ResRegEn_S = 1'b1;
        if (CntZero_S) begin
          State_SN = FINISH;
        end
      end
      FINISH: begin
        OutVld_SO = 1'b1;
        if (OutRdy_SI) begin
          State_SN = IDLE;
        end
      end
      default: ;
    endcase
  end

  assign RemSel_SN  = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN  = (LoadEn_S) ? ((~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] &
                                    (OpA_DI[31] ^ OpBSign_SI)) : ResInv_SP;

  assign AReg_DN   = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN   = (BRegEn_S) ? BMux_D : BReg_DP;
  assign ResReg_DN = (LoadEn_S) ? 32'b0 : ((ResRegEn_S) ? {ABComp_S, ResReg_DP[31:1]} : ResReg_DP);

  always_ff @(posedge Clk_CI or negedge Rst_RBI) begin
    if (~Rst_RBI) begin
      State_SP   <= IDLE;
      AReg_DP    <= 32'b0;
      BReg_DP    <= 32'b0;
      ResReg_DP  <= 32'b0;
      Cnt_DP     <= 6'b0;
      RemSel_SP  <= 1'b0;
      CompInv_SP <= 1'b0;
      ResInv_SP  <= 1'b0;
    end else begin
      State_SP   <= State_SN;
      AReg_DP    <= AReg_DN;
      BReg_DP    <= BReg_DN;
      ResReg_DP  <= ResReg_DN;
      Cnt_DP     <= Cnt_DN;
      RemSel_SP  <= RemSel_SN;
      CompInv_SP <= CompInv_SN;
      ResInv_SP  <= ResInv_SN;
    end
  end
endmodule

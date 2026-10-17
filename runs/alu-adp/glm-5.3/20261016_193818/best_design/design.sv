// Optimized CV32E40P ALU (from starting_point.sv, restructured).
module cv32e40p_alu
(
    input  logic        clk,
    input  logic        rst_n,
    input  logic        enable_i,
    input  logic [6:0]  operator_i,
    input  logic [31:0] operand_a_i,
    input  logic [31:0] operand_b_i,
    input  logic [31:0] operand_c_i,
    input  logic [1:0]  vector_mode_i,
    input  logic [4:0]  bmask_a_i,
    input  logic [4:0]  bmask_b_i,
    input  logic [1:0]  imm_vec_ext_i,
    input  logic        is_clpx_i,
    input  logic        is_subrot_i,
    input  logic [1:0]  clpx_shift_i,
    output logic [31:0] result_o,
    output logic        comparison_result_o,
    output logic        ready_o,
    input  logic        ex_ready_i
);

  // ---- decode (one-hot-ish group signals) ----
  logic [6:0] op;
  assign op = operator_i;

  logic f_add, f_sub, f_round, f_addu, f_subu;    // adder family
  logic f_logic;
  logic f_sll, f_srl, f_sra, f_ror;
  logic f_bext, f_bextu, f_bins, f_bclr, f_bset, f_brev;
  logic f_ff1, f_fl1, f_cnt, f_clb;
  logic f_exts, f_ext, f_ins;
  logic f_lts, f_ltu, f_les, f_leu, f_gts, f_gtu, f_ges, f_geu, f_eq, f_ne;
  logic f_slts, f_sltu, f_slets, f_sletu;
  logic f_min, f_minu, f_max, f_maxu, f_abs, f_clip, f_clipu;
  logic f_shuf, f_shuf2, f_pcklo, f_pckhi;
  logic f_divu, f_div, f_remu, f_rem;

  always_comb begin
    {f_add, f_sub, f_round, f_addu, f_subu} = 5'b0;
    {f_sll, f_srl, f_sra, f_ror} = 4'b0;
    {f_bext, f_bextu, f_bins, f_bclr, f_bset, f_brev} = 6'b0;
    {f_ff1, f_fl1, f_cnt, f_clb} = 4'b0;
    {f_exts, f_ext, f_ins} = 3'b0;
    {f_lts, f_ltu, f_les, f_leu, f_gts, f_gtu, f_ges, f_geu, f_eq, f_ne} = 10'b0;
    {f_slts, f_sltu, f_slets, f_sletu} = 4'b0;
    {f_min, f_minu, f_max, f_maxu, f_abs, f_clip, f_clipu} = 7'b0;
    {f_shuf, f_shuf2, f_pcklo, f_pckhi} = 4'b0;
    {f_divu, f_div, f_remu, f_rem} = 4'b0;
    f_logic = 1'b0;
    case (op)
      7'h18: f_add  = 1'b1;
      7'h19: f_sub  = 1'b1;
      7'h1a: f_addu = 1'b1;
      7'h1b: f_subu = 1'b1;
      7'h1c, 7'h1d: f_round = 1'b1;   // ADDR/SUBR  (rounded, signed)
      7'h1e, 7'h1f: f_round = 1'b1;   // ADDUR/SUBUR (rounded, unsigned)
      7'h15: f_logic = 1'b1;          // AND
      7'h2e: f_logic = 1'b1;          // OR
      7'h2f: f_logic = 1'b1;          // XOR
      7'h24: f_sra  = 1'b1;
      7'h25: f_srl  = 1'b1;
      7'h26: f_ror  = 1'b1;
      7'h27: f_sll  = 1'b1;
      7'h28: f_bext = 1'b1;
      7'h29: f_bextu= 1'b1;
      7'h2a: f_bins = 1'b1;
      7'h2b: f_bclr = 1'b1;
      7'h2c: f_bset = 1'b1;
      7'h49: f_brev = 1'b1;
      7'h36: f_ff1  = 1'b1;
      7'h37: f_fl1  = 1'b1;
      7'h34: f_cnt  = 1'b1;
      7'h35: f_clb  = 1'b1;
      7'h3e: f_exts = 1'b1;
      7'h3f: f_ext  = 1'b1;
      7'h2d: f_ins  = 1'b1;
      7'h00: f_lts  = 1'b1;
      7'h01: f_ltu  = 1'b1;
      7'h04: f_les  = 1'b1;
      7'h05: f_leu  = 1'b1;
      7'h08: f_gts  = 1'b1;
      7'h09: f_gtu  = 1'b1;
      7'h0a: f_ges  = 1'b1;
      7'h0b: f_geu  = 1'b1;
      7'h0c: f_eq   = 1'b1;
      7'h0d: f_ne   = 1'b1;
      7'h02: f_slts = 1'b1;
      7'h03: f_sltu = 1'b1;
      7'h06: f_slets= 1'b1;
      7'h07: f_sletu= 1'b1;
      7'h10: f_min  = 1'b1;
      7'h11: f_minu = 1'b1;
      7'h12: f_max  = 1'b1;
      7'h13: f_maxu = 1'b1;
      7'h14: f_abs  = 1'b1;
      7'h16: f_clip = 1'b1;
      7'h17: f_clipu= 1'b1;
      7'h3a: f_shuf = 1'b1;
      7'h3b: f_shuf2= 1'b1;
      7'h38: f_pcklo= 1'b1;
      7'h39: f_pckhi= 1'b1;
      7'h30: f_divu = 1'b1;
      7'h31: f_div  = 1'b1;
      7'h32: f_remu = 1'b1;
      7'h33: f_rem  = 1'b1;
      default: ;
    endcase
  end

  logic f_addfam, f_addsub, f_rounding, f_minmax, f_bitcnt, f_shufgrp, f_divgrp, f_veccmp, f_scalcmp;
  assign f_addfam   = f_add | f_sub | f_addu | f_subu | f_round;
  assign f_addsub   = f_add | f_sub;                 // ADD/SUB: arithmetic post-shift
  assign f_rounding = f_round;
  assign f_minmax   = f_min | f_minu | f_max | f_maxu;
  assign f_bitcnt   = f_ff1 | f_fl1 | f_cnt | f_clb;
  assign f_shufgrp  = f_shuf | f_shuf2 | f_pcklo | f_pckhi | f_exts | f_ext | f_ins;
  assign f_divgrp   = f_divu | f_div | f_remu | f_rem;
  assign f_veccmp   = f_lts | f_ltu | f_les | f_leu | f_gts | f_gtu | f_ges | f_geu | f_eq | f_ne;
  assign f_scalcmp  = f_slts | f_sltu | f_slets | f_sletu;

  logic f_bextins, f_bmaskops;
  assign f_bextins  = f_bext | f_bextu | f_bins;
  assign f_bmaskops = f_bextins | f_bclr | f_bset;

  // ---- vector mode ----
  logic vm16, vm8;
  assign vm16 = (vector_mode_i == 2'b10);
  assign vm8  = (vector_mode_i == 2'b11);

  // ---- operand preparation ----
  logic [31:0] operand_a_rev, operand_a_neg, operand_a_neg_rev, operand_b_neg;
  genvar k;
  generate
    for (k = 0; k < 32; k++) begin : gen_rev
      assign operand_a_rev[k]     = operand_a_i[31-k];
      assign operand_a_neg_rev[k] = operand_a_i[31-k] ? 1'b0 : 1'b1;
    end
  endgenerate
  assign operand_a_neg = ~operand_a_i;
  assign operand_b_neg = ~operand_b_i;

  //////////////////////////////////////////////////////////////////////////////////////
  // ADDER
  //////////////////////////////////////////////////////////////////////////////////////
  logic        adder_op_b_negate;
  logic [31:0] adder_op_a, adder_op_b;
  logic [35:0] adder_in_a, adder_in_b;
  logic [31:0] adder_result;
  logic [36:0] adder_result_expanded;

  assign adder_op_b_negate = (f_addfam & op[0]) | is_subrot_i;   // SUB/SUBU/SUBR/SUBUR

  assign adder_op_a = f_abs ? operand_a_neg : (is_subrot_i ? {operand_b_i[15:0], operand_a_i[31:16]} : operand_a_i);
  assign adder_op_b = adder_op_b_negate ? (is_subrot_i ? ~{operand_a_i[15:0], operand_b_i[31:16]} : operand_b_neg)
                                        : operand_b_i;

  always_comb begin
    adder_in_a[0]     = 1'b1;
    adder_in_a[8:1]   = adder_op_a[7:0];
    adder_in_a[9]     = 1'b1;
    adder_in_a[17:10] = adder_op_a[15:8];
    adder_in_a[18]    = 1'b1;
    adder_in_a[26:19] = adder_op_a[23:16];
    adder_in_a[27]    = 1'b1;
    adder_in_a[35:28] = adder_op_a[31:24];

    adder_in_b[0]     = 1'b0;
    adder_in_b[8:1]   = adder_op_b[7:0];
    adder_in_b[9]     = 1'b0;
    adder_in_b[17:10] = adder_op_b[15:8];
    adder_in_b[18]    = 1'b0;
    adder_in_b[26:19] = adder_op_b[23:16];
    adder_in_b[27]    = 1'b0;
    adder_in_b[35:28] = adder_op_b[31:24];

    if (adder_op_b_negate || f_abs || f_clip) begin
      adder_in_b[0] = 1'b1;
      case (vector_mode_i)
        2'b10: adder_in_b[18] = 1'b1;
        2'b11: begin adder_in_b[9] = 1'b1; adder_in_b[18] = 1'b1; adder_in_b[27] = 1'b1; end
        default: ;
      endcase
    end else begin
      case (vector_mode_i)
        2'b10: adder_in_a[18] = 1'b0;
        2'b11: begin adder_in_a[9] = 1'b0; adder_in_a[18] = 1'b0; adder_in_a[27] = 1'b0; end
        default: ;
      endcase
    end
  end

  assign adder_result_expanded = sk36(adder_in_a, adder_in_b);
  logic [36:0] adder_sum_expanded;
  assign adder_sum_expanded = adder_result_expanded;
  assign adder_result = {adder_result_expanded[35:28], adder_result_expanded[26:19],
                         adder_result_expanded[17:10], adder_result_expanded[8:1]};

  //////////////////////////////////////////////////////////////////////////////////////
  // normalization (round + post-shift) for the adder family
  //////////////////////////////////////////////////////////////////////////////////////
  logic [31:0] bmask;
  assign bmask = (~(32'hFFFFFFFE << bmask_a_i)) << bmask_b_i;

  // ---- fast parallel-prefix (Sklansky) adders ----
  function automatic [36:0] sk36(input logic [35:0] a, input logic [35:0] b);
    logic [35:0] pt, pp, g, c;
    begin
      pt = a ^ b;
      pp = a ^ b;
      g  = a & b;
      for (int L = 1; L < 36; L = L * 2) begin
        for (int i = 35; i >= L; i = i - 1) begin
          g[i] = g[i] | (pp[i] & g[i-L]);
          if (i >= 2*L) pp[i] = pp[i] & pp[i-L];
        end
      end
      c = {g[34:0], 1'b0};
      sk36 = {g[35] ^ a[35] ^ b[35], pt ^ c};
    end
  endfunction

  function automatic [31:0] sk32(input logic [31:0] a, input logic [31:0] b);
    logic [31:0] pt, pp, g, c;
    begin
      pt = a ^ b;
      pp = a ^ b;
      g  = a & b;
      for (int L = 1; L < 32; L = L * 2) begin
        for (int i = 31; i >= L; i = i - 1) begin
          g[i] = g[i] | (pp[i] & g[i-L]);
          if (i >= 2*L) pp[i] = pp[i] & pp[i-L];
        end
      end
      c = {g[30:0], 1'b0};
      sk32 = pt ^ c;
    end
  endfunction

  logic [31:0] adder_round_value, adder_round_result;
  assign adder_round_value  = f_rounding ? {1'b0, bmask[31:1]} : '0;
  assign adder_round_result = sk32(adder_result, adder_round_value);

  //////////////////////////////////////////////////////////////////////////////////////
  // SHIFTER
  //////////////////////////////////////////////////////////////////////////////////////
  logic        shift_left, shift_use_round, shift_arithmetic;
  logic [31:0] shift_amt, shift_amt_left, shift_amt_int, shift_amt_norm, shift_op_a;
  logic [31:0] shift_result, shift_right_result, shift_left_result;
  logic [15:0] clpx_shift_ex;
  logic [ 5:0] div_shift;
  logic        div_valid;

  assign shift_amt = div_valid ? {26'd0, div_shift} : operand_b_i;

  always_comb begin
    case (vector_mode_i)
      2'b10: begin
        shift_amt_left[15:0]  = shift_amt[31:16];
        shift_amt_left[31:16] = shift_amt[15:0];
      end
      2'b11: begin
        shift_amt_left[7:0]   = shift_amt[31:24];
        shift_amt_left[15:8]  = shift_amt[23:16];
        shift_amt_left[23:16] = shift_amt[15:8];
        shift_amt_left[31:24] = shift_amt[7:0];
      end
      default: shift_amt_left[31:0] = shift_amt[31:0];
    endcase
  end

  assign shift_left     = f_sll | f_bins | f_fl1 | f_clb | f_divgrp | f_brev;
  assign shift_use_round= f_addfam;
  assign shift_arithmetic = f_sra | f_bext | (f_addfam & ~op[1]);

  assign shift_op_a    = shift_left ? operand_a_rev : (shift_use_round ? adder_round_result : operand_a_i);
  assign shift_amt_int = shift_use_round ? shift_amt_norm : (shift_left ? shift_amt_left : shift_amt);
  assign shift_amt_norm = is_clpx_i ? {clpx_shift_ex, clpx_shift_ex} : {4{3'b000, bmask_b_i}};
  assign clpx_shift_ex = $unsigned(clpx_shift_i);

  logic [63:0] shift_op_a_32;
  assign shift_op_a_32 = f_ror ? {shift_op_a, shift_op_a}
                               : $signed({{32{shift_arithmetic & shift_op_a[31]}}, shift_op_a});

  always_comb begin
    case (vector_mode_i)
      2'b10: begin
        shift_right_result[31:16] = $signed({shift_arithmetic & shift_op_a[31], shift_op_a[31:16]}) >>> shift_amt_int[19:16];
        shift_right_result[15:0]  = $signed({shift_arithmetic & shift_op_a[15], shift_op_a[15:0]}) >>> shift_amt_int[3:0];
      end
      2'b11: begin
        shift_right_result[31:24] = $signed({shift_arithmetic & shift_op_a[31], shift_op_a[31:24]}) >>> shift_amt_int[26:24];
        shift_right_result[23:16] = $signed({shift_arithmetic & shift_op_a[23], shift_op_a[23:16]}) >>> shift_amt_int[18:16];
        shift_right_result[15:8]  = $signed({shift_arithmetic & shift_op_a[15], shift_op_a[15:8]})  >>> shift_amt_int[10:8];
        shift_right_result[7:0]   = $signed({shift_arithmetic & shift_op_a[7],  shift_op_a[7:0]})   >>> shift_amt_int[2:0];
      end
      default: shift_right_result = shift_op_a_32 >> shift_amt_int[4:0];
    endcase
  end

  generate
    for (k = 0; k < 32; k++) begin : gen_slr
      assign shift_left_result[k] = shift_right_result[31-k];
    end
  endgenerate

  assign shift_result = shift_left ? shift_left_result : shift_right_result;

  //////////////////////////////////////////////////////////////////////////////////////
  // COMPARISONS
  //////////////////////////////////////////////////////////////////////////////////////
  logic [3:0] is_equal, is_greater, cmp_result;
  logic [3:0] is_equal_vec, is_greater_vec;
  logic [3:0] cmp_signed;

  always_comb begin
    cmp_signed = 4'b0;
    if (f_gts | f_ges | f_lts | f_les | f_slts | f_slets | f_min | f_max | f_abs | f_clip | f_clipu) begin
      case (vector_mode_i)
        2'b11:   cmp_signed[3:0] = 4'b1111;
        2'b10:   cmp_signed[3:0] = 4'b1010;
        default: cmp_signed[3:0] = 4'b1000;
      endcase
    end
  end

  genvar i;
  generate
    for (i = 0; i < 4; i++) begin : gen_is_vec
      assign is_equal_vec[i]   = (operand_a_i[8*i+7:8*i] == operand_b_i[8*i+7:8*i]);
      assign is_greater_vec[i] = $signed({operand_a_i[8*i+7] & cmp_signed[i], operand_a_i[8*i+7:8*i]}) >
                                 $signed({operand_b_i[8*i+7] & cmp_signed[i], operand_b_i[8*i+7:8*i]});
    end
  endgenerate

  always_comb begin
    is_equal[3:0]   = {4{is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_equal_vec[0]}};
    is_greater[3:0] = {4{(is_greater_vec[3] | (is_equal_vec[3] & (is_greater_vec[2] |
                     (is_equal_vec[2] & (is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0]))))))}};
    case (vector_mode_i)
      2'b10: begin
        is_equal[1:0]   = {2{is_equal_vec[0] & is_equal_vec[1]}};
        is_equal[3:2]   = {2{is_equal_vec[2] & is_equal_vec[3]}};
        is_greater[1:0] = {2{is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0])}};
        is_greater[3:2] = {2{is_greater_vec[3] | (is_equal_vec[3] & is_greater_vec[2])}};
      end
      2'b11: begin
        is_equal[3:0]   = is_equal_vec[3:0];
        is_greater[3:0] = is_greater_vec[3:0];
      end
      default: ;
    endcase
  end

  always_comb begin
    cmp_result = is_equal;
    if (f_gts | f_gtu)                cmp_result = is_greater;
    else if (f_ges | f_geu)           cmp_result = is_greater | is_equal;
    else if (f_lts | f_ltu | f_slts | f_sltu) cmp_result = ~(is_greater | is_equal);
    else if (f_les | f_leu | f_slets | f_sletu) cmp_result = ~is_greater;
    else if (f_ne)                    cmp_result = ~is_equal;
  end

  assign comparison_result_o = cmp_result[3];

  // min/max/abs/clip
  logic [31:0] result_minmax;
  logic [3:0] sel_minmax;
  logic do_min;
  logic [31:0] minmax_b;

  assign minmax_b  = f_abs ? adder_result : operand_b_i;
  assign do_min    = f_min | f_minu | f_clip | f_clipu;
  assign sel_minmax[3:0] = is_greater ^ {4{do_min}};

  assign result_minmax[31:24] = sel_minmax[3] ? operand_a_i[31:24] : minmax_b[31:24];
  assign result_minmax[23:16] = sel_minmax[2] ? operand_a_i[23:16] : minmax_b[23:16];
  assign result_minmax[15:8]  = sel_minmax[1] ? operand_a_i[15:8]  : minmax_b[15:8];
  assign result_minmax[7:0]   = sel_minmax[0] ? operand_a_i[7:0]   : minmax_b[7:0];

  logic is_equal_clip;
  logic [31:0] operand_b_eq;
  assign operand_b_eq = f_clipu ? 32'b0 : operand_b_neg;
  assign is_equal_clip = (operand_a_i == operand_b_eq);

  logic [31:0] clip_result;
  always_comb begin
    clip_result = result_minmax;
    if (f_clipu) begin
      if (operand_a_i[31] || is_equal_clip) clip_result = '0;
    end else if (f_clip) begin
      if (adder_result_expanded[36] || is_equal_clip) clip_result = operand_b_neg;
    end
  end

  //////////////////////////////////////////////////////////////////////////////////////
  // SHUFFLE / PACK
  //////////////////////////////////////////////////////////////////////////////////////
  logic [3:0][1:0] shuffle_byte_sel;
  logic [3:0]      shuffle_reg_sel;
  logic [1:0]      shuffle_reg1_sel, shuffle_reg0_sel;
  logic [3:0]      shuffle_through;
  logic [31:0] shuffle_r1, shuffle_r0, shuffle_r1_in, shuffle_r0_in, shuffle_result, pack_result;

  always_comb begin
    shuffle_reg_sel  = '0;
    shuffle_reg1_sel = 2'b01;
    shuffle_reg0_sel = 2'b10;
    shuffle_through  = '1;
    if (f_ext | f_exts) begin
      if (f_exts) shuffle_reg1_sel = 2'b11;
      if (vm8) begin
        shuffle_reg_sel[3:1] = 3'b111;
        shuffle_reg_sel[0]   = 1'b0;
      end else begin
        shuffle_reg_sel[3:2] = 2'b11;
        shuffle_reg_sel[1:0] = 2'b00;
      end
    end else if (f_pcklo) begin
      shuffle_reg1_sel = 2'b00;
      if (vm8) begin
        shuffle_through = 4'b0011;
        shuffle_reg_sel = 4'b0001;
      end else begin
        shuffle_reg_sel = 4'b0011;
      end
    end else if (f_pckhi) begin
      shuffle_reg1_sel = 2'b00;
      if (vm8) begin
        shuffle_through = 4'b1100;
        shuffle_reg_sel = 4'b0100;
      end else begin
        shuffle_reg_sel = 4'b0011;
      end
    end else if (f_shuf2) begin
      if (vm8) begin
        shuffle_reg_sel[3] = ~operand_b_i[26];
        shuffle_reg_sel[2] = ~operand_b_i[18];
        shuffle_reg_sel[1] = ~operand_b_i[10];
        shuffle_reg_sel[0] = ~operand_b_i[2];
      end else begin
        shuffle_reg_sel[3] = ~operand_b_i[17];
        shuffle_reg_sel[2] = ~operand_b_i[17];
        shuffle_reg_sel[1] = ~operand_b_i[1];
        shuffle_reg_sel[0] = ~operand_b_i[1];
      end
    end else if (f_ins) begin
      if (vm8) begin
        shuffle_reg0_sel = 2'b00;
        case (imm_vec_ext_i)
          2'b00: shuffle_reg_sel[3:0] = 4'b1110;
          2'b01: shuffle_reg_sel[3:0] = 4'b1101;
          2'b10: shuffle_reg_sel[3:0] = 4'b1011;
          2'b11: shuffle_reg_sel[3:0] = 4'b0111;
        endcase
      end else begin
        shuffle_reg0_sel   = 2'b01;
        shuffle_reg_sel[3] = ~imm_vec_ext_i[0];
        shuffle_reg_sel[2] = ~imm_vec_ext_i[0];
        shuffle_reg_sel[1] = imm_vec_ext_i[0];
        shuffle_reg_sel[0] = imm_vec_ext_i[0];
      end
    end
  end

  always_comb begin
    shuffle_byte_sel = '0;
    if (f_ext | f_exts) begin
      if (vm8) begin
        shuffle_byte_sel[3] = imm_vec_ext_i[1:0];
        shuffle_byte_sel[2] = imm_vec_ext_i[1:0];
        shuffle_byte_sel[1] = imm_vec_ext_i[1:0];
        shuffle_byte_sel[0] = imm_vec_ext_i[1:0];
      end else begin
        shuffle_byte_sel[3] = {imm_vec_ext_i[0], 1'b1};
        shuffle_byte_sel[2] = {imm_vec_ext_i[0], 1'b1};
        shuffle_byte_sel[1] = {imm_vec_ext_i[0], 1'b1};
        shuffle_byte_sel[0] = {imm_vec_ext_i[0], 1'b0};
      end
    end else if (f_pcklo) begin
      if (vm8) begin
        shuffle_byte_sel[3] = 2'b00; shuffle_byte_sel[2] = 2'b00;
        shuffle_byte_sel[1] = 2'b00; shuffle_byte_sel[0] = 2'b00;
      end else begin
        shuffle_byte_sel[3] = 2'b01; shuffle_byte_sel[2] = 2'b00;
        shuffle_byte_sel[1] = 2'b01; shuffle_byte_sel[0] = 2'b00;
      end
    end else if (f_pckhi) begin
      if (vm8) begin
        shuffle_byte_sel[3] = 2'b00; shuffle_byte_sel[2] = 2'b00;
        shuffle_byte_sel[1] = 2'b00; shuffle_byte_sel[0] = 2'b00;
      end else begin
        shuffle_byte_sel[3] = 2'b11; shuffle_byte_sel[2] = 2'b10;
        shuffle_byte_sel[1] = 2'b11; shuffle_byte_sel[0] = 2'b10;
      end
    end else if (f_shuf | f_shuf2) begin
      if (vm8) begin
        shuffle_byte_sel[3] = operand_b_i[25:24];
        shuffle_byte_sel[2] = operand_b_i[17:16];
        shuffle_byte_sel[1] = operand_b_i[9:8];
        shuffle_byte_sel[0] = operand_b_i[1:0];
      end else begin
        shuffle_byte_sel[3] = {operand_b_i[16], 1'b1};
        shuffle_byte_sel[2] = {operand_b_i[16], 1'b0};
        shuffle_byte_sel[1] = {operand_b_i[0], 1'b1};
        shuffle_byte_sel[0] = {operand_b_i[0], 1'b0};
      end
    end else if (f_ins) begin
      shuffle_byte_sel[3] = 2'b11;
      shuffle_byte_sel[2] = 2'b10;
      shuffle_byte_sel[1] = 2'b01;
      shuffle_byte_sel[0] = 2'b00;
    end
  end

  assign shuffle_r0_in = shuffle_reg0_sel[1] ? operand_a_i :
                         (shuffle_reg0_sel[0] ? {2{operand_a_i[15:0]}} : {4{operand_a_i[7:0]}});
  assign shuffle_r1_in = shuffle_reg1_sel[1] ? {{8{operand_a_i[31]}}, {8{operand_a_i[23]}},
                                                 {8{operand_a_i[15]}}, {8{operand_a_i[7]}}} :
                         (shuffle_reg1_sel[0] ? operand_c_i : operand_b_i);

  generate
    for (k = 0; k < 4; k++) begin : gen_shfl
      assign shuffle_r0[8*k+:8] = shuffle_byte_sel[k][1] ?
            (shuffle_byte_sel[k][0] ? shuffle_r0_in[31:24] : shuffle_r0_in[23:16]) :
            (shuffle_byte_sel[k][0] ? shuffle_r0_in[15: 8] : shuffle_r0_in[ 7: 0]);
      assign shuffle_r1[8*k+:8] = shuffle_byte_sel[k][1] ?
            (shuffle_byte_sel[k][0] ? shuffle_r1_in[31:24] : shuffle_r1_in[23:16]) :
            (shuffle_byte_sel[k][0] ? shuffle_r1_in[15: 8] : shuffle_r1_in[ 7: 0]);
      assign shuffle_result[8*k+:8] = shuffle_reg_sel[k] ? shuffle_r1[8*k+:8] : shuffle_r0[8*k+:8];
      assign pack_result[8*k+:8]    = shuffle_through[k] ? shuffle_result[8*k+:8] : operand_c_i[8*k+:8];
    end
  endgenerate

  //////////////////////////////////////////////////////////////////////////////////////
  // BIT COUNTING
  //////////////////////////////////////////////////////////////////////////////////////
  logic [31:0] ff_input;
  logic [5:0] cnt_result, clb_result, bitop_result;
  logic [4:0] ff1_result, fl1_result;
  logic ff_no_one;

  cv32e40p_popcnt popcnt_i (.in_i(operand_a_i), .result_o(cnt_result));

  always_comb begin
    case (op)
      7'h36:   ff_input = operand_a_i;                                  // FF1
      7'h30, 7'h32, 7'h37: ff_input = operand_a_rev;                    // DIVU/REMU/FL1
      7'h31, 7'h33, 7'h35: ff_input = operand_a_i[31] ? operand_a_neg_rev : operand_a_rev; // DIV/REM/CLB
      default: ff_input = operand_a_i;
    endcase
  end

  cv32e40p_ff_one ff_one_i (.in_i(ff_input), .first_one_o(ff1_result), .no_ones_o(ff_no_one));

  assign fl1_result = 5'd31 - ff1_result;
  assign clb_result = ff1_result - 5'd1;

  always_comb begin
    bitop_result = '0;
    case (op)
      7'h36: bitop_result = ff_no_one ? 6'd32 : {1'b0, ff1_result};
      7'h37: bitop_result = ff_no_one ? 6'd32 : {1'b0, fl1_result};
      7'h34: bitop_result = cnt_result;
      7'h35: bitop_result = ff_no_one ? (operand_a_i[31] ? 6'd31 : 6'd0) : clb_result;
      default: ;
    endcase
  end

  //////////////////////////////////////////////////////////////////////////////////////
  // BIT MANIPULATION
  //////////////////////////////////////////////////////////////////////////////////////
  logic extract_sign;
  logic [31:0] bmask_first, bmask_inv, bextins_and, bextins_result, bclr_result, bset_result;

  assign bmask_first = 32'hFFFFFFFE << bmask_a_i;
  assign bmask_inv   = ~bmask;
  assign bextins_and = f_bins ? operand_c_i : {32{extract_sign}};
  assign extract_sign = f_bext & shift_result[bmask_a_i];
  assign bextins_result = (bmask & shift_result) | (bextins_and & bmask_inv);
  assign bclr_result = operand_a_i & bmask_inv;
  assign bset_result = operand_a_i | bmask;

  //////////////////////////////////////////////////////////////////////////////////////
  // BIT REVERSE (BREV)
  //////////////////////////////////////////////////////////////////////////////////////
  logic [31:0] radix_2_rev, radix_4_rev, radix_8_rev, reverse_result;
  logic [1:0] radix_mux_sel;
  assign radix_mux_sel = bmask_a_i[1:0];

  generate
    for (k = 0; k < 32; k++) begin : gen_r2
      assign radix_2_rev[k] = shift_result[31-k];
    end
    for (k = 0; k < 16; k++) begin : gen_r4
      assign radix_4_rev[2*k+1:2*k] = shift_result[31-2*k:31-2*k-1];
    end
    for (k = 0; k < 10; k++) begin : gen_r8
      assign radix_8_rev[3*k+2:3*k] = shift_result[31-3*k:31-3*k-2];
    end
    assign radix_8_rev[31:30] = 2'b0;
  endgenerate

  always_comb begin
    unique case (radix_mux_sel)
      2'b00:  reverse_result = radix_2_rev;
      2'b01:  reverse_result = radix_4_rev;
      2'b10:  reverse_result = radix_8_rev;
      default: reverse_result = radix_2_rev;
    endcase
  end

  //////////////////////////////////////////////////////////////////////////////////////
  // DIVIDER
  //////////////////////////////////////////////////////////////////////////////////////
  logic [31:0] result_div;
  logic        div_ready, div_signed, div_op_a_signed;
  logic [5:0]  div_shift_int;

  assign div_signed     = op[0];
  assign div_op_a_signed = operand_a_i[31] & div_signed;
  assign div_shift_int  = ff_no_one ? 6'd31 : clb_result;
  assign div_shift      = div_shift_int + (div_op_a_signed ? 6'd0 : 6'd1);
  assign div_valid      = enable_i & f_divgrp;

  cv32e40p_alu_div alu_div_i (
      .Clk_CI(clk), .Rst_RBI(rst_n),
      .OpA_DI(operand_b_i), .OpB_DI(shift_left_result), .OpBShift_DI(div_shift),
      .OpBIsZero_SI(cnt_result == 0), .OpBSign_SI(div_op_a_signed), .OpCode_SI(op[1:0]),
      .Res_DO(result_div), .InVld_SI(div_valid), .OutRdy_SI(ex_ready_i), .OutVld_SO(div_ready)
  );

  //////////////////////////////////////////////////////////////////////////////////////
  // RESULT MUX
  //////////////////////////////////////////////////////////////////////////////////////
  always_comb begin
    result_o = '0;
    case (op)
      7'h15:   result_o = operand_a_i & operand_b_i;
      7'h2e:   result_o = operand_a_i | operand_b_i;
      7'h2f:   result_o = operand_a_i ^ operand_b_i;
      7'h18, 7'h19, 7'h1a, 7'h1b, 7'h1c, 7'h1d, 7'h1e, 7'h1f, 7'h27, 7'h25, 7'h24, 7'h26:
               result_o = shift_result;
      7'h2a, 7'h28, 7'h29: result_o = bextins_result;
      7'h2b:   result_o = bclr_result;
      7'h2c:   result_o = bset_result;
      7'h49:   result_o = reverse_result;
      7'h3a, 7'h3b, 7'h38, 7'h39, 7'h3e, 7'h3f, 7'h2d: result_o = pack_result;
      7'h10, 7'h11, 7'h12, 7'h13: result_o = result_minmax;
      7'h14:   result_o = is_clpx_i ? {adder_result[31:16], operand_a_i[15:0]} : result_minmax;
      7'h16, 7'h17: result_o = clip_result;
      7'h00, 7'h01, 7'h04, 7'h05, 7'h08, 7'h09, 7'h0a, 7'h0b, 7'h0c, 7'h0d:
               result_o = {{8{cmp_result[3]}}, {8{cmp_result[2]}}, {8{cmp_result[1]}}, {8{cmp_result[0]}}};
      7'h02, 7'h03, 7'h06, 7'h07: result_o = {31'b0, cmp_result[3]};
      7'h36, 7'h37, 7'h34, 7'h35: result_o = {26'h0, bitop_result[5:0]};
      7'h30, 7'h31, 7'h32, 7'h33: result_o = result_div;
      default: ;
    endcase
  end

  assign ready_o = div_ready;

endmodule

// ---------------------------------------------------------------------------
// popcount (tree)
// ---------------------------------------------------------------------------
module cv32e40p_popcnt (
    input  logic [31:0] in_i,
    output logic [ 5:0] result_o
);
  logic [15:0][1:0] cnt_l1;
  logic [7:0][2:0]  cnt_l2;
  logic [3:0][3:0]  cnt_l3;
  logic [1:0][4:0]  cnt_l4;
  genvar l, m, n, p;
  generate
    for (l = 0; l < 16; l++) begin : gen_cnt_l1
      assign cnt_l1[l] = {1'b0, in_i[2*l]} + {1'b0, in_i[2*l+1]};
    end
    for (m = 0; m < 8; m++) begin : gen_cnt_l2
      assign cnt_l2[m] = {1'b0, cnt_l1[2*m]} + {1'b0, cnt_l1[2*m+1]};
    end
    for (n = 0; n < 4; n++) begin : gen_cnt_l3
      assign cnt_l3[n] = {1'b0, cnt_l2[2*n]} + {1'b0, cnt_l2[2*n+1]};
    end
    for (p = 0; p < 2; p++) begin : gen_cnt_l4
      assign cnt_l4[p] = {1'b0, cnt_l3[2*p]} + {1'b0, cnt_l3[2*p+1]};
    end
  endgenerate
  assign result_o = {1'b0, cnt_l4[0]} + {1'b0, cnt_l4[1]};
endmodule

// ---------------------------------------------------------------------------
// find first one
// ---------------------------------------------------------------------------
module cv32e40p_ff_one (
    input  logic [31:0] in_i,
    output logic [ 4:0] first_one_o,
    output logic        no_ones_o
);
  // isolate lowest set bit, then one-hot -> binary encode (OR tree)
  logic [31:0][4:0] idx;
  logic [31:0][4:0] nodes;
  logic [31:0] sel;
  genvar j;
  generate
    for (j = 0; j < 32; j++) begin : g_idx
      assign idx[j] = 5'(j);
    end
    for (j = 0; j < 16; j++) begin : g_l0
      assign sel[j] = in_i[2*j] | in_i[2*j+1];
      assign nodes[j] = in_i[2*j] ? idx[2*j] : idx[2*j+1];
    end
    for (j = 0; j < 8; j++) begin : g_l1
      assign sel[16+j] = sel[2*j] | sel[2*j+1];
      assign nodes[16+j] = sel[2*j] ? nodes[2*j] : nodes[2*j+1];
    end
    for (j = 0; j < 4; j++) begin : g_l2
      assign sel[24+j] = sel[16+2*j] | sel[16+2*j+1];
      assign nodes[24+j] = sel[16+2*j] ? nodes[16+2*j] : nodes[16+2*j+1];
    end
    for (j = 0; j < 2; j++) begin : g_l3
      assign sel[28+j] = sel[24+2*j] | sel[24+2*j+1];
      assign nodes[28+j] = sel[24+2*j] ? nodes[24+2*j] : nodes[24+2*j+1];
    end
    assign sel[30] = sel[28] | sel[29];
    assign nodes[30] = sel[28] ? nodes[28] : nodes[29];
  endgenerate
  assign first_one_o = nodes[30];
  assign no_ones_o = ~sel[30];
endmodule

// ---------------------------------------------------------------------------
// serial divider (unchanged)
// ---------------------------------------------------------------------------
module cv32e40p_alu_div #(
    parameter C_WIDTH     = 32,
    parameter C_LOG_WIDTH = 6
) (
    input  logic                   Clk_CI,
    input  logic                   Rst_RBI,
    input  logic [    C_WIDTH-1:0] OpA_DI,
    input  logic [    C_WIDTH-1:0] OpB_DI,
    input  logic [C_LOG_WIDTH-1:0] OpBShift_DI,
    input  logic                   OpBIsZero_SI,
    input  logic                   OpBSign_SI,
    input  logic [            1:0] OpCode_SI,
    input  logic                   InVld_SI,
    input  logic                   OutRdy_SI,
    output logic                   OutVld_SO,
    output logic [    C_WIDTH-1:0] Res_DO
);
  logic [C_WIDTH-1:0] ResReg_DP, ResReg_DN, ResReg_DP_rev, AReg_DP, AReg_DN, BReg_DP, BReg_DN;
  logic RemSel_SN, RemSel_SP, CompInv_SN, CompInv_SP, ResInv_SN, ResInv_SP;
  logic [C_WIDTH-1:0] AddMux_D, AddOut_D, AddTmp_D, BMux_D, OutMux_D;
  logic [C_LOG_WIDTH-1:0] Cnt_DP, Cnt_DN;
  logic CntZero_S, ARegEn_S, BRegEn_S, ResRegEn_S, ABComp_S, PmSel_S, LoadEn_S;

  enum logic [1:0] {IDLE, DIVIDE, FINISH} State_SN, State_SP;

  assign PmSel_S  = LoadEn_S & ~(OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI));
  assign AddMux_D = (LoadEn_S) ? OpA_DI : BReg_DP;
  assign BMux_D   = (LoadEn_S) ? OpB_DI : {CompInv_SP, (BReg_DP[$high(BReg_DP):1])};
  genvar index;
  generate
    for (index = 0; index < C_WIDTH; index++) begin : gen_bit_swapping
      assign ResReg_DP_rev[index] = ResReg_DP[C_WIDTH-1-index];
    end
  endgenerate
  assign OutMux_D = (RemSel_SP) ? AReg_DP : ResReg_DP_rev;
  assign Res_DO = (ResInv_SP) ? -$signed(OutMux_D) : OutMux_D;
  assign ABComp_S = ((AReg_DP == BReg_DP) | ((AReg_DP > BReg_DP) ^ CompInv_SP)) & ((|AReg_DP) | OpBIsZero_SI);
  assign AddTmp_D = (LoadEn_S) ? 0 : AReg_DP;
  assign AddOut_D = (PmSel_S) ? AddTmp_D + AddMux_D : AddTmp_D - $signed(AddMux_D);
  assign Cnt_DN = (LoadEn_S) ? OpBShift_DI : (~CntZero_S) ? Cnt_DP - 1 : Cnt_DP;
  assign CntZero_S = ~(|Cnt_DP);

  always_comb begin : p_fsm
    State_SN = State_SP;
    OutVld_SO = 1'b0;
    LoadEn_S = 1'b0;
    ARegEn_S = 1'b0;
    BRegEn_S = 1'b0;
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
        if (CntZero_S) State_SN = FINISH;
      end
      FINISH: begin
        OutVld_SO = 1'b1;
        if (OutRdy_SI) State_SN = IDLE;
      end
      default: ;
    endcase
  end

  assign RemSel_SN = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN = (LoadEn_S) ? ((~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI)) : ResInv_SP;
  assign AReg_DN = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN = (BRegEn_S) ? BMux_D : BReg_DP;
  assign ResReg_DN = (LoadEn_S) ? '0 : (ResRegEn_S) ? {ABComp_S, ResReg_DP[$high(ResReg_DP):1]} : ResReg_DP;

  always_ff @(posedge Clk_CI or negedge Rst_RBI) begin : p_regs
    if (~Rst_RBI) begin
      State_SP <= IDLE; AReg_DP <= '0; BReg_DP <= '0; ResReg_DP <= '0; Cnt_DP <= '0;
      RemSel_SP <= 1'b0; CompInv_SP <= 1'b0; ResInv_SP <= 1'b0;
    end else begin
      State_SP <= State_SN; AReg_DP <= AReg_DN; BReg_DP <= BReg_DN; ResReg_DP <= ResReg_DN;
      Cnt_DP <= Cnt_DN; RemSel_SP <= RemSel_SN; CompInv_SP <= CompInv_SN; ResInv_SP <= ResInv_SN;
    end
  end
endmodule

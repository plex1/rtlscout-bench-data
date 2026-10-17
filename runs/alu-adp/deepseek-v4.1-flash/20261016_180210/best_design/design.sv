// -----------------------------------------------------------------------------
// cv32e40p_alu - optimized re-implementation of the CV32E40P (RI5CY) ALU.
// Behaviourally/cycle identical to starting_point.sv.
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

  localparam logic [6:0] ALU_ADD   = 7'b0011000;
  localparam logic [6:0] ALU_SUB   = 7'b0011001;
  localparam logic [6:0] ALU_ADDU  = 7'b0011010;
  localparam logic [6:0] ALU_SUBU  = 7'b0011011;
  localparam logic [6:0] ALU_ADDR  = 7'b0011100;
  localparam logic [6:0] ALU_SUBR  = 7'b0011101;
  localparam logic [6:0] ALU_ADDUR = 7'b0011110;
  localparam logic [6:0] ALU_SUBUR = 7'b0011111;
  localparam logic [6:0] ALU_XOR   = 7'b0101111;
  localparam logic [6:0] ALU_OR    = 7'b0101110;
  localparam logic [6:0] ALU_AND   = 7'b0010101;
  localparam logic [6:0] ALU_SRA   = 7'b0100100;
  localparam logic [6:0] ALU_SRL   = 7'b0100101;
  localparam logic [6:0] ALU_ROR   = 7'b0100110;
  localparam logic [6:0] ALU_SLL   = 7'b0100111;
  localparam logic [6:0] ALU_BEXT  = 7'b0101000;
  localparam logic [6:0] ALU_BEXTU = 7'b0101001;
  localparam logic [6:0] ALU_BINS  = 7'b0101010;
  localparam logic [6:0] ALU_BCLR  = 7'b0101011;
  localparam logic [6:0] ALU_BSET  = 7'b0101100;
  localparam logic [6:0] ALU_BREV  = 7'b1001001;
  localparam logic [6:0] ALU_FF1   = 7'b0110110;
  localparam logic [6:0] ALU_FL1   = 7'b0110111;
  localparam logic [6:0] ALU_CNT   = 7'b0110100;
  localparam logic [6:0] ALU_CLB   = 7'b0110101;
  localparam logic [6:0] ALU_EXTS  = 7'b0111110;
  localparam logic [6:0] ALU_EXT   = 7'b0111111;
  localparam logic [6:0] ALU_LTS   = 7'b0000000;
  localparam logic [6:0] ALU_LTU   = 7'b0000001;
  localparam logic [6:0] ALU_LES   = 7'b0000100;
  localparam logic [6:0] ALU_LEU   = 7'b0000101;
  localparam logic [6:0] ALU_GTS   = 7'b0001000;
  localparam logic [6:0] ALU_GTU   = 7'b0001001;
  localparam logic [6:0] ALU_GES   = 7'b0001010;
  localparam logic [6:0] ALU_GEU   = 7'b0001011;
  localparam logic [6:0] ALU_EQ    = 7'b0001100;
  localparam logic [6:0] ALU_NE    = 7'b0001101;
  localparam logic [6:0] ALU_SLTS  = 7'b0000010;
  localparam logic [6:0] ALU_SLTU  = 7'b0000011;
  localparam logic [6:0] ALU_SLETS = 7'b0000110;
  localparam logic [6:0] ALU_SLETU = 7'b0000111;
  localparam logic [6:0] ALU_ABS   = 7'b0010100;
  localparam logic [6:0] ALU_CLIP  = 7'b0010110;
  localparam logic [6:0] ALU_CLIPU = 7'b0010111;
  localparam logic [6:0] ALU_INS   = 7'b0101101;
  localparam logic [6:0] ALU_MIN   = 7'b0010000;
  localparam logic [6:0] ALU_MINU  = 7'b0010001;
  localparam logic [6:0] ALU_MAX   = 7'b0010010;
  localparam logic [6:0] ALU_MAXU  = 7'b0010011;
  localparam logic [6:0] ALU_DIVU  = 7'b0110000;
  localparam logic [6:0] ALU_DIV   = 7'b0110001;
  localparam logic [6:0] ALU_REMU  = 7'b0110010;
  localparam logic [6:0] ALU_REM   = 7'b0110011;
  localparam logic [6:0] ALU_SHUF  = 7'b0111010;
  localparam logic [6:0] ALU_SHUF2 = 7'b0111011;
  localparam logic [6:0] ALU_PCKLO = 7'b0111000;
  localparam logic [6:0] ALU_PCKHI = 7'b0111001;
  localparam logic [1:0] VEC_MODE32 = 2'b00;
  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8  = 2'b11;

  logic [31:0] operand_a_rev;
  logic [31:0] operand_a_neg;
  logic [31:0] operand_a_neg_rev;

  assign operand_a_neg = ~operand_a_i;

  genvar k;
  generate
    for (k = 0; k < 32; k++) begin : gen_operand_a_rev
      assign operand_a_rev[k] = operand_a_i[31-k];
      assign operand_a_neg_rev[k] = operand_a_neg[31-k];
    end
  endgenerate

  logic [31:0] operand_b_neg;
  assign operand_b_neg = ~operand_b_i;

  logic [ 5:0] div_shift;
  logic        div_valid;
  logic [31:0] bmask;

  //////////////////////////////////////////////////////////////////////////
  // Adder
  //////////////////////////////////////////////////////////////////////////
  logic        adder_op_b_negate;
  logic [31:0] adder_op_a, adder_op_b;
  logic [35:0] adder_in_a, adder_in_b;
  logic [31:0] adder_result;
  logic [36:0] adder_result_expanded;

  assign adder_op_b_negate = (operator_i == ALU_SUB)  || (operator_i == ALU_SUBR) ||
                             (operator_i == ALU_SUBU) || (operator_i == ALU_SUBUR) ||
                             is_subrot_i;

  assign adder_op_a = (operator_i == ALU_ABS) ? operand_a_neg :
                      (is_subrot_i ? {operand_b_i[15:0], operand_a_i[31:16]} : operand_a_i);

  assign adder_op_b = adder_op_b_negate ?
                      (is_subrot_i ? ~{operand_a_i[15:0], operand_b_i[31:16]} : operand_b_neg) :
                      operand_b_i;

  logic [3:0] lane_cin;

  always_comb begin
    adder_in_a = 36'h0;
    adder_in_b = 36'h0;

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

    if (adder_op_b_negate || (operator_i == ALU_ABS) || (operator_i == ALU_CLIP)) begin
      adder_in_b[0] = 1'b1;
      case (vector_mode_i)
        VEC_MODE16: adder_in_b[18] = 1'b1;
        VEC_MODE8: begin
          adder_in_b[9]  = 1'b1;
          adder_in_b[18] = 1'b1;
          adder_in_b[27] = 1'b1;
        end
        default: ;
      endcase
    end else begin
      case (vector_mode_i)
        VEC_MODE16: adder_in_a[18] = 1'b0;
        VEC_MODE8: begin
          adder_in_a[9]  = 1'b0;
          adder_in_a[18] = 1'b0;
          adder_in_a[27] = 1'b0;
        end
        default: ;
      endcase
    end
  end

  assign adder_result_expanded = {adder_in_a[35], adder_in_a} + {adder_in_b[35], adder_in_b};
  assign adder_result = {
    adder_result_expanded[35:28],
    adder_result_expanded[26:19],
    adder_result_expanded[17:10],
    adder_result_expanded[8:1]
  };

  // normalization stage
  logic [31:0] adder_round_value;
  logic [31:0] adder_round_result;

  assign adder_round_value  = ((operator_i == ALU_ADDR) || (operator_i == ALU_SUBR) ||
                               (operator_i == ALU_ADDUR) || (operator_i == ALU_SUBUR)) ?
                                {1'b0, bmask[31:1]} : 32'h0;
  assign adder_round_result = adder_result + adder_round_value;

  //////////////////////////////////////////////////////////////////////////
  // Shifter
  //////////////////////////////////////////////////////////////////////////
  logic        shift_left;
  logic        shift_use_round;
  logic        shift_arithmetic;

  logic [31:0] shift_amt_left;
  logic [31:0] shift_amt;
  logic [31:0] shift_amt_int;
  logic [31:0] shift_amt_norm;
  logic [31:0] shift_op_a;
  logic [31:0] shift_result;
  logic [31:0] shift_right_result;
  logic [31:0] shift_left_result;
  logic [15:0] clpx_shift_ex;

  assign shift_amt = div_valid ? div_shift : operand_b_i;

  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_amt_left[15:0]  = shift_amt[31:16];
        shift_amt_left[31:16] = shift_amt[15:0];
      end
      VEC_MODE8: begin
        shift_amt_left[7:0]   = shift_amt[31:24];
        shift_amt_left[15:8]  = shift_amt[23:16];
        shift_amt_left[23:16] = shift_amt[15:8];
        shift_amt_left[31:24] = shift_amt[7:0];
      end
      default: shift_amt_left[31:0] = shift_amt[31:0];
    endcase
  end

  assign shift_left = (operator_i == ALU_SLL) || (operator_i == ALU_BINS) ||
                      (operator_i == ALU_FL1) || (operator_i == ALU_CLB)  ||
                      (operator_i == ALU_DIV) || (operator_i == ALU_DIVU) ||
                      (operator_i == ALU_REM) || (operator_i == ALU_REMU) ||
                      (operator_i == ALU_BREV);

  assign shift_use_round = (operator_i == ALU_ADD)   || (operator_i == ALU_SUB)   ||
                           (operator_i == ALU_ADDR)  || (operator_i == ALU_SUBR)  ||
                           (operator_i == ALU_ADDU)  || (operator_i == ALU_SUBU)  ||
                           (operator_i == ALU_ADDUR) || (operator_i == ALU_SUBUR);

  assign shift_arithmetic = (operator_i == ALU_SRA)  || (operator_i == ALU_BEXT) ||
                            (operator_i == ALU_ADD)  || (operator_i == ALU_SUB)  ||
                            (operator_i == ALU_ADDR) || (operator_i == ALU_SUBR);

  assign shift_op_a    = shift_left ? operand_a_rev :
                          (shift_use_round ? adder_round_result : operand_a_i);
  assign shift_amt_int = shift_use_round ? shift_amt_norm :
                          (shift_left ? shift_amt_left : shift_amt);

  assign shift_amt_norm = is_clpx_i ? {clpx_shift_ex, clpx_shift_ex} : {4{3'b000, bmask_b_i}};
  assign clpx_shift_ex  = $unsigned(clpx_shift_i);

  logic [63:0] shift_op_a_32;

  assign shift_op_a_32 = (operator_i == ALU_ROR) ? {shift_op_a, shift_op_a} :
                         $signed({{32{shift_arithmetic & shift_op_a[31]}}, shift_op_a});

  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_right_result[31:16] = $signed(
            {shift_arithmetic & shift_op_a[31], shift_op_a[31:16]}
        ) >>> shift_amt_int[19:16];
        shift_right_result[15:0] = $signed(
            {shift_arithmetic & shift_op_a[15], shift_op_a[15:0]}
        ) >>> shift_amt_int[3:0];
      end
      VEC_MODE8: begin
        shift_right_result[31:24] = $signed(
            {shift_arithmetic & shift_op_a[31], shift_op_a[31:24]}
        ) >>> shift_amt_int[26:24];
        shift_right_result[23:16] = $signed(
            {shift_arithmetic & shift_op_a[23], shift_op_a[23:16]}
        ) >>> shift_amt_int[18:16];
        shift_right_result[15:8] = $signed(
            {shift_arithmetic & shift_op_a[15], shift_op_a[15:8]}
        ) >>> shift_amt_int[10:8];
        shift_right_result[7:0] = $signed(
            {shift_arithmetic & shift_op_a[7], shift_op_a[7:0]}
        ) >>> shift_amt_int[2:0];
      end
      default: shift_right_result = shift_op_a_32 >> shift_amt_int[4:0];
    endcase
  end

  genvar j;
  generate
    for (j = 0; j < 32; j++) begin : gen_shift_left_result
      assign shift_left_result[j] = shift_right_result[31-j];
    end
  endgenerate

  assign shift_result = shift_left ? shift_left_result : shift_right_result;

  //////////////////////////////////////////////////////////////////////////
  // Comparators
  //////////////////////////////////////////////////////////////////////////
  logic [ 3:0] is_equal;
  logic [ 3:0] is_greater;

  logic [ 3:0] cmp_signed;
  logic [ 3:0] is_equal_vec;
  logic [ 3:0] is_greater_vec;
  logic [31:0] operand_b_eq;
  logic        is_equal_clip;

  always_comb begin
    operand_b_eq = operand_b_neg;
    if (operator_i == ALU_CLIPU) operand_b_eq = 32'h0;
  end
  assign is_equal_clip = (operand_a_i == operand_b_eq);

  always_comb begin
    cmp_signed = 4'b0;
    unique case (operator_i)
      ALU_GTS, ALU_GES, ALU_LTS, ALU_LES, ALU_SLTS, ALU_SLETS,
      ALU_MIN, ALU_MAX, ALU_ABS, ALU_CLIP, ALU_CLIPU: begin
        case (vector_mode_i)
          VEC_MODE8:  cmp_signed[3:0] = 4'b1111;
          VEC_MODE16: cmp_signed[3:0] = 4'b1010;
          default:    cmp_signed[3:0] = 4'b1000;
        endcase
      end
      default: ;
    endcase
  end

  genvar i;
  generate
    for (i = 0; i < 4; i++) begin : gen_is_vec
      assign is_equal_vec[i] = (operand_a_i[8*i+7:8*i] == operand_b_i[8*i+7:8*i]);
      assign is_greater_vec[i] = $signed(
          {operand_a_i[8*i+7] & cmp_signed[i], operand_a_i[8*i+7:8*i]}
      ) > $signed(
          {operand_b_i[8*i+7] & cmp_signed[i], operand_b_i[8*i+7:8*i]}
      );
    end
  endgenerate

  always_comb begin
    is_equal[3:0] = {4{is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_equal_vec[0]}};
    is_greater[3:0] = {4{is_greater_vec[3] | (is_equal_vec[3] & (is_greater_vec[2]
                                            | (is_equal_vec[2] & (is_greater_vec[1]
                                             | (is_equal_vec[1] & (is_greater_vec[0]))))))}};

    case (vector_mode_i)
      VEC_MODE16: begin
        is_equal[1:0]   = {2{is_equal_vec[0] & is_equal_vec[1]}};
        is_equal[3:2]   = {2{is_equal_vec[2] & is_equal_vec[3]}};
        is_greater[1:0] = {2{is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0])}};
        is_greater[3:2] = {2{is_greater_vec[3] | (is_equal_vec[3] & is_greater_vec[2])}};
      end
      VEC_MODE8: begin
        is_equal[3:0]   = is_equal_vec[3:0];
        is_greater[3:0] = is_greater_vec[3:0];
      end
      default: ;
    endcase
  end

  logic [3:0] cmp_result;

  always_comb begin
    cmp_result = is_equal;
    unique case (operator_i)
      ALU_EQ:                                 cmp_result = is_equal;
      ALU_NE:                                 cmp_result = ~is_equal;
      ALU_GTS, ALU_GTU:                       cmp_result = is_greater;
      ALU_GES, ALU_GEU:                       cmp_result = is_greater | is_equal;
      ALU_LTS, ALU_SLTS, ALU_LTU, ALU_SLTU:   cmp_result = ~(is_greater | is_equal);
      ALU_SLETS, ALU_SLETU, ALU_LES, ALU_LEU: cmp_result = ~is_greater;
      default:                                ;
    endcase
  end

  assign comparison_result_o = cmp_result[3];

  // min/max/abs handling
  logic [31:0] result_minmax;
  logic [ 3:0] sel_minmax;
  logic        do_min;
  logic [31:0] minmax_b;

  assign minmax_b = (operator_i == ALU_ABS) ? adder_result : operand_b_i;

  assign do_min   = (operator_i == ALU_MIN)  || (operator_i == ALU_MINU) ||
                    (operator_i == ALU_CLIP) || (operator_i == ALU_CLIPU);

  assign sel_minmax[3:0] = is_greater ^ {4{do_min}};

  assign result_minmax[31:24] = sel_minmax[3] ? operand_a_i[31:24] : minmax_b[31:24];
  assign result_minmax[23:16] = sel_minmax[2] ? operand_a_i[23:16] : minmax_b[23:16];
  assign result_minmax[15:8]  = sel_minmax[1] ? operand_a_i[15:8]  : minmax_b[15:8];
  assign result_minmax[7:0]   = sel_minmax[0] ? operand_a_i[7:0]   : minmax_b[7:0];

  //////////////////////////////////////////////////////////////////////////
  // Clip
  //////////////////////////////////////////////////////////////////////////
  logic [31:0] clip_result;

  always_comb begin
    clip_result = result_minmax;
    if (operator_i == ALU_CLIPU) begin
      if (operand_a_i[31] || is_equal_clip) clip_result = 32'h0;
    end else begin
      // CLIP
      if (adder_result_expanded[36] || is_equal_clip) clip_result = operand_b_neg;
    end
  end

  //////////////////////////////////////////////////////////////////////////
  // Shuffle
  //////////////////////////////////////////////////////////////////////////
  logic [3:0][1:0] shuffle_byte_sel;
  logic [3:0]      shuffle_reg_sel;
  logic [1:0]      shuffle_reg1_sel;
  logic [1:0]      shuffle_reg0_sel;
  logic [3:0]      shuffle_through;

  logic [31:0] shuffle_r1, shuffle_r0;
  logic [31:0] shuffle_r1_in, shuffle_r0_in;
  logic [31:0] shuffle_result;
  logic [31:0] pack_result;

  always_comb begin
    shuffle_reg_sel  = 4'h0;
    shuffle_reg1_sel = 2'b01;
    shuffle_reg0_sel = 2'b10;
    shuffle_through  = 4'hF;

    unique case (operator_i)
      ALU_EXT, ALU_EXTS: begin
        if (operator_i == ALU_EXTS) shuffle_reg1_sel = 2'b11;

        if (vector_mode_i == VEC_MODE8) begin
          shuffle_reg_sel[3:1] = 3'b111;
          shuffle_reg_sel[0]   = 1'b0;
        end else begin
          shuffle_reg_sel[3:2] = 2'b11;
          shuffle_reg_sel[1:0] = 2'b00;
        end
      end

      ALU_PCKLO: begin
        shuffle_reg1_sel = 2'b00;
        if (vector_mode_i == VEC_MODE8) begin
          shuffle_through = 4'b0011;
          shuffle_reg_sel = 4'b0001;
        end else begin
          shuffle_reg_sel = 4'b0011;
        end
      end

      ALU_PCKHI: begin
        shuffle_reg1_sel = 2'b00;
        if (vector_mode_i == VEC_MODE8) begin
          shuffle_through = 4'b1100;
          shuffle_reg_sel = 4'b0100;
        end else begin
          shuffle_reg_sel = 4'b0011;
        end
      end

      ALU_SHUF2: begin
        unique case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_reg_sel[3] = ~operand_b_i[26];
            shuffle_reg_sel[2] = ~operand_b_i[18];
            shuffle_reg_sel[1] = ~operand_b_i[10];
            shuffle_reg_sel[0] = ~operand_b_i[2];
          end
          VEC_MODE16: begin
            shuffle_reg_sel[3] = ~operand_b_i[17];
            shuffle_reg_sel[2] = ~operand_b_i[17];
            shuffle_reg_sel[1] = ~operand_b_i[1];
            shuffle_reg_sel[0] = ~operand_b_i[1];
          end
          default: ;
        endcase
      end

      ALU_INS: begin
        unique case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_reg0_sel = 2'b00;
            unique case (imm_vec_ext_i)
              2'b00:   shuffle_reg_sel[3:0] = 4'b1110;
              2'b01:   shuffle_reg_sel[3:0] = 4'b1101;
              2'b10:   shuffle_reg_sel[3:0] = 4'b1011;
              default: shuffle_reg_sel[3:0] = 4'b0111;
            endcase
          end
          VEC_MODE16: begin
            shuffle_reg0_sel   = 2'b01;
            shuffle_reg_sel[3] = ~imm_vec_ext_i[0];
            shuffle_reg_sel[2] = ~imm_vec_ext_i[0];
            shuffle_reg_sel[1] = imm_vec_ext_i[0];
            shuffle_reg_sel[0] = imm_vec_ext_i[0];
          end
          default: ;
        endcase
      end

      default: ;
    endcase
  end

  always_comb begin
    shuffle_byte_sel = 8'h0;

    unique case (operator_i)
      ALU_EXTS, ALU_EXT: begin
        unique case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_byte_sel[3] = imm_vec_ext_i[1:0];
            shuffle_byte_sel[2] = imm_vec_ext_i[1:0];
            shuffle_byte_sel[1] = imm_vec_ext_i[1:0];
            shuffle_byte_sel[0] = imm_vec_ext_i[1:0];
          end
          VEC_MODE16: begin
            shuffle_byte_sel[3] = {imm_vec_ext_i[0], 1'b1};
            shuffle_byte_sel[2] = {imm_vec_ext_i[0], 1'b1};
            shuffle_byte_sel[1] = {imm_vec_ext_i[0], 1'b1};
            shuffle_byte_sel[0] = {imm_vec_ext_i[0], 1'b0};
          end
          default: ;
        endcase
      end

      ALU_PCKLO: begin
        unique case (vector_mode_i)
          VEC_MODE8:  shuffle_byte_sel = 8'h00;
          VEC_MODE16: begin
            shuffle_byte_sel[3] = 2'b01;
            shuffle_byte_sel[2] = 2'b00;
            shuffle_byte_sel[1] = 2'b01;
            shuffle_byte_sel[0] = 2'b00;
          end
          default: ;
        endcase
      end

      ALU_PCKHI: begin
        unique case (vector_mode_i)
          VEC_MODE8:  shuffle_byte_sel = 8'h00;
          VEC_MODE16: begin
            shuffle_byte_sel[3] = 2'b11;
            shuffle_byte_sel[2] = 2'b10;
            shuffle_byte_sel[1] = 2'b11;
            shuffle_byte_sel[0] = 2'b10;
          end
          default: ;
        endcase
      end

      ALU_SHUF2, ALU_SHUF: begin
        unique case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_byte_sel[3] = operand_b_i[25:24];
            shuffle_byte_sel[2] = operand_b_i[17:16];
            shuffle_byte_sel[1] = operand_b_i[9:8];
            shuffle_byte_sel[0] = operand_b_i[1:0];
          end
          VEC_MODE16: begin
            shuffle_byte_sel[3] = {operand_b_i[16], 1'b1};
            shuffle_byte_sel[2] = {operand_b_i[16], 1'b0};
            shuffle_byte_sel[1] = {operand_b_i[0], 1'b1};
            shuffle_byte_sel[0] = {operand_b_i[0], 1'b0};
          end
          default: ;
        endcase
      end

      ALU_INS: begin
        shuffle_byte_sel[3] = 2'b11;
        shuffle_byte_sel[2] = 2'b10;
        shuffle_byte_sel[1] = 2'b01;
        shuffle_byte_sel[0] = 2'b00;
      end

      default: ;
    endcase
  end

  assign shuffle_r0_in = shuffle_reg0_sel[1] ? operand_a_i :
                          (shuffle_reg0_sel[0] ? {2{operand_a_i[15:0]}} : {4{operand_a_i[7:0]}});

  assign shuffle_r1_in = shuffle_reg1_sel[1] ? {
    {8{operand_a_i[31]}}, {8{operand_a_i[23]}}, {8{operand_a_i[15]}}, {8{operand_a_i[7]}}
  } : (shuffle_reg1_sel[0] ? operand_c_i : operand_b_i);

  assign shuffle_r0[31:24] = shuffle_byte_sel[3][1] ?
                              (shuffle_byte_sel[3][0] ? shuffle_r0_in[31:24] : shuffle_r0_in[23:16]) :
                              (shuffle_byte_sel[3][0] ? shuffle_r0_in[15: 8] : shuffle_r0_in[ 7: 0]);
  assign shuffle_r0[23:16] = shuffle_byte_sel[2][1] ?
                              (shuffle_byte_sel[2][0] ? shuffle_r0_in[31:24] : shuffle_r0_in[23:16]) :
                              (shuffle_byte_sel[2][0] ? shuffle_r0_in[15: 8] : shuffle_r0_in[ 7: 0]);
  assign shuffle_r0[15: 8] = shuffle_byte_sel[1][1] ?
                              (shuffle_byte_sel[1][0] ? shuffle_r0_in[31:24] : shuffle_r0_in[23:16]) :
                              (shuffle_byte_sel[1][0] ? shuffle_r0_in[15: 8] : shuffle_r0_in[ 7: 0]);
  assign shuffle_r0[ 7: 0] = shuffle_byte_sel[0][1] ?
                              (shuffle_byte_sel[0][0] ? shuffle_r0_in[31:24] : shuffle_r0_in[23:16]) :
                              (shuffle_byte_sel[0][0] ? shuffle_r0_in[15: 8] : shuffle_r0_in[ 7: 0]);

  assign shuffle_r1[31:24] = shuffle_byte_sel[3][1] ?
                              (shuffle_byte_sel[3][0] ? shuffle_r1_in[31:24] : shuffle_r1_in[23:16]) :
                              (shuffle_byte_sel[3][0] ? shuffle_r1_in[15: 8] : shuffle_r1_in[ 7: 0]);
  assign shuffle_r1[23:16] = shuffle_byte_sel[2][1] ?
                              (shuffle_byte_sel[2][0] ? shuffle_r1_in[31:24] : shuffle_r1_in[23:16]) :
                              (shuffle_byte_sel[2][0] ? shuffle_r1_in[15: 8] : shuffle_r1_in[ 7: 0]);
  assign shuffle_r1[15: 8] = shuffle_byte_sel[1][1] ?
                              (shuffle_byte_sel[1][0] ? shuffle_r1_in[31:24] : shuffle_r1_in[23:16]) :
                              (shuffle_byte_sel[1][0] ? shuffle_r1_in[15: 8] : shuffle_r1_in[ 7: 0]);
  assign shuffle_r1[ 7: 0] = shuffle_byte_sel[0][1] ?
                              (shuffle_byte_sel[0][0] ? shuffle_r1_in[31:24] : shuffle_r1_in[23:16]) :
                              (shuffle_byte_sel[0][0] ? shuffle_r1_in[15: 8] : shuffle_r1_in[ 7: 0]);

  assign shuffle_result[31:24] = shuffle_reg_sel[3] ? shuffle_r1[31:24] : shuffle_r0[31:24];
  assign shuffle_result[23:16] = shuffle_reg_sel[2] ? shuffle_r1[23:16] : shuffle_r0[23:16];
  assign shuffle_result[15:8]  = shuffle_reg_sel[1] ? shuffle_r1[15:8]  : shuffle_r0[15:8];
  assign shuffle_result[7:0]   = shuffle_reg_sel[0] ? shuffle_r1[7:0]   : shuffle_r0[7:0];

  assign pack_result[31:24] = shuffle_through[3] ? shuffle_result[31:24] : operand_c_i[31:24];
  assign pack_result[23:16] = shuffle_through[2] ? shuffle_result[23:16] : operand_c_i[23:16];
  assign pack_result[15:8]  = shuffle_through[1] ? shuffle_result[15:8]  : operand_c_i[15:8];
  assign pack_result[7:0]   = shuffle_through[0] ? shuffle_result[7:0]   : operand_c_i[7:0];

  //////////////////////////////////////////////////////////////////////////
  // Bit counting
  //////////////////////////////////////////////////////////////////////////
  logic [31:0] ff_input;
  logic [ 5:0] cnt_result;
  logic [ 5:0] clb_result;
  logic [ 4:0] ff1_result;
  logic        ff_no_one;
  logic [ 4:0] fl1_result;
  logic [ 5:0] bitop_result;

  cv32e40p_popcnt popcnt_i (
      .in_i     (operand_a_i),
      .result_o (cnt_result)
  );

  always_comb begin
    ff_input = 32'h0;
    case (operator_i)
      ALU_FF1:                     ff_input = operand_a_i;
      ALU_DIVU, ALU_REMU, ALU_FL1: ff_input = operand_a_rev;
      ALU_DIV, ALU_REM, ALU_CLB:   ff_input = operand_a_i[31] ? operand_a_neg_rev : operand_a_rev;
      default: ;
    endcase
  end

  cv32e40p_ff_one ff_one_i (
      .in_i        (ff_input),
      .first_one_o (ff1_result),
      .no_ones_o   (ff_no_one)
  );

  assign fl1_result = 5'd31 - ff1_result;
  // NOTE: context-width arithmetic on purpose: this is evaluated as a 6-bit
  // subtraction (matching the upstream RTL), so ff1_result == 0 yields 6'h3F.
  assign clb_result = ff1_result - 5'd1;

  always_comb begin
    bitop_result = 6'h0;
    case (operator_i)
      ALU_FF1: bitop_result = ff_no_one ? 6'd32 : {1'b0, ff1_result};
      ALU_FL1: bitop_result = ff_no_one ? 6'd32 : {1'b0, fl1_result};
      ALU_CNT: bitop_result = cnt_result;
      ALU_CLB: begin
        if (ff_no_one) bitop_result = operand_a_i[31] ? 6'd31 : 6'h0;
        else           bitop_result = clb_result;
      end
      default: ;
    endcase
  end

  //////////////////////////////////////////////////////////////////////////
  // Bit manipulation
  //////////////////////////////////////////////////////////////////////////
  logic        extract_is_signed;
  logic        extract_sign;
  logic [31:0] bmask_first, bmask_inv;
  logic [31:0] bextins_and;
  logic [31:0] bextins_result, bclr_result, bset_result;

  assign bmask_first = 32'hFFFFFFFE << bmask_a_i;
  assign bmask       = (~bmask_first) << bmask_b_i;
  assign bmask_inv   = ~bmask;

  assign extract_is_signed = (operator_i == ALU_BEXT);
  assign extract_sign      = extract_is_signed & shift_result[bmask_a_i];
  assign bextins_and       = (operator_i == ALU_BINS) ? operand_c_i : {32{extract_sign}};
  assign bextins_result    = (bmask & shift_result) | (bextins_and & bmask_inv);

  assign bclr_result = operand_a_i & bmask_inv;
  assign bset_result = operand_a_i | bmask;

  //////////////////////////////////////////////////////////////////////////
  // Bit reverse
  //////////////////////////////////////////////////////////////////////////
  logic [31:0] radix_2_rev;
  logic [31:0] radix_4_rev;
  logic [31:0] radix_8_rev;
  logic [31:0] reverse_result;
  logic [ 1:0] radix_mux_sel;

  assign radix_mux_sel = bmask_a_i[1:0];

  generate
    for (j = 0; j < 32; j++) begin : gen_radix_2_rev
      assign radix_2_rev[j] = shift_result[31-j];
    end
    for (j = 0; j < 16; j++) begin : gen_radix_4_rev
      assign radix_4_rev[2*j+1:2*j] = shift_result[31-j*2:31-j*2-1];
    end
    for (j = 0; j < 10; j++) begin : gen_radix_8_rev
      assign radix_8_rev[3*j+2:3*j] = shift_result[31-j*3:31-j*3-2];
    end
  endgenerate
  assign radix_8_rev[31:30] = 2'b0;

  always_comb begin
    unique case (radix_mux_sel)
      2'b00:   reverse_result = radix_2_rev;
      2'b01:   reverse_result = radix_4_rev;
      2'b10:   reverse_result = radix_8_rev;
      default: reverse_result = radix_2_rev;
    endcase
  end

  //////////////////////////////////////////////////////////////////////////
  // Divider
  //////////////////////////////////////////////////////////////////////////
  logic [31:0] result_div;
  logic        div_ready;
  logic        div_signed;
  logic        div_op_a_signed;
  logic [ 5:0] div_shift_int;

  assign div_signed      = operator_i[0];
  assign div_op_a_signed = operand_a_i[31] & div_signed;

  assign div_shift_int = ff_no_one ? 6'd31 : clb_result;
  assign div_shift     = div_shift_int + (div_op_a_signed ? 6'd0 : 6'd1);

  assign div_valid = enable_i & ((operator_i == ALU_DIV) || (operator_i == ALU_DIVU) ||
                                 (operator_i == ALU_REM) || (operator_i == ALU_REMU));

  // inputs A and B are swapped
  cv32e40p_alu_div alu_div_i (
      .Clk_CI       (clk),
      .Rst_RBI      (rst_n),
      .OpA_DI       (operand_b_i),
      .OpB_DI       (shift_left_result),
      .OpBShift_DI  (div_shift),
      .OpBIsZero_SI ((cnt_result == 6'd0)),
      .OpBSign_SI   (div_op_a_signed),
      .OpCode_SI    (operator_i[1:0]),
      .Res_DO       (result_div),
      .InVld_SI     (div_valid),
      .OutRdy_SI    (ex_ready_i),
      .OutVld_SO    (div_ready)
  );

  //////////////////////////////////////////////////////////////////////////
  // Result selection (flattened, one-hot masked OR)
  //////////////////////////////////////////////////////////////////////////
  logic [31:0] res_abs, res_cmp, res_sl, res_bitop;
  logic        sel_shift, sel_and, sel_or, sel_xor, sel_bextins, sel_bclr, sel_bset, sel_brev;
  logic        sel_pack, sel_mm, sel_abs, sel_clip, sel_cmp, sel_sl, sel_bitop, sel_div;

  assign res_abs   = is_clpx_i ? {adder_result[31:16], operand_a_i[15:0]} : result_minmax;
  assign res_cmp   = {{8{cmp_result[3]}}, {8{cmp_result[2]}}, {8{cmp_result[1]}}, {8{cmp_result[0]}}};
  assign res_sl    = {31'b0, comparison_result_o};
  assign res_bitop = {26'h0, bitop_result};

  assign sel_shift = (operator_i == ALU_ADD)   || (operator_i == ALU_SUB)   ||
                     (operator_i == ALU_ADDU)  || (operator_i == ALU_SUBU)  ||
                     (operator_i == ALU_ADDR)  || (operator_i == ALU_SUBR)  ||
                     (operator_i == ALU_ADDUR) || (operator_i == ALU_SUBUR) ||
                     (operator_i == ALU_SLL)   || (operator_i == ALU_SRL)   ||
                     (operator_i == ALU_SRA)   || (operator_i == ALU_ROR);

  assign sel_and = (operator_i == ALU_AND);
  assign sel_or  = (operator_i == ALU_OR);
  assign sel_xor = (operator_i == ALU_XOR);
  assign sel_bextins = (operator_i == ALU_BINS) || (operator_i == ALU_BEXT) || (operator_i == ALU_BEXTU);
  assign sel_bclr    = (operator_i == ALU_BCLR);
  assign sel_bset    = (operator_i == ALU_BSET);
  assign sel_brev    = (operator_i == ALU_BREV);
  assign sel_pack    = (operator_i == ALU_SHUF) || (operator_i == ALU_SHUF2) ||
                       (operator_i == ALU_PCKLO) || (operator_i == ALU_PCKHI) ||
                       (operator_i == ALU_EXT) || (operator_i == ALU_EXTS) || (operator_i == ALU_INS);
  assign sel_mm      = (operator_i == ALU_MIN) || (operator_i == ALU_MINU) ||
                       (operator_i == ALU_MAX) || (operator_i == ALU_MAXU);
  assign sel_abs     = (operator_i == ALU_ABS);
  assign sel_clip    = (operator_i == ALU_CLIP) || (operator_i == ALU_CLIPU);
  assign sel_sl      = (operator_i == ALU_SLTS) || (operator_i == ALU_SLTU) ||
                       (operator_i == ALU_SLETS) || (operator_i == ALU_SLETU);
  assign sel_cmp     = (operator_i == ALU_EQ)  || (operator_i == ALU_NE)  ||
                       (operator_i == ALU_GTS) || (operator_i == ALU_GTU) ||
                       (operator_i == ALU_GES) || (operator_i == ALU_GEU) ||
                       (operator_i == ALU_LTS) || (operator_i == ALU_LTU) ||
                       (operator_i == ALU_LES) || (operator_i == ALU_LEU);
  assign sel_bitop   = (operator_i == ALU_FF1) || (operator_i == ALU_FL1) ||
                       (operator_i == ALU_CNT) || (operator_i == ALU_CLB);
  assign sel_div     = (operator_i == ALU_DIV) || (operator_i == ALU_DIVU) ||
                       (operator_i == ALU_REM) || (operator_i == ALU_REMU);

  assign result_o = ({32{sel_shift}}   & shift_result)
                  | ({32{sel_and}}     & (operand_a_i & operand_b_i))
                  | ({32{sel_or}}      & (operand_a_i | operand_b_i))
                  | ({32{sel_xor}}     & (operand_a_i ^ operand_b_i))
                  | ({32{sel_bextins}} & bextins_result)
                  | ({32{sel_bclr}}    & bclr_result)
                  | ({32{sel_bset}}    & bset_result)
                  | ({32{sel_brev}}    & reverse_result)
                  | ({32{sel_pack}}    & pack_result)
                  | ({32{sel_mm}}      & result_minmax)
                  | ({32{sel_abs}}     & res_abs)
                  | ({32{sel_clip}}    & clip_result)
                  | ({32{sel_cmp}}     & res_cmp)
                  | ({32{sel_sl}}      & res_sl)
                  | ({32{sel_bitop}}   & res_bitop)
                  | ({32{sel_div}}     & result_div);

  assign ready_o = div_ready;

endmodule


// -----------------------------------------------------------------------------
// Population count
// -----------------------------------------------------------------------------
module cv32e40p_popcnt (
    input  logic [31:0] in_i,
    output logic [ 5:0] result_o
);

  logic [15:0][1:0] cnt_l1;
  logic [ 7:0][2:0] cnt_l2;
  logic [ 3:0][3:0] cnt_l3;
  logic [ 1:0][4:0] cnt_l4;

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


// -----------------------------------------------------------------------------
// Find first one
// -----------------------------------------------------------------------------
module cv32e40p_ff_one #(
    parameter LEN = 32
) (
    input  logic [LEN-1:0]           in_i,
    output logic [$clog2(LEN)-1:0]   first_one_o,
    output logic                     no_ones_o
);
  localparam NUM_LEVELS = $clog2(LEN);

  logic [          LEN-1:0][NUM_LEVELS-1:0] index_lut;
  logic [2**NUM_LEVELS-1:0]                 sel_nodes;
  logic [2**NUM_LEVELS-1:0][NUM_LEVELS-1:0] index_nodes;

  generate
    genvar j;
    for (j = 0; j < LEN; j++) begin : gen_index_lut
      assign index_lut[j] = $unsigned(j);
    end
  endgenerate

  generate
    genvar k, l, level;
    assign sel_nodes[2**NUM_LEVELS-1] = 1'b0;

    for (level = 0; level < NUM_LEVELS; level++) begin : gen_tree
      if (level < NUM_LEVELS - 1) begin : gen_non_root_level
        for (l = 0; l < 2 ** level; l++) begin : gen_node
          assign sel_nodes[2**level-1+l]   = sel_nodes[2**(level+1)-1+l*2] | sel_nodes[2**(level+1)-1+l*2+1];
          assign index_nodes[2**level-1+l] = (sel_nodes[2**(level+1)-1+l*2] == 1'b1) ?
                                           index_nodes[2**(level+1)-1+l*2] : index_nodes[2**(level+1)-1+l*2+1];
        end
      end
      if (level == NUM_LEVELS - 1) begin : gen_root_level
        for (k = 0; k < 2 ** level; k++) begin : gen_node
          if (k * 2 < LEN - 1) begin : gen_two
            assign sel_nodes[2**level-1+k] = in_i[k*2] | in_i[k*2+1];
            assign index_nodes[2**level-1+k] = (in_i[k*2] == 1'b1) ? index_lut[k*2] : index_lut[k*2+1];
          end
          if (k * 2 == LEN - 1) begin : gen_one
            assign sel_nodes[2**level-1+k]   = in_i[k*2];
            assign index_nodes[2**level-1+k] = index_lut[k*2];
          end
          if (k * 2 > LEN - 1) begin : gen_out_of_range
            assign sel_nodes[2**level-1+k]   = 1'b0;
            assign index_nodes[2**level-1+k] = '0;
          end
        end
      end
    end
  endgenerate

  assign first_one_o = index_nodes[0];
  assign no_ones_o   = ~sel_nodes[0];

endmodule


// -----------------------------------------------------------------------------
// Simple serial divider (cycle-exact copy of the upstream implementation)
// -----------------------------------------------------------------------------
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

  logic [C_WIDTH-1:0] ResReg_DP, ResReg_DN;
  logic [C_WIDTH-1:0] ResReg_DP_rev;
  logic [C_WIDTH-1:0] AReg_DP, AReg_DN;
  logic [C_WIDTH-1:0] BReg_DP, BReg_DN;

  logic RemSel_SN, RemSel_SP;
  logic CompInv_SN, CompInv_SP;
  logic ResInv_SN, ResInv_SP;

  logic [C_WIDTH-1:0] AddMux_D;
  logic [C_WIDTH-1:0] AddOut_D;
  logic [C_WIDTH-1:0] AddTmp_D;
  logic [C_WIDTH-1:0] BMux_D;
  logic [C_WIDTH-1:0] OutMux_D;

  logic [C_LOG_WIDTH-1:0] Cnt_DP, Cnt_DN;
  logic CntZero_S;

  logic ARegEn_S, BRegEn_S, ResRegEn_S, ABComp_S, PmSel_S, LoadEn_S;

  enum logic [1:0] {
    IDLE,
    DIVIDE,
    FINISH
  }
      State_SN, State_SP;

  assign PmSel_S  = LoadEn_S & ~(OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI));

  assign AddMux_D = (LoadEn_S) ? OpA_DI : BReg_DP;

  // attention: logical shift in case of negative operand B!
  assign BMux_D   = (LoadEn_S) ? OpB_DI : {CompInv_SP, (BReg_DP[$high(BReg_DP):1])};

  genvar index;
  generate
    for (index = 0; index < C_WIDTH; index++) begin : gen_bit_swapping
      assign ResReg_DP_rev[index] = ResReg_DP[C_WIDTH-1-index];
    end
  endgenerate

  assign OutMux_D = (RemSel_SP) ? AReg_DP : ResReg_DP_rev;

  assign Res_DO = (ResInv_SP) ? -$signed(OutMux_D) : OutMux_D;

  assign ABComp_S    = ((AReg_DP == BReg_DP) | ((AReg_DP > BReg_DP) ^ CompInv_SP)) & ((|AReg_DP) | OpBIsZero_SI);

  assign AddTmp_D = (LoadEn_S) ? 0 : AReg_DP;
  assign AddOut_D = (PmSel_S) ? AddTmp_D + AddMux_D : AddTmp_D - $signed(AddMux_D);

  assign Cnt_DN = (LoadEn_S) ? OpBShift_DI : (~CntZero_S) ? Cnt_DP - 1 : Cnt_DP;

  assign CntZero_S = ~(|Cnt_DP);

  always_comb begin : p_fsm
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

      default:  /* default */;
    endcase
  end

  assign RemSel_SN  = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN  = (LoadEn_S) ? (~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI) : ResInv_SP;

  assign AReg_DN   = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN   = (BRegEn_S) ? BMux_D : BReg_DP;
  assign ResReg_DN = (LoadEn_S) ? '0 : (ResRegEn_S) ? {ABComp_S, ResReg_DP[$high(ResReg_DP):1]} : ResReg_DP;

  always_ff @(posedge Clk_CI or negedge Rst_RBI) begin : p_regs
    if (~Rst_RBI) begin
      State_SP   <= IDLE;
      AReg_DP    <= '0;
      BReg_DP    <= '0;
      ResReg_DP  <= '0;
      Cnt_DP     <= '0;
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

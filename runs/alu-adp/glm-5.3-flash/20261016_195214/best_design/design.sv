// ============================================================================
// Optimized execution-stage ALU of the CV32E40P (RI5CY) RISC-V core.
//
// I/O behaviour is cycle-identical to the upstream implementation
// (starting_point.sv).  Restructured for area/delay:
//   * compact bit-slice opcode decode (no 7-bit equality comparators)
//   * per-lane 8-bit adder with explicit carry chaining (the upstream's
//     36-bit adder with interleaved carry bits, flattened per lane)
//   * one shared per-lane barrel shifter (left shifts / rotations via the
//     bit-reversal trick of the original)
//   * flat result selection
//   * the serial divider (cv32e40p_alu_div) is kept unchanged: it defines the
//     cycle-accurate div/rem handshake.
// ============================================================================

module cv32e40p_alu (
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

  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8  = 2'b11;

  // ==========================================================================
  // Opcode decode.  Only the 57 documented encodings are ever driven, so the
  // families can be decoded from short bit-slices of operator_i.
  // ==========================================================================
  logic op_addfam;   // ADD/SUB/ADDU/SUBU/ADDR/SUBR/ADDUR/SUBUR   0x18..0x1f
  logic op_mmblk;    // MIN/MINU/MAX/MAXU/ABS/AND/CLIP/CLIPU      0x10..0x17
  logic op_bmxblk;   // BEXT/BEXTU/BINS/BCLR/BSET/INS/OR/XOR      0x28..0x2f
  logic op_div;      // DIVU/DIV/REMU/REM                         0x30..0x33
  logic op_bitcnt;   // CNT/CLB/FF1/FL1                           0x34..0x37
  logic op_cmpblk;   // the 14 comparison opcodes                 0x00..0x0d

  logic op_sub, op_round;
  logic op_minmax, op_abs, op_clip, op_clipu;
  logic op_sra, op_srl, op_ror, op_sll;
  logic op_bext, op_bins;
  logic op_ff1, op_fl1;
  logic op_exts;

  assign op_addfam  = (operator_i[6:3] == 4'b0011);
  assign op_mmblk   = (operator_i[6:3] == 4'b0010);
  assign op_bmxblk  = (operator_i[6:3] == 4'b0101);
  assign op_div     = (operator_i[6:2] == 5'b01100);
  assign op_bitcnt  = (operator_i[6:2] == 5'b01101);
  assign op_cmpblk  = ~operator_i[6] & ~operator_i[5] & ~operator_i[4];

  assign op_sub    = op_addfam &  operator_i[0];
  assign op_round  = op_addfam &  operator_i[2];
  assign op_minmax = op_mmblk   & ~operator_i[2] & ~operator_i[1];  // MIN..MAXU
  assign op_abs    = op_mmblk   &  operator_i[2] & ~operator_i[1];  // ABS
  assign op_clip   = op_mmblk   &  operator_i[2] &  operator_i[1];  // CLIP/CLIPU
  assign op_clipu  = op_clip   &  operator_i[0];
  assign op_sra    = (operator_i[6:2] == 5'b01001) & ~operator_i[1] & ~operator_i[0];
  assign op_srl    = (operator_i[6:2] == 5'b01001) & ~operator_i[1] &  operator_i[0];
  assign op_ror    = (operator_i[6:2] == 5'b01001) &  operator_i[1] & ~operator_i[0];
  assign op_sll    = (operator_i[6:2] == 5'b01001) &  operator_i[1] &  operator_i[0];
  assign op_bext   = op_bmxblk  & (operator_i[2:0] == 3'b000);
  assign op_bins   = op_bmxblk  & (operator_i[2:0] == 3'b010);
  assign op_ff1    = op_bitcnt  &  operator_i[1] & ~operator_i[0];
  assign op_fl1    = op_bitcnt  &  operator_i[1] &  operator_i[0];
  assign op_exts   = (operator_i[6:2] == 5'b01111) & ~operator_i[0];
  assign op_brev   = (operator_i[6:3] == 4'b1001) & (operator_i[2:0] == 3'b001);

  logic vm16, vm8;
  assign vm16 =  vector_mode_i[1] & ~vector_mode_i[0];
  assign vm8  =  vector_mode_i[1] &  vector_mode_i[0];

  // ==========================================================================
  // Bit mask (bit-manipulation group and the rounding adder family)
  // ==========================================================================
  logic [31:0] bmask, bmask_inv;
  assign bmask     = (~{32'hFFFFFFFE << bmask_a_i}) << bmask_b_i;
  assign bmask_inv = ~bmask;

  // ==========================================================================
  // Per-lane adder.
  //
  // The upstream feeds a 36-bit adder whose layout is, per lane k,
  //   bit 9k   : carry-in of the lane
  //   bit 9k+8k: operand byte (adder_in[9k+1 +: 8] = byte k)
  // and extracts result byte k from bits 9k+1..9k+8.  The static carry-in
  // constants implement (byte k = (a_k + b_k + cin_k) mod 256):
  //   lane 0      : cin = 1 for subtract / abs / clip, else 0
  //   inner lanes : cin = carry of the previous lane
  //   boundary lanes (partitioned vector modes):
  //                 cin = 1 in the negate branch (the crossing carry is
  //                 swallowed there), else 0
  // The lane carry-out equals the 8-bit adder carry in every case.
  // ==========================================================================
  logic        nb;            // negate branch: sub / subrot / abs / clip
  logic [31:0] adder_a, adder_b;
  logic [8:0]  adder_s0, adder_s1, adder_s2, adder_s3;
  logic [31:0] adder_result;

  assign nb = op_sub | is_subrot_i | op_abs | op_clip;

  assign adder_a = op_abs      ? ~operand_a_i
                 : is_subrot_i ? {operand_b_i[15:0], operand_a_i[31:16]}
                               : operand_a_i;
  assign adder_b = (op_sub | is_subrot_i)
                     ? (is_subrot_i ? ~{operand_a_i[15:0], operand_b_i[31:16]}
                                    : ~operand_b_i)
                     : operand_b_i;

  assign adder_s0 = {1'b0, adder_a[ 7: 0]} + {1'b0, adder_b[ 7: 0]} + {8'b0, nb};
  assign adder_s1 = {1'b0, adder_a[15: 8]} + {1'b0, adder_b[15: 8]}
                    + ((vm8 ? nb : adder_s0[8]) ? 9'd1 : 9'd0);
  assign adder_s2 = {1'b0, adder_a[23:16]} + {1'b0, adder_b[23:16]}
                    + (((vm16 | vm8) ? nb : adder_s1[8]) ? 9'd1 : 9'd0);
  assign adder_s3 = {1'b0, adder_a[31:24]} + {1'b0, adder_b[31:24]}
                    + ((vm8 ? nb : adder_s2[8]) ? 9'd1 : 9'd0);

  assign adder_result = {adder_s3[7:0], adder_s2[7:0], adder_s1[7:0], adder_s0[7:0]};

  // ==========================================================================
  // Normalization stage (adder family only): optional rounding constant,
  // then a per-lane right shift.
  // ==========================================================================
  logic [31:0] adder_round_result;
  logic [31:0] shift_amt_norm;
  logic [31:0] sh_amt, sh_op;
  logic        sh_left, sh_arith;
  logic [31:0] shift_amt_left;
  logic [31:0] shift_amt;
  logic [31:0] shift_right_result, shift_left_result, shift_result;
  logic [15:0] clpx_shift_ex;
  logic        op_brev;

  assign adder_round_result = op_round ? (adder_result + {1'b0, bmask[31:1]})
                                       : adder_result;

  assign clpx_shift_ex   = $unsigned(clpx_shift_i);
  assign shift_amt_norm  = is_clpx_i ? {clpx_shift_ex, clpx_shift_ex}
                                     : {4{3'b000, bmask_b_i}};

  // amounts for the right shifter; div/rem use the leading-bit count of a
  logic [5:0]  div_shift;
  logic        div_signed, div_op_a_signed;
  logic        div_valid;
  logic [31:0] operand_a_rev, operand_a_neg_rev;
  genvar kr;
  generate
    for (kr = 0; kr < 32; kr++) begin : gen_a_rev
      assign operand_a_rev[kr]     = operand_a_i[31-kr];
      assign operand_a_neg_rev[kr] = ~operand_a_i[31-kr];
    end
  endgenerate
  logic [5:0]  cnt_result;
  logic [4:0]  ff1_result;
  logic        ff_no_one;
  logic [31:0] ff_input;

  assign shift_amt = div_valid ? {26'b0, div_shift} : operand_b_i;

  // by reversing the bits of the input, the order of the shift amounts is
  // reversed as well
  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_amt_left[15:0]  = shift_amt[31:16];
        shift_amt_left[31:16] = shift_amt[15:0];
      end
      VEC_MODE8: begin
        shift_amt_left[ 7: 0] = shift_amt[31:24];
        shift_amt_left[15: 8] = shift_amt[23:16];
        shift_amt_left[23:16] = shift_amt[15:8];
        shift_amt_left[31:24] = shift_amt[ 7:0];
      end
      default: begin
        shift_amt_left[31:0] = shift_amt[31:0];
      end
    endcase
  end

  assign sh_left  = op_sll | op_bins | op_brev | op_div;
  assign sh_arith = op_sra | op_bext | (op_addfam & ~operator_i[1]);

  assign sh_op  = sh_left ? operand_a_rev
                          : (op_addfam ? adder_round_result : operand_a_i);
  assign sh_amt = op_addfam ? shift_amt_norm
                            : (sh_left ? shift_amt_left : shift_amt);

  // per-lane right shifter
  logic [63:0] sh_op_32;
  assign sh_op_32 = op_ror ? {sh_op, sh_op}
                           : $signed({{32{sh_arith & sh_op[31]}}, sh_op});

  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: begin
        shift_right_result[31:16] = $signed({sh_arith & sh_op[31], sh_op[31:16]}) >>> sh_amt[19:16];
        shift_right_result[15:0]  = $signed({sh_arith & sh_op[15], sh_op[15:0]})  >>> sh_amt[ 3:0];
      end
      VEC_MODE8: begin
        shift_right_result[31:24] = $signed({sh_arith & sh_op[31], sh_op[31:24]}) >>> sh_amt[26:24];
        shift_right_result[23:16] = $signed({sh_arith & sh_op[23], sh_op[23:16]}) >>> sh_amt[18:16];
        shift_right_result[15: 8] = $signed({sh_arith & sh_op[15], sh_op[15: 8]}) >>> sh_amt[10: 8];
        shift_right_result[ 7: 0] = $signed({sh_arith & sh_op[ 7], sh_op[ 7: 0]}) >>> sh_amt[ 2: 0];
      end
      default: begin
        shift_right_result = sh_op_32 >> sh_amt[4:0];
      end
    endcase
  end

  genvar j;
  generate
    for (j = 0; j < 32; j++) begin : gen_shift_left_result
      assign shift_left_result[j] = shift_right_result[31-j];
    end
  endgenerate

  assign shift_result = sh_left ? shift_left_result : shift_right_result;

  // ==========================================================================
  // Vector comparisons
  // ==========================================================================
  logic [3:0] cmp_signed;
  logic [3:0] is_equal_vec, is_greater_vec;
  logic [3:0] is_equal, is_greater;
  logic [3:0] cmp_result;
  logic       is_equal_clip;
  logic [31:0] operand_b_eq;

  // signed comparison selected by the opcode
  logic cmp_signed_en;
  assign cmp_signed_en = (op_cmpblk | op_mmblk) & (~operator_i[0] | op_clipu);

  always_comb begin
    case (vector_mode_i)
      VEC_MODE8:  cmp_signed = {4{cmp_signed_en}};
      VEC_MODE16: cmp_signed = cmp_signed_en ? 4'b1010 : 4'b0000;
      default:    cmp_signed = cmp_signed_en ? 4'b1000 : 4'b0000;
    endcase
  end

  genvar i;
  generate
    for (i = 0; i < 4; i++) begin : gen_is_vec
      assign is_equal_vec[i] = (operand_a_i[8*i+7:8*i] == operand_b_i[8*i+7:8*i]);
      assign is_greater_vec[i] = $signed(
          {operand_a_i[8*i+7] & cmp_signed[i], operand_a_i[8*i+7:8*i]}
        ) > $signed(
          {operand_b_i[8*i+7] & cmp_signed[i], operand_b_i[8*i+7:i*8]}
        );
    end
  endgenerate

  always_comb begin
    is_equal[3:0]   = {4{is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_equal_vec[0]}};
    is_greater[3:0] = {4{is_greater_vec[3] | (is_equal_vec[3] & (is_greater_vec[2]
                      | (is_equal_vec[2] & (is_greater_vec[1]
                       | (is_equal_vec[1] & is_greater_vec[0])))))}};
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

  // LT/LTU/SLTS/SLTU (op[2]=0): a<b;  LE*/SLET* (op[2]=1): a<=b
  // GT*/GES (op[3]=1, op[2]=0/1...); EQ/NE (op[3]=1, op[1]=1/0 after EQ=100)
  always_comb begin
    if (!operator_i[3]) begin
      // 0x00..0x07: LTS LTU SLTS SLTU LES LEU SLETS SLETU
      cmp_result = operator_i[2] ? ~is_greater : ~(is_greater | is_equal);
    end else begin
      // 0x08..0x0d: GTS GTU GES GEU EQ NE
      case (operator_i[2:1])
        2'b00:   cmp_result = is_greater;             // GTS / GTU
        2'b01:   cmp_result = is_greater | is_equal;  // GES / GEU
        2'b10:   cmp_result = operator_i[0] ? ~is_equal : is_equal;  // NE / EQ
        default: cmp_result = ~is_equal;
      endcase
    end
  end

  assign comparison_result_o = cmp_result[3];

  // ==========================================================================
  // Min / max / abs / clip
  // ==========================================================================
  logic [31:0] minmax_b;
  logic [3:0]  sel_minmax;
  logic [31:0] result_minmax;
  logic [31:0] clip_result;

  assign minmax_b = op_abs ? adder_result : operand_b_i;
  assign sel_minmax = is_greater ^ {4{(op_minmax | op_clip)}};

  assign result_minmax[31:24] = sel_minmax[3] ? operand_a_i[31:24] : minmax_b[31:24];
  assign result_minmax[23:16] = sel_minmax[2] ? operand_a_i[23:16] : minmax_b[23:16];
  assign result_minmax[15: 8] = sel_minmax[1] ? operand_a_i[15: 8] : minmax_b[15: 8];
  assign result_minmax[ 7: 0] = sel_minmax[0] ? operand_a_i[ 7: 0] : minmax_b[ 7: 0];

  // second equality comparator for clip
  assign operand_b_eq  = op_clipu ? 32'b0 : ~operand_b_i;
  assign is_equal_clip = (operand_a_i == operand_b_eq);

  // CLIP: clamp(a, ~b, b) -> select ~b iff a <s ~b (borrow of a-(~b) = ~carry
  // of a+b+1, both signs equal) or a == ~b
  assign clip_result = op_clipu ? ((operand_a_i[31] | is_equal_clip) ? 32'b0
                                                                     : result_minmax)
                                : (((operand_a_i[31] & operand_b_i[31]) |
                                    ((operand_a_i[31] ^ operand_b_i[31]) & ~adder_s3[8]) |
                                    is_equal_clip) ? ~operand_b_i
                                                   : result_minmax);

  // ==========================================================================
  // Shuffle / pack / extract / insert
  // ==========================================================================
  logic [3:0][1:0] shuffle_byte_sel;
  logic [3:0]      shuffle_reg_sel;
  logic [1:0]      shuffle_reg1_sel;
  logic [1:0]      shuffle_reg0_sel;
  logic [3:0]      shuffle_through;
  logic [31:0]     shuffle_r1, shuffle_r0;
  logic [31:0]     shuffle_r1_in, shuffle_r0_in;
  logic [31:0]     shuffle_result;
  logic [31:0]     pack_result;

  always_comb begin
    shuffle_reg_sel  = '0;
    shuffle_reg1_sel = 2'b01;
    shuffle_reg0_sel = 2'b10;
    shuffle_through  = '1;

    case (operator_i)
      7'h3f, 7'h3e: begin  // EXT / EXTS
        if (op_exts) shuffle_reg1_sel = 2'b11;
        if (vm8) begin
          shuffle_reg_sel[3:1] = 3'b111;
          shuffle_reg_sel[0]   = 1'b0;
        end else begin
          shuffle_reg_sel[3:2] = 2'b11;
          shuffle_reg_sel[1:0] = 2'b00;
        end
      end

      7'h38: begin  // PCKLO
        shuffle_reg1_sel = 2'b00;
        if (vm8) begin
          shuffle_through = 4'b0011;
          shuffle_reg_sel = 4'b0001;
        end else begin
          shuffle_reg_sel = 4'b0011;
        end
      end

      7'h39: begin  // PCKHI
        shuffle_reg1_sel = 2'b00;
        if (vm8) begin
          shuffle_through = 4'b1100;
          shuffle_reg_sel = 4'b0100;
        end else begin
          shuffle_reg_sel = 4'b0011;
        end
      end

      7'h3b: begin  // SHUF2
        case (vector_mode_i)
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

      7'h2d: begin  // INS
        case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_reg0_sel = 2'b00;
            case (imm_vec_ext_i)
              2'b00: shuffle_reg_sel[3:0] = 4'b1110;
              2'b01: shuffle_reg_sel[3:0] = 4'b1101;
              2'b10: shuffle_reg_sel[3:0] = 4'b1011;
              2'b11: shuffle_reg_sel[3:0] = 4'b0111;
            endcase
          end
          VEC_MODE16: begin
            shuffle_reg0_sel   = 2'b01;
            shuffle_reg_sel[3] = ~imm_vec_ext_i[0];
            shuffle_reg_sel[2] = ~imm_vec_ext_i[0];
            shuffle_reg_sel[1] =  imm_vec_ext_i[0];
            shuffle_reg_sel[0] =  imm_vec_ext_i[0];
          end
          default: ;
        endcase
      end

      default: ;
    endcase
  end

  always_comb begin
    shuffle_byte_sel = '0;
    case (operator_i)
      7'h3f, 7'h3e: begin  // EXT / EXTS
        case (vector_mode_i)
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

      7'h38: begin  // PCKLO
        if (vm16) begin
          shuffle_byte_sel[3] = 2'b01;
          shuffle_byte_sel[2] = 2'b00;
          shuffle_byte_sel[1] = 2'b01;
          shuffle_byte_sel[0] = 2'b00;
        end
      end

      7'h39: begin  // PCKHI
        if (vm16) begin
          shuffle_byte_sel[3] = 2'b11;
          shuffle_byte_sel[2] = 2'b10;
          shuffle_byte_sel[1] = 2'b11;
          shuffle_byte_sel[0] = 2'b10;
        end
      end

      7'h3b, 7'h3a: begin  // SHUF2 / SHUF
        case (vector_mode_i)
          VEC_MODE8: begin
            shuffle_byte_sel[3] = operand_b_i[25:24];
            shuffle_byte_sel[2] = operand_b_i[17:16];
            shuffle_byte_sel[1] = operand_b_i[ 9: 8];
            shuffle_byte_sel[0] = operand_b_i[ 1: 0];
          end
          VEC_MODE16: begin
            shuffle_byte_sel[3] = {operand_b_i[16], 1'b1};
            shuffle_byte_sel[2] = {operand_b_i[16], 1'b0};
            shuffle_byte_sel[1] = {operand_b_i[ 0], 1'b1};
            shuffle_byte_sel[0] = {operand_b_i[ 0], 1'b0};
          end
          default: ;
        endcase
      end

      7'h2d: begin  // INS
        shuffle_byte_sel[3] = 2'b11;
        shuffle_byte_sel[2] = 2'b10;
        shuffle_byte_sel[1] = 2'b01;
        shuffle_byte_sel[0] = 2'b00;
      end

      default: ;
    endcase
  end

  assign shuffle_r0_in = shuffle_reg0_sel[1] ? operand_a_i
                        : (shuffle_reg0_sel[0] ? {2{operand_a_i[15:0]}} : {4{operand_a_i[7:0]}});

  assign shuffle_r1_in = shuffle_reg1_sel[1] ? {
                            {8{operand_a_i[31]}}, {8{operand_a_i[23]}},
                            {8{operand_a_i[15]}}, {8{operand_a_i[7]}}
                          }
                        : (shuffle_reg1_sel[0] ? operand_c_i : operand_b_i);

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
  assign shuffle_result[15: 8] = shuffle_reg_sel[1] ? shuffle_r1[15: 8] : shuffle_r0[15: 8];
  assign shuffle_result[ 7: 0] = shuffle_reg_sel[0] ? shuffle_r1[ 7: 0] : shuffle_r0[ 7: 0];

  assign pack_result[31:24] = shuffle_through[3] ? shuffle_result[31:24] : operand_c_i[31:24];
  assign pack_result[23:16] = shuffle_through[2] ? shuffle_result[23:16] : operand_c_i[23:16];
  assign pack_result[15: 8] = shuffle_through[1] ? shuffle_result[15: 8] : operand_c_i[15: 8];
  assign pack_result[ 7: 0] = shuffle_through[0] ? shuffle_result[ 7: 0] : operand_c_i[ 7: 0];

  // ==========================================================================
  // Bit counting (popcount, find-first-one)
  // ==========================================================================
  cv32e40p_popcnt popcnt_i (
      .in_i    (operand_a_i),
      .result_o(cnt_result)
  );

  // the find-first-one tree is shared with the divider normalization
  always_comb begin
    ff_input = '0;
    if (op_bitcnt | op_div) begin
      if (op_ff1) begin
        ff_input = operand_a_i;
      end else if (op_fl1 | (op_div & ~operator_i[0])) begin
        ff_input = operand_a_rev;
      end else begin
        // CLB and the signed div/rem ops
        ff_input = operand_a_i[31] ? operand_a_neg_rev : operand_a_rev;
      end
    end
  end

  cv32e40p_ff_one ff_one_i (
      .in_i       (ff_input),
      .first_one_o(ff1_result),
      .no_ones_o  (ff_no_one)
  );

  logic [5:0] clb_result;
  logic [5:0] bitop_result;

  assign clb_result = {1'b0, ff1_result - 5'd1};  // 5-bit wrap, as upstream

  always_comb begin
    bitop_result = '0;
    case (operator_i)
      7'h36:   bitop_result = ff_no_one ? 6'd32 : {1'b0, ff1_result};              // FF1
      7'h37:   bitop_result = ff_no_one ? 6'd32 : {1'b0, 5'd31 - ff1_result};      // FL1
      7'h34:   bitop_result = cnt_result;                                          // CNT
      7'h35:   bitop_result = ff_no_one ? (operand_a_i[31] ? 6'd31 : 6'd0)         // CLB
                                        : clb_result;
      default: ;
    endcase
  end

  // ==========================================================================
  // Bit manipulation (insert / extract / clear / set / reverse)
  // ==========================================================================
  logic        extract_sign;
  logic [31:0] bextins_and, bextins_result, bclr_result, bset_result;
  logic [31:0] radix_2_rev, radix_4_rev, radix_8_rev, reverse_result;

  assign extract_sign   = op_bext & shift_result[bmask_a_i];
  assign bextins_and    = op_bins ? operand_c_i : {32{extract_sign}};
  assign bextins_result = (bmask & shift_result) | (bextins_and & bmask_inv);
  assign bclr_result    = operand_a_i & bmask_inv;
  assign bset_result    = operand_a_i | bmask;

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
    assign radix_8_rev[31:30] = 2'b0;
  endgenerate

  always_comb begin
    case (bmask_a_i[1:0])
      2'b00:   reverse_result = radix_2_rev;
      2'b01:   reverse_result = radix_4_rev;
      2'b10:   reverse_result = radix_8_rev;
      default: reverse_result = radix_2_rev;
    endcase
  end

  // ==========================================================================
  // Serial divider (dividend = operand_b_i, divisor = operand_a_i)
  // ==========================================================================
  logic [31:0] result_div;
  logic        div_ready;

  assign div_op_a_signed = operand_a_i[31] & operator_i[0];
  assign div_valid       = enable_i & op_div;
  assign div_shift       = (ff_no_one ? 6'd31 : clb_result)
                           + (div_op_a_signed ? 6'd0 : 6'd1);

  cv32e40p_alu_div alu_div_i (
      .Clk_CI      (clk),
      .Rst_RBI     (rst_n),
      .OpA_DI      (operand_b_i),
      .OpB_DI      (shift_left_result),
      .OpBShift_DI (div_shift),
      .OpBIsZero_SI(cnt_result == 6'b0),
      .OpBSign_SI  (div_op_a_signed),
      .OpCode_SI   (operator_i[1:0]),
      .Res_DO      (result_div),
      .InVld_SI    (div_valid),
      .OutRdy_SI   (ex_ready_i),
      .OutVld_SO   (div_ready)
  );

  // ==========================================================================
  // Result selection
  // ==========================================================================
  always_comb begin
    result_o = '0;
    case (operator_i)
      // logic
      7'h15:   result_o = operand_a_i & operand_b_i;  // AND
      7'h2e:   result_o = operand_a_i | operand_b_i;  // OR
      7'h2f:   result_o = operand_a_i ^ operand_b_i;  // XOR

      // shifts and the adder family (normalized through the shifter)
      7'h18, 7'h19, 7'h1a, 7'h1b, 7'h1c, 7'h1d, 7'h1e, 7'h1f,
      7'h24, 7'h25, 7'h26, 7'h27:
               result_o = shift_result;

      // bit manipulation
      7'h28, 7'h29, 7'h2a:
               result_o = bextins_result;               // BEXT/BEXTU/BINS
      7'h2b:   result_o = bclr_result;                 // BCLR
      7'h2c:   result_o = bset_result;                 // BSET
      7'h49:   result_o = reverse_result;              // BREV

      // shuffle / pack / extract / insert
      7'h3a, 7'h3b, 7'h38, 7'h39, 7'h3e, 7'h3f, 7'h2d:
               result_o = pack_result;

      // min / max / abs / clip
      7'h10, 7'h11, 7'h12, 7'h13:
               result_o = result_minmax;
      7'h14:   result_o = is_clpx_i ? {adder_result[31:16], operand_a_i[15:0]}
                                    : result_minmax;    // ABS
      7'h16, 7'h17:
               result_o = clip_result;

      // vector comparisons: replicate the flag into every lane
      7'h00, 7'h01, 7'h04, 7'h05, 7'h08, 7'h09, 7'h0a, 7'h0b, 7'h0c, 7'h0d: begin
        result_o[31:24] = {8{cmp_result[3]}};
        result_o[23:16] = {8{cmp_result[2]}};
        result_o[15: 8] = {8{cmp_result[1]}};
        result_o[ 7: 0] = {8{cmp_result[0]}};
      end

      // scalar set-less
      7'h02, 7'h03, 7'h06, 7'h07:
               result_o = {31'b0, cmp_result[3]};

      // bit counting
      7'h34, 7'h35, 7'h36, 7'h37:
               result_o = {26'h0, bitop_result};

      // division
      7'h30, 7'h31, 7'h32, 7'h33:
               result_o = result_div;

      default: ;
    endcase
  end

  assign ready_o = div_ready;

endmodule

// ---------------------------------------------------------------------------
// population count
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// find first one
// ---------------------------------------------------------------------------
module cv32e40p_ff_one #(
    parameter LEN = 32
) (
    input  logic [LEN-1:0]          in_i,
    output logic [$clog2(LEN)-1:0]  first_one_o,
    output logic                    no_ones_o
);

  localparam NUM_LEVELS = $clog2(LEN);

  logic [LEN-1:0][NUM_LEVELS-1:0]           index_lut;
  logic [2**NUM_LEVELS-1:0]                 sel_nodes;
  logic [2**NUM_LEVELS-1:0][NUM_LEVELS-1:0] index_nodes;

  generate
    genvar j;
    for (j = 0; j < LEN; j++) begin : gen_index_lut
      assign index_lut[j] = $unsigned(j);
    end
  endgenerate

  generate
    genvar k, l2;
    assign sel_nodes[2**NUM_LEVELS-1] = 1'b0;

    for (genvar level = 0; level < NUM_LEVELS; level++) begin : gen_tree
      if (level < NUM_LEVELS - 1) begin : gen_non_root_level
        for (l2 = 0; l2 < 2 ** level; l2++) begin : gen_node
          assign sel_nodes[2**level-1+l2]   = sel_nodes[2**(level+1)-1+l2*2] | sel_nodes[2**(level+1)-1+l2*2+1];
          assign index_nodes[2**level-1+l2] = (sel_nodes[2**(level+1)-1+l2*2] == 1'b1) ?
                                              index_nodes[2**(level+1)-1+l2*2] :
                                              index_nodes[2**(level+1)-1+l2*2+1];
        end
      end
      if (level == NUM_LEVELS - 1) begin : gen_root_level
        for (k = 0; k < 2 ** level; k++) begin : gen_node
          if (k * 2 < LEN - 1) begin : gen_two
            assign sel_nodes[2**level-1+k]   = in_i[k*2] | in_i[k*2+1];
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

// ---------------------------------------------------------------------------
// simple serial divider for signed integers (int32) - unchanged
// ---------------------------------------------------------------------------
module cv32e40p_alu_div #(
    parameter C_WIDTH     = 32,
    parameter C_LOG_WIDTH = 6
) (
    input  logic                   Clk_CI,
    input  logic                   Rst_RBI,
    // input IF
    input  logic [    C_WIDTH-1:0] OpA_DI,
    input  logic [    C_WIDTH-1:0] OpB_DI,
    input  logic [C_LOG_WIDTH-1:0] OpBShift_DI,
    input  logic                   OpBIsZero_SI,
    //
    input  logic                   OpBSign_SI,  // gate this to 0 in case of unsigned ops
    input  logic [            1:0] OpCode_SI,   // 0: udiv, 2: urem, 1: div, 3: rem
    // handshake
    input  logic                   InVld_SI,
    // output IF
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

  // datapath
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

  assign ABComp_S = ((AReg_DP == BReg_DP) | ((AReg_DP > BReg_DP) ^ CompInv_SP)) & ((|AReg_DP) | OpBIsZero_SI);

  assign AddTmp_D = (LoadEn_S) ? 0 : AReg_DP;
  assign AddOut_D = (PmSel_S) ? AddTmp_D + AddMux_D : AddTmp_D - $signed(AddMux_D);

  // counter
  assign Cnt_DN = (LoadEn_S) ? OpBShift_DI : (~CntZero_S) ? Cnt_DP - 1 : Cnt_DP;

  assign CntZero_S = ~(|Cnt_DP);

  // FSM
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

      default: ;
    endcase
  end

  // regs
  assign RemSel_SN  = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN  = (LoadEn_S) ? (~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] & (OpA_DI[$high(
      OpA_DI
  )] ^ OpBSign_SI) : ResInv_SP;

  assign AReg_DN   = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN   = (BRegEn_S) ? BMux_D : BReg_DP;
  assign ResReg_DN = (LoadEn_S) ? '0 : (ResRegEn_S) ? {
    ABComp_S, ResReg_DP[$high(ResReg_DP):1]
  } : ResReg_DP;

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

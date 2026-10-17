// -----------------------------------------------------------------------------
// Self-contained (flattened) copy of the CV32E40P ALU.
//
// Source:  https://github.com/openhwgroup/cv32e40p  (rtl/cv32e40p_alu.sv,
//          rtl/cv32e40p_popcnt.sv, rtl/cv32e40p_ff_one.sv, rtl/cv32e40p_alu_div.sv,
//          rtl/include/cv32e40p_pkg.sv), master branch, fetched 2026-07-15.
// License: Solderpad Hardware License v0.51 (headers retained below).
//
// Mechanical changes only (no logic edits):
//   * cv32e40p_pkg import removed; ALU_* opcode and VEC_MODE* constants inlined
//     as localparams; the alu_opcode_e port is typed logic [6:0].
//   * The four modules are concatenated into this single file.
// Pristine originals: see _debug/orig/ in this benchmark's source directory.
// -----------------------------------------------------------------------------

// Copyright 2018 ETH Zurich and University of Bologna.
// Copyright and related rights are licensed under the Solderpad Hardware
// License, Version 0.51 (the "License"); you may not use this file except in
// compliance with the License.  You may obtain a copy of the License at
// http://solderpad.org/licenses/SHL-0.51. Unless required by applicable law
// or agreed to in writing, software, hardware and materials distributed under
// this License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

////////////////////////////////////////////////////////////////////////////////
// Engineer:       Matthias Baer - baermatt@student.ethz.ch                   //
//                                                                            //
// Additional contributions by:                                               //
//                 Igor Loi - igor.loi@unibo.it                               //
//                 Andreas Traber - atraber@student.ethz.ch                   //
//                 Michael Gautschi - gautschi@iis.ee.ethz.ch                 //
//                 Davide Schiavone - pschiavo@iis.ee.ethz.ch                 //
//                                                                            //
// Design Name:    ALU                                                        //
// Project Name:   RI5CY                                                      //
// Language:       SystemVerilog                                              //
//                                                                            //
// Description:    Arithmetic logic unit of the pipelined processor           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

module cv32e40p_alu
(
    input logic               clk,
    input logic               rst_n,
    input logic               enable_i,
    input logic [ 6:0] operator_i,
    input logic        [31:0] operand_a_i,
    input logic        [31:0] operand_b_i,
    input logic        [31:0] operand_c_i,

    input logic [1:0] vector_mode_i,
    input logic [4:0] bmask_a_i,
    input logic [4:0] bmask_b_i,
    input logic [1:0] imm_vec_ext_i,

    input logic       is_clpx_i,
    input logic       is_subrot_i,
    input logic [1:0] clpx_shift_i,

    output logic [31:0] result_o,
    output logic        comparison_result_o,

    output logic ready_o,
    input  logic ex_ready_i
);

  // ---- constants inlined from cv32e40p_pkg.sv (mechanical flatten) ----
  localparam logic [6:0] ALU_ADD = 7'b0011000;
  localparam logic [6:0] ALU_SUB = 7'b0011001;
  localparam logic [6:0] ALU_ADDU = 7'b0011010;
  localparam logic [6:0] ALU_SUBU = 7'b0011011;
  localparam logic [6:0] ALU_ADDR = 7'b0011100;
  localparam logic [6:0] ALU_SUBR = 7'b0011101;
  localparam logic [6:0] ALU_ADDUR = 7'b0011110;
  localparam logic [6:0] ALU_SUBUR = 7'b0011111;
  localparam logic [6:0] ALU_XOR = 7'b0101111;
  localparam logic [6:0] ALU_OR = 7'b0101110;
  localparam logic [6:0] ALU_AND = 7'b0010101;
  localparam logic [6:0] ALU_SRA = 7'b0100100;
  localparam logic [6:0] ALU_SRL = 7'b0100101;
  localparam logic [6:0] ALU_ROR = 7'b0100110;
  localparam logic [6:0] ALU_SLL = 7'b0100111;
  localparam logic [6:0] ALU_BEXT = 7'b0101000;
  localparam logic [6:0] ALU_BEXTU = 7'b0101001;
  localparam logic [6:0] ALU_BINS = 7'b0101010;
  localparam logic [6:0] ALU_BCLR = 7'b0101011;
  localparam logic [6:0] ALU_BSET = 7'b0101100;
  localparam logic [6:0] ALU_BREV = 7'b1001001;
  localparam logic [6:0] ALU_FF1 = 7'b0110110;
  localparam logic [6:0] ALU_FL1 = 7'b0110111;
  localparam logic [6:0] ALU_CNT = 7'b0110100;
  localparam logic [6:0] ALU_CLB = 7'b0110101;
  localparam logic [6:0] ALU_EXTS = 7'b0111110;
  localparam logic [6:0] ALU_EXT = 7'b0111111;
  localparam logic [6:0] ALU_LTS = 7'b0000000;
  localparam logic [6:0] ALU_LTU = 7'b0000001;
  localparam logic [6:0] ALU_LES = 7'b0000100;
  localparam logic [6:0] ALU_LEU = 7'b0000101;
  localparam logic [6:0] ALU_GTS = 7'b0001000;
  localparam logic [6:0] ALU_GTU = 7'b0001001;
  localparam logic [6:0] ALU_GES = 7'b0001010;
  localparam logic [6:0] ALU_GEU = 7'b0001011;
  localparam logic [6:0] ALU_EQ = 7'b0001100;
  localparam logic [6:0] ALU_NE = 7'b0001101;
  localparam logic [6:0] ALU_SLTS = 7'b0000010;
  localparam logic [6:0] ALU_SLTU = 7'b0000011;
  localparam logic [6:0] ALU_SLETS = 7'b0000110;
  localparam logic [6:0] ALU_SLETU = 7'b0000111;
  localparam logic [6:0] ALU_ABS = 7'b0010100;
  localparam logic [6:0] ALU_CLIP = 7'b0010110;
  localparam logic [6:0] ALU_CLIPU = 7'b0010111;
  localparam logic [6:0] ALU_INS = 7'b0101101;
  localparam logic [6:0] ALU_MIN = 7'b0010000;
  localparam logic [6:0] ALU_MINU = 7'b0010001;
  localparam logic [6:0] ALU_MAX = 7'b0010010;
  localparam logic [6:0] ALU_MAXU = 7'b0010011;
  localparam logic [6:0] ALU_DIVU = 7'b0110000;
  localparam logic [6:0] ALU_DIV = 7'b0110001;
  localparam logic [6:0] ALU_REMU = 7'b0110010;
  localparam logic [6:0] ALU_REM = 7'b0110011;
  localparam logic [6:0] ALU_SHUF = 7'b0111010;
  localparam logic [6:0] ALU_SHUF2 = 7'b0111011;
  localparam logic [6:0] ALU_PCKLO = 7'b0111000;
  localparam logic [6:0] ALU_PCKHI = 7'b0111001;
  localparam logic [1:0] VEC_MODE32 = 2'b00;
  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8 = 2'b11;

  logic [31:0] operand_a_rev;
  logic [31:0] operand_a_neg;
  logic [31:0] operand_a_neg_rev;

  assign operand_a_neg = ~operand_a_i;

  // bit reverse operand_a for left shifts and bit counting
  generate
    genvar k;
    for (k = 0; k < 32; k++) begin : gen_operand_a_rev
      assign operand_a_rev[k] = operand_a_i[31-k];
    end
  endgenerate

  // bit reverse operand_a_neg for left shifts and bit counting
  generate
    genvar m;
    for (m = 0; m < 32; m++) begin : gen_operand_a_neg_rev
      assign operand_a_neg_rev[m] = operand_a_neg[31-m];
    end
  endgenerate

  logic [31:0] operand_b_neg;

  assign operand_b_neg = ~operand_b_i;


  logic [ 5:0] div_shift;
  logic        div_valid;
  logic [31:0] bmask;

  //////////////////////////////////////////////////////////////////////////////////////////
  //   ____            _   _ _   _                      _      _       _     _            //
  //  |  _ \ __ _ _ __| |_(_) |_(_) ___  _ __   ___  __| |    / \   __| | __| | ___ _ __  //
  //  | |_) / _` | '__| __| | __| |/ _ \| '_ \ / _ \/ _` |   / _ \ / _` |/ _` |/ _ \ '__| //
  //  |  __/ (_| | |  | |_| | |_| | (_) | | | |  __/ (_| |  / ___ \ (_| | (_| |  __/ |    //
  //  |_|   \__,_|_|   \__|_|\__|_|\___/|_| |_|\___|\__,_| /_/   \_\__,_|\__,_|\___|_|    //
  //                                                                                      //
  //////////////////////////////////////////////////////////////////////////////////////////

  logic        adder_op_b_negate;
  logic [31:0] adder_op_a, adder_op_b;
  logic [31:0] adder_result;


  assign adder_op_b_negate = (operator_i == ALU_SUB) || (operator_i == ALU_SUBR) ||
                             (operator_i == ALU_SUBU) || (operator_i == ALU_SUBUR) || is_subrot_i;

  // prepare operand a
  assign adder_op_a = (operator_i == ALU_ABS) ? operand_a_neg : (is_subrot_i ? {
    operand_b_i[15:0], operand_a_i[31:16]
  } : operand_a_i);

  // prepare operand b
  assign adder_op_b = adder_op_b_negate ? (is_subrot_i ? ~{
    operand_a_i[15:0], operand_b_i[31:16]
  } : operand_b_neg) : operand_b_i;

  // per-lane adders with group carry lookahead.  Lane semantics match the
  // original separator adder exactly: carries chain across lanes in 32-bit
  // mode, are killed in vector modes, carry-in is adder_cin for sub-style ops.
  logic        adder_cin, clip_cout, clip_flag;
  logic [3:0]  lane_gg, lane_gp, lane_ci;
  logic [3:0]  g_c;
  logic        g_c3_m16;

  assign adder_cin = adder_op_b_negate || (operator_i == ALU_ABS) || (operator_i == ALU_CLIP);

  // full carry chain (32-bit mode)
  assign g_c[0] = adder_cin;
  assign g_c[1] = lane_gg[0] | (lane_gp[0] & g_c[0]);
  assign g_c[2] = lane_gg[1] | (lane_gp[1] & lane_gg[0]) | (lane_gp[1] & lane_gp[0] & g_c[0]);
  assign g_c[3] = lane_gg[2] | (lane_gp[2] & lane_gg[1]) | (lane_gp[2] & lane_gp[1] & lane_gg[0]) |
                  (lane_gp[2] & lane_gp[1] & lane_gp[0] & g_c[0]);
  assign g_c3_m16 = lane_gg[2] | (lane_gp[2] & adder_cin);

  assign lane_ci[0] = adder_cin;
  assign lane_ci[1] = (vector_mode_i == VEC_MODE8) ? adder_cin : g_c[1];
  assign lane_ci[2] = (vector_mode_i == VEC_MODE32) ? g_c[2] : adder_cin;
  assign lane_ci[3] = (vector_mode_i == VEC_MODE8) ? adder_cin :
                      (vector_mode_i == VEC_MODE32) ? g_c[3] : g_c3_m16;

  pp_add #(.W(8)) i_lane0 (
      .a(adder_op_a[7:0]), .b(adder_op_b[7:0]), .cin(lane_ci[0]),
      .sum(adder_result[7:0]), .cout(), .g_out(lane_gg[0]), .p_out(lane_gp[0]));
  pp_add #(.W(8)) i_lane1 (
      .a(adder_op_a[15:8]), .b(adder_op_b[15:8]), .cin(lane_ci[1]),
      .sum(adder_result[15:8]), .cout(), .g_out(lane_gg[1]), .p_out(lane_gp[1]));
  pp_add #(.W(8)) i_lane2 (
      .a(adder_op_a[23:16]), .b(adder_op_b[23:16]), .cin(lane_ci[2]),
      .sum(adder_result[23:16]), .cout(), .g_out(lane_gg[2]), .p_out(lane_gp[2]));
  pp_add #(.W(8)) i_lane3 (
      .a(adder_op_a[31:24]), .b(adder_op_b[31:24]), .cin(lane_ci[3]),
      .sum(adder_result[31:24]), .cout(), .g_out(lane_gg[3]), .p_out(lane_gp[3]));

  assign clip_cout = (vector_mode_i == VEC_MODE32) ? (lane_gg[3] | (lane_gp[3] & g_c[3])) : 1'b0;
  assign clip_flag = clip_cout ^ adder_op_a[31] ^ adder_op_b[31];


  // normalization stage
  logic [31:0] adder_round_value;
  logic [31:0] adder_round_result;

  assign adder_round_value  = ((operator_i == ALU_ADDR) || (operator_i == ALU_SUBR) ||
                               (operator_i == ALU_ADDUR) || (operator_i == ALU_SUBUR)) ?
                                {
    1'b0, bmask[31:1]
  } : '0;
  pp_add #(.W(32)) i_add_round (
      .a(adder_result),
      .b(adder_round_value),
      .cin(1'b0),
      .sum(adder_round_result),
      .cout()
  );


  ////////////////////////////////////////
  //  ____  _   _ ___ _____ _____       //
  // / ___|| | | |_ _|  ___|_   _|      //
  // \___ \| |_| || || |_    | |        //
  //  ___) |  _  || ||  _|   | |        //
  // |____/|_| |_|___|_|     |_|        //
  //                                    //
  ////////////////////////////////////////

  logic        shift_left;  // should we shift left
  logic        shift_use_round;
  logic        shift_arithmetic;

  logic [31:0] shift_amt_left;  // amount of shift, if to the left
  logic [31:0] shift_amt;  // amount of shift, to the right
  logic [31:0] shift_amt_int;  // amount of shift, used for the actual shifters
  logic [31:0] shift_amt_norm;  // amount of shift, used for normalization
  logic [31:0] shift_op_a;  // input of the shifter
  logic [31:0] shift_result;
  logic [31:0] shift_right_result;
  logic [31:0] shift_left_result;
  logic [15:0] clpx_shift_ex;

  // shifter is also used for preparing operand for division
  assign shift_amt = div_valid ? div_shift : operand_b_i;

  // by reversing the bits of the input, we also have to reverse the order of shift amounts
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

      default: // VEC_MODE32
      begin
        shift_amt_left[31:0] = shift_amt[31:0];
      end
    endcase
  end

  // ALU_FL1 and ALU_CBL are used for the bit counting ops later
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

  // choose the bit reversed or the normal input for shift operand a
  assign shift_op_a    = shift_left ? operand_a_rev :
                          (shift_use_round ? adder_round_result : operand_a_i);
  assign shift_amt_int = shift_use_round ? shift_amt_norm :
                          (shift_left ? shift_amt_left : shift_amt);

  assign shift_amt_norm = is_clpx_i ? {clpx_shift_ex, clpx_shift_ex} : {4{3'b000, bmask_b_i}};

  assign clpx_shift_ex = $unsigned(clpx_shift_i);

  // right shifts, we let the synthesizer optimize this
  logic [63:0] shift_op_a_32;

  assign shift_op_a_32 = (operator_i == ALU_ROR) ? {
        shift_op_a, shift_op_a
      } : $signed(
          {{32{shift_arithmetic & shift_op_a[31]}}, shift_op_a}
      );

  // unified per-lane barrel shifter: stages shift right by 1,2,4,8,16 within
  // the lanes of the active vector mode; the per-mode shifters and their
  // output mux are merged into these stages (select-aware stages).
  logic [5:0][31:0] sh_st;
  logic [31:0][5:0] sh_en;
  logic [31:0]      sh_fill;
  logic             is_rot32;

  assign is_rot32 = (operator_i == ALU_ROR) && (vector_mode_i == VEC_MODE32);

  genvar gs, gt;
  generate
    for (gs = 0; gs < 32; gs++) begin : gen_sh_bit
      // fill bit: lane sign for arithmetic shift, rotate source for ROR
      assign sh_fill[gs] = (shift_arithmetic &
                            ((vector_mode_i == VEC_MODE8) ? shift_op_a[{gs[4:3], 3'b111}] :
                             (vector_mode_i == VEC_MODE16) ? shift_op_a[{gs[4], 4'b1111}] :
                             shift_op_a[31])) | (is_rot32 & shift_op_a[gs]);
      // per-lane shift-amount bits
      for (gt = 0; gt < 5; gt++) begin : gen_sh_en
        assign sh_en[gs][gt] =
            (vector_mode_i == VEC_MODE8) ? (((gt < 3) ? shift_amt_int[{gs[4:3], 3'b000}+gt] : 1'b0)) :
            (vector_mode_i == VEC_MODE16) ? (((gt < 4) ? shift_amt_int[{gs[4], 4'b0000}+gt] : 1'b0)) :
            shift_amt_int[gt];
      end
    end

    assign sh_st[0] = shift_op_a;

    for (gs = 0; gs < 32; gs++) begin : gen_sh_b
      for (gt = 0; gt < 5; gt++) begin : gen_sh_s
        if (((1 << gt) + gs) < 32) begin : gen_has_src
          assign sh_st[gt+1][gs] = sh_en[gs][gt] ?
              ((((vector_mode_i == VEC_MODE8) && (((gs & 7) + (1 << gt)) > 7)) ||
                ((vector_mode_i == VEC_MODE16) && (((gs & 15) + (1 << gt)) > 15))) ?
               sh_fill[gs] : sh_st[gt][gs+(1<<gt)]) : sh_st[gt][gs];
        end else begin : gen_no_src
          assign sh_st[gt+1][gs] = sh_en[gs][gt] ?
              (is_rot32 ? sh_st[gt][gs+(1<<gt)-32] : sh_fill[gs]) : sh_st[gt][gs];
        end
      end
    end
  endgenerate

  assign shift_right_result = sh_st[5];

  // bit reverse the shift_right_result for left shifts
  genvar j;
  generate
    for (j = 0; j < 32; j++) begin : gen_shift_left_result
      assign shift_left_result[j] = shift_right_result[31-j];
    end
  endgenerate

  assign shift_result = shift_left ? shift_left_result : shift_right_result;


  //////////////////////////////////////////////////////////////////
  //   ____ ___  __  __ ____   _    ____  ___ ____   ___  _   _   //
  //  / ___/ _ \|  \/  |  _ \ / \  |  _ \|_ _/ ___| / _ \| \ | |  //
  // | |  | | | | |\/| | |_) / _ \ | |_) || |\___ \| | | |  \| |  //
  // | |__| |_| | |  | |  __/ ___ \|  _ < | | ___) | |_| | |\  |  //
  //  \____\___/|_|  |_|_| /_/   \_\_| \_\___|____/ \___/|_| \_|  //
  //                                                              //
  //////////////////////////////////////////////////////////////////

  logic [ 3:0] is_equal;
  logic [ 3:0] is_greater;  // handles both signed and unsigned forms

  // 8-bit vector comparisons, basic building blocks
  logic [ 3:0] cmp_signed;
  logic [ 3:0] is_equal_vec;
  logic [ 3:0] is_greater_vec;
  logic [31:0] operand_b_eq;
  logic        is_equal_clip;


  //second == comparator for CLIP instructions
  always_comb begin
    operand_b_eq = operand_b_neg;
    if (operator_i == ALU_CLIPU) operand_b_eq = '0;
    else operand_b_eq = operand_b_neg;
  end
  assign is_equal_clip = operand_a_i == operand_b_eq;

  always_comb begin
    cmp_signed = 4'b0;

    unique case (operator_i)
      ALU_GTS,
      ALU_GES,
      ALU_LTS,
      ALU_LES,
      ALU_SLTS,
      ALU_SLETS,
      ALU_MIN,
      ALU_MAX,
      ALU_ABS,
      ALU_CLIP,
      ALU_CLIPU: begin
        case (vector_mode_i)
          VEC_MODE8:  cmp_signed[3:0] = 4'b1111;
          VEC_MODE16: cmp_signed[3:0] = 4'b1010;
          default:    cmp_signed[3:0] = 4'b1000;
        endcase
      end

      default: ;
    endcase
  end

  // generate vector equal and greater than signals, cmp_signed decides if the
  // comparison is done signed or unsigned
  genvar i;
  generate
    for (i = 0; i < 4; i++) begin : gen_is_vec
      logic [8:0] cmp_xa, cmp_xb;
      logic       cmp_ge;
      assign is_equal_vec[i] = (operand_a_i[8*i+7:8*i] == operand_b_i[8*i+7:i*8]);
      assign cmp_xa = {1'b0, operand_a_i[8*i+7:8*i] ^ {cmp_signed[i], 7'b0}};
      assign cmp_xb = {1'b0, operand_b_i[8*i+7:i*8] ^ {cmp_signed[i], 7'b0}};
      pp_add #(.W(9)) i_cmp (
          .a(cmp_xa),
          .b(~cmp_xb),
          .cin(1'b1),
          .sum(),
          .cout(cmp_ge)
      );
      assign is_greater_vec[i] = cmp_ge & ~is_equal_vec[i];
    end
  endgenerate

  // generate the real equal and greater than signals that take the vector
  // mode into account (flat AND-OR combines, lanes computed in parallel)
  logic eq32, gt32, eq16h, gt16h, eq16l, gt16l;

  assign eq32  = is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_equal_vec[0];
  assign gt32  = is_greater_vec[3]
               | (is_equal_vec[3] & is_greater_vec[2])
               | (is_equal_vec[3] & is_equal_vec[2] & is_greater_vec[1])
               | (is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_greater_vec[0]);
  assign eq16h = is_equal_vec[3] & is_equal_vec[2];
  assign gt16h = is_greater_vec[3] | (is_equal_vec[3] & is_greater_vec[2]);
  assign eq16l = is_equal_vec[1] & is_equal_vec[0];
  assign gt16l = is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0]);

  always_comb begin
    case (vector_mode_i)
      VEC_MODE8: begin
        is_equal[3:0]   = is_equal_vec[3:0];
        is_greater[3:0] = is_greater_vec[3:0];
      end

      VEC_MODE16: begin
        is_equal[3:0]   = {{2{eq16h}}, {2{eq16l}}};
        is_greater[3:0] = {{2{gt16h}}, {2{gt16l}}};
      end

      default: begin
        is_equal[3:0]   = {4{eq32}};
        is_greater[3:0] = {4{gt32}};
      end
    endcase
  end

  // generate comparison result (balanced AND-OR select, not a priority chain)
  logic [3:0] cmp_result;

  assign cmp_result =
      ({4{(operator_i == ALU_EQ)}} & is_equal) |
      ({4{(operator_i == ALU_NE)}} & ~is_equal) |
      ({4{(operator_i == ALU_GTS) | (operator_i == ALU_GTU)}} & is_greater) |
      ({4{(operator_i == ALU_GES) | (operator_i == ALU_GEU)}} & (is_greater | is_equal)) |
      ({4{(operator_i == ALU_LTS) | (operator_i == ALU_SLTS) | (operator_i == ALU_LTU) |
          (operator_i == ALU_SLTU)}} & ~(is_greater | is_equal)) |
      ({4{(operator_i == ALU_LES) | (operator_i == ALU_SLETS) | (operator_i == ALU_LEU) |
          (operator_i == ALU_SLETU)}} & ~is_greater);

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

  assign result_minmax[31:24] = (sel_minmax[3] == 1'b1) ? operand_a_i[31:24] : minmax_b[31:24];
  assign result_minmax[23:16] = (sel_minmax[2] == 1'b1) ? operand_a_i[23:16] : minmax_b[23:16];
  assign result_minmax[15:8] = (sel_minmax[1] == 1'b1) ? operand_a_i[15:8] : minmax_b[15:8];
  assign result_minmax[7:0] = (sel_minmax[0] == 1'b1) ? operand_a_i[7:0] : minmax_b[7:0];

  //////////////////////////////////////////////////
  // Clip
  //////////////////////////////////////////////////
  logic [31:0] clip_result;  // result of clip and clip

  always_comb begin
    clip_result = result_minmax;
    if (operator_i == ALU_CLIPU) begin
      if (operand_a_i[31] || is_equal_clip) begin
        clip_result = '0;
      end else begin
        clip_result = result_minmax;
      end
    end else begin
      //CLIP
      if (clip_flag || is_equal_clip) begin
        clip_result = operand_b_neg;
      end else begin
        clip_result = result_minmax;
      end
    end

  end

  //////////////////////////////////////////////////
  //  ____  _   _ _   _ _____ _____ _     _____   //
  // / ___|| | | | | | |  ___|  ___| |   | ____|  //
  // \___ \| |_| | | | | |_  | |_  | |   |  _|    //
  //  ___) |  _  | |_| |  _| |  _| | |___| |___   //
  // |____/|_| |_|\___/|_|   |_|   |_____|_____|  //
  //                                              //
  //////////////////////////////////////////////////

  logic [3:0][1:0] shuffle_byte_sel;  // select byte in register: 31:24, 23:16, 15:8, 7:0
  logic [3:0]      shuffle_reg_sel;  // select register: rD/rS2 or rS1
  logic [1:0]      shuffle_reg1_sel;  // select register rD or rS2 for next stage
  logic [1:0]      shuffle_reg0_sel;
  logic [3:0]      shuffle_through;

  logic [31:0] shuffle_r1, shuffle_r0;
  logic [31:0] shuffle_r1_in, shuffle_r0_in;
  logic [31:0] shuffle_result;
  logic [31:0] pack_result;


  always_comb begin
    shuffle_reg_sel  = '0;
    shuffle_reg1_sel = 2'b01;
    shuffle_reg0_sel = 2'b10;
    shuffle_through  = '1;

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
              2'b00: begin
                shuffle_reg_sel[3:0] = 4'b1110;
              end
              2'b01: begin
                shuffle_reg_sel[3:0] = 4'b1101;
              end
              2'b10: begin
                shuffle_reg_sel[3:0] = 4'b1011;
              end
              2'b11: begin
                shuffle_reg_sel[3:0] = 4'b0111;
              end
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
    shuffle_byte_sel = '0;

    // byte selector
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
          VEC_MODE8: begin
            shuffle_byte_sel[3] = 2'b00;
            shuffle_byte_sel[2] = 2'b00;
            shuffle_byte_sel[1] = 2'b00;
            shuffle_byte_sel[0] = 2'b00;
          end

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
          VEC_MODE8: begin
            shuffle_byte_sel[3] = 2'b00;
            shuffle_byte_sel[2] = 2'b00;
            shuffle_byte_sel[1] = 2'b00;
            shuffle_byte_sel[0] = 2'b00;
          end

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
  assign shuffle_r0_in = shuffle_reg0_sel[1] ?
                          operand_a_i :
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
  assign shuffle_result[15:8] = shuffle_reg_sel[1] ? shuffle_r1[15:8] : shuffle_r0[15:8];
  assign shuffle_result[7:0] = shuffle_reg_sel[0] ? shuffle_r1[7:0] : shuffle_r0[7:0];

  assign pack_result[31:24] = shuffle_through[3] ? shuffle_result[31:24] : operand_c_i[31:24];
  assign pack_result[23:16] = shuffle_through[2] ? shuffle_result[23:16] : operand_c_i[23:16];
  assign pack_result[15:8] = shuffle_through[1] ? shuffle_result[15:8] : operand_c_i[15:8];
  assign pack_result[7:0] = shuffle_through[0] ? shuffle_result[7:0] : operand_c_i[7:0];


  /////////////////////////////////////////////////////////////////////
  //   ____  _ _      ____                  _      ___               //
  //  | __ )(_) |_   / ___|___  _   _ _ __ | |_   / _ \ _ __  ___    //
  //  |  _ \| | __| | |   / _ \| | | | '_ \| __| | | | | '_ \/ __|   //
  //  | |_) | | |_  | |__| (_) | |_| | | | | |_  | |_| | |_) \__ \_  //
  //  |____/|_|\__|  \____\___/ \__,_|_| |_|\__|  \___/| .__/|___(_) //
  //                                                   |_|           //
  /////////////////////////////////////////////////////////////////////

  logic [31:0] ff_input;  // either op_a_i or its bit reversed version
  logic [ 5:0] cnt_result;  // population count
  logic [ 5:0] clb_result;  // count leading bits
  logic [ 4:0] ff1_result;  // holds the index of the first '1'
  logic        ff_no_one;  // if no ones are found
  logic [ 4:0] fl1_result;  // holds the index of the last '1'
  logic [ 5:0] bitop_result;  // result of all bitop operations muxed together

  cv32e40p_popcnt popcnt_i (
      .in_i    (operand_a_i),
      .result_o(cnt_result)
  );

  logic ff_sel_ff1, ff_sel_rev, ff_sel_neg;

  assign ff_sel_ff1 = (operator_i == ALU_FF1);
  assign ff_sel_rev = (operator_i == ALU_DIVU) | (operator_i == ALU_REMU) | (operator_i == ALU_FL1);
  assign ff_sel_neg = (operator_i == ALU_DIV) | (operator_i == ALU_REM) | (operator_i == ALU_CLB);

  assign ff_input = ({32{ff_sel_ff1}} & operand_a_i) |
                    ({32{ff_sel_rev | ff_sel_neg}} &
                     (operand_a_rev ^ {32{ff_sel_neg & operand_a_i[31]}}));

  cv32e40p_ff_one ff_one_i (
      .in_i       (ff_input),
      .first_one_o(ff1_result),
      .no_ones_o  (ff_no_one)
  );

  // special case if ff1_res is 0 (no 1 found), then we keep the 0
  // this is done in the result mux
  assign fl1_result = 5'd31 - ff1_result;
  assign clb_result = ff1_result - 5'd1;

  assign bitop_result =
      ({6{(operator_i == ALU_FF1)}} & (ff_no_one ? 6'd32 : {1'b0, ff1_result})) |
      ({6{(operator_i == ALU_FL1)}} & (ff_no_one ? 6'd32 : {1'b0, fl1_result})) |
      ({6{(operator_i == ALU_CNT)}} & cnt_result) |
      ({6{(operator_i == ALU_CLB)}} &
       (ff_no_one ? ({6{operand_a_i[31]}} & 6'd31) : clb_result));


  ////////////////////////////////////////////////
  //  ____  _ _     __  __             _        //
  // | __ )(_) |_  |  \/  | __ _ _ __ (_)_ __   //
  // |  _ \| | __| | |\/| |/ _` | '_ \| | '_ \  //
  // | |_) | | |_  | |  | | (_| | | | | | |_) | //
  // |____/|_|\__| |_|  |_|\__,_|_| |_|_| .__/  //
  //                                    |_|     //
  ////////////////////////////////////////////////

  logic extract_is_signed;
  logic extract_sign;
  logic [31:0] bmask_first, bmask_inv;
  logic [31:0] bextins_and;
  logic [31:0] bextins_result, bclr_result, bset_result;


  // construct bit mask for insert/extract/bclr/bset
  // bmask looks like this 00..0011..1100..00
  //   bmask[k] = (k >= bmask_b_i) && (k < bmask_b_i + bmask_a_i + 1)
  logic [5:0] bmask_end;

  assign bmask_end = {1'b0, bmask_a_i} + {1'b0, bmask_b_i} + 6'd1;
  assign bmask_inv = ~bmask;

  genvar bm;
  generate
    for (bm = 0; bm < 32; bm++) begin : gen_bmask
      assign bmask[bm] = (bmask_b_i <= bm) && (bm < bmask_end);
    end
  endgenerate

  assign bextins_and       = (operator_i == ALU_BINS) ? operand_c_i : {32{extract_sign}};

  assign extract_is_signed = (operator_i == ALU_BEXT);
  assign extract_sign      = extract_is_signed & shift_result[bmask_a_i];

  assign bextins_result    = (bmask & shift_result) | (bextins_and & bmask_inv);

  assign bclr_result       = operand_a_i & bmask_inv;
  assign bset_result       = operand_a_i | bmask;

  /////////////////////////////////////////////////////////////////////////////////
  //  ____ _____ _______     _____  ________      ________ _____   _____ ______  //
  // |  _ \_   _|__   __|   |  __ \|  ____\ \    / /  ____|  __ \ / ____|  ____| //
  // | |_) || |    | |______| |__) | |__   \ \  / /| |__  | |__) | (___ | |__    //
  // |  _ < | |    | |______|  _  /|  __|   \ \/ / |  __| |  _  / \___ \|  __|   //
  // | |_) || |_   | |      | | \ \| |____   \  /  | |____| | \ \ ____) | |____  //
  // |____/_____|  |_|      |_|  \_\______|   \/   |______|_|  \_\_____/|______| //
  //                                                                             //
  /////////////////////////////////////////////////////////////////////////////////

  logic [31:0] radix_2_rev;
  logic [31:0] radix_4_rev;
  logic [31:0] radix_8_rev;
  logic [31:0] reverse_result;
  logic [ 1:0] radix_mux_sel;

  assign radix_mux_sel = bmask_a_i[1:0];

  generate
    // radix-2 bit reverse
    for (j = 0; j < 32; j++) begin : gen_radix_2_rev
      assign radix_2_rev[j] = shift_result[31-j];
    end
    // radix-4 bit reverse
    for (j = 0; j < 16; j++) begin : gen_radix_4_rev
      assign radix_4_rev[2*j+1:2*j] = shift_result[31-j*2:31-j*2-1];
    end
    // radix-8 bit reverse
    for (j = 0; j < 10; j++) begin : gen_radix_8_rev
      assign radix_8_rev[3*j+2:3*j] = shift_result[31-j*3:31-j*3-2];
    end
    assign radix_8_rev[31:30] = 2'b0;
  endgenerate

  always_comb begin
    reverse_result = '0;

    unique case (radix_mux_sel)
      2'b00: reverse_result = radix_2_rev;
      2'b01: reverse_result = radix_4_rev;
      2'b10: reverse_result = radix_8_rev;

      default: reverse_result = radix_2_rev;
    endcase
  end

  ////////////////////////////////////////////////////
  //  ____ _____     __     __  ____  _____ __  __  //
  // |  _ \_ _\ \   / /    / / |  _ \| ____|  \/  | //
  // | | | | | \ \ / /    / /  | |_) |  _| | |\/| | //
  // | |_| | |  \ V /    / /   |  _ <| |___| |  | | //
  // |____/___|  \_/    /_/    |_| \_\_____|_|  |_| //
  //                                                //
  ////////////////////////////////////////////////////

  logic [31:0] result_div;
  logic        div_ready;
  logic        div_signed;
  logic        div_op_a_signed;
  logic [ 5:0] div_shift_int;

  assign div_signed = operator_i[0];

  assign div_op_a_signed = operand_a_i[31] & div_signed;

  assign div_shift_int = ff_no_one ? 6'd31 : clb_result;
  assign div_shift = div_shift_int + (div_op_a_signed ? 6'd0 : 6'd1);

  assign div_valid = enable_i & ((operator_i == ALU_DIV) || (operator_i == ALU_DIVU) ||
                     (operator_i == ALU_REM) || (operator_i == ALU_REMU));

  // inputs A and B are swapped
  cv32e40p_alu_div alu_div_i (
      .Clk_CI (clk),
      .Rst_RBI(rst_n),

      // input IF
      .OpA_DI      (operand_b_i),
      .OpB_DI      (shift_left_result),
      .OpBShift_DI (div_shift),
      .OpBIsZero_SI((cnt_result == 0)),

      .OpBSign_SI(div_op_a_signed),
      .OpCode_SI (operator_i[1:0]),

      .Res_DO(result_div),

      // Hand-Shake
      .InVld_SI (div_valid),
      .OutRdy_SI(ex_ready_i),
      .OutVld_SO(div_ready)
  );

  ////////////////////////////////////////////////////////
  //   ____                 _ _     __  __              //
  //  |  _ \ ___  ___ _   _| | |_  |  \/  |_   ___  __  //
  //  | |_) / _ \/ __| | | | | __| | |\/| | | | \ \/ /  //
  //  |  _ <  __/\__ \ |_| | | |_  | |  | | |_| |>  <   //
  //  |_| \_\___||___/\__,_|_|\__| |_|  |_|\__,_/_/\_\  //
  //                                                    //
  ////////////////////////////////////////////////////////

  // source-group enables for the result mux (balanced AND-OR tree)
  logic [15:0] rsel;

  assign rsel[0]  = (operator_i == ALU_AND);
  assign rsel[1]  = (operator_i == ALU_OR);
  assign rsel[2]  = (operator_i == ALU_XOR);
  assign rsel[3]  = (operator_i == ALU_ADD) | (operator_i == ALU_ADDR) | (operator_i == ALU_ADDU) |
                    (operator_i == ALU_ADDUR) | (operator_i == ALU_SUB) | (operator_i == ALU_SUBR) |
                    (operator_i == ALU_SUBU) | (operator_i == ALU_SUBUR) | (operator_i == ALU_SLL) |
                    (operator_i == ALU_SRL) | (operator_i == ALU_SRA) | (operator_i == ALU_ROR);
  assign rsel[4]  = (operator_i == ALU_BINS) | (operator_i == ALU_BEXT) | (operator_i == ALU_BEXTU);
  assign rsel[5]  = (operator_i == ALU_BCLR);
  assign rsel[6]  = (operator_i == ALU_BSET);
  assign rsel[7]  = (operator_i == ALU_BREV);
  assign rsel[8]  = (operator_i == ALU_SHUF) | (operator_i == ALU_SHUF2) | (operator_i == ALU_PCKLO) |
                    (operator_i == ALU_PCKHI) | (operator_i == ALU_EXT) | (operator_i == ALU_EXTS) |
                    (operator_i == ALU_INS);
  assign rsel[9]  = (operator_i == ALU_MIN) | (operator_i == ALU_MINU) | (operator_i == ALU_MAX) |
                    (operator_i == ALU_MAXU) | ((operator_i == ALU_ABS) & ~is_clpx_i);
  assign rsel[10] = (operator_i == ALU_ABS) & is_clpx_i;
  assign rsel[11] = (operator_i == ALU_CLIP) | (operator_i == ALU_CLIPU);
  assign rsel[12] = (operator_i == ALU_EQ) | (operator_i == ALU_NE) | (operator_i == ALU_GTS) |
                    (operator_i == ALU_GES) | (operator_i == ALU_GTU) | (operator_i == ALU_GEU) |
                    (operator_i == ALU_LTS) | (operator_i == ALU_LES) | (operator_i == ALU_LTU) |
                    (operator_i == ALU_LEU);
  assign rsel[13] = (operator_i == ALU_SLTS) | (operator_i == ALU_SLTU) | (operator_i == ALU_SLETS) |
                    (operator_i == ALU_SLETU);
  assign rsel[14] = (operator_i == ALU_FF1) | (operator_i == ALU_FL1) | (operator_i == ALU_CLB) |
                    (operator_i == ALU_CNT);
  assign rsel[15] = (operator_i == ALU_DIV) | (operator_i == ALU_DIVU) | (operator_i == ALU_REM) |
                    (operator_i == ALU_REMU);

  always_comb begin
    result_o = ({32{rsel[0]}} & (operand_a_i & operand_b_i)) |
               ({32{rsel[1]}} & (operand_a_i | operand_b_i)) |
               ({32{rsel[2]}} & (operand_a_i ^ operand_b_i)) |
               ({32{rsel[3]}} & shift_result) |
               ({32{rsel[4]}} & bextins_result) |
               ({32{rsel[5]}} & bclr_result) |
               ({32{rsel[6]}} & bset_result) |
               ({32{rsel[7]}} & reverse_result) |
               ({32{rsel[8]}} & pack_result) |
               ({32{rsel[9]}} & result_minmax) |
               ({32{rsel[10]}} & {adder_result[31:16], operand_a_i[15:0]}) |
               ({32{rsel[11]}} & clip_result) |
               ({32{rsel[12]}} & {{8{cmp_result[3]}}, {8{cmp_result[2]}},
                                  {8{cmp_result[1]}}, {8{cmp_result[0]}}}) |
               ({32{rsel[13]}} & {31'b0, comparison_result_o}) |
               ({32{rsel[14]}} & {26'h0, bitop_result[5:0]}) |
               ({32{rsel[15]}} & result_div);
  end

  assign ready_o = div_ready;

endmodule

// Copyright 2018 ETH Zurich and University of Bologna.
// Copyright and related rights are licensed under the Solderpad Hardware
// License, Version 0.51 (the "License"); you may not use this file except in
// compliance with the License.  You may obtain a copy of the License at
// http://solderpad.org/licenses/SHL-0.51. Unless required by applicable law
// or agreed to in writing, software, hardware and materials distributed under
// this License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

////////////////////////////////////////////////////////////////////////////////
// Engineer:       Andreas Traber - atraber@student.ethz.ch                   //
//                                                                            //
// Additional contributions by:                                               //
//                 Davide Schiavone - pschiavo@iis.ee.ethz.ch                 //
//                                                                            //
// Design Name:    cv32e40p_popcnt                                               //
// Project Name:   RI5CY                                                      //
// Language:       SystemVerilog                                              //
//                                                                            //
// Description:    Count the number of '1's in a word                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

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
  endgenerate

  generate
    for (m = 0; m < 8; m++) begin : gen_cnt_l2
      assign cnt_l2[m] = {1'b0, cnt_l1[2*m]} + {1'b0, cnt_l1[2*m+1]};
    end
  endgenerate

  generate
    for (n = 0; n < 4; n++) begin : gen_cnt_l3
      assign cnt_l3[n] = {1'b0, cnt_l2[2*n]} + {1'b0, cnt_l2[2*n+1]};
    end
  endgenerate

  generate
    for (p = 0; p < 2; p++) begin : gen_cnt_l4
      assign cnt_l4[p] = {1'b0, cnt_l3[2*p]} + {1'b0, cnt_l3[2*p+1]};
    end
  endgenerate

  assign result_o = {1'b0, cnt_l4[0]} + {1'b0, cnt_l4[1]};

endmodule

// Copyright 2018 ETH Zurich and University of Bologna.
// Copyright and related rights are licensed under the Solderpad Hardware
// License, Version 0.51 (the "License"); you may not use this file except in
// compliance with the License.  You may obtain a copy of the License at
// http://solderpad.org/licenses/SHL-0.51. Unless required by applicable law
// or agreed to in writing, software, hardware and materials distributed under
// this License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

////////////////////////////////////////////////////////////////////////////////
// Engineer:       Andreas Traber - atraber@student.ethz.ch                   //
//                                                                            //
// Additional contributions by:                                               //
//                 Davide Schiavone - pschiavo@iis.ee.ethz.ch                 //
//                                                                            //
// Design Name:    cv32e40p_ff_one                                               //
// Project Name:   RI5CY                                                      //
// Language:       SystemVerilog                                              //
//                                                                            //
// Description:    Find First One                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

module cv32e40p_ff_one #(
    parameter LEN = 32
) (
    input logic [LEN-1:0] in_i,

    output logic [$clog2(LEN)-1:0] first_one_o,
    output logic                   no_ones_o
);

  localparam NUM_LEVELS = $clog2(LEN);

  logic [          LEN-1:0][NUM_LEVELS-1:0] index_lut;
  logic [2**NUM_LEVELS-1:0]                 sel_nodes;
  logic [2**NUM_LEVELS-1:0][NUM_LEVELS-1:0] index_nodes;


  //////////////////////////////////////////////////////////////////////////////
  // generate tree structure
  //////////////////////////////////////////////////////////////////////////////

  generate
    genvar j;
    for (j = 0; j < LEN; j++) begin : gen_index_lut
      assign index_lut[j] = $unsigned(j);
    end
  endgenerate

  generate
    genvar k;
    genvar l;
    genvar level;

    assign sel_nodes[2**NUM_LEVELS-1] = 1'b0;

    for (level = 0; level < NUM_LEVELS; level++) begin : gen_tree
      //------------------------------------------------------------
      if (level < NUM_LEVELS - 1) begin : gen_non_root_level
        for (l = 0; l < 2 ** level; l++) begin : gen_node
          assign sel_nodes[2**level-1+l]   = sel_nodes[2**(level+1)-1+l*2] | sel_nodes[2**(level+1)-1+l*2+1];
          assign index_nodes[2**level-1+l] = (sel_nodes[2**(level+1)-1+l*2] == 1'b1) ?
                                           index_nodes[2**(level+1)-1+l*2] : index_nodes[2**(level+1)-1+l*2+1];
        end
      end
      //------------------------------------------------------------
      if (level == NUM_LEVELS - 1) begin : gen_root_level
        for (k = 0; k < 2 ** level; k++) begin : gen_node
          // if two successive indices are still in the vector...
          if (k * 2 < LEN - 1) begin : gen_two
            assign sel_nodes[2**level-1+k] = in_i[k*2] | in_i[k*2+1];
            assign index_nodes[2**level-1+k] = (in_i[k*2] == 1'b1) ? index_lut[k*2] : index_lut[k*2+1];
          end
          // if only the first index is still in the vector...
          if (k * 2 == LEN - 1) begin : gen_one
            assign sel_nodes[2**level-1+k]   = in_i[k*2];
            assign index_nodes[2**level-1+k] = index_lut[k*2];
          end
          // if index is out of range
          if (k * 2 > LEN - 1) begin : gen_out_of_range
            assign sel_nodes[2**level-1+k]   = 1'b0;
            assign index_nodes[2**level-1+k] = '0;
          end
        end
      end
      //------------------------------------------------------------
    end
  endgenerate

  //////////////////////////////////////////////////////////////////////////////
  // connect output
  //////////////////////////////////////////////////////////////////////////////

  assign first_one_o = index_nodes[0];
  assign no_ones_o   = ~sel_nodes[0];

endmodule

// Copyright 2018 ETH Zurich and University of Bologna.
// Copyright and related rights are licensed under the Solderpad Hardware
// License, Version 0.51 (the "License"); you may not use this file except in
// compliance with the License.  You may obtain a copy of the License at
// http://solderpad.org/licenses/SHL-0.51. Unless required by applicable law
// or agreed to in writing, software, hardware and materials distributed under
// this License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

///////////////////////////////////////////////////////////////////////////////
// File       : Simple Serial Divider
// Ver        : 1.0
// Date       : 15.03.2016
///////////////////////////////////////////////////////////////////////////////
//
// Description: this is a simple serial divider for signed integers (int32).
//
///////////////////////////////////////////////////////////////////////////////
//
// Authors    : Michael Schaffner (schaffner@iis.ee.ethz.ch)
//              Andreas Traber    (atraber@iis.ee.ethz.ch)
//
///////////////////////////////////////////////////////////////////////////////

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
    input  logic [            1:0] OpCode_SI,  // 0: udiv, 2: urem, 1: div, 3: rem
    // handshake
    input  logic                   InVld_SI,
    // output IF
    input  logic                   OutRdy_SI,
    output logic                   OutVld_SO,
    output logic [    C_WIDTH-1:0] Res_DO
);

  ///////////////////////////////////////////////////////////////////////////////
  // signal declarations
  ///////////////////////////////////////////////////////////////////////////////

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
  logic AddCout_S, AddInv_S;

  enum logic [1:0] {
    IDLE,
    DIVIDE,
    FINISH
  }
      State_SN, State_SP;


  ///////////////////////////////////////////////////////////////////////////////
  // datapath
  ///////////////////////////////////////////////////////////////////////////////

  assign PmSel_S  = LoadEn_S & ~(OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI));

  // muxes
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

  // invert if necessary
  assign Res_DO = (ResInv_SP) ? -$signed(OutMux_D) : OutMux_D;

  // main adder (parallel-prefix), also provides the AReg>=BReg compare used in DIVIDE
  assign AddTmp_D = (LoadEn_S) ? 0 : AReg_DP;
  assign AddInv_S = ~PmSel_S;

  pp_add #(.W(C_WIDTH)) i_div_add (
      .a(AddTmp_D),
      .b(AddInv_S ? ~AddMux_D : AddMux_D),
      .cin(AddInv_S),
      .sum(AddOut_D),
      .cout(AddCout_S)
  );

  // main comparator: in DIVIDE, AddOut_D = AReg-BReg, so
  //   A==B  <=>  ~|AddOut_D,  A>B  <=>  AddCout_S & |AddOut_D
  assign ABComp_S = ((CompInv_SP ? ~AddCout_S : AddCout_S) | ~(|AddOut_D)) & ((|AReg_DP) | OpBIsZero_SI);

  ///////////////////////////////////////////////////////////////////////////////
  // counter
  ///////////////////////////////////////////////////////////////////////////////

  assign Cnt_DN = (LoadEn_S) ? OpBShift_DI : (~CntZero_S) ? Cnt_DP - 1 : Cnt_DP;

  assign CntZero_S = ~(|Cnt_DP);

  ///////////////////////////////////////////////////////////////////////////////
  // FSM
  ///////////////////////////////////////////////////////////////////////////////

  always_comb begin : p_fsm
    // default
    State_SN   = State_SP;

    OutVld_SO  = 1'b0;

    LoadEn_S   = 1'b0;

    ARegEn_S   = 1'b0;
    BRegEn_S   = 1'b0;
    ResRegEn_S = 1'b0;

    case (State_SP)
      /////////////////////////////////
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
      /////////////////////////////////
      DIVIDE: begin

        ARegEn_S   = ABComp_S;
        BRegEn_S   = 1'b1;
        ResRegEn_S = 1'b1;

        // calculation finished
        // one more divide cycle (32nd divide cycle)
        if (CntZero_S) begin
          State_SN = FINISH;
        end
      end
      /////////////////////////////////
      FINISH: begin
        OutVld_SO = 1'b1;

        if (OutRdy_SI) begin
          State_SN = IDLE;
        end
      end
      /////////////////////////////////
      default:  /* default */;
      /////////////////////////////////
    endcase
  end


  ///////////////////////////////////////////////////////////////////////////////
  // regs
  ///////////////////////////////////////////////////////////////////////////////

  // get flags
  assign RemSel_SN = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN = (LoadEn_S) ? (~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] & (OpA_DI[$high(
      OpA_DI
  )] ^ OpBSign_SI) : ResInv_SP;

  assign AReg_DN = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN = (BRegEn_S) ? BMux_D : BReg_DP;
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

  ///////////////////////////////////////////////////////////////////////////////
  // assertions
  ///////////////////////////////////////////////////////////////////////////////

`ifdef CV32E40P_ASSERT_ON
  initial begin : p_assertions
    assert (C_LOG_WIDTH == $clog2(C_WIDTH + 1))
    else $error("C_LOG_WIDTH must be $clog2(C_WIDTH+1)");
  end
`endif

endmodule  // serDiv

////////////////////////////////////////////////////////////////////////////////
// Sklansky parallel-prefix adder.  Written explicitly because technology
// mapping of the inferred "+" operator yields ripple-carry chains here.
//   sum = a + b + cin   (mod 2^W),  cout = carry out
////////////////////////////////////////////////////////////////////////////////

module pp_add #(
    parameter W = 32
) (
    input  logic [W-1:0] a,
    input  logic [W-1:0] b,
    input  logic         cin,
    output logic [W-1:0] sum,
    output logic         cout,
    output logic         g_out,
    output logic         p_out
);

  localparam L = $clog2(W);

  logic [W-1:0] g0, p0, t;
  logic [L:0][W-1:0] G, P;

  assign g0   = a & b;
  assign p0   = a | b;
  assign t    = a ^ b;
  assign G[0] = g0;
  assign P[0] = p0;

  genvar l, i;
  generate
    for (l = 1; l <= L; l++) begin : gen_lvl
      for (i = 0; i < W; i++) begin : gen_bit
        if ((i % (2 << (l - 1))) < (1 << (l - 1))) begin : gen_pass
          assign G[l][i] = G[l-1][i];
          assign P[l][i] = P[l-1][i];
        end else begin : gen_comb
          assign G[l][i] = G[l-1][i] | (P[l-1][i] &
                           G[l-1][i - (i % (2 << (l-1))) + (1 << (l-1)) - 1]);
          assign P[l][i] = P[l-1][i] &
                           P[l-1][i - (i % (2 << (l-1))) + (1 << (l-1)) - 1];
        end
      end
    end
  endgenerate

  assign sum[0] = t[0] ^ cin;

  genvar j;
  generate
    for (j = 1; j < W; j++) begin : gen_sum
      assign sum[j] = t[j] ^ (G[L][j-1] | (P[L][j-1] & cin));
    end
  endgenerate

  assign cout = G[L][W-1] | (P[L][W-1] & cin);
  assign g_out = G[L][W-1];
  assign p_out = P[L][W-1];

endmodule

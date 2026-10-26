// Optimized, self-contained CV32E40P execution-stage ALU.
module cv32e40p_alu (
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
  localparam logic [6:0]
    ALU_LTS=7'h00, ALU_LTU=7'h01, ALU_SLTS=7'h02, ALU_SLTU=7'h03,
    ALU_LES=7'h04, ALU_LEU=7'h05, ALU_SLETS=7'h06, ALU_SLETU=7'h07,
    ALU_GTS=7'h08, ALU_GTU=7'h09, ALU_GES=7'h0a, ALU_GEU=7'h0b,
    ALU_EQ=7'h0c, ALU_NE=7'h0d,
    ALU_MIN=7'h10, ALU_MINU=7'h11, ALU_MAX=7'h12, ALU_MAXU=7'h13,
    ALU_ABS=7'h14, ALU_AND=7'h15, ALU_CLIP=7'h16, ALU_CLIPU=7'h17,
    ALU_ADD=7'h18, ALU_SUB=7'h19, ALU_ADDU=7'h1a, ALU_SUBU=7'h1b,
    ALU_ADDR=7'h1c, ALU_SUBR=7'h1d, ALU_ADDUR=7'h1e, ALU_SUBUR=7'h1f,
    ALU_SRA=7'h24, ALU_SRL=7'h25, ALU_ROR=7'h26, ALU_SLL=7'h27,
    ALU_BEXT=7'h28, ALU_BEXTU=7'h29, ALU_BINS=7'h2a, ALU_BCLR=7'h2b,
    ALU_BSET=7'h2c, ALU_INS=7'h2d, ALU_OR=7'h2e, ALU_XOR=7'h2f,
    ALU_DIVU=7'h30, ALU_DIV=7'h31, ALU_REMU=7'h32, ALU_REM=7'h33,
    ALU_CNT=7'h34, ALU_CLB=7'h35, ALU_FF1=7'h36, ALU_FL1=7'h37,
    ALU_PCKLO=7'h38, ALU_PCKHI=7'h39, ALU_SHUF=7'h3a, ALU_SHUF2=7'h3b,
    ALU_EXTS=7'h3e, ALU_EXT=7'h3f, ALU_BREV=7'h49;

  logic [31:0] operand_b_neg;
  logic [31:0] bmask, bmask_inv;
  assign operand_b_neg = ~operand_b_i;
  assign bmask = (~(32'hffff_fffe << bmask_a_i)) << bmask_b_i;
  assign bmask_inv = ~bmask;

  // -----------------------------------------------------------------------
  // Segmented lane adder.  It is shared by the add/sub family and ABS.
  // -----------------------------------------------------------------------
  logic add_negate, add_abs, add_cin;
  logic [31:0] adder_op_a, adder_op_b, adder_result;
  logic [31:0] adder_round_result;
  logic [4:0] add_shift_amt;

  assign add_negate = operator_i[0] | is_subrot_i;
  assign add_abs = operator_i == ALU_ABS;
  assign add_cin = add_negate | add_abs;

  always_comb begin
    if (add_abs) begin
      adder_op_a = ~operand_a_i;
      adder_op_b = 32'b0;
    end else begin
      adder_op_a = is_subrot_i ? {operand_b_i[15:0], operand_a_i[31:16]} : operand_a_i;
      adder_op_b = add_negate ?
                   (is_subrot_i ? ~{operand_a_i[15:0], operand_b_i[31:16]} : operand_b_neg) :
                   operand_b_i;
    end
  end

  always_comb begin
    case (vector_mode_i)
      2'b10: begin
        adder_result[31:16] = adder_op_a[31:16] + adder_op_b[31:16] + add_cin;
        adder_result[15:0]  = adder_op_a[15:0]  + adder_op_b[15:0]  + add_cin;
      end
      2'b11: begin
        adder_result[31:24] = adder_op_a[31:24] + adder_op_b[31:24] + add_cin;
        adder_result[23:16] = adder_op_a[23:16] + adder_op_b[23:16] + add_cin;
        adder_result[15:8]  = adder_op_a[15:8]  + adder_op_b[15:8]  + add_cin;
        adder_result[7:0]   = adder_op_a[7:0]   + adder_op_b[7:0]   + add_cin;
      end
      default: adder_result = adder_op_a + adder_op_b + add_cin;
    endcase
  end

  assign adder_round_result = adder_result + ((operator_i[2] && operator_i[6:3] == 4'b0011) ? {1'b0,bmask[31:1]} : 32'b0);
  assign add_shift_amt = is_clpx_i ? {3'b0,clpx_shift_i} : bmask_b_i;

  logic [31:0] add_shift_result;
  always_comb begin
    case (vector_mode_i)
      2'b10: begin
        if (!operator_i[1]) begin
          add_shift_result[31:16] = $signed(adder_round_result[31:16]) >>> add_shift_amt[3:0];
          add_shift_result[15:0]  = $signed(adder_round_result[15:0])  >>> add_shift_amt[3:0];
        end else begin
          add_shift_result[31:16] = adder_round_result[31:16] >> add_shift_amt[3:0];
          add_shift_result[15:0]  = adder_round_result[15:0]  >> add_shift_amt[3:0];
        end
      end
      2'b11: begin
        if (!operator_i[1]) begin
          add_shift_result[31:24] = $signed(adder_round_result[31:24]) >>> add_shift_amt[2:0];
          add_shift_result[23:16] = $signed(adder_round_result[23:16]) >>> add_shift_amt[2:0];
          add_shift_result[15:8]  = $signed(adder_round_result[15:8])  >>> add_shift_amt[2:0];
          add_shift_result[7:0]   = $signed(adder_round_result[7:0])   >>> add_shift_amt[2:0];
        end else begin
          add_shift_result[31:24] = adder_round_result[31:24] >> add_shift_amt[2:0];
          add_shift_result[23:16] = adder_round_result[23:16] >> add_shift_amt[2:0];
          add_shift_result[15:8]  = adder_round_result[15:8]  >> add_shift_amt[2:0];
          add_shift_result[7:0]   = adder_round_result[7:0]   >> add_shift_amt[2:0];
        end
      end
      default: begin
        if (!operator_i[1]) add_shift_result = $signed(adder_round_result) >>> add_shift_amt;
        else                add_shift_result = adder_round_result >> add_shift_amt;
      end
    endcase
  end

  // -----------------------------------------------------------------------
  // Ordinary and rotate shifts.  ROR only rotates in 32-bit mode.
  // -----------------------------------------------------------------------
  logic [31:0] shift_result;
  logic [63:0] rotate_wide;
  assign rotate_wide = {operand_a_i,operand_a_i} >> operand_b_i[4:0];
  always_comb begin
    shift_result = 32'b0;
    case (vector_mode_i)
      2'b10: begin
        case (operator_i)
          ALU_SLL: begin
            shift_result[31:16] = operand_a_i[31:16] << operand_b_i[19:16];
            shift_result[15:0]  = operand_a_i[15:0]  << operand_b_i[3:0];
          end
          ALU_SRA: begin
            shift_result[31:16] = $signed(operand_a_i[31:16]) >>> operand_b_i[19:16];
            shift_result[15:0]  = $signed(operand_a_i[15:0])  >>> operand_b_i[3:0];
          end
          default: begin
            shift_result[31:16] = operand_a_i[31:16] >> operand_b_i[19:16];
            shift_result[15:0]  = operand_a_i[15:0]  >> operand_b_i[3:0];
          end
        endcase
      end
      2'b11: begin
        case (operator_i)
          ALU_SLL: begin
            shift_result[31:24] = operand_a_i[31:24] << operand_b_i[26:24];
            shift_result[23:16] = operand_a_i[23:16] << operand_b_i[18:16];
            shift_result[15:8]  = operand_a_i[15:8]  << operand_b_i[10:8];
            shift_result[7:0]   = operand_a_i[7:0]   << operand_b_i[2:0];
          end
          ALU_SRA: begin
            shift_result[31:24] = $signed(operand_a_i[31:24]) >>> operand_b_i[26:24];
            shift_result[23:16] = $signed(operand_a_i[23:16]) >>> operand_b_i[18:16];
            shift_result[15:8]  = $signed(operand_a_i[15:8])  >>> operand_b_i[10:8];
            shift_result[7:0]   = $signed(operand_a_i[7:0])   >>> operand_b_i[2:0];
          end
          default: begin
            shift_result[31:24] = operand_a_i[31:24] >> operand_b_i[26:24];
            shift_result[23:16] = operand_a_i[23:16] >> operand_b_i[18:16];
            shift_result[15:8]  = operand_a_i[15:8]  >> operand_b_i[10:8];
            shift_result[7:0]   = operand_a_i[7:0]   >> operand_b_i[2:0];
          end
        endcase
      end
      default: begin
        case (operator_i)
          ALU_SLL: shift_result = operand_a_i << operand_b_i[4:0];
          ALU_SRA: shift_result = $signed(operand_a_i) >>> operand_b_i[4:0];
          ALU_ROR: shift_result = rotate_wide[31:0];
          default: shift_result = operand_a_i >> operand_b_i[4:0];
        endcase
      end
    endcase
  end

  // -----------------------------------------------------------------------
  // Comparisons, min/max, absolute value and clip.
  // -----------------------------------------------------------------------
  logic [3:0] cmp_signed, is_equal_vec, is_greater_vec;
  logic [3:0] is_equal, is_greater, cmp_result;
  genvar gi;
  generate
    for (gi=0; gi<4; gi++) begin : g_cmp
      assign is_equal_vec[gi] = operand_a_i[8*gi +: 8] == operand_b_i[8*gi +: 8];
      assign is_greater_vec[gi] = $signed({operand_a_i[8*gi+7] & cmp_signed[gi],operand_a_i[8*gi +: 8]}) >
                                  $signed({operand_b_i[8*gi+7] & cmp_signed[gi],operand_b_i[8*gi +: 8]});
    end
  endgenerate
  always_comb begin
    case (operator_i)
      ALU_GTS,ALU_GES,ALU_LTS,ALU_LES,ALU_SLTS,ALU_SLETS,ALU_MIN,ALU_MAX: begin
        case (vector_mode_i)
          2'b11: cmp_signed = 4'b1111;
          2'b10: cmp_signed = 4'b1010;
          default: cmp_signed = 4'b1000;
        endcase
      end
      default: cmp_signed = 4'b0000;
    endcase
  end
  always_comb begin
    is_equal = {4{&is_equal_vec}};
    is_greater = {4{is_greater_vec[3] |
                    (is_equal_vec[3] & (is_greater_vec[2] |
                    (is_equal_vec[2] & (is_greater_vec[1] |
                    (is_equal_vec[1] & is_greater_vec[0])))))}};
    case (vector_mode_i)
      2'b10: begin
        is_equal[1:0] = {2{is_equal_vec[1] & is_equal_vec[0]}};
        is_equal[3:2] = {2{is_equal_vec[3] & is_equal_vec[2]}};
        is_greater[1:0] = {2{is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0])}};
        is_greater[3:2] = {2{is_greater_vec[3] | (is_equal_vec[3] & is_greater_vec[2])}};
      end
      2'b11: begin
        is_equal = is_equal_vec;
        is_greater = is_greater_vec;
      end
      default: ;
    endcase
  end
  always_comb begin
    case (operator_i)
      ALU_NE: cmp_result = ~is_equal;
      ALU_GTS,ALU_GTU: cmp_result = is_greater;
      ALU_GES,ALU_GEU: cmp_result = is_greater | is_equal;
      ALU_LTS,ALU_LTU,ALU_SLTS,ALU_SLTU: cmp_result = ~(is_greater | is_equal);
      ALU_LES,ALU_LEU,ALU_SLETS,ALU_SLETU: cmp_result = ~is_greater;
      default: cmp_result = is_equal;
    endcase
  end
  assign comparison_result_o = cmp_result[3];

  logic do_min;
  logic [3:0] sel_minmax;
  logic [31:0] result_minmax, abs_result, clip_result;
  assign do_min = operator_i == ALU_MIN || operator_i == ALU_MINU;
  assign sel_minmax = is_greater ^ {4{do_min}};
  always_comb begin
    result_minmax[31:24] = sel_minmax[3] ? operand_a_i[31:24] : operand_b_i[31:24];
    result_minmax[23:16] = sel_minmax[2] ? operand_a_i[23:16] : operand_b_i[23:16];
    result_minmax[15:8]  = sel_minmax[1] ? operand_a_i[15:8]  : operand_b_i[15:8];
    result_minmax[7:0]   = sel_minmax[0] ? operand_a_i[7:0]   : operand_b_i[7:0];
  end
  always_comb begin
    case (vector_mode_i)
      2'b10: begin
        if (is_clpx_i) begin
          abs_result[31:16] = adder_result[31:16];
          abs_result[15:0]  = operand_a_i[15:0];
        end else begin
          abs_result[31:16] = operand_a_i[31] ? adder_result[31:16] : operand_a_i[31:16];
          abs_result[15:0]  = operand_a_i[15] ? adder_result[15:0]  : operand_a_i[15:0];
        end
      end
      2'b11: begin
        abs_result[31:24] = operand_a_i[31] ? adder_result[31:24] : operand_a_i[31:24];
        abs_result[23:16] = operand_a_i[23] ? adder_result[23:16] : operand_a_i[23:16];
        abs_result[15:8]  = operand_a_i[15] ? adder_result[15:8]  : operand_a_i[15:8];
        abs_result[7:0]   = operand_a_i[7]  ? adder_result[7:0]   : operand_a_i[7:0];
      end
      default: abs_result = operand_a_i[31] ? adder_result : operand_a_i;
    endcase
  end
  always_comb begin
    if (operator_i == ALU_CLIPU) begin
      if ($signed(operand_a_i) <= 0) clip_result = 32'b0;
      else if ($signed(operand_a_i) > $signed(operand_b_i)) clip_result = operand_b_i;
      else clip_result = operand_a_i;
    end else begin
      if ($signed(operand_a_i) <= $signed(operand_b_neg)) clip_result = operand_b_neg;
      else if ($signed(operand_a_i) > $signed(operand_b_i)) clip_result = operand_b_i;
      else clip_result = operand_a_i;
    end
  end

  // -----------------------------------------------------------------------
  // Bit manipulation and bit reverse.  Left shifts reuse the right shifter
  // through bit reversal, as in the upstream datapath.
  // -----------------------------------------------------------------------
  logic bit_shift_left, bit_shift_arith;
  logic [31:0] bit_shift_input, bit_shift_result, bit_reverse_input;
  logic [32:0] bit_shift_wide;
  logic [31:0] bextins_result, reverse_result;
  logic extract_sign;
  genvar gr;
  generate
    for (gr=0; gr<32; gr++) begin : g_brev_in
      assign bit_reverse_input[gr] = bit_shift_left ? operand_a_i[31-gr] : operand_a_i[gr];
    end
  endgenerate
  assign bit_shift_left = operator_i == ALU_BINS || operator_i == ALU_BREV;
  assign bit_shift_arith = operator_i == ALU_BEXT;
  assign bit_shift_input = bit_reverse_input;
  assign bit_shift_wide = $signed({bit_shift_arith & bit_shift_input[31],bit_shift_input}) >>> operand_b_i[4:0];
  generate
    for (gr=0; gr<32; gr++) begin : g_brev_out
      assign bit_shift_result[gr] = bit_shift_left ? bit_shift_wide[31-gr] : bit_shift_wide[gr];
    end
  endgenerate
  assign extract_sign = (operator_i == ALU_BEXT) & bit_shift_result[bmask_a_i];
  assign bextins_result = (bmask & bit_shift_result) |
                          ((operator_i == ALU_BINS ? operand_c_i : {32{extract_sign}}) & bmask_inv);
  integer br;
  always_comb begin
    reverse_result = 32'b0;
    case (bmask_a_i[1:0])
      2'b01: begin
        for (br=0;br<16;br++) reverse_result[2*br +: 2] = bit_shift_result[31-2*br -: 2];
      end
      2'b10: begin
        for (br=0;br<10;br++) reverse_result[3*br +: 3] = bit_shift_result[31-3*br -: 3];
      end
      default: begin
        for (br=0;br<32;br++) reverse_result[br] = bit_shift_result[31-br];
      end
    endcase
  end

  // -----------------------------------------------------------------------
  // Byte/halfword extract, insert, shuffle and pack.
  // -----------------------------------------------------------------------
  function automatic logic [7:0] pick_byte(input logic [31:0] x,input logic [1:0] s);
    case (s)
      2'd0: pick_byte=x[7:0]; 2'd1: pick_byte=x[15:8];
      2'd2: pick_byte=x[23:16]; default: pick_byte=x[31:24];
    endcase
  endfunction
  function automatic logic [15:0] pick_half(input logic [31:0] x,input logic s);
    pick_half = s ? x[31:16] : x[15:0];
  endfunction
  logic [31:0] shuffle_result;
  logic [7:0] ext_byte;
  logic [15:0] ext_half;
  always_comb begin
    shuffle_result = operand_c_i;
    ext_byte = pick_byte(operand_a_i,imm_vec_ext_i);
    ext_half = pick_half(operand_a_i,imm_vec_ext_i[0]);
    case (operator_i)
      ALU_EXT,ALU_EXTS: begin
        if (vector_mode_i == 2'b11)
          shuffle_result = (operator_i == ALU_EXTS) ? {{24{ext_byte[7]}},ext_byte} : {24'b0,ext_byte};
        else
          shuffle_result = (operator_i == ALU_EXTS) ? {{16{ext_half[15]}},ext_half} : {16'b0,ext_half};
      end
      ALU_INS: begin
        if (vector_mode_i == 2'b11) begin
          case (imm_vec_ext_i)
            2'd0: shuffle_result[7:0]=operand_a_i[7:0];
            2'd1: shuffle_result[15:8]=operand_a_i[7:0];
            2'd2: shuffle_result[23:16]=operand_a_i[7:0];
            default: shuffle_result[31:24]=operand_a_i[7:0];
          endcase
        end else if (imm_vec_ext_i[0]) shuffle_result[31:16]=operand_a_i[15:0];
        else shuffle_result[15:0]=operand_a_i[15:0];
      end
      ALU_SHUF,ALU_SHUF2: begin
        if (vector_mode_i == 2'b11) begin
          for (int sk=0;sk<4;sk++)
            shuffle_result[8*sk +: 8] = pick_byte((operator_i == ALU_SHUF2 && !operand_b_i[8*sk+2]) ? operand_c_i : operand_a_i,
                                                   operand_b_i[8*sk +: 2]);
        end else begin
          for (int sk=0;sk<2;sk++)
            shuffle_result[16*sk +: 16] = pick_half((operator_i == ALU_SHUF2 && !operand_b_i[16*sk+1]) ? operand_c_i : operand_a_i,
                                                     operand_b_i[16*sk]);
        end
      end
      ALU_PCKLO: begin
        if (vector_mode_i == 2'b11) shuffle_result={operand_c_i[31:16],operand_a_i[7:0],operand_b_i[7:0]};
        else shuffle_result={operand_a_i[15:0],operand_b_i[15:0]};
      end
      ALU_PCKHI: begin
        if (vector_mode_i == 2'b11) shuffle_result={operand_a_i[7:0],operand_b_i[7:0],operand_c_i[15:0]};
        else shuffle_result={operand_a_i[31:16],operand_b_i[31:16]};
      end
      default: ;
    endcase
  end

  // -----------------------------------------------------------------------
  // Bit counting and divider normalization.  The divisor is operand A and
  // the dividend is operand B, matching the CV32E40P interface convention.
  // -----------------------------------------------------------------------
  logic [31:0] operand_a_rev, operand_a_neg_rev, ff_input;
  logic [5:0] cnt_result, clb_result, bitop_result, div_shift;
  logic [4:0] ff1_result;
  logic ff_no_one, div_op_a_signed, div_valid, div_ready;
  logic [31:0] div_operand_b, result_div;
  generate
    for (gr=0; gr<32; gr++) begin : g_arev
      assign operand_a_rev[gr]=operand_a_i[31-gr];
      assign operand_a_neg_rev[gr]=~operand_a_i[31-gr];
    end
  endgenerate
  cv32e40p_popcnt u_popcnt(.in_i(operand_a_i),.result_o(cnt_result));
  always_comb begin
    case (operator_i)
      ALU_FF1: ff_input=operand_a_i;
      ALU_DIVU,ALU_REMU,ALU_FL1: ff_input=operand_a_rev;
      ALU_DIV,ALU_REM,ALU_CLB: ff_input=operand_a_i[31] ? operand_a_neg_rev : operand_a_rev;
      default: ff_input=32'b0;
    endcase
  end
  cv32e40p_ff_one u_ff1(.in_i(ff_input),.first_one_o(ff1_result),.no_ones_o(ff_no_one));
  assign clb_result = ff1_result - 1'b1;
  always_comb begin
    case (operator_i)
      ALU_FF1: bitop_result=ff_no_one ? 6'd32 : {1'b0,ff1_result};
      ALU_FL1: bitop_result=ff_no_one ? 6'd32 : {1'b0,5'd31-ff1_result};
      ALU_CNT: bitop_result=cnt_result;
      ALU_CLB: begin
        if (ff_no_one) bitop_result=operand_a_i[31] ? 6'd31 : 6'd0;
        else bitop_result=clb_result;
      end
      default: bitop_result=6'b0;
    endcase
  end
  assign div_op_a_signed = operand_a_i[31] & operator_i[0];
  assign div_shift = (ff_no_one ? 6'd31 : clb_result) + (div_op_a_signed ? 6'd0 : 6'd1);
  assign div_operand_b = operand_a_i << div_shift;
  assign div_valid = enable_i && operator_i[6:2] == 5'b01100;
  cv32e40p_alu_div u_div (
    .Clk_CI(clk),.Rst_RBI(rst_n),.OpA_DI(operand_b_i),.OpB_DI(div_operand_b),
    .OpBShift_DI(div_shift),.OpBIsZero_SI(operand_a_i==0),.OpBSign_SI(div_op_a_signed),
    .OpCode_SI(operator_i[1:0]),.Res_DO(result_div),.InVld_SI(div_valid),
    .OutRdy_SI(ex_ready_i),.OutVld_SO(div_ready)
  );
  assign ready_o=div_ready;

  // -----------------------------------------------------------------------
  // Final operation selection.
  // -----------------------------------------------------------------------
  always_comb begin
    result_o=32'b0;
    case (operator_i)
      ALU_AND: result_o=operand_a_i & operand_b_i;
      ALU_OR:  result_o=operand_a_i | operand_b_i;
      ALU_XOR: result_o=operand_a_i ^ operand_b_i;
      ALU_ADD,ALU_SUB,ALU_ADDU,ALU_SUBU,ALU_ADDR,ALU_SUBR,ALU_ADDUR,ALU_SUBUR:
        result_o=add_shift_result;
      ALU_SLL,ALU_SRL,ALU_SRA,ALU_ROR: result_o=shift_result;
      ALU_BINS,ALU_BEXT,ALU_BEXTU: result_o=bextins_result;
      ALU_BCLR: result_o=operand_a_i & bmask_inv;
      ALU_BSET: result_o=operand_a_i | bmask;
      ALU_BREV: result_o=reverse_result;
      ALU_SHUF,ALU_SHUF2,ALU_PCKLO,ALU_PCKHI,ALU_EXT,ALU_EXTS,ALU_INS: result_o=shuffle_result;
      ALU_MIN,ALU_MINU,ALU_MAX,ALU_MAXU: result_o=result_minmax;
      ALU_ABS: result_o=abs_result;
      ALU_CLIP,ALU_CLIPU: result_o=clip_result;
      ALU_EQ,ALU_NE,ALU_GTU,ALU_GEU,ALU_LTU,ALU_LEU,ALU_GTS,ALU_GES,ALU_LTS,ALU_LES:
        result_o={{8{cmp_result[3]}},{8{cmp_result[2]}},{8{cmp_result[1]}},{8{cmp_result[0]}}};
      ALU_SLTS,ALU_SLTU,ALU_SLETS,ALU_SLETU: result_o={31'b0,comparison_result_o};
      ALU_FF1,ALU_FL1,ALU_CLB,ALU_CNT: result_o={26'b0,bitop_result};
      ALU_DIV,ALU_DIVU,ALU_REM,ALU_REMU: result_o=result_div;
      default: ;
    endcase
  end
endmodule

// Balanced population count tree.
module cv32e40p_popcnt(input logic [31:0] in_i,output logic [5:0] result_o);
  logic [15:0][1:0] c1;
  logic [7:0][2:0] c2;
  logic [3:0][3:0] c3;
  logic [1:0][4:0] c4;
  genvar i;
  generate
    for(i=0;i<16;i++) begin:g1 assign c1[i]={1'b0,in_i[2*i]}+{1'b0,in_i[2*i+1]}; end
    for(i=0;i<8;i++) begin:g2 assign c2[i]={1'b0,c1[2*i]}+{1'b0,c1[2*i+1]}; end
    for(i=0;i<4;i++) begin:g3 assign c3[i]={1'b0,c2[2*i]}+{1'b0,c2[2*i+1]}; end
    for(i=0;i<2;i++) begin:g4 assign c4[i]={1'b0,c3[2*i]}+{1'b0,c3[2*i+1]}; end
  endgenerate
  assign result_o={1'b0,c4[0]}+{1'b0,c4[1]};
endmodule

// Balanced find-first-one tree (least-significant set bit).
module cv32e40p_ff_one #(parameter LEN=32)(
  input logic [LEN-1:0] in_i,
  output logic [$clog2(LEN)-1:0] first_one_o,
  output logic no_ones_o
);
  localparam NUM_LEVELS=$clog2(LEN);
  logic [LEN-1:0][NUM_LEVELS-1:0] index_lut;
  logic [2**NUM_LEVELS-1:0] sel_nodes;
  logic [2**NUM_LEVELS-1:0][NUM_LEVELS-1:0] index_nodes;
  genvar j,k,l,level;
  generate
    for(j=0;j<LEN;j++) begin:g_lut assign index_lut[j]=$unsigned(j); end
    assign sel_nodes[2**NUM_LEVELS-1]=1'b0;
    for(level=0;level<NUM_LEVELS;level++) begin:g_tree
      if(level<NUM_LEVELS-1) begin:g_nr
        for(l=0;l<2**level;l++) begin:g_node
          assign sel_nodes[2**level-1+l]=sel_nodes[2**(level+1)-1+l*2]|sel_nodes[2**(level+1)-1+l*2+1];
          assign index_nodes[2**level-1+l]=sel_nodes[2**(level+1)-1+l*2] ? index_nodes[2**(level+1)-1+l*2] : index_nodes[2**(level+1)-1+l*2+1];
        end
      end
      if(level==NUM_LEVELS-1) begin:g_root
        for(k=0;k<2**level;k++) begin:g_node
          if(k*2<LEN-1) begin:g_two
            assign sel_nodes[2**level-1+k]=in_i[k*2]|in_i[k*2+1];
            assign index_nodes[2**level-1+k]=in_i[k*2] ? index_lut[k*2] : index_lut[k*2+1];
          end
          if(k*2==LEN-1) begin:g_one
            assign sel_nodes[2**level-1+k]=in_i[k*2];
            assign index_nodes[2**level-1+k]=index_lut[k*2];
          end
          if(k*2>LEN-1) begin:g_out
            assign sel_nodes[2**level-1+k]=1'b0;
            assign index_nodes[2**level-1+k]='0;
          end
        end
      end
    end
  endgenerate
  assign first_one_o=index_nodes[0];
  assign no_ones_o=~sel_nodes[0];
endmodule

// Cycle-identical CV32E40P serial divider.  Its FSM, counter and recurrence
// are deliberately kept in its original form.
module cv32e40p_alu_div #(parameter C_WIDTH=32,parameter C_LOG_WIDTH=6)(
  input logic Clk_CI,input logic Rst_RBI,
  input logic [C_WIDTH-1:0] OpA_DI,input logic [C_WIDTH-1:0] OpB_DI,
  input logic [C_LOG_WIDTH-1:0] OpBShift_DI,input logic OpBIsZero_SI,
  input logic OpBSign_SI,input logic [1:0] OpCode_SI,input logic InVld_SI,
  input logic OutRdy_SI,output logic OutVld_SO,output logic [C_WIDTH-1:0] Res_DO
);
  logic [C_WIDTH-1:0] ResReg_DP,ResReg_DN,ResReg_DP_rev,AReg_DP,AReg_DN,BReg_DP,BReg_DN;
  logic RemSel_SN,RemSel_SP,CompInv_SN,CompInv_SP,ResInv_SN,ResInv_SP;
  logic [C_WIDTH-1:0] AddMux_D,AddOut_D,AddTmp_D,BMux_D,OutMux_D;
  logic [C_LOG_WIDTH-1:0] Cnt_DP,Cnt_DN;
  logic CntZero_S,ARegEn_S,BRegEn_S,ResRegEn_S,ABComp_S,PmSel_S,LoadEn_S;
  enum logic [1:0] {IDLE,DIVIDE,FINISH} State_SN,State_SP;
  assign PmSel_S=LoadEn_S & ~(OpCode_SI[0] & (OpA_DI[$high(OpA_DI)]^OpBSign_SI));
  assign AddMux_D=LoadEn_S ? OpA_DI : BReg_DP;
  assign BMux_D=LoadEn_S ? OpB_DI : {CompInv_SP,BReg_DP[$high(BReg_DP):1]};
  genvar index;
  generate for(index=0;index<C_WIDTH;index++) begin:g_rev assign ResReg_DP_rev[index]=ResReg_DP[C_WIDTH-1-index]; end endgenerate
  assign OutMux_D=RemSel_SP ? AReg_DP : ResReg_DP_rev;
  assign Res_DO=ResInv_SP ? -$signed(OutMux_D) : OutMux_D;
  assign ABComp_S=((AReg_DP==BReg_DP)|((AReg_DP>BReg_DP)^CompInv_SP))&((|AReg_DP)|OpBIsZero_SI);
  assign AddTmp_D=LoadEn_S ? 0 : AReg_DP;
  assign AddOut_D=PmSel_S ? AddTmp_D+AddMux_D : AddTmp_D-$signed(AddMux_D);
  assign Cnt_DN=LoadEn_S ? OpBShift_DI : (~CntZero_S) ? Cnt_DP-1'b1 : Cnt_DP;
  assign CntZero_S=~(|Cnt_DP);
  always_comb begin
    State_SN=State_SP; OutVld_SO=1'b0; LoadEn_S=1'b0;
    ARegEn_S=1'b0; BRegEn_S=1'b0; ResRegEn_S=1'b0;
    case(State_SP)
      IDLE: begin
        OutVld_SO=1'b1;
        if(InVld_SI) begin OutVld_SO=1'b0; ARegEn_S=1'b1; BRegEn_S=1'b1; LoadEn_S=1'b1; State_SN=DIVIDE; end
      end
      DIVIDE: begin
        ARegEn_S=ABComp_S; BRegEn_S=1'b1; ResRegEn_S=1'b1;
        if(CntZero_S) State_SN=FINISH;
      end
      FINISH: begin OutVld_SO=1'b1; if(OutRdy_SI) State_SN=IDLE; end
      default: ;
    endcase
  end
  assign RemSel_SN=LoadEn_S ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN=LoadEn_S ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN=LoadEn_S ? ((~OpBIsZero_SI)|OpCode_SI[1])&OpCode_SI[0]&(OpA_DI[$high(OpA_DI)]^OpBSign_SI) : ResInv_SP;
  assign AReg_DN=ARegEn_S ? AddOut_D : AReg_DP;
  assign BReg_DN=BRegEn_S ? BMux_D : BReg_DP;
  assign ResReg_DN=LoadEn_S ? '0 : ResRegEn_S ? {ABComp_S,ResReg_DP[$high(ResReg_DP):1]} : ResReg_DP;
  always_ff @(posedge Clk_CI or negedge Rst_RBI) begin
    if(~Rst_RBI) begin
      State_SP<=IDLE; AReg_DP<='0; BReg_DP<='0; ResReg_DP<='0; Cnt_DP<='0;
      RemSel_SP<=1'b0; CompInv_SP<=1'b0; ResInv_SP<=1'b0;
    end else begin
      State_SP<=State_SN; AReg_DP<=AReg_DN; BReg_DP<=BReg_DN; ResReg_DP<=ResReg_DN; Cnt_DP<=Cnt_DN;
      RemSel_SP<=RemSel_SN; CompInv_SP<=CompInv_SN; ResInv_SP<=ResInv_SN;
    end
  end
endmodule

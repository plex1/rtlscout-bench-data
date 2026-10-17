module cv32e40p_alu
(
    input logic               clk,
    input logic               rst_n,
    input logic               enable_i,
    input logic [ 6:0] operator_i,
    input logic [31:0] operand_a_i,
    input logic [31:0] operand_b_i,
    input logic [31:0] operand_c_i,

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

  localparam logic [6:0] ALU_ADD  = 7'b0011000;
  localparam logic [6:0] ALU_SUB  = 7'b0011001;
  localparam logic [6:0] ALU_ADDU = 7'b0011010;
  localparam logic [6:0] ALU_SUBU = 7'b0011011;
  localparam logic [6:0] ALU_ADDR = 7'b0011100;
  localparam logic [6:0] ALU_SUBR = 7'b0011101;
  localparam logic [6:0] ALU_ADDUR= 7'b0011110;
  localparam logic [6:0] ALU_SUBUR= 7'b0011111;
  localparam logic [6:0] ALU_XOR  = 7'b0101111;
  localparam logic [6:0] ALU_OR   = 7'b0101110;
  localparam logic [6:0] ALU_AND  = 7'b0010101;
  localparam logic [6:0] ALU_SRA  = 7'b0100100;
  localparam logic [6:0] ALU_SRL  = 7'b0100101;
  localparam logic [6:0] ALU_ROR  = 7'b0100110;
  localparam logic [6:0] ALU_SLL  = 7'b0100111;
  localparam logic [6:0] ALU_BEXT = 7'b0101000;
  localparam logic [6:0] ALU_BEXTU= 7'b0101001;
  localparam logic [6:0] ALU_BINS = 7'b0101010;
  localparam logic [6:0] ALU_BCLR = 7'b0101011;
  localparam logic [6:0] ALU_BSET = 7'b0101100;
  localparam logic [6:0] ALU_BREV = 7'b1001001;
  localparam logic [6:0] ALU_FF1  = 7'b0110110;
  localparam logic [6:0] ALU_FL1  = 7'b0110111;
  localparam logic [6:0] ALU_CNT  = 7'b0110100;
  localparam logic [6:0] ALU_CLB  = 7'b0110101;
  localparam logic [6:0] ALU_EXTS = 7'b0111110;
  localparam logic [6:0] ALU_EXT  = 7'b0111111;
  localparam logic [6:0] ALU_LTS  = 7'b0000000;
  localparam logic [6:0] ALU_LTU  = 7'b0000001;
  localparam logic [6:0] ALU_LES  = 7'b0000100;
  localparam logic [6:0] ALU_LEU  = 7'b0000101;
  localparam logic [6:0] ALU_GTS  = 7'b0001000;
  localparam logic [6:0] ALU_GTU  = 7'b0001001;
  localparam logic [6:0] ALU_GES  = 7'b0001010;
  localparam logic [6:0] ALU_GEU  = 7'b0001011;
  localparam logic [6:0] ALU_EQ   = 7'b0001100;
  localparam logic [6:0] ALU_NE   = 7'b0001101;
  localparam logic [6:0] ALU_SLTS = 7'b0000010;
  localparam logic [6:0] ALU_SLTU = 7'b0000011;
  localparam logic [6:0] ALU_SLETS= 7'b0000110;
  localparam logic [6:0] ALU_SLETU= 7'b0000111;
  localparam logic [6:0] ALU_ABS  = 7'b0010100;
  localparam logic [6:0] ALU_CLIP = 7'b0010110;
  localparam logic [6:0] ALU_CLIPU= 7'b0010111;
  localparam logic [6:0] ALU_INS  = 7'b0101101;
  localparam logic [6:0] ALU_MIN  = 7'b0010000;
  localparam logic [6:0] ALU_MINU = 7'b0010001;
  localparam logic [6:0] ALU_MAX  = 7'b0010010;
  localparam logic [6:0] ALU_MAXU = 7'b0010011;
  localparam logic [6:0] ALU_DIVU = 7'b0110000;
  localparam logic [6:0] ALU_DIV  = 7'b0110001;
  localparam logic [6:0] ALU_REMU = 7'b0110010;
  localparam logic [6:0] ALU_REM  = 7'b0110011;
  localparam logic [6:0] ALU_SHUF = 7'b0111010;
  localparam logic [6:0] ALU_SHUF2= 7'b0111011;
  localparam logic [6:0] ALU_PCKLO= 7'b0111000;
  localparam logic [6:0] ALU_PCKHI= 7'b0111001;

  localparam logic [1:0] VEC_MODE32 = 2'b00;
  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8  = 2'b11;

  logic [31:0] bmask;
  assign bmask = (~(32'hFFFFFFFE << bmask_a_i)) << bmask_b_i;

  function automatic [4:0] ff1_lsb(input [31:0] x);
    for (int j = 0; j < 32; j = j + 1)
      if (x[j]) return j[4:0];
    return 5'd0;
  endfunction

  function automatic [4:0] ff1_msb(input [31:0] x);
    for (int j = 31; j >= 0; j = j - 1)
      if (x[j]) return j[4:0];
    return 5'd0;
  endfunction

  function automatic [31:0] byte_mux(input [31:0] w, input [1:0] s);
    case (s)
      2'd0:    byte_mux = w[7:0];
      2'd1:    byte_mux = w[15:8];
      2'd2:    byte_mux = w[23:16];
      default: byte_mux = w[31:24];
    endcase
  endfunction

  function automatic [15:0] half_mux(input [31:0] w, input s);
    half_mux = s ? w[31:16] : w[15:0];
  endfunction

  function automatic [31:0] lane_addsub(input [31:0] a, input [31:0] b,
                                        input [1:0] vm, input sub);
    logic [31:0] r;
    case (vm)
      VEC_MODE16: begin
        if (sub) r = {a[31:16] - b[31:16], a[15:0] - b[15:0]};
        else     r = {a[31:16] + b[31:16], a[15:0] + b[15:0]};
      end
      VEC_MODE8: begin
        if (sub) r = {a[31:24] - b[31:24], a[23:16] - b[23:16],
                      a[15:8]  - b[15:8],  a[7:0]   - b[7:0]};
        else     r = {a[31:24] + b[31:24], a[23:16] + b[23:16],
                      a[15:8]  + b[15:8],  a[7:0]   + b[7:0]};
      end
      default: begin
        if (sub) r = a - b;
        else     r = a + b;
      end
    endcase
    lane_addsub = r;
  endfunction

  function automatic [31:0] lane_shr(input [31:0] x, input [1:0] vm,
                                     input [31:0] amounts, input arith);
    case (vm)
      VEC_MODE16: begin
        if (arith)
          lane_shr = {$signed(x[31:16]) >>> amounts[19:16],
                      $signed(x[15:0])  >>> amounts[3:0]};
        else
          lane_shr = {x[31:16] >> amounts[19:16],
                      x[15:0]  >> amounts[3:0]};
      end
      VEC_MODE8: begin
        if (arith)
          lane_shr = {$signed(x[31:24]) >>> amounts[26:24],
                      $signed(x[23:16]) >>> amounts[18:16],
                      $signed(x[15:8])  >>> amounts[10:8],
                      $signed(x[7:0])   >>> amounts[2:0]};
        else
          lane_shr = {x[31:24] >> amounts[26:24],
                      x[23:16] >> amounts[18:16],
                      x[15:8]  >> amounts[10:8],
                      x[7:0]   >> amounts[2:0]};
      end
      default: begin
        if (arith) lane_shr = $signed(x) >>> amounts[4:0];
        else       lane_shr = x >> amounts[4:0];
      end
    endcase
  endfunction

  function automatic [31:0] lane_shl(input [31:0] x, input [1:0] vm,
                                     input [31:0] amounts);
    case (vm)
      VEC_MODE16: lane_shl = {x[31:16] << amounts[19:16], x[15:0] << amounts[3:0]};
      VEC_MODE8:  lane_shl = {x[31:24] << amounts[26:24],
                              x[23:16] << amounts[18:16],
                              x[15:8]  << amounts[10:8],
                              x[7:0]   << amounts[2:0]};
      default:     lane_shl = x << amounts[4:0];
    endcase
  endfunction

  function automatic [31:0] ror32(input [31:0] x, input [4:0] sh);
    ror32 = (sh == 5'd0) ? x : ((x >> sh) | (x << (6'd32 - {1'b0, sh})));
  endfunction

  logic [3:0] cmp_signed;
  always_comb begin
    cmp_signed = 4'b0000;
    unique case (operator_i)
      ALU_GTS, ALU_GES, ALU_LTS, ALU_LES,
      ALU_SLTS, ALU_SLETS, ALU_MIN, ALU_MAX,
      ALU_ABS, ALU_CLIP, ALU_CLIPU: begin
        case (vector_mode_i)
          VEC_MODE8:  cmp_signed = 4'b1111;
          VEC_MODE16: cmp_signed = 4'b1010;
          default:    cmp_signed = 4'b1000;
        endcase
      end
      default: ;
    endcase
  end

  logic [3:0] beq;
  logic [3:0] bgt;
  genvar gi;
  generate
    for (gi = 0; gi < 4; gi = gi + 1) begin : gen_bytecmp
      assign beq[gi] = (operand_a_i[8*gi +: 8] == operand_b_i[8*gi +: 8]);
      assign bgt[gi] = $signed({operand_a_i[8*gi+7] & cmp_signed[gi],
                                operand_a_i[8*gi +: 8]}) >
                       $signed({operand_b_i[8*gi+7] & cmp_signed[gi],
                                operand_b_i[8*gi +: 8]});
    end
  endgenerate

  logic [3:0] is_equal;
  logic [3:0] is_greater;
  always_comb begin
    is_equal   = {4{beq[3] & beq[2] & beq[1] & beq[0]}};
    is_greater = {4{bgt[3] | (beq[3] & (bgt[2] |
                  (beq[2] & (bgt[1] | (beq[1] & bgt[0])))))}};

    case (vector_mode_i)
      VEC_MODE16: begin
        is_equal[1:0]   = {2{beq[0] & beq[1]}};
        is_equal[3:2]   = {2{beq[2] & beq[3]}};
        is_greater[1:0] = {2{bgt[1] | (beq[1] & bgt[0])}};
        is_greater[3:2] = {2{bgt[3] | (beq[3] & bgt[2])}};
      end
      VEC_MODE8: begin
        is_equal   = beq;
        is_greater = bgt;
      end
      default: ;
    endcase
  end

  logic [3:0] cmp_result;
  always_comb begin
    cmp_result = is_equal;
    unique case (operator_i)
      ALU_EQ:                                  cmp_result = is_equal;
      ALU_NE:                                  cmp_result = ~is_equal;
      ALU_GTS, ALU_GTU:                        cmp_result = is_greater;
      ALU_GES, ALU_GEU:                        cmp_result = is_greater | is_equal;
      ALU_LTS, ALU_SLTS, ALU_LTU, ALU_SLTU:    cmp_result = ~(is_greater | is_equal);
      ALU_SLETS, ALU_SLETU, ALU_LES, ALU_LEU:  cmp_result = ~is_greater;
      default: ;
    endcase
  end

  assign comparison_result_o = cmp_result[3];

  logic        adder_sub;
  logic [31:0] adder_a;
  logic [31:0] adder_b;
  logic [31:0] adder_raw;
  logic [31:0] adder_rounded;
  logic [31:0] adder_shift_amounts;
  logic        adder_shift_arith;
  logic        adder_shift_use;

  assign adder_sub = (operator_i == ALU_SUB)  || (operator_i == ALU_SUBR) ||
                     (operator_i == ALU_SUBU) || (operator_i == ALU_SUBUR) ||
                     is_subrot_i;

  always_comb begin
    adder_a = operand_a_i;
    adder_b = operand_b_i;
    if (is_subrot_i) begin
      adder_a = {operand_b_i[15:0], operand_a_i[31:16]};
      adder_b = {operand_a_i[15:0], operand_b_i[31:16]};
    end
  end

  always_comb begin
    adder_raw = '0;
    if (operator_i == ALU_ABS)
      adder_raw = lane_addsub('0, operand_a_i, vector_mode_i, 1'b1);
    else
      adder_raw = lane_addsub(adder_a, adder_b, vector_mode_i, adder_sub);
  end

  always_comb begin
    if ((operator_i == ALU_ADDR) || (operator_i == ALU_SUBR) ||
        (operator_i == ALU_ADDUR) || (operator_i == ALU_SUBUR))
      adder_rounded = adder_raw + (bmask >> 1);
    else
      adder_rounded = adder_raw;
  end

  assign adder_shift_use = (operator_i == ALU_ADD)  || (operator_i == ALU_SUB)  ||
                           (operator_i == ALU_ADDR) || (operator_i == ALU_SUBR) ||
                           (operator_i == ALU_ADDU) || (operator_i == ALU_SUBU) ||
                           (operator_i == ALU_ADDUR)|| (operator_i == ALU_SUBUR);
  assign adder_shift_arith = (operator_i == ALU_ADD) || (operator_i == ALU_SUB) ||
                             (operator_i == ALU_ADDR)|| (operator_i == ALU_SUBR);

  always_comb begin
    if (adder_shift_use && is_clpx_i)
      adder_shift_amounts = {{14'b0, clpx_shift_i}, {14'b0, clpx_shift_i}};
    else
      adder_shift_amounts = {4{3'b000, bmask_b_i}};
  end

  logic [31:0] minmax_b;
  logic [3:0]  minmax_sel;
  logic        do_min;
  logic [31:0] minmax_result;

  always_comb begin
    minmax_b = operand_b_i;
    if (operator_i == ALU_ABS)
      minmax_b = lane_addsub('0, operand_a_i, vector_mode_i, 1'b1);

    do_min = (operator_i == ALU_MIN)  || (operator_i == ALU_MINU) ||
             (operator_i == ALU_CLIP) || (operator_i == ALU_CLIPU);

    minmax_sel = is_greater ^ {4{do_min}};

    minmax_result[31:24] = minmax_sel[3] ? operand_a_i[31:24] : minmax_b[31:24];
    minmax_result[23:16] = minmax_sel[2] ? operand_a_i[23:16] : minmax_b[23:16];
    minmax_result[15:8]  = minmax_sel[1] ? operand_a_i[15:8]  : minmax_b[15:8];
    minmax_result[7:0]   = minmax_sel[0] ? operand_a_i[7:0]   : minmax_b[7:0];
  end

  logic [31:0] clip_result;
  always_comb begin
    clip_result = minmax_result;
    if (operator_i == ALU_CLIPU) begin
      if ($signed(operand_a_i) <= 32'sd0)
        clip_result = '0;
    end else if (operator_i == ALU_CLIP) begin
      if ($signed(operand_a_i) <= $signed(~operand_b_i))
        clip_result = ~operand_b_i;
    end
  end

  logic [5:0] cnt_result;
  assign cnt_result = $countones(operand_a_i);

  logic [5:0] bitcount_result;
  always_comb begin
    bitcount_result = '0;
    unique case (operator_i)
      ALU_CNT: bitcount_result = cnt_result;
      ALU_FF1: bitcount_result = (operand_a_i == 32'h0) ? 6'd32 : {1'b0, ff1_lsb(operand_a_i)};
      ALU_FL1: bitcount_result = (operand_a_i == 32'h0) ? 6'd32 : {1'b0, ff1_msb(operand_a_i)};
      ALU_CLB: begin
        if (operand_a_i == 32'h0)
          bitcount_result = 6'd0;
        else if (operand_a_i[31])
          bitcount_result = ((~operand_a_i) == 32'h0) ? 6'd31 :
                            6'd30 - {1'b0, ff1_lsb(~operand_a_i)};
        else
          bitcount_result = 6'd30 - {1'b0, ff1_msb(operand_a_i)};
      end
      default: ;
    endcase
  end

  logic [31:0] reverse_result;
  always_comb begin
    reverse_result = '0;
    if (operator_i == ALU_BREV) begin
      logic [31:0] t;
      t = operand_a_i << operand_b_i[4:0];
      case (bmask_a_i[1:0])
        2'b00, 2'b11: begin
          for (int j = 0; j < 32; j = j + 1) reverse_result[j] = t[31-j];
        end
        2'b01: begin
          for (int j = 0; j < 16; j = j + 1)
            reverse_result[2*j +: 2] = t[31-2*j -: 2];
        end
        2'b10: begin
          for (int j = 0; j < 10; j = j + 1)
            reverse_result[3*j +: 3] = t[31-3*j -: 3];
          reverse_result[31:30] = 2'b00;
        end
        default: ;
      endcase
    end
  end

  logic [31:0] shuf_result;
  always_comb begin
    shuf_result = '0;
    unique case (operator_i)
      ALU_EXT, ALU_EXTS: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_result = {24'b0, operand_a_i[8*imm_vec_ext_i[1:0] +: 8]};
          if (operator_i == ALU_EXTS)
            shuf_result = {{24{operand_a_i[8*imm_vec_ext_i[1:0]+7]}},
                           operand_a_i[8*imm_vec_ext_i[1:0] +: 8]};
        end else begin
          shuf_result = {16'b0, operand_a_i[16*imm_vec_ext_i[0] +: 16]};
          if (operator_i == ALU_EXTS)
            shuf_result = {{16{operand_a_i[16*imm_vec_ext_i[0]+15]}},
                           operand_a_i[16*imm_vec_ext_i[0] +: 16]};
        end
      end
      ALU_INS: begin
        shuf_result = operand_c_i;
        if (vector_mode_i == VEC_MODE8)
          shuf_result[8*imm_vec_ext_i[1:0] +: 8] = operand_a_i[7:0];
        else
          shuf_result[16*imm_vec_ext_i[0] +: 16] = operand_a_i[15:0];
      end
      ALU_SHUF: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_result = {byte_mux(operand_a_i, operand_b_i[25:24]),
                         byte_mux(operand_a_i, operand_b_i[17:16]),
                         byte_mux(operand_a_i, operand_b_i[9:8]),
                         byte_mux(operand_a_i, operand_b_i[1:0])};
        else
          shuf_result = {half_mux(operand_a_i, operand_b_i[16]),
                         half_mux(operand_a_i, operand_b_i[0])};
      end
      ALU_SHUF2: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_result = {byte_mux(operand_b_i[26] ? operand_a_i : operand_c_i, operand_b_i[25:24]),
                         byte_mux(operand_b_i[18] ? operand_a_i : operand_c_i, operand_b_i[17:16]),
                         byte_mux(operand_b_i[10] ? operand_a_i : operand_c_i, operand_b_i[9:8]),
                         byte_mux(operand_b_i[2]  ? operand_a_i : operand_c_i, operand_b_i[1:0])};
        else
          shuf_result = {half_mux(operand_b_i[17] ? operand_a_i : operand_c_i, operand_b_i[16]),
                         half_mux(operand_b_i[1]  ? operand_a_i : operand_c_i, operand_b_i[0])};
      end
      ALU_PCKLO: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_result = {operand_c_i[31:16], operand_a_i[7:0], operand_b_i[7:0]};
        else
          shuf_result = {operand_a_i[15:0], operand_b_i[15:0]};
      end
      ALU_PCKHI: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_result = {operand_a_i[7:0], operand_b_i[7:0], operand_c_i[15:0]};
        else
          shuf_result = {operand_a_i[31:16], operand_b_i[31:16]};
      end
      default: ;
    endcase
  end

  logic        div_signed;
  logic [5:0]  div_shift;
  logic [31:0] div_shifted_a;
  logic [31:0] div_result;
  logic        div_ready;
  logic        div_valid;

  assign div_signed = operator_i[0];
  assign div_valid  = enable_i & ((operator_i == ALU_DIV)  || (operator_i == ALU_DIVU) ||
                                  (operator_i == ALU_REM)  || (operator_i == ALU_REMU));

  always_comb begin
    if (operand_a_i == 32'h0)
      div_shift = 6'd32;
    else if (div_signed & operand_a_i[31]) begin
      if ((~operand_a_i) == 32'h0)
        div_shift = 6'd31;
      else
        div_shift = 6'd30 - {1'b0, ff1_lsb(~operand_a_i)};
    end else begin
      div_shift = 6'd31 - {1'b0, ff1_msb(operand_a_i)};
    end
  end

  assign div_shifted_a = operand_a_i << div_shift[4:0];

  cv32e40p_alu_div alu_div_i (
      .Clk_CI      (clk),
      .Rst_RBI     (rst_n),
      .OpA_DI      (operand_b_i),
      .OpB_DI      (div_shifted_a),
      .OpBShift_DI (div_shift),
      .OpBIsZero_SI(operand_a_i == 32'h0),
      .OpBSign_SI  (div_signed & operand_a_i[31]),
      .OpCode_SI   (operator_i[1:0]),
      .Res_DO      (div_result),
      .InVld_SI    (div_valid),
      .OutRdy_SI   (ex_ready_i),
      .OutVld_SO   (div_ready)
  );

  always_comb begin
    result_o = '0;
    unique case (operator_i)
      ALU_AND: result_o = operand_a_i & operand_b_i;
      ALU_OR:  result_o = operand_a_i | operand_b_i;
      ALU_XOR: result_o = operand_a_i ^ operand_b_i;

      ALU_SLL: result_o = lane_shl(operand_a_i, vector_mode_i, operand_b_i);
      ALU_SRL: result_o = lane_shr(operand_a_i, vector_mode_i, operand_b_i, 1'b0);
      ALU_SRA: result_o = lane_shr(operand_a_i, vector_mode_i, operand_b_i, 1'b1);
      ALU_ROR: begin
        if (vector_mode_i == VEC_MODE32)
          result_o = ror32(operand_a_i, operand_b_i[4:0]);
        else
          result_o = lane_shr(operand_a_i, vector_mode_i, operand_b_i, 1'b0);
      end

      ALU_ADD, ALU_ADDR, ALU_ADDU, ALU_ADDUR,
      ALU_SUB, ALU_SUBR, ALU_SUBU, ALU_SUBUR:
        result_o = lane_shr(adder_rounded, vector_mode_i, adder_shift_amounts, adder_shift_arith);

      ALU_BCLR: result_o = operand_a_i & ~bmask;
      ALU_BSET: result_o = operand_a_i | bmask;
      ALU_BEXTU: result_o = (operand_a_i >> operand_b_i[4:0]) & bmask;
      ALU_BEXT: begin
        logic [31:0] t;
        t = $signed(operand_a_i) >>> operand_b_i[4:0];
        result_o = (t & bmask) | ({32{t[bmask_a_i]}} & ~bmask);
      end
      ALU_BINS: result_o = ((operand_a_i << operand_b_i[4:0]) & bmask) |
                           (operand_c_i & ~bmask);
      ALU_BREV: result_o = reverse_result;

      ALU_FF1, ALU_FL1, ALU_CNT, ALU_CLB:
        result_o = {26'b0, bitcount_result};

      ALU_SHUF, ALU_SHUF2, ALU_PCKLO, ALU_PCKHI, ALU_EXT, ALU_EXTS, ALU_INS:
        result_o = shuf_result;

      ALU_MIN, ALU_MINU, ALU_MAX, ALU_MAXU:
        result_o = minmax_result;

      ALU_ABS: begin
        if (is_clpx_i)
          result_o = {adder_raw[31:16], operand_a_i[15:0]};
        else
          result_o = minmax_result;
      end

      ALU_CLIP, ALU_CLIPU:
        result_o = clip_result;

      ALU_EQ, ALU_NE, ALU_GTU, ALU_GEU, ALU_LTU, ALU_LEU,
      ALU_GTS, ALU_GES, ALU_LTS, ALU_LES: begin
        result_o[31:24] = {8{cmp_result[3]}};
        result_o[23:16] = {8{cmp_result[2]}};
        result_o[15:8]  = {8{cmp_result[1]}};
        result_o[7:0]   = {8{cmp_result[0]}};
      end

      ALU_SLTS, ALU_SLTU, ALU_SLETS, ALU_SLETU:
        result_o = {31'b0, comparison_result_o};

      ALU_DIV, ALU_DIVU, ALU_REM, ALU_REMU:
        result_o = div_result;

      default: ;
    endcase
  end

  assign ready_o = div_ready;

endmodule

// Same serial divider as the reference.
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

  enum logic [1:0] { IDLE, DIVIDE, FINISH } State_SN, State_SP;

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

  assign RemSel_SN = (LoadEn_S) ? OpCode_SI[1] : RemSel_SP;
  assign CompInv_SN = (LoadEn_S) ? OpBSign_SI : CompInv_SP;
  assign ResInv_SN = (LoadEn_S) ? (~OpBIsZero_SI | OpCode_SI[1]) & OpCode_SI[0] & (OpA_DI[$high(OpA_DI)] ^ OpBSign_SI) : ResInv_SP;

  assign AReg_DN   = (ARegEn_S) ? AddOut_D : AReg_DP;
  assign BReg_DN   = (BRegEn_S) ? BMux_D : BReg_DP;
  assign ResReg_DN = (LoadEn_S) ? '0 : (ResRegEn_S) ? { ABComp_S, ResReg_DP[$high(ResReg_DP):1] } : ResReg_DP;

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
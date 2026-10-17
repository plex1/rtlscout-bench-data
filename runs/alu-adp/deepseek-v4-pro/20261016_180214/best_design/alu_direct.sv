// -----------------------------------------------------------------------------
// Direct behavioral implementation of cv32e40p_alu
// -----------------------------------------------------------------------------

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

  // mask shared by adder rounding and bit-manipulation
  logic [31:0] bmask;
  assign bmask = (~(32'hFFFFFFFE << bmask_a_i)) << bmask_b_i;

  // ------------------------------------------------------------------
  // small helper functions
  // ------------------------------------------------------------------
  function automatic [4:0] fl1_func(input [31:0] x);
    integer j;
    begin
      fl1_func = 5'd0;
      for (j = 0; j < 32; j = j + 1) begin
        if (x[j]) fl1_func = j[4:0];
      end
    end
  endfunction

  function automatic [4:0] ff1_func(input [31:0] x);
    integer j;
    begin
      ff1_func = 5'd0;
      for (j = 31; j >= 0; j = j - 1) begin
        if (x[j]) ff1_func = j[4:0];
      end
    end
  endfunction

  function automatic [7:0] byte_mux(input [31:0] word, input [1:0] sel);
    case (sel)
      2'd0:    byte_mux = word[7:0];
      2'd1:    byte_mux = word[15:8];
      2'd2:    byte_mux = word[23:16];
      default: byte_mux = word[31:24];
    endcase
  endfunction

  function automatic [15:0] half_mux(input [31:0] word, input sel);
    half_mux = sel ? word[31:16] : word[15:0];
  endfunction

  // lane add/sub with no cross-lane carries inside the raw result
  function automatic [31:0] lane_addsub(input [31:0] a, input [31:0] b,
                                        input [1:0] mode, input sub);
    case (mode)
      VEC_MODE16: begin
        if (sub)
          lane_addsub = {a[31:16] - b[31:16], a[15:0] - b[15:0]};
        else
          lane_addsub = {a[31:16] + b[31:16], a[15:0] + b[15:0]};
      end
      VEC_MODE8: begin
        if (sub)
          lane_addsub = {a[31:24] - b[31:24],
                         a[23:16] - b[23:16],
                         a[15:8]  - b[15:8],
                         a[7:0]   - b[7:0]};
        else
          lane_addsub = {a[31:24] + b[31:24],
                         a[23:16] + b[23:16],
                         a[15:8]  + b[15:8],
                         a[7:0]   + b[7:0]};
      end
      default: begin
        lane_addsub = sub ? (a - b) : (a + b);
      end
    endcase
  endfunction

  function automatic [31:0] lane_shr(input [31:0] x, input [1:0] mode,
                                     input [31:0] amt_packed, input arith);
    logic [31:0] out;
    begin
      out = '0;
      case (mode)
        VEC_MODE16: begin
          if (arith) begin
            out[15:0]  = $signed(x[15:0])  >>> amt_packed[3:0];
            out[31:16] = $signed(x[31:16]) >>> amt_packed[19:16];
          end else begin
            out[15:0]  = x[15:0]  >> amt_packed[3:0];
            out[31:16] = x[31:16] >> amt_packed[19:16];
          end
        end
        VEC_MODE8: begin
          if (arith) begin
            out[7:0]   = $signed(x[7:0])   >>> amt_packed[2:0];
            out[15:8]  = $signed(x[15:8])  >>> amt_packed[10:8];
            out[23:16] = $signed(x[23:16]) >>> amt_packed[18:16];
            out[31:24] = $signed(x[31:24]) >>> amt_packed[26:24];
          end else begin
            out[7:0]   = x[7:0]   >> amt_packed[2:0];
            out[15:8]  = x[15:8]  >> amt_packed[10:8];
            out[23:16] = x[23:16] >> amt_packed[18:16];
            out[31:24] = x[31:24] >> amt_packed[26:24];
          end
        end
        default: begin
          if (arith) out = $signed(x) >>> amt_packed[4:0];
          else       out = x >> amt_packed[4:0];
        end
      endcase
      lane_shr = out;
    end
  endfunction

  function automatic [31:0] lane_shl(input [31:0] x, input [1:0] mode,
                                     input [31:0] amt_packed);
    case (mode)
      VEC_MODE16: lane_shl = {x[31:16] << amt_packed[19:16],
                              x[15:0]  << amt_packed[3:0]};
      VEC_MODE8:  lane_shl = {x[31:24] << amt_packed[26:24],
                              x[23:16] << amt_packed[18:16],
                              x[15:8]  << amt_packed[10:8],
                              x[7:0]   << amt_packed[2:0]};
      default:     lane_shl = x << amt_packed[4:0];
    endcase
  endfunction

  function automatic [31:0] ror32_func(input [31:0] x, input [4:0] amt);
    ror32_func = (amt == 5'd0) ? x : ((x >> amt) | (x << (6'd32 - {1'b0, amt})));
  endfunction

  // ------------------------------------------------------------------
  // byte-level compare primitives, aggregated like the reference
  // ------------------------------------------------------------------
  logic [3:0] cmp_signed_bytes;
  always_comb begin
    cmp_signed_bytes = 4'b0000;

    unique case (operator_i)
      ALU_GTS, ALU_GES, ALU_LTS, ALU_LES,
      ALU_SLTS, ALU_SLETS, ALU_MIN, ALU_MAX,
      ALU_ABS, ALU_CLIP, ALU_CLIPU: begin
        case (vector_mode_i)
          VEC_MODE8:  cmp_signed_bytes = 4'b1111;
          VEC_MODE16: cmp_signed_bytes = 4'b1010;
          default:    cmp_signed_bytes = 4'b1000;
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
      assign bgt[gi] = $signed({operand_a_i[8*gi+7] & cmp_signed_bytes[gi],
                                operand_a_i[8*gi +: 8]}) >
                       $signed({operand_b_i[8*gi+7] & cmp_signed_bytes[gi],
                                operand_b_i[8*gi +: 8]});
    end
  endgenerate

  logic [3:0] is_equal;
  logic [3:0] is_greater;

  always_comb begin
    // 32-bit default
    is_equal[3:0] = {4{beq[3] & beq[2] & beq[1] & beq[0]}};
    is_greater[3:0] = {4{bgt[3] | (beq[3] & (bgt[2] |
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

  // ------------------------------------------------------------------
  // adder / shifter datapath
  // ------------------------------------------------------------------
  logic        adder_sub;
  logic [31:0] adder_a;
  logic [31:0] adder_b;
  logic [31:0] adder_raw;
  logic [31:0] adder_rounded;
  logic [31:0] adder_shift_amt_packed;
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
    adder_raw = lane_addsub(adder_a, adder_b, vector_mode_i, adder_sub);
    if (operator_i == ALU_ABS)
      adder_raw = lane_addsub('0, operand_a_i, vector_mode_i, 1'b1);
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
    if (is_clpx_i && adder_shift_use)
      adder_shift_amt_packed = {{14'b0, clpx_shift_i}, {14'b0, clpx_shift_i}};
    else
      adder_shift_amt_packed = {4{3'b000, bmask_b_i}};
  end

  // ------------------------------------------------------------------
  // main result mux
  // ------------------------------------------------------------------
  logic [31:0] minmax_b;
  logic [ 3:0] minmax_sel;
  logic        do_min;
  logic [31:0] minmax_result;
  logic [31:0] clip_result;
  logic [31:0] pack_result;
  logic [ 5:0] cnt_result;
  logic [ 5:0] bitcount_result;
  logic [31:0] reverse_result;
  logic [31:0] shuf_word;
  logic [31:0] shuf_res;
  logic [31:0] div_shifted_a;

  // popcount submodule (small, used by CNT and divider zero check)
  cv32e40p_popcnt popcnt_i (
      .in_i    (operand_a_i),
      .result_o(cnt_result)
  );

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

  always_comb begin
    if (operator_i == ALU_CLIPU) begin
      if ($signed(operand_a_i) <= 32'sd0)
        clip_result = '0;
      else if ($signed(operand_a_i) > $signed(operand_b_i))
        clip_result = operand_b_i;
      else
        clip_result = operand_a_i;
    end else if (operator_i == ALU_CLIP) begin
      if ($signed(operand_a_i) <= $signed(~operand_b_i))
        clip_result = ~operand_b_i;
      else if ($signed(operand_a_i) > $signed(operand_b_i))
        clip_result = operand_b_i;
      else
        clip_result = operand_a_i;
    end else begin
      clip_result = minmax_result;
    end
  end

  // bit-count results
  always_comb begin
    bitcount_result = '0;
    unique case (operator_i)
      ALU_CNT: bitcount_result = cnt_result;
      ALU_FF1: bitcount_result = (operand_a_i == 32'h0) ? 6'd32 : {1'b0, ff1_func(operand_a_i)};
      ALU_FL1: bitcount_result = (operand_a_i == 32'h0) ? 6'd32 : {1'b0, fl1_func(operand_a_i)};
      ALU_CLB: begin
        if (operand_a_i == 32'h0)
          bitcount_result = 6'd0;
        else if (operand_a_i[31])
          bitcount_result = ((~operand_a_i) == 32'h0) ? 6'd31 :
                            {2'b00, 30 - fl1_func(~operand_a_i)};
        else
          bitcount_result = {2'b00, 30 - fl1_func(operand_a_i)};
      end
      default: ;
    endcase
  end

  // BREV result
  always_comb begin
    reverse_result = '0;
    if (operator_i == ALU_BREV) begin
      logic [31:0] t;
      t = operand_a_i << operand_b_i[4:0];

      case (bmask_a_i[1:0])
        2'b00, 2'b11: begin
          integer j;
          for (j = 0; j < 32; j = j + 1)
            reverse_result[j] = t[31-j];
        end
        2'b01: begin
          integer j;
          for (j = 0; j < 16; j = j + 1)
            reverse_result[2*j +: 2] = {t[31-2*j], t[30-2*j]};
        end
        2'b10: begin
          integer j;
          for (j = 0; j < 10; j = j + 1)
            reverse_result[3*j +: 3] = {t[31-3*j], t[30-3*j], t[29-3*j]};
          reverse_result[31:30] = 2'b00;
        end
        default: ;
      endcase
    end
  end

  // shuffle/pack/extract path
  always_comb begin
    shuf_word = '0;
    shuf_res  = '0;

    unique case (operator_i)
      ALU_EXT, ALU_EXTS: begin
        if (vector_mode_i == VEC_MODE8) begin
          if (operator_i == ALU_EXTS)
            shuf_res = {{24{operand_a_i[8*imm_vec_ext_i[1:0]+7]}},
                        operand_a_i[8*imm_vec_ext_i[1:0] +: 8]};
          else
            shuf_res = {24'b0, operand_a_i[8*imm_vec_ext_i[1:0] +: 8]};
        end else begin
          if (operator_i == ALU_EXTS)
            shuf_res = {{16{operand_a_i[16*imm_vec_ext_i[0]+15]}},
                        operand_a_i[16*imm_vec_ext_i[0] +: 16]};
          else
            shuf_res = {16'b0, operand_a_i[16*imm_vec_ext_i[0] +: 16]};
        end
      end

      ALU_INS: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_res = operand_c_i;
          shuf_res[8*imm_vec_ext_i[1:0] +: 8] = operand_a_i[7:0];
        end else begin
          shuf_res = operand_c_i;
          shuf_res[16*imm_vec_ext_i[0] +: 16] = operand_a_i[15:0];
        end
      end

      ALU_SHUF: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_res[7:0]   = byte_mux(operand_a_i, operand_b_i[1:0]);
          shuf_res[15:8]  = byte_mux(operand_a_i, operand_b_i[9:8]);
          shuf_res[23:16] = byte_mux(operand_a_i, operand_b_i[17:16]);
          shuf_res[31:24] = byte_mux(operand_a_i, operand_b_i[25:24]);
        end else begin
          shuf_res[15:0]  = half_mux(operand_a_i, operand_b_i[0]);
          shuf_res[31:16] = half_mux(operand_a_i, operand_b_i[16]);
        end
      end

      ALU_SHUF2: begin
        if (vector_mode_i == VEC_MODE8) begin
          shuf_res[7:0]   = byte_mux(operand_b_i[2]  ? operand_a_i : operand_c_i, operand_b_i[1:0]);
          shuf_res[15:8]  = byte_mux(operand_b_i[10] ? operand_a_i : operand_c_i, operand_b_i[9:8]);
          shuf_res[23:16] = byte_mux(operand_b_i[18] ? operand_a_i : operand_c_i, operand_b_i[17:16]);
          shuf_res[31:24] = byte_mux(operand_b_i[26] ? operand_a_i : operand_c_i, operand_b_i[25:24]);
        end else begin
          shuf_res[15:0]  = half_mux(operand_b_i[1] ? operand_a_i : operand_c_i, operand_b_i[0]);
          shuf_res[31:16] = half_mux(operand_b_i[17] ? operand_a_i : operand_c_i, operand_b_i[16]);
        end
      end

      ALU_PCKLO: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_res = {operand_c_i[31:16], operand_a_i[7:0], operand_b_i[7:0]};
        else
          shuf_res = {operand_a_i[15:0], operand_b_i[15:0]};
      end

      ALU_PCKHI: begin
        if (vector_mode_i == VEC_MODE8)
          shuf_res = {operand_a_i[7:0], operand_b_i[7:0], operand_c_i[15:0]};
        else
          shuf_res = {operand_a_i[31:16], operand_b_i[31:16]};
      end

      default: ;
    endcase
  end

  always_comb begin
    pack_result = shuf_res;
  end

  // ------------------------------------------------------------------
  // serial divider interface
  // ------------------------------------------------------------------
  logic        div_signed;
  logic [31:0] div_result;
  logic        div_ready;
  logic [ 5:0] div_shift;
  logic        div_valid;

  assign div_signed = operator_i[0];
  assign div_valid  = enable_i & ((operator_i == ALU_DIV)  || (operator_i == ALU_DIVU) ||
                                  (operator_i == ALU_REM)  || (operator_i == ALU_REMU));

  // Replicate the reference divider's shift amount computation.
  always_comb begin
    if (operand_a_i == 32'h0)
      div_shift = 6'd32;
    else if (div_signed & operand_a_i[31]) begin
      if ((~operand_a_i) == 32'h0)
        div_shift = 6'd31;
      else
        div_shift = {2'b00, 30 - fl1_func(~operand_a_i)};
    end else begin
      div_shift = {2'b00, 31 - fl1_func(operand_a_i)};
    end
  end

  // The reference ALU feeds the divider with operand_a left-shifted by div_shift.
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

  // ------------------------------------------------------------------
  // final result mux
  // ------------------------------------------------------------------
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
          result_o = ror32_func(operand_a_i, operand_b_i[4:0]);
        else
          result_o = lane_shr(operand_a_i, vector_mode_i, operand_b_i, 1'b0);
      end

      ALU_ADD, ALU_ADDR, ALU_ADDU, ALU_ADDUR,
      ALU_SUB, ALU_SUBR, ALU_SUBU, ALU_SUBUR:
        result_o = lane_shr(adder_rounded, vector_mode_i,
                            adder_shift_amt_packed, adder_shift_arith);

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
        result_o = pack_result;

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
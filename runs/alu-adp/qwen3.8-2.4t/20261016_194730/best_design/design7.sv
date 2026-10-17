// cv32e40p_alu - v7: design4 + divider comparator derived from the subtractor
// (removes the independent == and > 32-bit compare trees from the divider).

module cv32e40p_alu
(
    input logic               clk,
    input logic               rst_n,
    input logic               enable_i,
    input logic [ 6:0]        operator_i,
    input logic [31:0]        operand_a_i,
    input logic [31:0]        operand_b_i,
    input logic [31:0]        operand_c_i,

    input logic [1:0]         vector_mode_i,
    input logic [4:0]         bmask_a_i,
    input logic [4:0]         bmask_b_i,
    input logic [1:0]         imm_vec_ext_i,

    input logic               is_clpx_i,
    input logic               is_subrot_i,
    input logic [1:0]         clpx_shift_i,

    output logic [31:0]       result_o,
    output logic              comparison_result_o,

    output logic              ready_o,
    input  logic              ex_ready_i
);

  /* verilator lint_off WIDTHTRUNC */
  /* verilator lint_off CASEINCOMPLETE */
  /* verilator lint_off UNUSEDSIGNAL */

  localparam logic [6:0] ALU_ADD   = 7'b0011000;
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
  localparam logic [6:0] ALU_ABS   = 7'b0010100;
  localparam logic [6:0] ALU_CLIP  = 7'b0010110;
  localparam logic [6:0] ALU_CLIPU = 7'b0010111;
  localparam logic [6:0] ALU_INS   = 7'b0101101;
  localparam logic [6:0] ALU_MIN   = 7'b0010000;
  localparam logic [6:0] ALU_MINU  = 7'b0010001;
  localparam logic [6:0] ALU_MAX   = 7'b0010010;
  localparam logic [6:0] ALU_MAXU  = 7'b0010011;
  localparam logic [6:0] ALU_SHUF  = 7'b0111010;
  localparam logic [6:0] ALU_SHUF2 = 7'b0111011;
  localparam logic [6:0] ALU_PCKLO = 7'b0111000;
  localparam logic [6:0] ALU_PCKHI = 7'b0111001;
  localparam logic [1:0] VEC_MODE16 = 2'b10;
  localparam logic [1:0] VEC_MODE8  = 2'b11;

  // ----------------------------------------------------------------
  // opcode group decodes
  // ----------------------------------------------------------------
  wire mode32 = (vector_mode_i == 2'b00);
  wire mode16 = (vector_mode_i == VEC_MODE16);
  wire mode8  = (vector_mode_i == VEC_MODE8);

  wire g_addfam = (operator_i[6:3] == 4'b0011);
  wire g_shift  = (operator_i[6:2] == 5'b01001);         // SRA SRL ROR SLL
  wire g_div    = (operator_i[6:2] == 5'b01100);         // DIVU DIV REMU REM
  wire g_bitcnt = (operator_i[6:2] == 5'b01101);         // CNT CLB FF1 FL1
  wire g_brev   = (operator_i == ALU_BREV);
  wire g_bmanip = (operator_i[6:2] == 5'b01010) ||       // BEXT BEXTU BINS BCLR
                  (operator_i == ALU_BSET) || g_brev;
  wire g_bext   = (operator_i == ALU_BEXT);
  wire g_bins   = (operator_i == ALU_BINS);
  wire g_cmp    = (operator_i[6:4] == 3'b000);
  wire g_minmax = (operator_i[6:2] == 5'b00100);
  wire g_abs    = (operator_i == ALU_ABS);
  wire g_clip   = (operator_i == ALU_CLIP) || (operator_i == ALU_CLIPU);
  wire g_shuf   = (operator_i[6:2] == 5'b01110);
  wire g_extins = (operator_i == ALU_EXT) || (operator_i == ALU_EXTS) || (operator_i == ALU_INS);
  wire g_logic  = (operator_i == ALU_AND) || (operator_i == ALU_OR) || (operator_i == ALU_XOR);

  wire sll_op = (operator_i == ALU_SLL);
  wire sra_op = (operator_i == ALU_SRA);
  wire ror_op = (operator_i == ALU_ROR);

  // ----------------------------------------------------------------
  // segmented adder (adder family, ABS, CLIP)
  // ----------------------------------------------------------------
  wire addfam_sub = g_addfam & operator_i[0];
  wire b_invert   = addfam_sub | is_subrot_i;
  wire abs_op     = g_abs;
  wire clip_op    = (operator_i == ALU_CLIP);
  wire neg        = b_invert | abs_op | clip_op;

  wire [31:0] add_op_a = abs_op ? ~operand_a_i :
                         is_subrot_i ? {operand_b_i[15:0], operand_a_i[31:16]} : operand_a_i;
  wire [31:0] add_op_b = b_invert ? (is_subrot_i ? ~{operand_a_i[15:0], operand_b_i[31:16]}
                                                  : ~operand_b_i)
                                   : operand_b_i;

  wire sep9_a  = neg | ~mode8;
  wire sep18_a = neg | mode32;
  wire sep27_a = neg | ~mode8;
  wire sep0_b  = neg;
  wire sep9_b  = neg & mode8;
  wire sep18_b = neg & (mode8 | mode16);
  wire sep27_b = neg & mode8;

  wire [35:0] adder_in_a = {add_op_a[31:24], sep27_a,
                            add_op_a[23:16], sep18_a,
                            add_op_a[15:8],  sep9_a,
                            add_op_a[7:0],   1'b1};
  wire [35:0] adder_in_b = {add_op_b[31:24], sep27_b,
                            add_op_b[23:16], sep18_b,
                            add_op_b[15:8],  sep9_b,
                            add_op_b[7:0],   sep0_b};

  wire [36:0] adder_result_expanded = $signed(adder_in_a) + $signed(adder_in_b);
  wire [31:0] adder_result = {adder_result_expanded[35:28],
                              adder_result_expanded[26:19],
                              adder_result_expanded[17:10],
                              adder_result_expanded[8:1]};

  // ----------------------------------------------------------------
  // rounding add (adder family)
  // ----------------------------------------------------------------
  wire [31:0] bmask = (~((32'hFFFFFFFE) << bmask_a_i)) << bmask_b_i;
  wire addfam_round = operator_i[2];   // meaningful within g_addfam
  wire [31:0] round_sum = adder_result + (addfam_round ? {1'b0, bmask[31:1]} : 32'h0);

  // ----------------------------------------------------------------
  // RS_A: adder-family per-lane right shifter (data = round_sum)
  // ----------------------------------------------------------------
  wire addfam_arith = ~operator_i[1];
  wire [4:0] fam_sh32 = is_clpx_i ? {3'b000, clpx_shift_i} : bmask_b_i;
  wire [3:0] fam_sh16 = is_clpx_i ? {2'b00, clpx_shift_i}  : bmask_b_i[3:0];
  wire [2:0] fam_sh8  = is_clpx_i ? {1'b0, clpx_shift_i}   : bmask_b_i[2:0];

  logic [31:0] rs_a_result;
  always_comb begin
    if (mode16) begin
      rs_a_result[31:16] = $signed({addfam_arith & round_sum[31], round_sum[31:16]}) >>> fam_sh16;
      rs_a_result[15:0]  = $signed({addfam_arith & round_sum[15], round_sum[15:0]})  >>> fam_sh16;
    end else if (mode8) begin
      rs_a_result[31:24] = $signed({addfam_arith & round_sum[31], round_sum[31:24]}) >>> fam_sh8;
      rs_a_result[23:16] = $signed({addfam_arith & round_sum[23], round_sum[23:16]}) >>> fam_sh8;
      rs_a_result[15:8]  = $signed({addfam_arith & round_sum[15], round_sum[15:8]})  >>> fam_sh8;
      rs_a_result[7:0]   = $signed({addfam_arith & round_sum[7],  round_sum[7:0]})   >>> fam_sh8;
    end else begin
      rs_a_result = $signed({{32{addfam_arith & round_sum[31]}}, round_sum}) >>> fam_sh32;
    end
  end

  // ----------------------------------------------------------------
  // RS_B: shift ops + BEXT/BEXTU per-lane right shifter (data = a)
  // ----------------------------------------------------------------
  wire rs_b_arith = sra_op | g_bext;
  wire ror32 = ror_op & mode32;
  wire [63:0] rs_b_in = ror32 ? {operand_a_i, operand_a_i} :
                                {{32{rs_b_arith & operand_a_i[31]}}, operand_a_i};

  logic [31:0] rs_b_result;
  always_comb begin
    if (mode16) begin
      rs_b_result[31:16] = $signed({rs_b_arith & operand_a_i[31], operand_a_i[31:16]}) >>> operand_b_i[19:16];
      rs_b_result[15:0]  = $signed({rs_b_arith & operand_a_i[15], operand_a_i[15:0]})  >>> operand_b_i[3:0];
    end else if (mode8) begin
      rs_b_result[31:24] = $signed({rs_b_arith & operand_a_i[31], operand_a_i[31:24]}) >>> operand_b_i[26:24];
      rs_b_result[23:16] = $signed({rs_b_arith & operand_a_i[23], operand_a_i[23:16]}) >>> operand_b_i[18:16];
      rs_b_result[15:8]  = $signed({rs_b_arith & operand_a_i[15], operand_a_i[15:8]})  >>> operand_b_i[10:8];
      rs_b_result[7:0]   = $signed({rs_b_arith & operand_a_i[7],  operand_a_i[7:0]})   >>> operand_b_i[2:0];
    end else begin
      rs_b_result = rs_b_in >> operand_b_i[4:0];
    end
  end

  // ----------------------------------------------------------------
  // RS_L: right shifter over bit-reversed data => left shift
  //   users: SLL (per lane), BINS/BREV (32-bit), div divisor prep
  // ----------------------------------------------------------------
  wire [31:0] operand_a_rev;
  genvar ri;
  generate
    for (ri = 0; ri < 32; ri++) begin : gen_a_rev
      assign operand_a_rev[ri] = operand_a_i[31-ri];
    end
  endgenerate

  // ff1 machinery (bit-count ops + div_shift)
  wire ff1_op = (operator_i == ALU_FF1);
  wire sign_flip = (g_div & operator_i[0]) | (operator_i == ALU_CLB);
  wire [31:0] ff_input = ff1_op ? operand_a_i : (operand_a_rev ^ {32{sign_flip & operand_a_i[31]}});
  wire [31:0] ff_onehot = ff_input & (~ff_input + 32'd1);
  wire        ff_no_one = ~|ff_input;
  logic [4:0] ff1_result;
  always_comb begin
    ff1_result = 5'd0;
    for (int i = 1; i < 32; i++)
      if (ff_onehot[i]) ff1_result = i[4:0];
  end
  wire [5:0] clb_result = ff1_result - 5'd1;

  wire         div_op_a_signed = operand_a_i[31] & operator_i[0];
  wire [5:0]   div_shift_int   = ff_no_one ? 6'd31 : clb_result;
  wire [5:0]   div_shift       = div_shift_int + (div_op_a_signed ? 6'd0 : 6'd1);

  // reversed shift amount for per-lane left shifts
  logic [31:0] ls_amt_rev;
  wire [31:0] ls_amt_src = g_div ? {26'b0, div_shift[4:0]} : operand_b_i;
  always_comb begin
    case (vector_mode_i)
      VEC_MODE16: ls_amt_rev = {ls_amt_src[15:0], ls_amt_src[31:16]};
      VEC_MODE8:  ls_amt_rev = {ls_amt_src[7:0], ls_amt_src[15:8], ls_amt_src[23:16], ls_amt_src[31:24]};
      default:    ls_amt_rev = ls_amt_src;
    endcase
  end

  logic [31:0] rs_l_out;   // right-shift of reversed data
  always_comb begin
    if (mode16 & sll_op) begin
      rs_l_out[31:16] = operand_a_rev[31:16] >> ls_amt_rev[19:16];
      rs_l_out[15:0]  = operand_a_rev[15:0]  >> ls_amt_rev[3:0];
    end else if (mode8 & sll_op) begin
      rs_l_out[31:24] = operand_a_rev[31:24] >> ls_amt_rev[26:24];
      rs_l_out[23:16] = operand_a_rev[23:16] >> ls_amt_rev[18:16];
      rs_l_out[15:8]  = operand_a_rev[15:8]  >> ls_amt_rev[10:8];
      rs_l_out[7:0]   = operand_a_rev[7:0]   >> ls_amt_rev[2:0];
    end else begin
      rs_l_out = operand_a_rev >> ls_amt_rev[4:0];
    end
  end

  wire [31:0] shift_left_result;
  generate
    for (ri = 0; ri < 32; ri++) begin : gen_left_result
      assign shift_left_result[ri] = rs_l_out[31-ri];
    end
  endgenerate

  // ----------------------------------------------------------------
  // bit manipulation (32-bit mode)
  // ----------------------------------------------------------------
  wire        extract_sign   = g_bext & rs_b_result[bmask_a_i];
  wire [31:0] bextins_and    = g_bins ? operand_c_i : {32{extract_sign}};
  wire [31:0] bextins_result = (bmask & (g_bins ? shift_left_result : rs_b_result))
                             | (bextins_and & ~bmask);
  wire [31:0] bclr_result    = operand_a_i & ~bmask;
  wire [31:0] bset_result    = operand_a_i | bmask;

  logic [31:0] radix_2_rev, radix_4_rev, radix_8_rev, brev_result;
  genvar bj;
  generate
    for (bj = 0; bj < 32; bj++) begin : gen_radix_2_rev
      assign radix_2_rev[bj] = shift_left_result[31-bj];
    end
    for (bj = 0; bj < 16; bj++) begin : gen_radix_4_rev
      assign radix_4_rev[2*bj+1:2*bj] = shift_left_result[31-bj*2:31-bj*2-1];
    end
    for (bj = 0; bj < 10; bj++) begin : gen_radix_8_rev
      assign radix_8_rev[3*bj+2:3*bj] = shift_left_result[31-bj*3:31-bj*3-2];
    end
    assign radix_8_rev[31:30] = 2'b0;
  endgenerate
  always_comb begin
    unique case (bmask_a_i[1:0])
      2'b01:   brev_result = radix_4_rev;
      2'b10:   brev_result = radix_8_rev;
      default: brev_result = radix_2_rev;
    endcase
  end

  logic [31:0] bmanip_result;
  always_comb begin
    if (operator_i == ALU_BCLR)       bmanip_result = bclr_result;
    else if (operator_i == ALU_BSET)  bmanip_result = bset_result;
    else if (g_brev)                  bmanip_result = brev_result;
    else                              bmanip_result = bextins_result;
  end

  // ----------------------------------------------------------------
  // comparisons + min/max/abs/clip
  // ----------------------------------------------------------------
  logic [3:0] cmp_signed;
  wire [3:0] is_equal_vec;
  wire [3:0] is_greater_vec;

  wire cmp_grp_signed = g_cmp & ~operator_i[0];
  wire minmax_signed  = (operator_i == ALU_MIN) || (operator_i == ALU_MAX);
  wire signed_family  = cmp_grp_signed | minmax_signed | g_abs | g_clip;

  always_comb begin
    if (signed_family) begin
      if (mode8)       cmp_signed = 4'b1111;
      else if (mode16) cmp_signed = 4'b1010;
      else             cmp_signed = 4'b1000;
    end else cmp_signed = 4'b0000;
  end

  genvar ci;
  generate
    for (ci = 0; ci < 4; ci++) begin : gen_cmp_vec
      assign is_equal_vec[ci] = (operand_a_i[8*ci+7:8*ci] == operand_b_i[8*ci+7:8*ci]);
      assign is_greater_vec[ci] = $signed({operand_a_i[8*ci+7] & cmp_signed[ci], operand_a_i[8*ci+7:8*ci]})
                                > $signed({operand_b_i[8*ci+7] & cmp_signed[ci], operand_b_i[8*ci+7:8*ci]});
    end
  endgenerate

  logic [3:0] is_equal, is_greater;
  always_comb begin
    is_equal   = {4{is_equal_vec[3] & is_equal_vec[2] & is_equal_vec[1] & is_equal_vec[0]}};
    is_greater = {4{is_greater_vec[3] | (is_equal_vec[3] & (is_greater_vec[2]
                                       | (is_equal_vec[2] & (is_greater_vec[1]
                                       | (is_equal_vec[1] & is_greater_vec[0])))))}};
    if (mode16) begin
      is_equal[1:0]   = {2{is_equal_vec[0] & is_equal_vec[1]}};
      is_equal[3:2]   = {2{is_equal_vec[2] & is_equal_vec[3]}};
      is_greater[1:0] = {2{is_greater_vec[1] | (is_equal_vec[1] & is_greater_vec[0])}};
      is_greater[3:2] = {2{is_greater_vec[3] | (is_equal_vec[3] & is_greater_vec[2])}};
    end else if (mode8) begin
      is_equal   = is_equal_vec;
      is_greater = is_greater_vec;
    end
  end

  logic [3:0] cmp_result;
  always_comb begin
    unique case (operator_i[6:2])
      5'b00000: cmp_result = ~(is_greater | is_equal);
      5'b00001: cmp_result = ~is_greater;
      5'b00010: cmp_result = operator_i[1] ? (is_greater | is_equal) : is_greater;
      5'b00011: cmp_result = operator_i[0] ? ~is_equal : is_equal;
      default:  cmp_result = is_equal;
    endcase
  end

  assign comparison_result_o = cmp_result[3];

  wire do_min = (operator_i == ALU_MIN) || (operator_i == ALU_MINU) || g_clip;
  wire [3:0] sel_minmax = is_greater ^ {4{do_min}};
  wire [31:0] minmax_b = g_abs ? adder_result : operand_b_i;
  wire [31:0] result_minmax = {sel_minmax[3] ? operand_a_i[31:24] : minmax_b[31:24],
                               sel_minmax[2] ? operand_a_i[23:16] : minmax_b[23:16],
                               sel_minmax[1] ? operand_a_i[15:8]  : minmax_b[15:8],
                               sel_minmax[0] ? operand_a_i[7:0]   : minmax_b[7:0]};

  wire [31:0] operand_b_neg = ~operand_b_i;
  wire is_equal_clip = (operator_i == ALU_CLIPU) ? (operand_a_i == 32'h0)
                                                 : (operand_a_i == operand_b_neg);
  logic [31:0] clip_result;
  always_comb begin
    clip_result = result_minmax;
    if (operator_i == ALU_CLIPU) begin
      if (operand_a_i[31] || is_equal_clip) clip_result = '0;
    end else begin
      if (adder_result_expanded[36] || is_equal_clip) clip_result = operand_b_neg;
    end
  end

  // ----------------------------------------------------------------
  // shuffle / pack / ext / ins
  // ----------------------------------------------------------------
  logic [31:0] shuf_result;
  always_comb begin
    shuf_result = '0;
    unique case (operator_i)
      ALU_PCKLO: shuf_result = mode8 ? {operand_c_i[31:16], operand_a_i[7:0], operand_b_i[7:0]}
                                      : {operand_a_i[15:0], operand_b_i[15:0]};
      ALU_PCKHI: shuf_result = mode8 ? {operand_a_i[7:0], operand_b_i[7:0], operand_c_i[15:0]}
                                      : {operand_a_i[31:16], operand_b_i[31:16]};
      ALU_SHUF: begin
        if (mode8) begin
          shuf_result = {operand_a_i[8*operand_b_i[25:24] +: 8],
                         operand_a_i[8*operand_b_i[17:16] +: 8],
                         operand_a_i[8*operand_b_i[9:8]   +: 8],
                         operand_a_i[8*operand_b_i[1:0]   +: 8]};
        end else begin
          shuf_result = {operand_b_i[16] ? operand_a_i[31:16] : operand_a_i[15:0],
                         operand_b_i[0]  ? operand_a_i[31:16] : operand_a_i[15:0]};
        end
      end
      ALU_SHUF2: begin
        if (mode8) begin
          shuf_result[31:24] = operand_b_i[26] ? operand_a_i[8*operand_b_i[25:24] +: 8]
                                               : operand_c_i[8*operand_b_i[25:24] +: 8];
          shuf_result[23:16] = operand_b_i[18] ? operand_a_i[8*operand_b_i[17:16] +: 8]
                                               : operand_c_i[8*operand_b_i[17:16] +: 8];
          shuf_result[15:8]  = operand_b_i[10] ? operand_a_i[8*operand_b_i[9:8] +: 8]
                                               : operand_c_i[8*operand_b_i[9:8] +: 8];
          shuf_result[7:0]   = operand_b_i[2]  ? operand_a_i[8*operand_b_i[1:0] +: 8]
                                               : operand_c_i[8*operand_b_i[1:0] +: 8];
        end else begin
          shuf_result[31:16] = operand_b_i[17] ? (operand_b_i[16] ? operand_a_i[31:16] : operand_a_i[15:0])
                                               : (operand_b_i[16] ? operand_c_i[31:16] : operand_c_i[15:0]);
          shuf_result[15:0]  = operand_b_i[1]  ? (operand_b_i[0] ? operand_a_i[31:16] : operand_a_i[15:0])
                                               : (operand_b_i[0] ? operand_c_i[31:16] : operand_c_i[15:0]);
        end
      end
      ALU_EXT: begin
        if (mode8)      shuf_result = {24'h0, operand_a_i[8*imm_vec_ext_i +: 8]};
        else            shuf_result = {16'h0, operand_a_i[16*imm_vec_ext_i[0] +: 16]};
      end
      ALU_EXTS: begin
        if (mode8)      shuf_result = {{24{operand_a_i[8*imm_vec_ext_i + 7]}}, operand_a_i[8*imm_vec_ext_i +: 8]};
        else            shuf_result = {{16{operand_a_i[16*imm_vec_ext_i[0] + 15]}}, operand_a_i[16*imm_vec_ext_i[0] +: 16]};
      end
      ALU_INS: begin
        if (mode8) begin
          shuf_result = operand_c_i;
          shuf_result[8*imm_vec_ext_i +: 8] = operand_a_i[7:0];
        end else begin
          shuf_result = imm_vec_ext_i[0] ? {operand_a_i[15:0], operand_c_i[15:0]}
                                         : {operand_c_i[31:16], operand_a_i[15:0]};
        end
      end
      default: ;
    endcase
  end

  // ----------------------------------------------------------------
  // bit counting
  // ----------------------------------------------------------------
  wire [5:0] cnt_result;
  cv32e40p_popcnt popcnt_i (.in_i(operand_a_i), .result_o(cnt_result));

  wire [4:0] fl1_result = 5'd31 - ff1_result;

  logic [5:0] bitop_result;
  always_comb begin
    bitop_result = '0;
    unique case (operator_i)
      ALU_FF1: bitop_result = ff_no_one ? 6'd32 : {1'b0, ff1_result};
      ALU_FL1: bitop_result = ff_no_one ? 6'd32 : {1'b0, fl1_result};
      ALU_CNT: bitop_result = cnt_result;
      ALU_CLB: bitop_result = ff_no_one ? (operand_a_i[31] ? 6'd31 : 6'd0) : clb_result;
      default: ;
    endcase
  end

  // ----------------------------------------------------------------
  // divider (cycle-exact FSM; comparator from subtractor)
  // ----------------------------------------------------------------
  logic [31:0] result_div;
  logic        div_ready;
  wire         div_valid = enable_i & g_div;

  cv32e40p_alu_div alu_div_i (
      .Clk_CI (clk),
      .Rst_RBI(rst_n),
      .OpA_DI      (operand_b_i),
      .OpB_DI      (shift_left_result),
      .OpBShift_DI (div_shift),
      .OpBIsZero_SI(~|operand_a_i),
      .OpBSign_SI(div_op_a_signed),
      .OpCode_SI (operator_i[1:0]),
      .Res_DO(result_div),
      .InVld_SI (div_valid),
      .OutRdy_SI(ex_ready_i),
      .OutVld_SO(div_ready)
  );

  assign ready_o = div_ready;

  // ----------------------------------------------------------------
  // result mux (grouped OR)
  // ----------------------------------------------------------------
  wire [31:0] logic_result = (operator_i == ALU_AND) ? (operand_a_i & operand_b_i) :
                             (operator_i == ALU_OR)  ? (operand_a_i | operand_b_i) :
                                                       (operand_a_i ^ operand_b_i);

  wire [31:0] cmp_fill = {{8{cmp_result[3]}}, {8{cmp_result[2]}}, {8{cmp_result[1]}}, {8{cmp_result[0]}}};
  wire scalar_setless = g_cmp & operator_i[1] & ~operator_i[3];

  wire [31:0] minmax_result = g_abs ? (is_clpx_i ? {adder_result[31:16], operand_a_i[15:0]}
                                                  : result_minmax)
                                    : (g_clip ? clip_result : result_minmax);

  wire g_rs_b_use = g_shift & ~sll_op;   // SRA SRL ROR

  assign result_o = ({32{g_logic}}   & logic_result)
                  | ({32{g_addfam}}  & rs_a_result)
                  | ({32{g_rs_b_use}} & rs_b_result)
                  | ({32{sll_op}}    & shift_left_result)
                  | ({32{g_bmanip}}  & bmanip_result)
                  | ({32{g_bitcnt}}  & {26'h0, bitop_result})
                  | ({32{g_extins | g_shuf}} & shuf_result)
                  | ({32{g_cmp & ~scalar_setless}} & cmp_fill)
                  | ({32{scalar_setless}} & {31'b0, cmp_result[3]})
                  | ({32{g_minmax | g_abs | g_clip}} & minmax_result)
                  | ({32{g_div}}     & result_div);

endmodule

// ---------------- popcount ----------------
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

// ---------------- serial divider (cycle-exact; comparator from subtractor) ----------------
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
  } State_SN, State_SP;

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

  // main adder. During DIVIDE PmSel_S is 0, so SubOut_D is AReg - BReg:
  // derive the upstream comparator ((A==B)|((A>B)^CompInv)) from it.
  assign AddTmp_D = (LoadEn_S) ? 0 : AReg_DP;
  logic [C_WIDTH:0] SubOut_D;
  assign SubOut_D = {1'b0, AddTmp_D} + {1'b0, ~AddMux_D} + {{(C_WIDTH-1){1'b0}}, 2'b01};
  assign AddOut_D = PmSel_S ? (AddTmp_D + AddMux_D) : SubOut_D[C_WIDTH-1:0];
  wire sub_ge  = SubOut_D[C_WIDTH];              // unsigned AReg >= BReg
  wire sub_eq  = ~(|SubOut_D[C_WIDTH-1:0]);      // AReg == BReg
  // A>B = ge & ~eq;  upstream = eq | (gt ^ CompInv):
  //   CompInv=0 -> eq | gt = ge ;  CompInv=1 -> eq | ~gt = ~ge | eq
  wire ab_cmp_raw = CompInv_SP ? (~sub_ge | sub_eq) : sub_ge;
  assign ABComp_S = ab_cmp_raw & ((|AReg_DP) | OpBIsZero_SI);

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

endmodule

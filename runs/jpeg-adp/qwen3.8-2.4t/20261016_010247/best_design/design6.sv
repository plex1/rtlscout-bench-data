//-----------------------------------------------------------------
// jpeg_idct_2d — optimized 8x8 2-D IDCT (bit-exact with core_jpeg original)
// design6: ALL register-to-register paths limited to <= 2 adder levels.
//  - even-set data repeats on beats {0,2,5,7}, odd on {1,3,4,6}: products
//    for later beats precomputed and stored in reused registers
//  - 181* split: partial181(t5)/partial181(t6) at beat 5 (q6/q2 reused),
//    combine at beat 6 (m1a/m7b reused), fix181 at beat 7
//  - direct-write transpose from x pass (no x block_out/streamer)
//-----------------------------------------------------------------

module idct_pass6
#(
     parameter OUT_SHIFT = 11
    ,parameter IN_BITS   = 16
    ,parameter BO_BITS   = 20
    ,parameter TRANSPOSE_IDX = 0
)
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [31:0]   inport_data0_i
    ,input  [31:0]   inport_data1_i
    ,input  [31:0]   inport_data2_i
    ,input  [31:0]   inport_data3_i
    ,input  [ 2:0]   inport_idx_i
    // direct-write outputs (x pass)
    ,output          dw_stro_a_o
    ,output          dw_stro_c_o
    ,output [ 2:0]   dw_row_o
    ,output [19:0]   dw_d0_o, dw_d3_o, dw_d4_o, dw_d7_o
    ,output [19:0]   dw_d1_o, dw_d2_o, dw_d5_o, dw_d6_o
    ,input           dw_accept_i
    // streamed outputs (y pass)
    ,output          outport_valid_o
    ,output [31:0]   outport_data_o
    ,output [ 5:0]   outport_idx_o
);

// depth-2 CSD building blocks (32-bit wrapping):
// C1*x = ((x<<12)+(x)) - ((x<<6)+(x<<4))
// C2*x = ((x<<12)+(x<<3)) - ((x<<8)+(x<<6))
// C5*x = ((x<<11)+(x<<2)) + ((x<<8)-(x<<5))
// C6*x = ((x<<10)-(x)) + ((x<<9)+(x<<5))
// C7*x = ((x<<10)+(x<<5)) - ((x<<8)+(x))
// C4*x = ((x<<12)-(x<<10)) + ((x<<4)-(x<<6)) - (x<<7)
// C3*x = ((x<<12)+(x<<6)) + ((x<<4)-(x<<1)) - ((x<<9)+(x<<8))
// partial181(x) = ((x<<6)+(x<<3)) + ((x<<1)+(x)); 181*x = (x<<8) - partial181
function automatic signed [31:0] mc1(input signed [31:0] x);
begin mc1 = ((x << 12) + x) - ((x << 6) + (x << 4)); end endfunction
function automatic signed [31:0] mc2(input signed [31:0] x);
begin mc2 = ((x << 12) + (x << 3)) - ((x << 8) + (x << 6)); end endfunction
function automatic signed [31:0] mc4(input signed [31:0] x);
begin mc4 = (((x << 12) - (x << 10)) + ((x << 4) - (x << 6))) - (x << 7); end endfunction
function automatic signed [31:0] mc5(input signed [31:0] x);
begin mc5 = ((x << 11) + (x << 2)) + ((x << 8) - (x << 5)); end endfunction
function automatic signed [31:0] mc6(input signed [31:0] x);
begin mc6 = ((x << 10) - x) + ((x << 9) + (x << 5)); end endfunction
function automatic signed [31:0] mc7(input signed [31:0] x);
begin mc7 = ((x << 10) + (x << 5)) - ((x << 8) + x); end endfunction
function automatic signed [31:0] mc3(input signed [31:0] x);
begin mc3 = (((x << 12) + (x << 6)) + ((x << 4) - (x << 1))) - ((x << 9) + (x << 8)); end endfunction
function automatic signed [31:0] partial181(input signed [31:0] x);
begin partial181 = ((x << 6) + (x << 3)) + ((x << 1) + x); end endfunction

// (p = x*181 mod 2^32) -> p/256 signed, truncating toward zero (depth 2)
function automatic signed [31:0] fix181(input [31:0] p);
reg [31:0] q;
begin
    q = p[31] ? (p + 32'd255) : p;
    fix181 = $signed(q) >>> 8;
end endfunction

wire signed [31:0] X0 = {{(32-IN_BITS){inport_data0_i[IN_BITS-1]}}, inport_data0_i[IN_BITS-1:0]};
wire signed [31:0] X1 = {{(32-IN_BITS){inport_data1_i[IN_BITS-1]}}, inport_data1_i[IN_BITS-1:0]};
wire signed [31:0] X2 = {{(32-IN_BITS){inport_data2_i[IN_BITS-1]}}, inport_data2_i[IN_BITS-1:0]};
wire signed [31:0] X3 = {{(32-IN_BITS){inport_data3_i[IN_BITS-1]}}, inport_data3_i[IN_BITS-1:0]};

reg signed [31:0] i0;
reg signed [31:0] q6, q2;   // beat0: mc6(X1),mc2(X3) -> beat5 s2; beat5: partial181(t6/t5)
reg signed [31:0] m1a, m7b; // beat0: mc2(X1),mc6(X3) -> beat1 t3; beat1: C1*c1,C7*c7;
                            // beat6: 181*t65, 181*t56
reg signed [31:0] m7a, m1b; // beat1: C7*c1, C1*c7 -> beat3 s4
reg signed [31:0] p_s6a, p_s6b, p_s5a, p_s5b; // beat1 -> beat2
reg signed [31:0] t0, t3, s7, s6, s5q, t0f, t7, s4q;
reg signed [31:0] t4, t5, t6, t65, t56, t1, t2, s5;

always @(posedge clk)
begin
    if (rst)
    begin
        i0 <= 0; q6 <= 0; q2 <= 0;
        m1a <= 0; m7b <= 0; m7a <= 0; m1b <= 0;
        p_s6a <= 0; p_s6b <= 0; p_s5a <= 0; p_s5b <= 0;
        t0 <= 0; t3 <= 0; s7 <= 0; s6 <= 0; s5q <= 0; t0f <= 0; t7 <= 0; s4q <= 0;
        t4 <= 0; t5 <= 0; t6 <= 0; t65 <= 0; t56 <= 0; t1 <= 0; t2 <= 0; s5 <= 0;
    end
    else
    begin
        case (inport_idx_i)
        3'd0:  // even set
        begin
            i0  <= X0 + X2;
            q6  <= mc6(X1);    // for s2 at beat 5
            q2  <= mc2(X3);
            m1a <= mc2(X1);    // for t3 (s3) at beat 1
            m7b <= mc6(X3);
        end
        3'd1:  // odd set
        begin
            t3    <= m1a + m7b;        // s3 = C2*c2 + C6*c6 (reads beat0 values)
            m1a   <= mc1(X0);          // C1*c1
            m7b   <= mc7(X3);          // C7*c7
            m7a   <= mc7(X0);          // C7*c1
            m1b   <= mc1(X3);          // C1*c7
            p_s6a <= mc5(X2);          // C5*c5
            p_s6b <= mc3(X1);          // C3*c3
            p_s5a <= mc3(X2);          // C3*c5
            p_s5b <= mc5(X1);          // C5*c3
            t0    <= mc4(i0);          // s0
        end
        3'd2:  // even set
        begin
            i0  <= X0 - X2;            // i1
            t0f <= t0 + t3;
            t3  <= t0 - t3;
            s7  <= m1a + m7b;
            s6  <= p_s6a + p_s6b;      // butterfly s6
            s5q <= p_s5a - p_s5b;      // butterfly s5
        end
        3'd3:  // odd set
        begin
            t7  <= s6 + s7;
            s4q <= m7a - m1b;          // s4
        end
        3'd4:  // odd set
        begin
            t4 <= s4q + s5q;           // s4 + s5
            t5 <= s4q - s5q;           // s4 - s5
            t6 <= s7 - s6;
        end
        3'd5:  // even set
        begin
            t3  <= q6 - q2;            // s2 = C6*c2 - C2*c6 (beat0 products)
            t4  <= mc4(i0);            // s1
            t65 <= t6 - t5;
            t56 <= t5 + t6;
            q6  <= partial181(t6);     // reused for 181* split
            q2  <= partial181(t5);
        end
        3'd6:  // odd set
        begin
            t1  <= t4 + t3;
            t2  <= t4 - t3;
            m1a <= (t65 << 8) - (q6 - q2);  // 181*t65 mod 2^32
            m7b <= (t56 << 8) - (q6 + q2);  // 181*t56 mod 2^32
        end
        default: // 3'd7
        begin
            s5 <= fix181(m1a);         // (t65*181)/256
            s6 <= fix181(m7b);         // (t56*181)/256
        end
        endcase
    end
end

// ---- strobes ----
wire b3 = inport_valid_i && (inport_idx_i == 3'd3);
wire b7 = inport_valid_i && (inport_idx_i == 3'd7);
reg b3_q, b3_qq, b7_q;
reg [2:0] row_q;
always @(posedge clk)
begin
    if (rst)
    begin
        b3_q <= 1'b0; b3_qq <= 1'b0; b7_q <= 1'b0; row_q <= 3'd0;
    end
    else
    begin
        b3_q  <= b3;
        b3_qq <= b3_q;
        b7_q  <= b7;
        if (img_start_i)
            row_q <= 3'd0;
        else if (b7_q)
            row_q <= row_q + 3'd1;
    end
end

wire signed [31:0] r0 = (t0f + t7) >>> OUT_SHIFT;
wire signed [31:0] r7 = (t0f - t7) >>> OUT_SHIFT;
wire signed [31:0] r3 = (t3 + t4) >>> OUT_SHIFT;
wire signed [31:0] r4 = (t3 - t4) >>> OUT_SHIFT;
wire signed [31:0] r2 = (t2 + s5) >>> OUT_SHIFT;
wire signed [31:0] r5 = (t2 - s5) >>> OUT_SHIFT;
wire signed [31:0] r1 = (t1 + s6) >>> OUT_SHIFT;
wire signed [31:0] r6 = (t1 - s6) >>> OUT_SHIFT;

// ---- direct write (x pass) ----
assign dw_stro_a_o = b3_qq && dw_accept_i;
assign dw_stro_c_o = b7_q && dw_accept_i;
assign dw_row_o    = row_q;
assign dw_d0_o = r0[19:0];
assign dw_d3_o = r3[19:0];
assign dw_d4_o = r4[19:0];
assign dw_d7_o = r7[19:0];
assign dw_d1_o = r1[19:0];
assign dw_d2_o = r2[19:0];
assign dw_d5_o = r5[19:0];
assign dw_d6_o = r6[19:0];

// ---- streamed output (y pass) ----
reg signed [BO_BITS-1:0] bo0, bo1, bo2, bo3, bo4, bo5, bo6, bo7, bo_tmp;
always @(posedge clk)
begin
    if (rst)
    begin
        bo0 <= 0; bo1 <= 0; bo2 <= 0; bo3 <= 0;
        bo4 <= 0; bo5 <= 0; bo6 <= 0; bo7 <= 0;
        bo_tmp <= 0;
    end
    else
    begin
        if (b3_qq)
        begin
            bo0 <= r0[BO_BITS-1:0];
            bo_tmp <= r7[BO_BITS-1:0];
            bo3 <= r3[BO_BITS-1:0];
            bo4 <= r4[BO_BITS-1:0];
        end
        if (b7_q)
        begin
            bo7 <= bo_tmp;
            bo2 <= r2[BO_BITS-1:0];
            bo5 <= r5[BO_BITS-1:0];
            bo1 <= r1[BO_BITS-1:0];
            bo6 <= r6[BO_BITS-1:0];
        end
    end
end

reg [7:0] vdly;
always @(posedge clk)
begin
    if (rst || img_start_i)
        vdly <= 8'b0;
    else
        vdly <= {vdly[6:0], inport_valid_i};
end

reg [5:0] ptr;
always @(posedge clk)
begin
    if (rst || img_start_i)
        ptr <= 6'd0;
    else if (outport_valid_o)
        ptr <= ptr + 6'd1;
end

reg [BO_BITS-1:0] bo_rd;
always @*
begin
    case (ptr[2:0])
    3'd0: bo_rd = bo0;
    3'd1: bo_rd = bo1;
    3'd2: bo_rd = bo2;
    3'd3: bo_rd = bo3;
    3'd4: bo_rd = bo4;
    3'd5: bo_rd = bo5;
    3'd6: bo_rd = bo6;
    default: bo_rd = bo7;
    endcase
end

assign outport_valid_o = vdly[7];
generate
if (BO_BITS < 32)
    assign outport_data_o = {{(32-BO_BITS){bo_rd[BO_BITS-1]}}, bo_rd};
else
    assign outport_data_o = bo_rd;
endgenerate
generate
if (TRANSPOSE_IDX)
    assign outport_idx_o = {ptr[2:0], ptr[5:3]};
else
    assign outport_idx_o = ptr;
endgenerate

endmodule

//-----------------------------------------------------------------
// Transpose: double-buffered corner turn with direct writes from x pass.
//-----------------------------------------------------------------
module idct_transpose_direct
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           wstro_a_i
    ,input           wstro_c_i
    ,input  [ 2:0]   wrow_i
    ,input  [19:0]   wd0_i, wd3_i, wd4_i, wd7_i
    ,input  [19:0]   wd1_i, wd2_i, wd5_i, wd6_i
    ,output          accept_o
    ,output          outport_valid_o
    ,output [31:0]   outport_data0_o
    ,output [31:0]   outport_data1_o
    ,output [31:0]   outport_data2_o
    ,output [31:0]   outport_data3_o
    ,output [ 2:0]   outport_idx_o
);

localparam DW = 20;

reg         block_wr_q;
reg         block_rd_q;
reg [1:0]   block_ready_q;

assign accept_o = ~block_ready_q[block_wr_q];

localparam ST_IDLE   = 2'd0;
localparam ST_SETUP  = 2'd1;
localparam ST_ACTIVE = 2'd2;

reg [1:0] state_q;
reg [5:0] rcnt_q;

wire active_w  = (state_q == ST_ACTIVE);
wire rd_last_w = active_w && (rcnt_q == 6'd63);
wire blk_done_w = wstro_c_i && (wrow_i == 3'd7);

always @(posedge clk)
begin
    if (rst || img_start_i)
    begin
        block_ready_q <= 2'b0;
        block_wr_q    <= 1'b0;
        block_rd_q    <= 1'b0;
    end
    else
    begin
        if (blk_done_w)
        begin
            block_ready_q[block_wr_q] <= 1'b1;
            block_wr_q                <= ~block_wr_q;
        end
        if (rd_last_w)
        begin
            block_ready_q[block_rd_q] <= 1'b0;
            block_rd_q                <= ~block_rd_q;
        end
    end
end

always @(posedge clk)
begin
    if (rst || img_start_i)
    begin
        state_q <= ST_IDLE;
        rcnt_q  <= 6'd0;
    end
    else
    begin
        case (state_q)
        ST_IDLE:
            if (block_ready_q[block_rd_q])
                state_q <= ST_SETUP;
        ST_SETUP:
            state_q <= ST_ACTIVE;
        ST_ACTIVE:
        begin
            rcnt_q <= rcnt_q + 6'd1;
            if (rcnt_q == 6'd63)
                state_q <= ST_IDLE;
        end
        default: ;
        endcase
    end
end

// data sequence per column (y beats 0..7): lo, hi, lo, hi, hi, lo, hi, lo[next]
wire [2:0] col_w   = rcnt_q[5:3];
wire [2:0] phase_w = rcnt_q[2:0];

reg        rhi_w;
reg [2:0]  rcol_w;
always @*
begin
    rcol_w = col_w;
    case (phase_w)
    3'd0: rhi_w = 1'b1;
    3'd1: rhi_w = 1'b0;
    3'd2: rhi_w = 1'b1;
    3'd3: rhi_w = 1'b1;
    3'd4: rhi_w = 1'b0;
    3'd5: rhi_w = 1'b1;
    3'd6:
    begin
        rhi_w  = 1'b0;
        rcol_w = col_w + 3'd1;
    end
    default:
    begin
        rhi_w  = 1'b0;
        rcol_w = col_w + 3'd1;
    end
    endcase
end

wire [4:0] raddr_w = active_w ? {block_rd_q, rhi_w, rcol_w}
                              : {block_rd_q, 1'b0, 3'd0};

reg [DW-1:0] mem0 [0:31];
reg [DW-1:0] mem1 [0:31];
reg [DW-1:0] mem2 [0:31];
reg [DW-1:0] mem3 [0:31];

reg [DW-1:0] rdq0, rdq1, rdq2, rdq3;

wire bank_sel0 = (wrow_i[2:1] == 2'd0);
wire bank_sel1 = (wrow_i[2:1] == 2'd1);
wire bank_sel2 = (wrow_i[2:1] == 2'd2);
wire bank_sel3 = (wrow_i[2:1] == 2'd3);
wire parity_w  = wrow_i[0];

genvar gi;
generate
for (gi = 0; gi < 32; gi = gi + 1)
begin : gen_wr
    localparam [2:0] COL = gi[2:0];
    wire is_a_col = (COL == 3'd0) || (COL == 3'd3) || (COL == 3'd4) || (COL == 3'd7);
    wire blk_match = (gi[4] == block_wr_q);
    wire we0 = ((is_a_col ? wstro_a_i : wstro_c_i)) && bank_sel0 && (gi[3] == parity_w) && blk_match;
    wire we1 = ((is_a_col ? wstro_a_i : wstro_c_i)) && bank_sel1 && (gi[3] == parity_w) && blk_match;
    wire we2 = ((is_a_col ? wstro_a_i : wstro_c_i)) && bank_sel2 && (gi[3] == parity_w) && blk_match;
    wire we3 = ((is_a_col ? wstro_a_i : wstro_c_i)) && bank_sel3 && (gi[3] == parity_w) && blk_match;

    wire [DW-1:0] wdata =
        (COL == 3'd0) ? wd0_i :
        (COL == 3'd1) ? wd1_i :
        (COL == 3'd2) ? wd2_i :
        (COL == 3'd3) ? wd3_i :
        (COL == 3'd4) ? wd4_i :
        (COL == 3'd5) ? wd5_i :
        (COL == 3'd6) ? wd6_i : wd7_i;

    always @(posedge clk)
    begin
        if (we0) mem0[gi] <= wdata;
        if (we1) mem1[gi] <= wdata;
        if (we2) mem2[gi] <= wdata;
        if (we3) mem3[gi] <= wdata;
    end
end
endgenerate

always @(posedge clk)
begin
    rdq0 <= mem0[raddr_w];
    rdq1 <= mem1[raddr_w];
    rdq2 <= mem2[raddr_w];
    rdq3 <= mem3[raddr_w];
end

assign outport_valid_o = active_w;
assign outport_idx_o   = phase_w;
assign outport_data0_o = {{(32-DW){rdq0[DW-1]}}, rdq0};
assign outport_data1_o = {{(32-DW){rdq1[DW-1]}}, rdq1};
assign outport_data2_o = {{(32-DW){rdq2[DW-1]}}, rdq2};
assign outport_data3_o = {{(32-DW){rdq3[DW-1]}}, rdq3};

endmodule

//-----------------------------------------------------------------
// Top
//-----------------------------------------------------------------
module jpeg_idct_2d
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           img_end_i
    ,input           inport_valid_i
    ,input  [15:0]   inport_data0_i
    ,input  [15:0]   inport_data1_i
    ,input  [15:0]   inport_data2_i
    ,input  [15:0]   inport_data3_i
    ,input  [ 2:0]   inport_idx_i
    ,output          outport_valid_o
    ,output [31:0]   outport_data_o
    ,output [ 5:0]   outport_idx_o
);

/* verilator lint_off UNUSED */
wire unused_ok = img_end_i;
/* verilator lint_on UNUSED */

wire        accept_w;
wire        wstro_a_w;
wire        wstro_c_w;
wire [ 2:0] wrow_w;
wire [19:0] wd0_w, wd1_w, wd2_w, wd3_w, wd4_w, wd5_w, wd6_w, wd7_w;

wire        x_dummy_valid;
wire [31:0] x_dummy_data;
wire [ 5:0] x_dummy_idx;

idct_pass6 #(.OUT_SHIFT(11), .IN_BITS(16), .BO_BITS(20), .TRANSPOSE_IDX(0))
u_idct_x
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(inport_valid_i)
    ,.inport_data0_i({{16{inport_data0_i[15]}}, inport_data0_i})
    ,.inport_data1_i({{16{inport_data1_i[15]}}, inport_data1_i})
    ,.inport_data2_i({{16{inport_data2_i[15]}}, inport_data2_i})
    ,.inport_data3_i({{16{inport_data3_i[15]}}, inport_data3_i})
    ,.inport_idx_i(inport_idx_i)
    ,.dw_stro_a_o(wstro_a_w)
    ,.dw_stro_c_o(wstro_c_w)
    ,.dw_row_o(wrow_w)
    ,.dw_d0_o(wd0_w)
    ,.dw_d3_o(wd3_w)
    ,.dw_d4_o(wd4_w)
    ,.dw_d7_o(wd7_w)
    ,.dw_d1_o(wd1_w)
    ,.dw_d2_o(wd2_w)
    ,.dw_d5_o(wd5_w)
    ,.dw_d6_o(wd6_w)
    ,.dw_accept_i(accept_w)
    ,.outport_valid_o(x_dummy_valid)
    ,.outport_data_o(x_dummy_data)
    ,.outport_idx_o(x_dummy_idx)
);

wire        t_valid_w;
wire [31:0] t_data0_w;
wire [31:0] t_data1_w;
wire [31:0] t_data2_w;
wire [31:0] t_data3_w;
wire [ 2:0] t_idx_w;

idct_transpose_direct u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.wstro_a_i(wstro_a_w)
    ,.wstro_c_i(wstro_c_w)
    ,.wrow_i(wrow_w)
    ,.wd0_i(wd0_w)
    ,.wd3_i(wd3_w)
    ,.wd4_i(wd4_w)
    ,.wd7_i(wd7_w)
    ,.wd1_i(wd1_w)
    ,.wd2_i(wd2_w)
    ,.wd5_i(wd5_w)
    ,.wd6_i(wd6_w)
    ,.accept_o(accept_w)
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
);

wire        y_dummy_a, y_dummy_c;
wire [ 2:0] y_dummy_row;
wire [19:0] y_dummy_d0, y_dummy_d1, y_dummy_d2, y_dummy_d3;
wire [19:0] y_dummy_d4, y_dummy_d5, y_dummy_d6, y_dummy_d7;

idct_pass6 #(.OUT_SHIFT(15), .IN_BITS(20), .BO_BITS(17), .TRANSPOSE_IDX(1))
u_idct_y
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(t_valid_w)
    ,.inport_data0_i(t_data0_w)
    ,.inport_data1_i(t_data1_w)
    ,.inport_data2_i(t_data2_w)
    ,.inport_data3_i(t_data3_w)
    ,.inport_idx_i(t_idx_w)
    ,.dw_stro_a_o(y_dummy_a)
    ,.dw_stro_c_o(y_dummy_c)
    ,.dw_row_o(y_dummy_row)
    ,.dw_d0_o(y_dummy_d0)
    ,.dw_d3_o(y_dummy_d3)
    ,.dw_d4_o(y_dummy_d4)
    ,.dw_d7_o(y_dummy_d7)
    ,.dw_d1_o(y_dummy_d1)
    ,.dw_d2_o(y_dummy_d2)
    ,.dw_d5_o(y_dummy_d5)
    ,.dw_d6_o(y_dummy_d6)
    ,.dw_accept_i(1'b1)
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

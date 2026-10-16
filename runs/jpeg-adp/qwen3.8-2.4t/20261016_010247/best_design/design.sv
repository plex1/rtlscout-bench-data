//-----------------------------------------------------------------
// jpeg_idct_2d — optimized 8x8 2-D IDCT (bit-exact with core_jpeg original)
//  - AAN butterfly with CSD constant multipliers (balanced adder trees)
//  - compressed pipeline: products formed combinationally at beat time
//  - transpose: 20-bit storage (x-pass outputs bounded by +/-346400 < 2^19)
//-----------------------------------------------------------------

module idct_pass
#(
     parameter OUT_SHIFT = 11
    ,parameter IN_BITS   = 16   // width of sign-extended input values
    ,parameter BO_BITS   = 20   // storage width of shifted outputs
    ,parameter TRANSPOSE_IDX = 0 // 1 for y pass: idx = {ptr[2:0],ptr[5:3]}
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
    ,output          outport_valid_o
    ,output [31:0]   outport_data_o
    ,output [ 5:0]   outport_idx_o
);

// Constant-multiply building blocks, all 32-bit wrapping, max depth 2:
// C1*x = ((x<<12)+(x)) - ((x<<6)+(x<<4))
// C2*x = ((x<<12)+(x<<3)) - ((x<<8)+(x<<6))
// C5*x = ((x<<11)+(x<<2)) + ((x<<8)-(x<<5))
// C6*x = ((x<<10)-(x)) + ((x<<9)+(x<<5))
// C7*x = ((x<<10)+(x<<5)) - ((x<<8)+(x))
// partial3(x) = ((x<<12)+(x<<6)) + ((x<<4)-(x<<1)) ; C3*x = partial3(x) - ((x<<9)+(x<<8))
// partial4(x) = ((x<<12)-(x<<10)) + ((x<<4)-(x<<6)); C4*x = partial4(x) - (x<<7)
// partial181(x) = ((x<<6)+(x<<3)) + ((x<<1)+(x)); 181*x = (x<<8) - partial181(x)
function automatic signed [31:0] mc1(input signed [31:0] x);
begin mc1 = ((x << 12) + x) - ((x << 6) + (x << 4)); end endfunction
function automatic signed [31:0] mc2(input signed [31:0] x);
begin mc2 = ((x << 12) + (x << 3)) - ((x << 8) + (x << 6)); end endfunction
function automatic signed [31:0] mc5(input signed [31:0] x);
begin mc5 = ((x << 11) + (x << 2)) + ((x << 8) - (x << 5)); end endfunction
function automatic signed [31:0] mc6(input signed [31:0] x);
begin mc6 = ((x << 10) - x) + ((x << 9) + (x << 5)); end endfunction
function automatic signed [31:0] mc7(input signed [31:0] x);
begin mc7 = ((x << 10) + (x << 5)) - ((x << 8) + x); end endfunction
function automatic signed [31:0] partial3(input signed [31:0] x);
begin partial3 = ((x << 12) + (x << 6)) + ((x << 4) - (x << 1)); end endfunction
function automatic signed [31:0] partial4(input signed [31:0] x);
begin partial4 = ((x << 12) - (x << 10)) + ((x << 4) - (x << 6)); end endfunction
function automatic signed [31:0] partial181(input signed [31:0] x);
begin partial181 = ((x << 6) + (x << 3)) + ((x << 1) + x); end endfunction

// (p = x*181 mod 2^32) -> p/256, signed division truncating toward zero
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

// ---- state registers (shared across the 8-beat schedule) ----
reg signed [31:0] i0, i1;      // c0+c4, c0-c4
reg signed [31:0] m2, m6;      // C2*c2, C6*c6 (later C6*c2, C2*c6)
reg signed [31:0] w0, w1;      // partial4(i0), partial4(i1)
reg signed [31:0] s0, s3, s7;  // even sums
reg signed [31:0] p1a, p7b, p7a, p1b; // C1*c1, C7*c7, C7*c1, C1*c7
reg signed [31:0] p5, p5b, p3, p3b;   // C5*c5, C5*c3, partial3(c3), partial3(c5)
reg signed [31:0] s4, s5, s6;
reg signed [31:0] t0f, t3f, t7, t4, t5, t65, t56;
reg signed [31:0] s1, n6, n2, q5, q6;
reg signed [31:0] t1, t2, m5p, m6p;
reg signed [31:0] g5, g6;

always @(posedge clk)
begin
    if (rst)
    begin
        i0 <= 0; i1 <= 0; m2 <= 0; m6 <= 0; w0 <= 0; w1 <= 0;
        s0 <= 0; s3 <= 0; s7 <= 0;
        p1a <= 0; p7b <= 0; p7a <= 0; p1b <= 0;
        p5 <= 0; p5b <= 0; p3 <= 0; p3b <= 0;
        s4 <= 0; s5 <= 0; s6 <= 0;
        t0f <= 0; t3f <= 0; t7 <= 0; t4 <= 0; t5 <= 0; t65 <= 0; t56 <= 0;
        s1 <= 0; n6 <= 0; n2 <= 0; q5 <= 0; q6 <= 0;
        t1 <= 0; t2 <= 0; m5p <= 0; m6p <= 0;
        g5 <= 0; g6 <= 0;
    end
    else
    begin
        case (inport_idx_i)
        3'd0:  // even set (c0,c2,c4,c6)
        begin
            i0 <= X0 + X2;
            m2 <= mc2(X1);
            m6 <= mc6(X3);
        end
        3'd1:  // odd set (c1,c3,c5,c7): all depth-2 partials/products
        begin
            p1a <= mc1(X0);        // C1*c1
            p7b <= mc7(X3);        // C7*c7
            p7a <= mc7(X0);        // C7*c1
            p1b <= mc1(X3);        // C1*c7
            p5  <= mc5(X2);        // C5*c5
            p5b <= mc5(X1);        // C5*c3
            p3  <= partial3(X1);   // C3*c3 + pair(c3)
            p3b <= partial3(X2);   // C3*c5 + pair(c5)
            w0  <= partial4(i0);   // s0 + (i0<<7)
        end
        3'd2:  // even set
        begin
            i1 <= X0 - X2;
            s0 <= w0 - (i0 << 7);  // C4*(c0+c4)
            s3 <= m2 + m6;
            s7 <= p1a + p7b;
        end
        3'd3:  // odd set
        begin
            t0f <= s0 + s3;
            t3f <= s0 - s3;
            s6  <= (p5 + p3) - ((X1 << 9) + (X1 << 8));      // C5*c5 + C3*c3
            s4  <= p7a - p1b;                                 // C7*c1 - C1*c7
            s5  <= (p3b - p5b) - ((X2 << 9) + (X2 << 8));     // C3*c5 - C5*c3
            w1  <= partial4(i1);   // s1 + (i1<<7)
        end
        3'd4:  // odd set (repeated)
        begin
            t7  <= s6 + s7;
            t4  <= s4 + s5;
            t5  <= s4 - s5;
            t65 <= (s7 - s6) - (s4 - s5);
            t56 <= (s7 - s6) + (s4 - s5);
        end
        3'd5:  // even set
        begin
            s1 <= w1 - (i1 << 7);  // C4*(c0-c4)
            n6 <= mc6(X1);         // C6*c2
            n2 <= mc2(X3);         // C2*c6
            q5 <= partial181(t65);
            q6 <= partial181(t56);
        end
        3'd6:  // odd set
        begin
            t1  <= s1 + n6 - n2;   // s1 + s2
            t2  <= s1 - n6 + n2;   // s1 - s2
            m5p <= (t65 << 8) - q5; // 181*t65 mod 2^32
            m6p <= (t56 << 8) - q6; // 181*t56 mod 2^32
        end
        default: // 3'd7
        begin
            g5 <= fix181(m5p);     // (t65*181)/256
            g6 <= fix181(m6p);     // (t56*181)/256
        end
        endcase
    end
end

// ---- block_out write strobes (fixed schedule timing) ----
wire b4 = inport_valid_i && (inport_idx_i == 3'd4);
wire b7 = inport_valid_i && (inport_idx_i == 3'd7);
reg b4_q, b7_q;
always @(posedge clk)
begin
    if (rst)
    begin
        b4_q <= 1'b0; b7_q <= 1'b0;
    end
    else
    begin
        b4_q <= b4;   // capture group A (t0f,t7,t3f,t4 ready) during next cycle
        b7_q <= b7;   // capture group B (t1,t2,g5,g6 ready) during next cycle
    end
end

reg signed [BO_BITS-1:0] bo0, bo1, bo2, bo3, bo4, bo5, bo6, bo7, bo_tmp;

wire signed [31:0] r0 = (t0f + t7) >>> OUT_SHIFT;
wire signed [31:0] r7 = (t0f - t7) >>> OUT_SHIFT;
wire signed [31:0] r3 = (t3f + t4) >>> OUT_SHIFT;
wire signed [31:0] r4 = (t3f - t4) >>> OUT_SHIFT;
wire signed [31:0] r2 = (t2 + g5) >>> OUT_SHIFT;
wire signed [31:0] r5 = (t2 - g5) >>> OUT_SHIFT;
wire signed [31:0] r1 = (t1 + g6) >>> OUT_SHIFT;
wire signed [31:0] r6 = (t1 - g6) >>> OUT_SHIFT;

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
        if (b4_q)
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

// ---- output streamer: valid = input valid delayed 8 cycles ----
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
    assign outport_idx_o = {ptr[2:0], ptr[5:3]};   // y-pass un-transpose
else
    assign outport_idx_o = ptr;
endgenerate

endmodule

//-----------------------------------------------------------------
// Transpose: double-buffered corner turn.
// 4 banks x 32 x 20 bits. Write addr {blk, idx[3:0]}; each RAM holds rows
// {2k, 2k+1} at lo/hi bit of the address. Read schedule fixed per column.
//-----------------------------------------------------------------
module idct_transpose
(
     input           clk
    ,input           rst
    ,input           img_start_i
    ,input           inport_valid_i
    ,input  [31:0]   inport_data_i
    ,input  [ 5:0]   inport_idx_i
    ,output          outport_valid_o
    ,output [31:0]   outport_data0_o
    ,output [31:0]   outport_data1_o
    ,output [31:0]   outport_data2_o
    ,output [31:0]   outport_data3_o
    ,output [ 2:0]   outport_idx_o
);

localparam DW = 20; // x-pass outputs bounded by +/-346400 (< 2^19)

reg         block_wr_q;
reg         block_rd_q;
reg [1:0]   block_ready_q;

wire inport_accept_w = ~block_ready_q[block_wr_q];

localparam ST_IDLE   = 2'd0;
localparam ST_SETUP  = 2'd1;
localparam ST_ACTIVE = 2'd2;

reg [1:0] state_q;
reg [5:0] rcnt_q;

wire active_w  = (state_q == ST_ACTIVE);
wire rd_last_w = active_w && (rcnt_q == 6'd63);

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
        if (inport_valid_i && (inport_idx_i == 6'd63) && inport_accept_w)
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

// ---- read address generation ----
// data sequence per column (y beats 0..7): lo, hi, lo, hi, hi, lo, hi, lo[next]
// raddr presented during cycle k is registered -> provides data for cycle k+1
wire [2:0] col_w   = rcnt_q[5:3];
wire [2:0] phase_w = rcnt_q[2:0];

reg        rhi_w;
reg [2:0]  rcol_w;
always @*
begin
    rcol_w = col_w;
    case (phase_w)
    3'd0: rhi_w = 1'b1;                    // next need: hi[col]
    3'd1: rhi_w = 1'b0;                    // next need: lo[col]
    3'd2: rhi_w = 1'b1;                    // next need: hi[col]
    3'd3: rhi_w = 1'b1;                    // next need: hi[col]
    3'd4: rhi_w = 1'b0;                    // next need: lo[col]
    3'd5: rhi_w = 1'b1;                    // next need: hi[col]
    3'd6:                                  // next need: lo[col+1]
    begin
        rhi_w  = 1'b0;
        rcol_w = col_w + 3'd1;
    end
    default:                               // next need: lo[col+1] (first beat of next col)
    begin
        rhi_w  = 1'b0;
        rcol_w = col_w + 3'd1;
    end
    endcase
end

wire [4:0] raddr_w = active_w ? {block_rd_q, rhi_w, rcol_w}
                              : {block_rd_q, 1'b0, 3'd0};

// ---- storage ----
wire [4:0] waddr_w = {block_wr_q, inport_idx_i[3:0]};
wire       wr_en_w = inport_valid_i && inport_accept_w;
wire       wr0_w   = wr_en_w && (inport_idx_i[5:4] == 2'd0);
wire       wr1_w   = wr_en_w && (inport_idx_i[5:4] == 2'd1);
wire       wr2_w   = wr_en_w && (inport_idx_i[5:4] == 2'd2);
wire       wr3_w   = wr_en_w && (inport_idx_i[5:4] == 2'd3);

reg [DW-1:0] mem0 [0:31];
reg [DW-1:0] mem1 [0:31];
reg [DW-1:0] mem2 [0:31];
reg [DW-1:0] mem3 [0:31];

reg [DW-1:0] rdq0, rdq1, rdq2, rdq3;

always @(posedge clk)
begin
    if (wr0_w) mem0[waddr_w] <= inport_data_i[DW-1:0];
    if (wr1_w) mem1[waddr_w] <= inport_data_i[DW-1:0];
    if (wr2_w) mem2[waddr_w] <= inport_data_i[DW-1:0];
    if (wr3_w) mem3[waddr_w] <= inport_data_i[DW-1:0];

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

wire        x_valid_w;
wire [31:0] x_data_w;
wire [ 5:0] x_idx_w;

idct_pass #(.OUT_SHIFT(11), .IN_BITS(16), .BO_BITS(20), .TRANSPOSE_IDX(0))
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
    ,.outport_valid_o(x_valid_w)
    ,.outport_data_o(x_data_w)
    ,.outport_idx_o(x_idx_w)
);

wire        t_valid_w;
wire [31:0] t_data0_w;
wire [31:0] t_data1_w;
wire [31:0] t_data2_w;
wire [31:0] t_data3_w;
wire [ 2:0] t_idx_w;

idct_transpose u_transpose
(
     .clk(clk)
    ,.rst(rst)
    ,.img_start_i(img_start_i)
    ,.inport_valid_i(x_valid_w)
    ,.inport_data_i(x_data_w)
    ,.inport_idx_i(x_idx_w)
    ,.outport_valid_o(t_valid_w)
    ,.outport_data0_o(t_data0_w)
    ,.outport_data1_o(t_data1_w)
    ,.outport_data2_o(t_data2_w)
    ,.outport_data3_o(t_data3_w)
    ,.outport_idx_o(t_idx_w)
);

idct_pass #(.OUT_SHIFT(15), .IN_BITS(20), .BO_BITS(32), .TRANSPOSE_IDX(1))
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
    ,.outport_valid_o(outport_valid_o)
    ,.outport_data_o(outport_data_o)
    ,.outport_idx_o(outport_idx_o)
);

endmodule

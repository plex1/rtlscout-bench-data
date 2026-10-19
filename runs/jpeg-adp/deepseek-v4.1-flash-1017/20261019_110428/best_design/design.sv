//-----------------------------------------------------------------
// jpeg_idct_2d  —  8x8 2-D IDCT (baseline JPEG)
//
// Bit-exact re-implementation of the ultraembedded/core_jpeg
// jpeg_idct_x -> jpeg_idct_transpose -> jpeg_idct_y composition.
//
// 1-D core (identical in both passes, only the widths/shift differ):
//   s0=(x0+x4)*C4    s1=(x0-x4)*C4
//   s2=x2*C6-x6*C2   s3=x2*C2+x6*C6
//   s4=x1*C7-x7*C1   s5=x5*C3-x3*C5
//   s6=x5*C5+x3*C3   s7=x1*C1+x7*C7
//   t0=s0+s3 t3=s0-s3 t1=s1+s2 t2=s1-s2
//   t7=s6+s7 t6=s7-s6 t4=s4+s5 t5=s4-s5
//   q5=(t6-t5)*181/256   q6=(t5+t6)*181/256
//   w0=t0+t7 w7=t0-t7 w3=t3+t4 w4=t3-t4
//   w2=t2+q5 w5=t2-q5 w1=t1+q6 w6=t1-q6
//   out = w >>> SH
//
// Pipeline (one register stage per cycle, so each stage is at most a
// constant multiply or a single adder level):
//   ph1 latch round-A operands        ph2 round-A products, latch B ops
//   ph3 round-B products              ph4 t values  (2-level product tree)
//   ph5 dv=t6-t5, t5+t6               ph6 q values
//   ph7 w values (written straight into the destination register file)
//
// Optimisations vs. the reference RTL:
//  * every multiply is by a *constant* (C1..C7, plus 181 for the 1/sqrt(2)
//    scale) instead of five fully general 32x32 multipliers per pass
//  * the round select uses register enables, never a datapath mux
//  * the 4x32x32 double-buffered corner-turn RAM is replaced by a 16x8x20
//    register file; a column is read with a registered one-hot AND-OR
//    select (a binary mux makes one flop drive hundreds of gates)
//-----------------------------------------------------------------
module idct_eng #(
    parameter integer IW = 16,   // signed input width
    parameter integer DW = 31,   // datapath width (signed, wrapping)
    parameter integer SH = 11    // output arithmetic shift
)(
    input  wire                 clk,
    input  wire [2:0]           ph,
    input  wire signed [IW-1:0] x0, x1, x2, x3,
    input  wire signed [IW-1:0] x4, x5, x6, x7,
    output wire signed [DW-1:0] w0, w1, w2, w3,
    output wire signed [DW-1:0] w4, w5, w6, w7
);

localparam signed [12:0] C1 = 13'sd4017;
localparam signed [12:0] C2 = 13'sd3784;
localparam signed [12:0] C3 = 13'sd3406;
localparam signed [12:0] C4 = 13'sd2896;
localparam signed [12:0] C5 = 13'sd2276;
localparam signed [12:0] C6 = 13'sd1567;
localparam signed [12:0] C7 = 13'sd 799;

localparam integer PW = IW + 14;     // exact product width

//---- round operands : one register set, reloaded per round -------------
reg signed [IW:0]   om4;
reg signed [IW-1:0] om6, om2, om7, om1, om3, om5;

wire signed [IW:0] sum4 = {x0[IW-1], x0} + {x4[IW-1], x4};
wire signed [IW:0] dif4 = {x0[IW-1], x0} - {x4[IW-1], x4};

wire signed [PW-1:0] p4 = om4 * C4;
wire signed [PW-1:0] p6 = om6 * C6;
wire signed [PW-1:0] p2 = om2 * C2;
wire signed [PW-1:0] p7 = om7 * C7;
wire signed [PW-1:0] p1 = om1 * C1;
wire signed [PW-1:0] p3 = om3 * C3;
wire signed [PW-1:0] p5 = om5 * C5;

reg signed [PW-1:0] A4, A6, A2, A7, A1, A3, A5;
reg signed [PW-1:0] B4, B6, B2, B7, B1, B3, B5;

reg signed [DW-1:0] t0, t1, t2, t3, t4, t5, t6, t7;
reg signed [DW-1:0] dv5, dv6;
reg signed [31:0]   pv5, pv6;
reg signed [DW-1:0] r0, r1, r2, r3, r4, r5, r6, r7;

function signed [DW-1:0] ext;
    input signed [PW-1:0] v;
    reg   signed [DW-1:0] r;
    begin
        r = v;   // implicit resize keeps the low DW bits (32-bit wrap)
        ext = r;
    end
endfunction

// t values : two-level tree straight from the product registers
wire signed [DW-1:0] a3b6 = ext(B2) + ext(B6);   // s3
wire signed [DW-1:0] a2b6 = ext(A6) - ext(A2);   // s2
wire signed [DW-1:0] a7a1 = ext(A7) - ext(A1);   // s4
wire signed [DW-1:0] a3a5 = ext(A3) - ext(A5);   // s5
wire signed [DW-1:0] b5b3 = ext(B5) + ext(B3);   // s6
wire signed [DW-1:0] b1b7 = ext(B1) + ext(B7);   // s7

wire signed [DW-1:0] nt0 = ext(A4) + a3b6;
wire signed [DW-1:0] nt3 = ext(A4) - a3b6;
wire signed [DW-1:0] nt1 = ext(B4) + a2b6;
wire signed [DW-1:0] nt2 = ext(B4) - a2b6;
wire signed [DW-1:0] nt7 = b5b3 + b1b7;
wire signed [DW-1:0] nt6 = b1b7 - b5b3;
wire signed [DW-1:0] nt4 = a7a1 + a3a5;
wire signed [DW-1:0] nt5 = a7a1 - a3a5;

wire signed [DW-1:0] nd5 = t6 - t5;
wire signed [DW-1:0] nd6 = t5 + t6;

// (p*1)>>8 with the reference's truncate-toward-zero division
wire signed [DW-1:0] nq5 = (pv5 + (pv5[31] ? 32'sd255 : 32'sd0)) >>> 8;
wire signed [DW-1:0] nq6 = (pv6 + (pv6[31] ? 32'sd255 : 32'sd0)) >>> 8;

always @(posedge clk)
begin
    case (ph)
    3'd1: begin
        om4 <= sum4; om6 <= x2; om2 <= x6;
        om7 <= x1;   om1 <= x7; om3 <= x5; om5 <= x3;
    end
    3'd2: begin
        om4 <= dif4; om6 <= x6; om2 <= x2;
        om7 <= x7;   om1 <= x1; om3 <= x3; om5 <= x5;
        A4 <= p4; A6 <= p6; A2 <= p2; A7 <= p7;
        A1 <= p1; A3 <= p3; A5 <= p5;
    end
    3'd3: begin
        B4 <= p4; B6 <= p6; B2 <= p2; B7 <= p7;
        B1 <= p1; B3 <= p3; B5 <= p5;
    end
    3'd4: begin
        t0 <= nt0; t1 <= nt1; t2 <= nt2; t3 <= nt3;
        t4 <= nt4; t5 <= nt5; t6 <= nt6; t7 <= nt7;
    end
    3'd5: begin
        dv5 <= nd5;
        dv6 <= nd6;
    end
    3'd6: begin
        pv5 <= dv5 * 32'sd181;
        pv6 <= dv6 * 32'sd181;
    end
    default: ;
    endcase
end

assign w0 = (t0 + t7) >>> SH;
assign w7 = (t0 - t7) >>> SH;
assign w3 = (t3 + t4) >>> SH;
assign w4 = (t3 - t4) >>> SH;
assign w2 = (t2 + nq5) >>> SH;
assign w5 = (t2 - nq5) >>> SH;
assign w1 = (t1 + nq6) >>> SH;
assign w6 = (t1 - nq6) >>> SH;

endmodule


//=================================================================
// top level
//=================================================================
module jpeg_idct_2d
(
     input  wire        clk
    ,input  wire        rst
    ,input  wire        img_start_i
    ,input  wire        img_end_i
    ,input  wire        inport_valid_i
    ,input  wire [15:0] inport_data0_i
    ,input  wire [15:0] inport_data1_i
    ,input  wire [15:0] inport_data2_i
    ,input  wire [15:0] inport_data3_i
    ,input  wire [ 2:0] inport_idx_i
    ,output wire        outport_valid_o
    ,output wire [31:0] outport_data_o
    ,output wire [ 5:0] outport_idx_o
);

//-----------------------------------------------------------------
// input capture: idx==0 -> (c0,c2,c4,c6), idx==1 -> (c1,c3,c5,c7)
//-----------------------------------------------------------------
reg signed [15:0] rx0, rx1, rx2, rx3, rx4, rx5, rx6, rx7;

wire beat0 = inport_valid_i && (inport_idx_i == 3'd0);
wire beat1 = inport_valid_i && (inport_idx_i == 3'd1);

always @(posedge clk)
begin
    if (beat0) begin
        rx0 <= inport_data0_i;
        rx2 <= inport_data1_i;
        rx4 <= inport_data2_i;
        rx6 <= inport_data3_i;
    end
    if (beat1) begin
        rx1 <= inport_data0_i;
        rx3 <= inport_data1_i;
        rx5 <= inport_data2_i;
        rx7 <= inport_data3_i;
    end
end

//-----------------------------------------------------------------
// row pass : 8-cycle period per input row, w stored at phase 7
//-----------------------------------------------------------------
reg [2:0] rph;
reg [2:0] rcnt;
reg [2:0] wrow;
reg       wb;
reg       ready0, ready1;
reg       xact;

wire xwr = (rph == 3'd7) && xact;

wire signed [30:0] rw0, rw1, rw2, rw3, rw4, rw5, rw6, rw7;

idct_eng #(.IW(16), .DW(31), .SH(11)) u_row
(
    .clk(clk), .ph(rph),
    .x0(rx0), .x1(rx1), .x2(rx2), .x3(rx3),
    .x4(rx4), .x5(rx5), .x6(rx6), .x7(rx7),
    .w0(rw0), .w1(rw1), .w2(rw2), .w3(rw3),
    .w4(rw4), .w5(rw5), .w6(rw6), .w7(rw7)
);

//-----------------------------------------------------------------
// corner turn : tr[{bank,col}][row] = element (row,col) after row pass
//-----------------------------------------------------------------
reg [19:0] tr [0:15][0:7];

always @(posedge clk)
begin
    if (xwr) begin
        tr[{wb,3'd0}][wrow] <= rw0[19:0];
        tr[{wb,3'd1}][wrow] <= rw1[19:0];
        tr[{wb,3'd2}][wrow] <= rw2[19:0];
        tr[{wb,3'd3}][wrow] <= rw3[19:0];
        tr[{wb,3'd4}][wrow] <= rw4[19:0];
        tr[{wb,3'd5}][wrow] <= rw5[19:0];
        tr[{wb,3'd6}][wrow] <= rw6[19:0];
        tr[{wb,3'd7}][wrow] <= rw7[19:0];
    end
end

//-----------------------------------------------------------------
// column pass : 8-cycle period per column
//   cn=0 capture column, cn=1..6 pipeline, cn=7 store result + advance
//-----------------------------------------------------------------
reg [2:0]  cn;
reg        yrun;
reg [2:0]  lcol;
reg [2:0]  wcol;
reg        rb2;
reg        yv;
reg [5:0]  ocnt;
reg [15:0] soh;

wire rdy = rb2 ? ready1 : ready0;
wire sy  = ~yrun && rdy;
wire [2:0] lnext = (lcol == 3'd7) ? 3'd0 : (lcol + 3'd1);

always @(posedge clk)
begin
    if (rst || img_start_i)      soh <= 16'd1;
    else if (sy)                 soh <= 16'd1 << {rb2, 3'd0};
    else if (yrun && cn == 3'd7) soh <= 16'd1 << {rb2, lnext};
end

wire [19:0] rdv [0:7];

genvar gk, gi;
generate
for (gk = 0; gk < 8; gk = gk + 1) begin : g_rd
    wire [19:0] m0 [0:7];
    wire [19:0] m1 [0:3];
    for (gi = 0; gi < 8; gi = gi + 1) begin : g1
        assign m0[gi] = ({20{soh[2*gi]}}   & tr[2*gi][gk])
                      | ({20{soh[2*gi+1]}} & tr[2*gi+1][gk]);
    end
    for (gi = 0; gi < 4; gi = gi + 1) begin : g2
        assign m1[gi] = m0[2*gi] | m0[2*gi+1];
    end
    assign rdv[gk] = (m1[0] | m1[1]) | (m1[2] | m1[3]);
end
endgenerate

reg signed [19:0] cx0, cx1, cx2, cx3, cx4, cx5, cx6, cx7;

wire signed [31:0] yw0, yw1, yw2, yw3, yw4, yw5, yw6, yw7;

idct_eng #(.IW(20), .DW(32), .SH(15)) u_col
(
    .clk(clk), .ph(cn),
    .x0(cx0), .x1(cx1), .x2(cx2), .x3(cx3),
    .x4(cx4), .x5(cx5), .x6(cx6), .x7(cx7),
    .w0(yw0), .w1(yw1), .w2(yw2), .w3(yw3),
    .w4(yw4), .w5(yw5), .w6(yw6), .w7(yw7)
);

//-----------------------------------------------------------------
// output buffer (one column) + streamer
//-----------------------------------------------------------------
reg signed [31:0] obuf [0:7];

always @(posedge clk)
begin
    if (yrun && (cn == 3'd7)) begin
        obuf[0] <= yw0;  obuf[1] <= yw1;  obuf[2] <= yw2;  obuf[3] <= yw3;
        obuf[4] <= yw4;  obuf[5] <= yw5;  obuf[6] <= yw6;  obuf[7] <= yw7;
    end
end

always @(posedge clk)
begin
    if (rst || img_start_i) begin
        rph     <= 3'd0;
        rcnt    <= 3'd0;
        wrow    <= 3'd0;
        wb      <= 1'b0;
        ready0  <= 1'b0;
        ready1  <= 1'b0;
        xact    <= 1'b0;
        cn      <= 3'd0;
        yrun    <= 1'b0;
        lcol    <= 3'd0;
        wcol    <= 3'd0;
        rb2     <= 1'b0;
        yv      <= 1'b0;
        ocnt    <= 6'd0;
        cx0 <= 20'sd0; cx1 <= 20'sd0; cx2 <= 20'sd0; cx3 <= 20'sd0;
        cx4 <= 20'sd0; cx5 <= 20'sd0; cx6 <= 20'sd0; cx7 <= 20'sd0;
    end
    else begin
        //---- row pass ------------------------------------------
        if (beat1) rph <= 3'd1;
        else       rph <= rph + 3'd1;

        if (beat1) begin
            wrow <= rcnt;
            rcnt <= (rcnt == 3'd7) ? 3'd0 : rcnt + 3'd1;
            xact <= 1'b1;
        end

        if (xwr) begin
            xact <= 1'b0;
            if (wrow == 3'd7) begin
                if (wb) ready1 <= 1'b1;
                else    ready0 <= 1'b1;
                wb <= ~wb;
            end
        end

        //---- column pass ---------------------------------------
        if (!yrun) begin
            if (rdy) begin
                if (rb2) ready1 <= 1'b0;
                else     ready0 <= 1'b0;
                yrun <= 1'b1;
                cn   <= 3'd0;
                lcol <= 3'd0;
                wcol <= 3'd0;
            end
        end
        else begin
            cn <= cn + 3'd1;

            if (cn == 3'd0) begin
                cx0 <= rdv[0]; cx1 <= rdv[1]; cx2 <= rdv[2]; cx3 <= rdv[3];
                cx4 <= rdv[4]; cx5 <= rdv[5]; cx6 <= rdv[6]; cx7 <= rdv[7];
            end
            if (cn == 3'd7) begin
                if (wcol == 3'd0) begin
                    yv   <= 1'b1;
                    ocnt <= 6'd0;
                end
                if (wcol == 3'd7) begin
                    yrun <= 1'b0;
                    rb2  <= ~rb2;
                end
                else begin
                    wcol <= wcol + 3'd1;
                    lcol <= lnext;
                end
            end
        end

        //---- output streamer -----------------------------------
        if (yv) begin
            if (ocnt == 6'd63) yv <= 1'b0;
            else               ocnt <= ocnt + 6'd1;
        end
    end
end

assign outport_valid_o = yv;
assign outport_data_o  = obuf[ocnt[2:0]];
assign outport_idx_o   = {ocnt[2:0], ocnt[5:3]};

endmodule

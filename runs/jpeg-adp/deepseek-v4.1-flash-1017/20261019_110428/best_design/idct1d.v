//-----------------------------------------------------------------
// 1-D IDCT (factorized AAN-style, cos*4096) — bit-exact re-implementation
// of ultraembedded/core_jpeg jpeg_idct_x / jpeg_idct_y arithmetic.
//
// s0 = (x0+x4)*C4        s1 = (x0-x4)*C4
// s2 = x2*C6 - x6*C2     s3 = x2*C2 + x6*C6
// s4 = x0*C7 - x6*C1     s5 = x4*C3 - x2*C5
// s6 = x4*C5 + x2*C3     s7 = x0*C1 + x6*C7
// t0 = s0+s3   t3 = s0-s3   t1 = s1+s2   t2 = s1-s2
// t7 = s6+s7   t6 = s7-s6   t4 = s4+s5   t5 = s4-s5
// d1 = (t6-t5)*181/256   d2 = (t5+t6)*181/256
// w0 = t0+t7   w7 = t0-t7   w3 = t3+t4   w4 = t3-t4
// w2 = t2+d1   w5 = t2-d1   w1 = t1+d2   w6 = t1-d2
// outputs w_k >>> SH
//-----------------------------------------------------------------
module idct_1d #(
    parameter integer IW = 16,   // input width  (signed)
    parameter integer DW = 32,   // datapath width
    parameter integer SH = 11    // output shift
)(
    input  wire                clk,
    input  wire                start,
    input  wire signed [IW-1:0] x0, x1, x2, x3,
    input  wire signed [IW-1:0] x4, x5, x6, x7,
    output wire signed [DW-1:0] w0, w1, w2, w3,
    output wire signed [DW-1:0] w4, w5, w6, w7,
    output wire                ovalid
);

localparam integer AW = IW + 1;      // width of x0+x4
localparam integer CW = 13;          // coefficient width

localparam signed [CW-1:0] C1 = 4017;
localparam signed [CW-1:0] C2 = 3784;
localparam signed [CW-1:0] C3 = 3406;
localparam signed [CW-1:0] C4 = 2896;
localparam signed [CW-1:0] C5 = 2276;
localparam signed [CW-1:0] C6 = 1567;
localparam signed [CW-1:0] C7 =  799;

reg  signed [AW-1:0] r0, r1, r2, r3, r4, r5, r6, r7;
reg  signed [AW-1:0] i0, i1;
reg  [2:0]           ph;

reg  signed [CW-1:0] mb0, mb1, mb2, mb3;
reg  signed [AW-1:0] ma0, ma1, ma2, ma3;

reg  signed [DW-1:0] s0, s1, s2, s3, s4, s5, s6, s7;
reg  signed [DW-1:0] t0, t1, t2, t3, t4, t5, t6, t7;
reg  signed [DW-1:0] w0r, w3r, w4r, w7r, d1r, d2r;

// multiplier
wire signed [AW+CW-1:0] p0 = ma0 * mb0;
wire signed [AW+CW-1:0] p1 = ma1 * mb1;
wire signed [AW+CW-1:0] p2 = ma2 * mb2;
wire signed [AW+CW-1:0] p3 = ma3 * mb3;

wire signed [DW-1:0] p0e = {{(DW-AW-CW){p0[AW+CW-1]}}, p0};
wire signed [DW-1:0] p1e = {{(DW-AW-CW){p1[AW+CW-1]}}, p1};
wire signed [DW-1:0] p2e = {{(DW-AW-CW){p2[AW+CW-1]}}, p2};
wire signed [DW-1:0] p3e = {{(DW-AW-CW){p3[AW+CW-1]}}, p3};

// (t*181)/256  ==  (t>>>8) + (t>>>5) + (t>>>3) + (t>>>1)
function signed [DW-1:0] sc181;
    input signed [DW-1:0] v;
    begin
        sc181 = (v >>> 8) + (v >>> 5) + (v >>> 3) + (v >>> 1);
    end
endfunction

wire signed [DW-1:0] q1 = sc181(t6 - t5);
wire signed [DW-1:0] q2 = sc181(t5 + t6);

always @(posedge clk)
begin
    if (start) begin
        r0 <= {{(AW-IW){x0[IW-1]}}, x0};
        r1 <= {{(AW-IW){x1[IW-1]}}, x1};
        r2 <= {{(AW-IW){x2[IW-1]}}, x2};
        r3 <= {{(AW-IW){x3[IW-1]}}, x3};
        r4 <= {{(AW-IW){x4[IW-1]}}, x4};
        r5 <= {{(AW-IW){x5[IW-1]}}, x5};
        r6 <= {{(AW-IW){x6[IW-1]}}, x6};
        r7 <= {{(AW-IW){x7[IW-1]}}, x7};
        i0 <= {{(AW-IW){x0[IW-1]}}, x0} + {{(AW-IW){x4[IW-1]}}, x4};
        i1 <= {{(AW-IW){x0[IW-1]}}, x0} - {{(AW-IW){x4[IW-1]}}, x4};
        mb0 <= C4; mb1 <= C4; mb2 <= C2; mb3 <= C6;
        ma0 <= {{(AW-IW){x0[IW-1]}}, x0} + {{(AW-IW){x4[IW-1]}}, x4};
        ma1 <= {{(AW-IW){x0[IW-1]}}, x0} - {{(AW-IW){x4[IW-1]}}, x4};
        ma2 <= {{(AW-IW){x2[IW-1]}}, x2};
        ma3 <= {{(AW-IW){x6[IW-1]}}, x6};
        ph  <= 3'd0;
    end
    else begin
        ph <= ph + 3'd1;
        case (ph)
        3'd0: begin // capture round0 products
            s0 <= p0e; s1 <= p1e; s3 <= p2e + p3e;
            ma0 <= r2; ma1 <= r6; ma2 <= r4; ma3 <= r2;
            mb0 <= C2; mb1 <= C6; mb2 <= C5; mb3 <= C3;
        end
        3'd1: begin // capture round1 products
            s2 <= p0e - p1e; s6 <= p2e + p3e;
            t0 <= s0 + s3; t3 <= s0 - s3;
            ma0 <= r4; ma1 <= r2; ma2 <= r0; ma3 <= r6;
            mb0 <= C3; mb1 <= C5; mb2 <= C1; mb3 <= C7;
        end
        3'd2: begin // capture round2 products
            s5 <= p0e - p1e; s7 <= p2e + p3e;
            t1 <= s1 + s2; t2 <= s1 - s2;
            ma0 <= r0; ma1 <= r6; ma2 <= r0; ma3 <= r6;
            mb0 <= C7; mb1 <= C1; mb2 <= C7; mb3 <= C1;
        end
        3'd3: begin // capture round3 products
            s4 <= p0e - p1e;
            t7 <= s6 + s7; t6 <= s7 - s6;
        end
        3'd4: begin
            t5 <= s4 - s5; t4 <= s4 + s5;
        end
        3'd5: begin
            w0r <= (t0 + t7) >>> SH;
            w7r <= (t0 - t7) >>> SH;
            w3r <= (t3 + t4) >>> SH;
            w4r <= (t3 - t4) >>> SH;
            d1r <= q1;
            d2r <= q2;
        end
        default: ;
        endcase
    end
end

assign w0 = w0r;
assign w1 = (t1 + d2r) >>> SH;
assign w2 = (t2 + d1r) >>> SH;
assign w3 = w3r;
assign w4 = w4r;
assign w5 = (t2 - d1r) >>> SH;
assign w6 = (t1 - d2r) >>> SH;
assign w7 = w7r;
assign ovalid = (ph == 3'd6);

endmodule

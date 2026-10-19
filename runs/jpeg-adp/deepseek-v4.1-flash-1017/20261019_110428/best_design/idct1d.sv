//-----------------------------------------------------------------
// 1-D IDCT core (factorized JPEG / AAN-style, cos*4096) — bit exact
// re-implementation of ultraembedded/core_jpeg jpeg_idct_x / jpeg_idct_y
//
//   s0 = (x0+x4)*C4    s1 = (x0-x4)*C4
//   s2 = x2*C6 - x6*C2 s3 = x2*C2 + x6*C6
//   u  = x1*C7 - x7*C1 v  = x5*C3 - x3*C5
//   p  = x1*C1 + x7*C7 q  = x5*C5 + x3*C3
//   t0 = s0+s3  t3 = s0-s3  t1 = s1+s2  t2 = s1-s2
//   t7 = p+q    t6 = p-q    t4 = u+v    t5 = u-v
//   q5 = (t6-t5)*181/256   q6 = (t5+t6)*181/256     (32-bit wrapping)
//   w0 = t0+t7  w7 = t0-t7  w3 = t3+t4  w4 = t3-t4
//   w2 = t2+q5  w5 = t2-q5  w1 = t1+q6  w6 = t1-q6
//   out = w >>> SH
//
// x0..x7 must be stable while ph runs 1..7 (they are captured externally
// at phase 0).  Results are combinationally valid at ph == 7.
//-----------------------------------------------------------------
module idct1d #(
    parameter integer IW = 16,
    parameter integer SH = 11
)(
    input  wire               clk,
    input  wire [2:0]         ph,
    input  wire signed [IW-1:0] x0, x1, x2, x3,
    input  wire signed [IW-1:0] x4, x5, x6, x7,
    output wire signed [31:0] w0, w1, w2, w3,
    output wire signed [31:0] w4, w5, w6, w7
);

localparam signed [12:0] C1 = 13'sd4017;
localparam signed [12:0] C2 = 13'sd3784;
localparam signed [12:0] C3 = 13'sd3406;
localparam signed [12:0] C4 = 13'sd2896;
localparam signed [12:0] C5 = 13'sd2276;
localparam signed [12:0] C6 = 13'sd1567;
localparam signed [12:0] C7 = 13'sd 799;

wire signed [31:0] X0 = {{(32-IW){x0[IW-1]}}, x0};
wire signed [31:0] X1 = {{(32-IW){x1[IW-1]}}, x1};
wire signed [31:0] X2 = {{(32-IW){x2[IW-1]}}, x2};
wire signed [31:0] X3 = {{(32-IW){x3[IW-1]}}, x3};
wire signed [31:0] X4 = {{(32-IW){x4[IW-1]}}, x4};
wire signed [31:0] X5 = {{(32-IW){x5[IW-1]}}, x5};
wire signed [31:0] X6 = {{(32-IW){x6[IW-1]}}, x6};
wire signed [31:0] X7 = {{(32-IW){x7[IW-1]}}, x7};

// ---- multiplier operand selection (three rounds over phases 1,2,3) -----
wire rA = (ph == 3'd1);
wire rC = (ph == 3'd3);

wire signed [IW-1:0] ma0 = rA ? x0 : x1;
wire signed [IW-1:0] ma1 = rA ? x4 : x7;
wire signed [IW-1:0] ma2 = rA ? x2 : x5;
wire signed [IW-1:0] ma3 = rA ? x6 : x3;

wire signed [12:0] mb0 = rA ? C4 : rC ? C1 : C7;
wire signed [12:0] mb1 = rA ? C4 : rC ? C7 : C1;
wire signed [12:0] mb2 = rA ? C2 : rC ? C5 : C3;
wire signed [12:0] mb3 = rA ? C6 : rC ? C3 : C5;

wire signed [IW+12:0] p0 = ma0 * mb0;
wire signed [IW+12:0] p1 = ma1 * mb1;
wire signed [IW+12:0] p2 = ma2 * mb2;
wire signed [IW+12:0] p3 = ma3 * mb3;
wire signed [IW+12:0] p4 = x2 * C6;
wire signed [IW+12:0] p5 = x6 * C2;

reg signed [IW+12:0] a0, a1, a2, a3, a4, a5;   // round A products
reg signed [IW+12:0] b0, b1, b2, b3;           // round B products
reg signed [IW+12:0] c0, c1, c2, c3;           // round C products

reg signed [31:0] t0, t1, t2, t3, t4, t5, t6, t7;
reg signed [31:0] d5, d6;
reg signed [31:0] q5, q6;

// t values, combinationally derived from the product registers
wire signed [31:0] s0 = {{(32-IW-13){a0[IW+12]}}, a0} + {{(32-IW-13){a1[IW+12]}}, a1};
wire signed [31:0] s1 = {{(32-IW-13){a0[IW+12]}}, a0} - {{(32-IW-13){a1[IW+12]}}, a1};
wire signed [31:0] s3 = {{(32-IW-13){a2[IW+12]}}, a2} + {{(32-IW-13){a3[IW+12]}}, a3};
wire signed [31:0] s2 = {{(32-IW-13){a4[IW+12]}}, a4} - {{(32-IW-13){a5[IW+12]}}, a5};
wire signed [31:0] uv = {{(32-IW-13){b0[IW+12]}}, b0} - {{(32-IW-13){b1[IW+12]}}, b1};
wire signed [31:0] vv = {{(32-IW-13){b2[IW+12]}}, b2} - {{(32-IW-13){b3[IW+12]}}, b3};
wire signed [31:0] pp = {{(32-IW-13){c0[IW+12]}}, c0} + {{(32-IW-13){c1[IW+12]}}, c1};
wire signed [31:0] qq = {{(32-IW-13){c2[IW+12]}}, c2} + {{(32-IW-13){c3[IW+12]}}, c3};

// (v*181) >>> 8 with 32-bit wrapping multiply and truncation toward zero
function signed [31:0] scale181;
    input signed [31:0] v;
    reg   signed [31:0] p;
    begin
        p = (v <<< 7) + (v <<< 5) + (v <<< 4) + (v <<< 2) + v;
        scale181 = (p + (p[31] ? 32'sd255 : 32'sd0)) >>> 8;
    end
endfunction

wire signed [31:0] nq5 = scale181(d6 - d5);
wire signed [31:0] nq6 = scale181(d5 + d6);

always @(posedge clk)
begin
    case (ph)
    3'd2: begin
        a0 <= p0; a1 <= p1; a2 <= p2; a3 <= p3; a4 <= p4; a5 <= p5;
    end
    3'd3: begin
        b0 <= p0; b1 <= p1; b2 <= p2; b3 <= p3;
    end
    3'd4: begin
        c0 <= p0; c1 <= p1; c2 <= p2; c3 <= p3;
        t0 <= s0 + s3;  t3 <= s0 - s3;
        t1 <= s1 + s2;  t2 <= s1 - s2;
        t7 <= pp + qq;  t6 <= pp - qq;
        t4 <= uv + vv;  t5 <= uv - vv;
    end
    3'd5: begin
        d5 <= t6 - t5;
        d6 <= t5 + t6;
    end
    3'd6: begin
        q5 <= nq5;
        q6 <= nq6;
    end
    default: ;
    endcase
end

assign w0 = (t0 + t7) >>> SH;
assign w7 = (t0 - t7) >>> SH;
assign w3 = (t3 + t4) >>> SH;
assign w4 = (t3 - t4) >>> SH;
assign w2 = (t2 + q5) >>> SH;
assign w5 = (t2 - q5) >>> SH;
assign w1 = (t1 + q6) >>> SH;
assign w6 = (t1 - q6) >>> SH;

endmodule

//-----------------------------------------------------------------
// 1-D 8-point IDCT (factorized, cos*4096) — bit exact re-implementation
// of ultraembedded/core_jpeg jpeg_idct_x / jpeg_idct_y arithmetic.
//
//   s0 = (x0+x4)*C4     s1 = (x0-x4)*C4
//   s2 = x2*C6 - x6*C2  s3 = x2*C2 + x6*C6
//   s4 = x1*C7 - x7*C1  s5 = x5*C3 - x3*C5
//   s6 = x5*C5 + x3*C3  s7 = x1*C1 + x7*C7
//   t0 = s0+s3  t3 = s0-s3  t1 = s1+s2  t2 = s1-s2
//   t7 = s6+s7  t6 = s7-s6  t4 = s4+s5  t5 = s4-s5
//   q5 = (t6-t5)*181/256    q6 = (t5+t6)*181/256
//   w0 = t0+t7  w7 = t0-t7  w3 = t3+t4  w4 = t3-t4
//   w2 = t2+q5  w5 = t2-q5  w1 = t1+q6  w6 = t1-q6
//   out = w >>> SH
//
// Seven *constant* multipliers, each used twice (two operand rounds), so
// 14 products are produced in 2 cycles.  go is asserted one cycle before
// the x inputs are required; results are registered 7 cycles later.
//-----------------------------------------------------------------
module idct1d #(
    parameter integer IW = 16,   // input width  (signed)
    parameter integer PW = 28,   // product width (signed, exact)
    parameter integer DW = 32,   // datapath width (signed, wrapping)
    parameter integer SH = 11    // output arithmetic shift
)(
    input  wire                 clk,
    input  wire                 go,
    input  wire signed [IW-1:0] x0, x1, x2, x3,
    input  wire signed [IW-1:0] x4, x5, x6, x7,
    output wire signed [DW-1:0] w0, w1, w2, w3,
    output wire signed [DW-1:0] w4, w5, w6, w7,
    output wire                 ovalid
);

localparam signed [12:0] K1 = 13'sd4017;
localparam signed [12:0] K2 = 13'sd3784;
localparam signed [12:0] K3 = 13'sd3406;
localparam signed [12:0] K4 = 13'sd2896;
localparam signed [12:0] K5 = 13'sd2276;
localparam signed [12:0] K6 = 13'sd1567;
localparam signed [12:0] K7 = 13'sd 799;

//-----------------------------------------------------------------
// Phase counter: 0 = idle, 1..7 = pipeline stages
//-----------------------------------------------------------------
reg [2:0] ph;

always @(posedge clk)
begin
    if (go)          ph <= 3'd1;
    else if (ph != 3'd0) ph <= ph + 3'd1;
end

//-----------------------------------------------------------------
// Multiplier operand selects (round A at ph==1, round B at ph==2)
//-----------------------------------------------------------------
wire       ra = (ph == 3'd1);

wire signed [IW-1:0] o0 = ra ? x0 : x4;   // * K4
wire signed [IW-1:0] o1 = ra ? x2 : x6;   // * K6
wire signed [IW-1:0] o2 = ra ? x2 : x6;   // * K2
wire signed [IW-1:0] o3 = ra ? x1 : x7;   // * K7
wire signed [IW-1:0] o4 = ra ? x5 : x3;   // * K3
wire signed [IW-1:0] o5 = ra ? x5 : x3;   // * K5
wire signed [IW-1:0] o6 = ra ? x1 : x7;   // * K1

wire signed [PW-1:0] p0 = o0 * K4;
wire signed [PW-1:0] p1 = o1 * K6;
wire signed [PW-1:0] p2 = o2 * K2;
wire signed [PW-1:0] p3 = o3 * K7;
wire signed [PW-1:0] p4 = o4 * K3;
wire signed [PW-1:0] p5 = o5 * K5;
wire signed [PW-1:0] p6 = o6 * K1;

reg signed [PW-1:0] A0,A1,A2,A3,A4,A5,A6;
reg signed [PW-1:0] B0,B1,B2,B3,B4,B5,B6;

reg signed [DW-1:0] s0,s1,s2,s3,s4,s5,s6,s7;
reg signed [DW-1:0] t0,t1,t2,t3,t4,t5,t6,t7;
reg signed [DW-1:0] q5,q6;
reg signed [DW-1:0] y0,y1,y2,y3,y4,y5,y6,y7;

function signed [DW-1:0] ext;
    input signed [PW-1:0] v;
    begin ext = {{(DW-PW){v[PW-1]}}, v}; end
endfunction

// (v*181)/256 with 32-bit wrapping multiply, truncation toward zero
wire signed [DW-1:0] dv5 = t6 - t5;
wire signed [DW-1:0] dv6 = t5 + t6;
wire signed [DW-1:0] m5  = (dv5 <<< 7) + (dv5 <<< 5) + (dv5 <<< 4) + (dv5 <<< 2) + dv5;
wire signed [DW-1:0] m6  = (dv6 <<< 7) + (dv6 <<< 5) + (dv6 <<< 4) + (dv6 <<< 2) + dv6;
wire signed [DW-1:0] n5  = (m5 + (m5[DW-1] ? {{(DW-8){1'b0}}, 8'd255} : {DW{1'b0}})) >>> 8;
wire signed [DW-1:0] n6  = (m6 + (m6[DW-1] ? {{(DW-8){1'b0}}, 8'd255} : {DW{1'b0}})) >>> 8;

always @(posedge clk)
begin
    case (ph)
    3'd1: begin
        A0 <= p0; A1 <= p1; A2 <= p2; A3 <= p3;
        A4 <= p4; A5 <= p5; A6 <= p6;
    end
    3'd2: begin
        B0 <= p0; B1 <= p1; B2 <= p2; B3 <= p3;
        B4 <= p4; B5 <= p5; B6 <= p6;
    end
    3'd3: begin
        s0 <= ext(A0) + ext(B0);
        s1 <= ext(A0) - ext(B0);
        s2 <= ext(A1) - ext(B2);
        s3 <= ext(A2) + ext(B1);
        s4 <= ext(A3) - ext(B6);
        s5 <= ext(A4) - ext(B5);
        s6 <= ext(A5) + ext(B4);
        s7 <= ext(A6) + ext(B3);
    end
    3'd4: begin
        t0 <= s0 + s3; t3 <= s0 - s3;
        t1 <= s1 + s2; t2 <= s1 - s2;
        t7 <= s6 + s7; t6 <= s7 - s6;
        t4 <= s4 + s5; t5 <= s4 - s5;
    end
    3'd5: begin
        q5 <= n5;
        q6 <= n6;
        y0 <= (t0 + t7) >>> SH;
        y7 <= (t0 - t7) >>> SH;
        y3 <= (t3 + t4) >>> SH;
        y4 <= (t3 - t4) >>> SH;
    end
    3'd6: begin
        y2 <= (t2 + q5) >>> SH;
        y5 <= (t2 - q5) >>> SH;
        y1 <= (t1 + q6) >>> SH;
        y6 <= (t1 - q6) >>> SH;
    end
    default: ;
    endcase
end

assign w0 = y0;
assign w1 = y1;
assign w2 = y2;
assign w3 = y3;
assign w4 = y4;
assign w5 = y5;
assign w6 = y6;
assign w7 = y7;
assign ovalid = (ph == 3'd7);

endmodule

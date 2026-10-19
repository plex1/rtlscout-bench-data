//-----------------------------------------------------------------
// jpeg_idct_2d  —  8x8 2-D IDCT (baseline JPEG)
//
// Bit-exact re-implementation of the ultraembedded/core_jpeg
// jpeg_idct_x -> jpeg_idct_transpose -> jpeg_idct_y composition.
//
//   * all multiplies are by *constants* (cos*4096), 4 shared
//     constant multipliers per 1-D engine
//   * the double-buffered 4x32x32 corner-turn RAM is replaced by
//     two 64x20 transpose banks written a whole row at a time by the
//     row pass and read as a whole column by the column pass
//   * 32-bit signed wrap-around semantics preserved exactly
//   * (x*181)/256 done as a shift/add tree
//-----------------------------------------------------------------
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

localparam signed [12:0] C1 = 13'sd4017;
localparam signed [12:0] C2 = 13'sd3784;
localparam signed [12:0] C3 = 13'sd3406;
localparam signed [12:0] C4 = 13'sd2896;
localparam signed [12:0] C5 = 13'sd2276;
localparam signed [12:0] C6 = 13'sd1567;
localparam signed [12:0] C7 = 13'sd 799;

//=================================================================
// Input capture : two coefficient register sets (ping-pong)
//=================================================================
reg  signed [15:0] xr [0:15];      // set0 = 0..7 , set1 = 8..15
reg  [2:0]         xrow_q [0:1];
reg                xrdy_q [0:1];
reg                wrset_q;

wire [5:0] beat64 = inport_valid_i ? 6'd1 : 6'd0;

always @(posedge clk)
begin
    if (rst) begin
        wrset_q   <= 1'b0;
        xrdy_q[0] <= 1'b0;
        xrdy_q[1] <= 1'b0;
        xrow_q[0] <= 3'd0;
        xrow_q[1] <= 3'd0;
    end
    else begin
        if (img_start_i) begin
            wrset_q   <= 1'b0;
            xrdy_q[0] <= 1'b0;
            xrdy_q[1] <= 1'b0;
        end
        else if (inport_valid_i) begin
            if (inport_idx_i == 3'd0) begin
                xr[{wrset_q,3'b000}] <= inport_data0_i;
                xr[{wrset_q,3'b000}+3'd1] <= inport_data1_i;
                xr[{wrset_q,3'b000}+3'd2] <= inport_data2_i;
                xr[{wrset_q,3'b000}+3'd3] <= inport_data3_i;
            end
            else if (inport_idx_i == 3'd1) begin
                xr[{wrset_q,3'b000}+3'd4] <= inport_data0_i;
                xr[{wrset_q,3'b000}+3'd5] <= inport_data1_i;
                xr[{wrset_q,3'b000}+3'd6] <= inport_data2_i;
                xr[{wrset_q,3'b000}+3'd7] <= inport_data3_i;
                xrdy_q[wrset_q] <= 1'b1;
                wrset_q         <= ~wrset_q;
            end
        end
    end
end

//=================================================================
// 1-D core (shared source for both passes)
//=================================================================
module idct1d_core #(parameter integer IW = 16) (
    input  wire                 clk,
    input  wire                 start,
    input  wire [2:0]           ph,
    input  wire signed [IW-1:0] x0, x1, x2, x3,
    input  wire signed [IW-1:0] x4, x5, x6, x7,
    output wire signed [31:0]   w0, w1, w2, w3,
    output wire signed [31:0]   w4, w5, w6, w7
);
localparam signed [12:0] C1 = 13'sd4017;
localparam signed [12:0] C2 = 13'sd3784;
localparam signed [12:0] C3 = 13'sd3406;
localparam signed [12:0] C4 = 13'sd2896;
localparam signed [12:0] C5 = 13'sd2276;
localparam signed [12:0] C6 = 13'sd1567;
localparam signed [12:0] C7 = 13'sd 799;

// --- multiplier operand / constant select (4 rounds over ph 0..3) -----
wire signed [IW-1:0] ma0 = (ph == 3'd0) ? x0 : (ph == 3'd1) ? x2 : (ph == 3'd2) ? x5 : x5;
wire signed [IW-1:0] ma1 = (ph == 3'd0) ? x4 : (ph == 3'd1) ? x6 : (ph == 3'd2) ? x3 : x3;
wire signed [IW-1:0] ma2 = (ph == 3'd0) ? x2 : (ph == 3'd1) ? x1 : (ph == 3'd2) ? x1 : 0;
wire signed [IW-1:0] ma3 = (ph == 3'd0) ? x6 : (ph == 3'd1) ? x7 : (ph == 3'd2) ? x7 : 0;

wire signed [12:0] mb0 = (ph == 3'd0) ? C4 : (ph == 3'd1) ? C6 : (ph == 3'd2) ? C5 : C3;
wire signed [12:0] mb1 = (ph == 3'd0) ? C4 : (ph == 3'd1) ? C2 : (ph == 3'd2) ? C3 : C5;
wire signed [12:0] mb2 = (ph == 3'd0) ? C2 : (ph == 3'd1) ? C1 : (ph == 3'd2) ? C7 : 13'sd0;
wire signed [12:0] mb3 = (ph == 3'd0) ? C6 : (ph == 3'd1) ? C7 : (ph == 3'd2) ? C1 : 13'sd0;

wire signed [IW+12:0] pr0 = ma0 * mb0;
wire signed [IW+12:0] pr1 = ma1 * mb1;
wire signed [IW+12:0] pr2 = ma2 * mb2;
wire signed [IW+12:0] pr3 = ma3 * mb3;

localparam integer PW = IW + 13;

reg signed [PW-1:0] q1_0, q1_1, q1_2, q1_3;
reg signed [PW-1:0] q2_0, q2_1, q2_2, q2_3;
reg signed [PW-1:0] q3_0, q3_1, q3_2, q3_3;
reg signed [PW-1:0] q4_0, q4_1;

reg signed [31:0] t0,t1,t2,t3,t4,t5,t6,t7;
reg signed [31:0] qq5, qq6;
reg signed [31:0] r0,r1,r2,r3,r4,r5,r6,r7;

function signed [31:0] se;
    input signed [PW-1:0] v;
    begin se = {{(32-PW){v[PW-1]}}, v}; end
endfunction

// (v*181)>>>8 with 32-bit wrapping multiply, truncation toward zero
function signed [31:0] sc181;
    input signed [31:0] v;
    reg   signed [31:0] p;
    begin
        p = (v <<< 7) + (v <<< 5) + (v <<< 4) + (v <<< 2) + v;
        sc181 = (p + (p[31] ? 32'sd255 : 32'sd0)) >>> 8;
    end
endfunction

wire signed [31:0] s0 = se(q1_0) + se(q1_1);
wire signed [31:0] s1 = se(q1_0) - se(q1_1);
wire signed [31:0] s3 = se(q1_2) + se(q1_3);
wire signed [31:0] s2 = se(q2_0) - se(q2_1);
wire signed [31:0] s7 = se(q2_2) + se(q2_3);
wire signed [31:0] s6 = se(q3_0) + se(q3_1);
wire signed [31:0] s4 = se(q3_2) - se(q3_3);
wire signed [31:0] s5 = se(q4_0) - se(q4_1);

wire signed [31:0] nt0 = s0 + s3;
wire signed [31:0] nt3 = s0 - s3;
wire signed [31:0] nt1 = s1 + s2;
wire signed [31:0] nt2 = s1 - s2;
wire signed [31:0] nt7 = s6 + s7;
wire signed [31:0] nt6 = s7 - s6;
wire signed [31:0] nt4 = s4 + s5;
wire signed [31:0] nt5 = s4 - s5;

always @(posedge clk)
begin
    case (ph)
    3'd0: begin
        q1_0 <= pr0; q1_1 <= pr1; q1_2 <= pr2; q1_3 <= pr3;
    end
    3'd1: begin
        q2_0 <= pr0; q2_1 <= pr1; q2_2 <= pr2; q2_3 <= pr3;
    end
    3'd2: begin
        q3_0 <= pr0; q3_1 <= pr1; q3_2 <= pr2; q3_3 <= pr3;
    end
    3'd3: begin
        q4_0 <= pr0; q4_1 <= pr1;
        t0 <= nt0; t3 <= nt3;
    end
    3'd4: begin
        t1 <= nt1; t2 <= nt2;
        t7 <= nt7; t6 <= nt6;
    end
    3'd5: begin
        t4 <= nt4; t5 <= nt5;
    end
    3'd6: begin
        qq5 <= sc181(t6 - t5);
        qq6 <= sc181(t5 + t6);
    end
    3'd7: begin
        r0 <= (t0 + t7) >>> 11;
        r7 <= (t0 - t7) >>> 11;
        r3 <= (t3 + t4) >>> 11;
        r4 <= (t3 - t4) >>> 11;
        r2 <= (t2 + qq5) >>> 11;
        r5 <= (t2 - qq5) >>> 11;
        r1 <= (t1 + qq6) >>> 11;
        r6 <= (t1 - qq6) >>> 11;
    end
    default: ;
    endcase
end

assign w0 = r0; assign w1 = r1; assign w2 = r2; assign w3 = r3;
assign w4 = r4; assign w5 = r5; assign w6 = r6; assign w7 = r7;

endmodule

endmodule

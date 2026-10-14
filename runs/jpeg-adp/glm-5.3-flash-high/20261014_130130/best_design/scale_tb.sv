// isolated testbench for idct_scale181 (from design7.v)
module idct_scale181 (
    input  wire signed [31:0] x,
    output wire        [31:0] m,
    output wire               incr
);
    wire [31:0] r0 = x;
    wire [31:0] r1 = x << 2;
    wire [31:0] r2 = x << 4;
    wire [31:0] r3 = x << 5;
    wire [31:0] r4 = x << 7;

    wire [31:0] s1 = r0 ^ r1 ^ r2;
    wire [31:0] c1 = ((r0 & r1) | (r0 & r2) | (r1 & r2)) << 1;
    wire [31:0] s2 = r3 ^ r4;
    wire [31:0] c2 = (r3 & r4) << 1;

    wire [31:0] s3 = s1 ^ c1 ^ s2;
    wire [31:0] c3 = ((s1 & c1) | (s1 & s2) | (c1 & s2)) << 1;

    wire [31:0] s4 = s3 ^ c3 ^ c2;
    wire [31:0] c4 = ((s3 & c3) | (s3 & c2) | (c3 & c2)) << 1;

    assign m = s4 + c4;
    assign incr = m[31] & (m[7:0] != 8'b0);
endmodule

module jpeg_idct_2d;
  reg signed [31:0] x;
  wire [31:0] m;
  wire incr;
  idct_scale181 dut (.x(x), .m(m), .incr(incr));

  function [31:0] refmul;
    input signed [31:0] xx;
    begin
      refmul = (xx * 181);
    end
  endfunction

  function signed [31:0] refdiv;
    input signed [31:0] xx;
    begin
      refdiv = (xx * 181) / 256;
    end
  endfunction

  integer i, errs;
  reg signed [31:0] v;
  reg signed [31:0] q;
  reg signed [31:0] qref;
  initial begin
    errs = 0;
    for (i = 0; i < 200000; i = i + 1) begin
      case (i % 6)
        0: v = 0;
        1: v = 32'sd2147483647;
        2: v = -32'sd2147483648;
        3: v = i * 37;
        4: v = -(i * 91);
        5: v = $random;
      endcase
      x = v;
      #1;
      if (m !== refmul(v)) begin
        if (errs < 5) $display("MULBAD x=%0d m=%h ref=%h", v, m, refmul(v));
        errs = errs + 1;
      end
      q = $signed(m) >>> 8;
      if (incr) q = q + 1;
      qref = refdiv(v);
      if (q !== qref) begin
        if (errs < 10) $display("DIVBAD x=%0d q=%0d ref=%0d m=%h", v, q, qref, m);
        errs = errs + 1;
      end
    end
    $display("DONE errs=%0d", errs);
    $finish;
  end
endmodule

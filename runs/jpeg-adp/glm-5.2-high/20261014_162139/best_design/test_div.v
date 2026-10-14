module test_div(input [31:0] a, output [31:0] q1, output [31:0] q2);
  wire signed [31:0] prod = a * 181;
  assign q1 = prod / 256;
  assign q2 = (prod >>> 8) + (prod[31] & |prod[7:0]);
endmodule

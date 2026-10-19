// probe: significand multiply only, to find the delay floor of the 11x11 product
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a, b;
  output [15:0] y;
  wire [10:0] ma = {|a[14:10], a[9:0]};
  wire [10:0] mb = {|b[14:10], b[9:0]};
  wire [21:0] P  = ma * mb;
  assign y = P[21:6];
endmodule

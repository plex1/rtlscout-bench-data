module tpu
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [63:0]	w_data,
    input		a_valid,
    input [63:0]	a_data,
    output		c_valid,
    output [255:0]	c_data
    );
   wire [15:0] prod_b;
   wire [15:0] prod_m;
   assign prod_b = $signed(a_data[7:0]) * $signed(w_data[7:0]);
   mul8s u_mul (.x(a_data[7:0]), .y(w_data[7:0]), .p(prod_m));
   assign c_valid = 1'b1;
   assign c_data = {255'b0, (prod_m != prod_b)};
endmodule

// Radix-4 Booth signed 8x8 -> 16 multiplier, partial products sign-extended
// to 16 bits and summed. Digits computed from multiplicand x, applied to y.
module mul8s
   (
    input  [7:0]  x,
    input  [7:0]  y,
    output [15:0] p
    );

   // For a signed multiplicand a (x), append implied 0 at bit -1 and sign
   // bit: a8 = {a[7], a, 1'b0}. Digit k = a8[2k+2 : 2k] for k=0..3.
   wire [9:0] a8 = {x[7], x, 1'b0};
   wire [2:0] d0 = a8[2:0];
   wire [2:0] d1 = a8[4:2];
   wire [2:0] d2 = a8[6:4];
   wire [2:0] d3 = a8[8:6];

   // Booth value in {-2,-1,0,1,2} times y, sign-extended to 16 bits.
   function [15:0] ppext;
      input [2:0] dig;
      input [7:0] m;
      reg [15:0] se;	// sign-extended m
      reg [15:0] t2;	// 2*m
      begin
	 se = {{8{m[7]}}, m};
	 t2 = { {7{m[7]}}, m, 1'b0 };
	 case (dig)
	   3'b000, 3'b111: ppext = 16'b0;
	   3'b001, 3'b010: ppext = se;
	   3'b011:	    ppext = t2;
	   3'b100:	    ppext = ~t2 + 16'd1;
	   default:	    ppext = ~se + 16'd1;	// 101,110 -> -1
	 endcase
      end
   endfunction

   wire [15:0] pp0 = ppext(d0, y);
   wire [15:0] pp1 = ppext(d1, y) << 2;
   wire [15:0] pp2 = ppext(d2, y) << 4;
   wire [15:0] pp3 = ppext(d3, y) << 6;

   assign p = pp0 + pp1 + pp2 + pp3;

endmodule

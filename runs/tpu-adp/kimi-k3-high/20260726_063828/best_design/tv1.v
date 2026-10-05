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
   wire [15:0] prod;
   mul8s u_mul (.x(a_data[7:0]), .y(w_data[7:0]), .p(prod));
   assign c_valid = 1'b0;
   assign c_data = {240'b0, prod};
endmodule

module mul8s
   (
    input  [7:0]  x,
    input  [7:0]  y,
    output [15:0] p
    );

   function [8:0] booth_pp;
      input [8:0] wsel;
      input [7:0] mult;
      reg   [2:0] t;
      reg   [8:0] u;
      begin
	 t = {wsel[4], wsel[2], wsel[0]};
	 case (t)
	   3'b000, 3'b111: u = 9'b0;
	   3'b001, 3'b010: u = {mult[7], mult};
	   3'b011:	    u = {mult, 1'b0};
	   3'b100:	    u = ~{mult, 1'b0};
	   default:	    u = ~{mult[7], mult};
	 endcase
	 booth_pp = {u, wsel[6] ^ wsel[4]};
      end
   endfunction

   wire y7 = y[7];
   wire [8:0] w1 = {x[0], 1'b0}  ^ {y7, y7};
   wire [8:0] w2 = {x[2], x[1]}  ^ {y7, y7};
   wire [8:0] w3 = {x[4], x[3]}  ^ {y7, y7};
   wire [8:0] w4 = {x[6], x[5]}  ^ {y7, y7};

   wire [8:0] pp0 = booth_pp(w1, y);
   wire [8:0] pp1 = booth_pp(w2, y);
   wire [8:0] pp2 = booth_pp(w3, y);
   wire [8:0] pp3 = booth_pp(w4, y);

   wire [15:0] c0 = {pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[8], pp0[7:0]};
   wire [15:0] c1 = {pp1[8], pp1[8], pp1[8], pp1[8], pp1[8], pp1[8], pp1[7:0], 2'b00};
   wire [15:0] c2 = {pp2[8], pp2[8], pp2[8], pp2[8], pp2[7:0], 4'b0000};
   wire [15:0] c3 = {pp3[8], pp3[8], pp3[7:0], 6'b000000};
   wire [15:0] bw = 16'hA002;

   assign p = c0 + c1 + c2 + c3 + bw;

endmodule

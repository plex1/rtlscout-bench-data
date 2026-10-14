module mul13_tb;
  reg clk = 0;
  reg rst = 1;
  reg signed [20:0] a;
  reg [12:0] b;
  wire [31:0] p;

  idct_mul13 #(.AW(21)) dut (.clk(clk), .rst(rst), .a_i(a), .b_i(b), .p_o(p));

  always #5 clk = ~clk;

  // reference: 32-bit wrap multiply (a signed, b positive)
  function [31:0] refmul;
    input signed [20:0] aa;
    input [12:0] bb;
    reg [63:0] prod;
    begin
      prod = $signed(aa) * $signed({1'b0, bb});
      refmul = prod[31:0];
    end
  endfunction

  reg [31:0] err;
  integer i, j;
  reg signed [20:0] av;
  reg [12:0] bv;

  initial begin
    err = 0;
    a = 0; b = 0;
    repeat (4) @(posedge clk);
    rst = 0;
    // sweep constants and a bunch of operands
    for (i = 0; i < 4096; i = i + 1) begin
      bv = i;
      for (j = 0; j < 64; j = j + 1) begin
        case (j % 8)
          0: av = 21'sd0;
          1: av = 21'sd1048575;
          2: av = -21'sd1048576;
          3: av = 21'sd1;
          4: av = -21'sd1;
          5: av = $random;
          6: av = -($random % 1048576);
          7: av = ($random % 1048576);
        endcase
        @(negedge clk);
        a = av; b = bv;
        @(negedge clk);
        if (p !== refmul(av, bv)) begin
          if (err < 10)
            $display("MISMATCH a=%0d b=%0d p=%h ref=%h", av, bv, p, refmul(av, bv));
          err = err + 1;
        end
      end
    end
    $display("DONE errors=%0d", err);
    $finish;
  end
endmodule

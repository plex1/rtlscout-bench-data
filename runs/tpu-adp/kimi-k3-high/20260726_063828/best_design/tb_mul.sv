module tb;
  logic clk, rst, w_valid, a_valid, c_valid;
  logic [63:0] w_data, a_data;
  logic [255:0] c_data;

  tpu dut (.clk(clk), .rst(rst), .w_valid(w_valid), .w_data(w_data),
           .a_valid(a_valid), .a_data(a_data), .c_valid(c_valid), .c_data(c_data));

  initial clk = 0;
  always #5 clk = ~clk;

  integer aa, bb, errors;
  integer got, exp;

  initial begin
    errors = 0;
    rst = 1;
    repeat (2) @(negedge clk);
    rst = 0;
    for (aa = 0; aa < 256; aa = aa + 1) begin
      for (bb = 0; bb < 256; bb = bb + 1) begin
        a_data = {56'b0, aa[7:0]};
        w_data = {56'b0, bb[7:0]};
        #1;
        got = $signed({1'b0, c_data[15:0]});
        exp = $signed(aa[7:0]) * $signed(bb[7:0]);
        if (got !== exp) begin
          errors = errors + 1;
          if (errors < 20)
            $display("MUL_ERR a=%0d b=%0d got=%0d exp=%0d", $signed(aa[7:0]), $signed(bb[7:0]), got, exp);
        end
      end
    end
    $display("TB_SUMMARY mul_errors=%0d", errors);
    $finish;
  end
endmodule

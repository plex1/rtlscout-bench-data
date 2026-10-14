module t;
  localparam int M = 7;
  logic [M-1:0] z;
  logic [M-1:0] o0, o1, o2, o3;
  logic [M-1:0] n0, n1, n2, n3;
  logic [M-1:0] pv;
  assign pv = {z[M-2:0], 1'b0};
  // construction A: fold right, then XOR with shifted original
  assign o0 = pv;
  assign o1 = pv | {1'b0, o0[M-2:0]};
  assign o2 = o1 | {{2{1'b0}}, o1[M-3:0]};
  assign o3 = o2 | {{4{1'b0}}, o2[M-5:0]};
  // construction B: propagate-minus then XOR with original
  assign n0 = z;
  assign n1 = z | {1'b0, n0[M-2:0]};
  assign n2 = n1 | {{2{1'b0}}, n1[M-3:0]};
  assign n3 = n2 | {{4{1'b0}}, n2[M-5:0]};
  initial begin
    for (int i = 0; i < 128; i++) begin
      z = i;
      #1;
      if (o3 != (-z)) $display("A %0d %04b %04b", i, o3, -z);
      if (n3 != (-z)) $display("B %0d %04b %04b", i, n3, -z);
    end
    $finish;
  end
endmodule

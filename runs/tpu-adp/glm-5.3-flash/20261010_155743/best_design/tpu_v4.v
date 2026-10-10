//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU-style MXU), 8x8.
// Cycle-accurate re-implementation of the logikbench `tpu` benchmark.
//
// Schedule is IDENTICAL to the original PE mesh: product (i,j) of activation
// row m is formed at compute index m+i+j, sampling the stationary weight as
// of that same compute index (row skew of i stages, east chain of j stages,
// weights shifted south while w_valid).  Output de-skew aligns the columns.
//
// ARITHMETIC (the optimization): each accumulator stage keeps its partial
// sum in CARRY-SAVE form (S, C) with S + C == true value (mod 2^19).  The
// 8x8 signed product is injected as raw partial-product rows (Baugh-Wooley
// style sign corrections), so the stage is a pure carry-save Wallace
// reduction -- no carry-propagate adder anywhere in the 8-stage chain.  The
// final S+C resolution (one 19-bit CPA) happens once, at the output.
//
//   stage: (S',C') = CSreduce( S, C, pp_rows(A,B) )
//   pp rows: (A & {8{b_k}}) << k  for k=0..7          (A_u*B_u)
//            -256*b7*A_u  as ~Y+1, Y=(A&{8{b7}})<<8
//            -256*a7*B_u  as ~Y+1, Y=(B&{8{a7}})<<8
//            +2^16*a7*b7  (single bit), +2 constant (the two +1s)
//   => S'+C' == S + C + A*B (mod 2^19)  exactly, signed-safe.
//#############################################################################

module csa #(parameter W = 19) (
    input  wire [W-1:0] x, y, z,
    output wire [W-1:0] s, c
);
  wire [W-1:0] maj = (x & y) | (x & z) | (y & z);
  assign s = x ^ y ^ z;
  assign c = {maj[W-2:0], 1'b0};
endmodule

// 14-operand Wallace reduction: 14 -> 10 -> 7 -> 5 -> 4 -> 3 -> 2
module csa_tree14 #(parameter W = 19) (
    input  wire [14*W-1:0] din,
    output wire [W-1:0]    s, c
);
  wire [W-1:0] l0 [0:13];
  wire [W-1:0] l1 [0:9];
  wire [W-1:0] l2 [0:6];
  wire [W-1:0] l3 [0:4];
  wire [W-1:0] l4 [0:3];
  wire [W-1:0] l5 [0:2];
  wire [W-1:0] l6 [0:1];
  genvar g, k;
  generate
    for (g = 0; g < 14; g = g + 1) begin : g_in
      assign l0[g] = din[g*W +: W];
    end
    // 14 -> 10 : four CSAs (12 operands) + two pass-through
    for (k = 0; k < 4; k = k + 1) begin : g_a
      csa #(.W(W)) u (.x(l0[3*k]), .y(l0[3*k+1]), .z(l0[3*k+2]),
                      .s(l1[2*k]), .c(l1[2*k+1]));
    end
    for (k = 0; k < 2; k = k + 1) begin : g_ap
      assign l1[8+k] = l0[12+k];
    end
    // 10 -> 7 : three CSAs + one pass-through
    for (k = 0; k < 3; k = k + 1) begin : g_b
      csa #(.W(W)) u (.x(l1[3*k]), .y(l1[3*k+1]), .z(l1[3*k+2]),
                      .s(l2[2*k]), .c(l2[2*k+1]));
    end
    assign l2[6] = l1[9];
    // 7 -> 5 : two CSAs + one pass-through
    for (k = 0; k < 2; k = k + 1) begin : g_c
      csa #(.W(W)) u (.x(l2[3*k]), .y(l2[3*k+1]), .z(l2[3*k+2]),
                      .s(l3[2*k]), .c(l3[2*k+1]));
    end
    assign l3[4] = l2[6];
    // 5 -> 4 : one CSA + two pass-through
    for (k = 0; k < 1; k = k + 1) begin : g_d
      csa #(.W(W)) u (.x(l3[0]), .y(l3[1]), .z(l3[2]), .s(l4[0]), .c(l4[1]));
    end
    assign l4[2] = l3[3];
    assign l4[3] = l3[4];
    // 4 -> 3 : one CSA + one pass-through
    for (k = 0; k < 1; k = k + 1) begin : g_e
      csa #(.W(W)) u (.x(l4[0]), .y(l4[1]), .z(l4[2]), .s(l5[0]), .c(l5[1]));
    end
    assign l5[2] = l4[3];
    // 3 -> 2
    for (k = 0; k < 1; k = k + 1) begin : g_f
      csa #(.W(W)) u (.x(l5[0]), .y(l5[1]), .z(l5[2]), .s(l6[0]), .c(l6[1]));
    end
    assign s = l6[0];
    assign c = l6[1];
  endgenerate
endmodule

// One accumulate stage: (S',C') = CSreduce(S, C, partial products of A*B).
// All operands are W-bit vectors; the identity is exact mod 2^W.
module tpu_mac #(parameter W = 19) (
    input  wire          clk, rst, en,
    input  wire [W-1:0]  s_in, c_in,
    input  wire [7:0]    a, b,
    output reg  [W-1:0]  s_out, c_out
);
  wire [14*W-1:0] ops;
  assign ops[0*W +: W] = s_in;
  assign ops[1*W +: W] = c_in;
  genvar k;
  generate
    for (k = 0; k < 8; k = k + 1) begin : g_row
      // widen to W bits BEFORE shifting (an 8-bit shift would truncate)
      assign ops[(2+k)*W +: W] = {{(W-8){1'b0}}, a & {8{b[k]}}} << k;
    end
  endgenerate
  // -( (a & {8{b7}}) << 8 )  as ~Y (the +1 lives in the shared +2 constant)
  assign ops[10*W +: W] = {3'b111, ~(a & {8{b[7]}}), 8'hFF};
  // -( (b & {8{a7}}) << 8 )
  assign ops[11*W +: W] = {3'b111, ~(b & {8{a[7]}}), 8'hFF};
  // +2 : the two end-around +1s of the complemented rows
  assign ops[12*W +: W] = {{(W-2){1'b0}}, 2'b10};
  // +2^16 * a7 * b7
  assign ops[13*W +: W] = {{(W-17){1'b0}}, a[7] & b[7], {16{1'b0}}};

  wire [W-1:0] sn, cn;
  csa_tree14 #(.W(W)) u_tree (.din(ops), .s(sn), .c(cn));

  always @(posedge clk) begin
    if (rst) begin
      s_out <= {W{1'b0}};
      c_out <= {W{1'b0}};
    end else if (en) begin
      s_out <= sn;
      c_out <= cn;
    end
  end
endmodule

module tpu #(parameter N = 8, DW = 8, ACCW = 32) (
    input  wire            clk,
    input  wire            rst,
    input  wire            w_valid,
    input  wire [N*DW-1:0] w_data,
    input  wire            a_valid,
    input  wire [N*DW-1:0] a_data,
    output wire            c_valid,
    output wire [N*ACCW-1:0] c_data
);

  localparam LAT = 2*N-1;    // 15
  localparam W   = 2*DW+3;   // 19
  wire en = ~w_valid;

  genvar gi, gj, gk;

  // ---------------- activation distribution ----------------
  wire [DW-1:0] a_h [0:N-1][0:N-1];

  generate
    for (gi = 0; gi < N; gi = gi + 1) begin : g_arow
      wire [DW-1:0] din = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      wire [DW-1:0] s0;
      if (gi == 0) begin : g_sk0
        assign s0 = din;
      end else begin : g_sk
        reg [DW-1:0] sr [0:gi-1];
        integer k;
        always @(posedge clk) begin
          if (rst) begin
            for (k = 0; k < gi; k = k + 1) sr[k] <= {DW{1'b0}};
          end else if (en) begin
            sr[0] <= din;
            for (k = 1; k < gi; k = k + 1) sr[k] <= sr[k-1];
          end
        end
        assign s0 = sr[gi-1];
      end
      assign a_h[gi][0] = s0;

      reg [DW-1:0] ech [0:N-2];
      integer e;
      always @(posedge clk) begin
        if (rst) begin
          for (e = 0; e < N-1; e = e + 1) ech[e] <= {DW{1'b0}};
        end else if (en) begin
          ech[0] <= s0;
          for (e = 1; e < N-1; e = e + 1) ech[e] <= ech[e-1];
        end
      end
      for (gk = 1; gk < N; gk = gk + 1) begin : g_ej
        assign a_h[gi][gk] = ech[gk-1];
      end
    end

    // ---------------- per-column weights + carry-save MAC chains ----------
    for (gj = 0; gj < N; gj = gj + 1) begin : g_col

      reg [DW-1:0] wj [0:N-1];
      integer s;
      always @(posedge clk) begin
        if (rst) begin
          for (s = 0; s < N; s = s + 1) wj[s] <= {DW{1'b0}};
        end else if (w_valid) begin
          wj[0] <= w_data[gj*DW +: DW];
          for (s = 1; s < N; s = s + 1) wj[s] <= wj[s-1];
        end
      end

      wire [W-1:0] Sv [0:N-1];
      wire [W-1:0] Cv [0:N-1];

      tpu_mac #(.W(W)) u_st0 (.clk(clk), .rst(rst), .en(en),
          .s_in({W{1'b0}}), .c_in({W{1'b0}}),
          .a(a_h[0][gj]), .b(wj[0]), .s_out(Sv[0]), .c_out(Cv[0]));

      for (gi = 1; gi < N; gi = gi + 1) begin : g_st
        tpu_mac #(.W(W)) u_st (.clk(clk), .rst(rst), .en(en),
            .s_in(Sv[gi-1]), .c_in(Cv[gi-1]),
            .a(a_h[gi][gj]), .b(wj[gi]), .s_out(Sv[gi]), .c_out(Cv[gi]));
      end

      // resolve the carry-save pair once (19-bit CPA)
      wire [W-1:0] res = Sv[N-1] + Cv[N-1];

      // output de-skew: align all columns of one result row
      wire [W-1:0] cout;
      if (N-1-gj == 0) begin : g_d0
        assign cout = res;
      end else begin : g_dn
        reg [W-1:0] dr [0:N-2-gj];
        integer d;
        always @(posedge clk) begin
          if (rst) begin
            for (d = 0; d < N-1-gj; d = d + 1) dr[d] <= {W{1'b0}};
          end else if (en) begin
            dr[0] <= res;
            for (d = 1; d < N-1-gj; d = d + 1) dr[d] <= dr[d-1];
          end
        end
        assign cout = dr[N-2-gj];
      end
      assign c_data[gj*ACCW +: ACCW] = {{(ACCW-W){cout[W-1]}}, cout};
    end
  endgenerate

  // ---------------- valid pipeline ----------------
  reg [LAT-1:0] vpipe;
  always @(posedge clk) begin
    if (rst) vpipe <= {LAT{1'b0}};
    else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
  end
  assign c_valid = vpipe[LAT-1];

endmodule

// v11: carry-save systolic array.
//  * each PE: partial products of a*w (8 shifted copies of a / ~a selected by
//    weight bits, plus one correction bit) compressed TOGETHER with the two
//    incoming carry-save rows in a 5-level CSA tree -> (qs,qc) rows (19b).
//    No CPA anywhere in the per-cycle path.
//  * one 19-bit CPA per column at the end (with a full cycle of slack for
//    columns 0..6, combinational at the pin for column 7).
//  * all arithmetic is modulo 2^19 (the result fits in 19 bits), so the
//    shifted-carry trick {c[17:0],1'b0} is exact.
module tpu_pe_cs (
  input gclk,
  input [18:0] ps, pc,   // incoming carry-save rows (partial = ps + 2*pc)
  input [7:0] a, w,      // activation (signed), weight (signed)
  output [18:0] qs, qc
);
  wire signed [18:0] sa = {{11{a[7]}}, a};
  wire signed [18:0] na = ~sa;                 // = sext(~a) = -a-1
  wire [18:0] b  = {pc[17:0], 1'b0};           // 2*pc (mod 2^19)
  // partial products: a*w = sum_{k=0..6} w_k*(a<<k)  +  w_7*((~a)<<7) + w_7<<7
  wire [18:0] t0 = {19{w[0]}} & sa;
  wire [18:0] t1 = {19{w[1]}} & (sa <<< 1);
  wire [18:0] t2 = {19{w[2]}} & (sa <<< 2);
  wire [18:0] t3 = {19{w[3]}} & (sa <<< 3);
  wire [18:0] t4 = {19{w[4]}} & (sa <<< 4);
  wire [18:0] t5 = {19{w[5]}} & (sa <<< 5);
  wire [18:0] t6 = {19{w[6]}} & (sa <<< 6);
  wire [18:0] t7 = {19{w[7]}} & (na <<< 7);
  wire [18:0] tb = {19{w[7]}} & 19'h0080;      // correction bit at position 7

  // level 1
  wire [18:0] a1 = ps ^ b ^ t0;   wire [18:0] c1 = (ps & b) | (t0 & (ps | b));
  wire [18:0] a2 = t1 ^ t2 ^ t3;  wire [18:0] c2 = (t1 & t2) | (t3 & (t1 | t2));
  wire [18:0] a3 = t4 ^ t5 ^ t6;  wire [18:0] c3 = (t4 & t5) | (t6 & (t4 | t5));
  wire [18:0] b1 = {c1[17:0], 1'b0};
  wire [18:0] b2 = {c2[17:0], 1'b0};
  wire [18:0] b3 = {c3[17:0], 1'b0};
  // level 2
  wire [18:0] d1 = a1 ^ b1 ^ a2;  wire [18:0] e1 = (a1 & b1) | (a2 & (a1 | b1));
  wire [18:0] d2 = b2 ^ a3 ^ b3;  wire [18:0] e2 = (b2 & a3) | (b3 & (b2 | a3));
  wire [18:0] f1 = {e1[17:0], 1'b0};
  wire [18:0] f2 = {e2[17:0], 1'b0};
  // level 3
  wire [18:0] g1 = d1 ^ f1 ^ d2;  wire [18:0] h1 = (d1 & f1) | (d2 & (d1 | f1));
  wire [18:0] g2 = f2 ^ t7 ^ tb;  wire [18:0] h2 = (f2 & t7) | (tb & (f2 | t7));
  wire [18:0] i1 = {h1[17:0], 1'b0};
  wire [18:0] i2 = {h2[17:0], 1'b0};
  // level 4
  wire [18:0] k1 = g1 ^ i1 ^ g2;  wire [18:0] l1 = (g1 & i1) | (g2 & (g1 | i1));
  wire [18:0] m1 = {l1[17:0], 1'b0};
  // level 5
  wire [18:0] sn = k1 ^ m1 ^ i2;
  wire [18:0] cn = (k1 & m1) | (i2 & (k1 | m1));

  reg [18:0] qs, qc;
  always @(posedge gclk) begin
    qs <= sn; qc <= cn;
  end
endmodule

module tpu (
  input clk, rst, w_valid, a_valid,
  input [63:0] w_data, a_data,
  output c_valid,
  output [255:0] c_data
);
  wire gclk = clk & ~w_valid;
  wire wclk = clk &  w_valid;
  wire en   = ~w_valid;

  // ---------------- weights: per-column south-shifting registers ----------
  reg [7:0] wq [0:7][0:7];
  integer ii, jj;
  always @(posedge wclk) begin
    for (jj = 0; jj < 8; jj = jj+1) begin
      wq[0][jj] <= w_data[jj*8 +: 8];
      for (ii = 1; ii < 8; ii = ii+1)
        wq[ii][jj] <= wq[ii-1][jj];
    end
  end

  // ---------------- activation skew --------------------------------------
  wire [7:0] act [0:7][0:7];
  genvar gi, gj;
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_arow
      wire [7:0] din = a_valid ? a_data[gi*8 +: 8] : 8'h00;
      reg [7:0] c [0:gi+6];
      integer q;
      always @(posedge gclk) begin
        c[0] <= din;
        for (q = 1; q <= gi+6; q = q+1) c[q] <= c[q-1];
      end
      for (gj = 0; gj < 8; gj = gj+1) begin : g_tap
        if (gi + gj == 0) begin : g_t0
          assign act[gi][gj] = din;
        end else begin : g_tn
          assign act[gi][gj] = c[gi+gj-1];
        end
      end
    end
  endgenerate

  // ---------------- PE mesh (carry-save partial sums flow south) ----------
  wire [18:0] ps [0:7][0:7];
  wire [18:0] pc [0:7][0:7];
  generate
    for (gi = 0; gi < 8; gi = gi+1) begin : g_pe_row
      for (gj = 0; gj < 8; gj = gj+1) begin : g_pe_col
        tpu_pe_cs u_pe (
          .gclk (gclk),
          .ps   (gi == 0 ? 19'd0 : ps[gi-1][gj]),
          .pc   (gi == 0 ? 19'd0 : pc[gi-1][gj]),
          .a    (act[gi][gj]),
          .w    (wq[gi][gj]),
          .qs   (ps[gi][gj]),
          .qc   (pc[gi][gj]));
      end
    end
  endgenerate

  // ---------------- final CPA + de-skew -----------------------------------
  generate
    for (gj = 0; gj < 8; gj = gj+1) begin : g_out
      wire [18:0] r = ps[7][gj] + {pc[7][gj][17:0], 1'b0};
      if (gj == 7) begin : g_comb
        assign c_data[255 -: 32] = {{13{r[18]}}, r};
      end else begin : g_reg
        reg [18:0] r0;
        always @(posedge gclk) r0 <= r;
        if (gj == 6) begin : g_d0
          assign c_data[gj*32 +: 32] = {{13{r0[18]}}, r0};
        end else begin : g_dn
          reg [18:0] dr [0:6-gj-1];
          integer k;
          always @(posedge gclk) begin
            dr[0] <= r0;
            for (k = 1; k < 6-gj; k = k+1) dr[k] <= dr[k-1];
          end
          wire [18:0] last = dr[6-gj-1];
          assign c_data[gj*32 +: 32] = {{13{last[18]}}, last};
        end
      end
    end
  endgenerate

  // ---------------- valid pipeline ----------------------------------------
  reg [14:0] vpipe;
  always @(posedge clk or posedge rst) begin
    if (rst) vpipe <= 15'h0;
    else if (en) vpipe <= {vpipe[13:0], a_valid};
  end
  assign c_valid = vpipe[14];

endmodule

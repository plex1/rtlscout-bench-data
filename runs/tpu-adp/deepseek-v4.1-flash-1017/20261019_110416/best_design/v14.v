// v14: exact golden staging, redundant accumulator, Booth radix-4 fused MAC.
//
// Same structure as v12 (see below) but the multiply uses radix-4 Booth
// encoding: 4 partial products instead of 8, so the per-stage compressor tree
// is only 3 levels deep (6 addends -> 4 -> 3 -> 2) instead of 5.
//
// Golden per-PE staging (a weight reload freezes every partial sum, so the
// 8-deep linear chain is mandatory):
//   acc(i,j) = acc(i-1,j) + A(T-i-j)[i]*B[i][j]
// Intermediate values are carried in redundant (sum,carry) form, so each stage
// is one compressor tree and only ONE carry-propagate adder per column is
// needed at the very end.

module csa3 #(parameter W = 19)
   (input [W-1:0] x, y, z,
    output [W-1:0] s,
    output [W-1:0] c);
   wire [W-1:0] m = (x & y) | (x & z) | (y & z);
   assign s = x ^ y ^ z;
   assign c = {m[W-2:0], 1'b0};
endmodule

// (s0,c0) + a*w with radix-4 Booth  ->  (so,co), all redundant.
module tpu_booth #(parameter W = 19)
   (input [W-1:0]	s0,
    input [W-1:0]	c0,
    input [7:0]		a,
    input [7:0]		w,
    output [W-1:0]	so,
    output [W-1:0]	co);

   wire [W-1:0] sext = {{(W-8){w[7]}}, w};

   wire [W-1:0] pp [0:3];
   genvar	k;
   generate
      for (k = 0; k < 4; k = k + 1) begin : g_b
         wire b1 = a[2*k+1];
         wire b0 = a[2*k];
         wire bm = (k == 0) ? 1'b0 : a[2*k-1];
         wire is1 = b0 ^ bm;
         wire is2 = (b0 & bm & ~b1) | (~b0 & ~bm & b1);
         wire zero = ~(b1 | b0 | bm) | (b1 & b0 & bm);
         wire neg = b1;
         wire [W-1:0] base = is2 ? (sext <<< 1) : sext;
         wire [W-1:0] val  = neg ? (~base + {{(W-1){1'b0}}, 1'b1}) : base;
         assign pp[k] = zero ? {W{1'b0}} : (val <<< (2*k));
      end
   endgenerate

   wire [W-1:0] A0s, A0c, A1s, A1c;
   wire [W-1:0] B0s, B0c;
   wire [W-1:0] C0s, C0c;

   // 6 addends -> 4
   csa3 #(.W(W)) u0 (.x(s0),    .y(c0),    .z(pp[0]), .s(A0s), .c(A0c));
   csa3 #(.W(W)) u1 (.x(pp[1]), .y(pp[2]), .z(pp[3]), .s(A1s), .c(A1c));
   // 4 -> 3
   csa3 #(.W(W)) u2 (.x(A0s),   .y(A0c),   .z(A1s),   .s(B0s), .c(B0c));
   // 3 -> 2
   csa3 #(.W(W)) u3 (.x(B0s),   .y(B0c),   .z(A1c),   .s(C0s), .c(C0c));

   assign so = C0s;
   assign co = C0c;
endmodule

module tpu
  #(parameter N = 8, DW = 8, ACCW = 32)
   (
    input		 clk,
    input		 rst,
    input		 w_valid,
    input [N*DW-1:0]	 w_data,
    input		 a_valid,
    input [N*DW-1:0]	 a_data,
    output		 c_valid,
    output [N*ACCW-1:0]	 c_data
    );

   localparam LAT = 2*N - 1;	// 15
   localparam AW  = 19;

   wire		      ld = w_valid;
   wire		      en = ~w_valid;

   genvar	      i, j;
   integer	      kk;

   //--------------------------------------------------------------------------
   // Activation skew: a_tap[i][j] = A row byte i delayed by (i+j) compute cycles.
   //--------------------------------------------------------------------------
   wire [DW-1:0] a_tap [0:N-1][0:N-1];
   generate
      for (i = 0; i < N; i = i + 1) begin : g_adly
         wire [DW-1:0] din = a_valid ? a_data[i*DW +: DW] : {DW{1'b0}};
         if (i + N - 1 == 0) begin : g_none
            assign a_tap[i][0] = din;
         end
         else begin : g_sr
            reg [DW-1:0] sr [0:i+N-2];
            always @(posedge clk) begin
               if (rst) begin
                  for (kk = 0; kk < i + N - 1; kk = kk + 1)
                    sr[kk] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (kk = 1; kk < i + N - 1; kk = kk + 1)
                    sr[kk] <= sr[kk-1];
               end
            end
            for (j = 0; j < N; j = j + 1) begin : g_tap
               if (i + j == 0) begin : g_t0
                  assign a_tap[i][j] = din;
               end
               else begin : g_tn
                  assign a_tap[i][j] = sr[i+j-1];
               end
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Stationary weights (shift south, bottom matrix row first).
   //--------------------------------------------------------------------------
   reg [DW-1:0] wreg [0:N-1][0:N-1];
   generate
      for (j = 0; j < N; j = j + 1) begin : g_wload
         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1)
                 wreg[kk][j] <= {DW{1'b0}};
            end
            else if (ld) begin
               wreg[0][j] <= w_data[j*DW +: DW];
               for (kk = 1; kk < N; kk = kk + 1)
                 wreg[kk][j] <= wreg[kk-1][j];
            end
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Per-column fused MAC cascade + output de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire [AW-1:0] ns [0:N-1];
         wire [AW-1:0] nc [0:N-1];
         reg  [AW-1:0] sr [0:N-1];
         reg  [AW-1:0] cr [0:N-1];

         for (i = 0; i < N; i = i + 1) begin : g_st
            tpu_booth #(.W(AW)) u_s
              (.s0(i == 0 ? {AW{1'b0}} : sr[i-1]),
               .c0(i == 0 ? {AW{1'b0}} : cr[i-1]),
               .a(a_tap[i][j]), .w(wreg[i][j]), .so(ns[i]), .co(nc[i]));
            always @(posedge clk) begin
               if (rst) begin
                  sr[i] <= {AW{1'b0}};
                  cr[i] <= {AW{1'b0}};
               end
               else if (en) begin
                  sr[i] <= ns[i];
                  cr[i] <= nc[i];
               end
            end
         end

         wire [AW-1:0] fin = sr[N-1] + cr[N-1];

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){fin[AW-1]}}, fin};
         end
         else begin : g_dsn
            reg [AW-1:0] ds [0:D-1];
            integer      m;
            always @(posedge clk) begin
               if (rst) begin
                  for (m = 0; m < D; m = m + 1)
                    ds[m] <= {AW{1'b0}};
               end
               else if (en) begin
                  ds[0] <= fin;
                  for (m = 1; m < D; m = m + 1)
                    ds[m] <= ds[m-1];
               end
            end
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){ds[D-1][AW-1]}}, ds[D-1]};
         end
      end
   endgenerate

   //--------------------------------------------------------------------------
   // Valid pipeline.
   //--------------------------------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

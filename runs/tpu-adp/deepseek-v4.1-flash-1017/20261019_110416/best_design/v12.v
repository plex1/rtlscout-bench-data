// v12: exact golden staging, every accumulate stage is ONE fused carry-save
// reduction (partial products of a*w plus the previous redundant partial sum).
//
// Golden per-PE staging (a weight reload freezes every partial sum, so the
// 8-deep linear chain is mandatory):
//   acc(i,j) = acc(i-1,j) + A(T-i-j)[i]*B[i][j]
// The intermediate values are internal, so they are carried in redundant
// (sum,carry) form: each stage is a 10-input Wallace reduction and only ONE
// carry-propagate adder per column is needed at the very end.

module csa3 #(parameter W = 19)
   (input [W-1:0] x, y, z,
    output [W-1:0] s,
    output [W-1:0] c);
   wire [W-1:0] m = (x & y) | (x & z) | (y & z);
   assign s = x ^ y ^ z;
   assign c = {m[W-2:0], 1'b0};
endmodule

// (s0,c0) + a*w  ->  (so,co), all in redundant form.
module tpu_csa #(parameter W = 19)
   (input [W-1:0]	s0,
    input [W-1:0]	c0,
    input [7:0]		a,
    input [7:0]		w,
    output [W-1:0]	so,
    output [W-1:0]	co);

   wire [W-1:0] sext = {{(W-8){w[7]}}, w};
   wire [W-1:0] pp [0:7];
   genvar	b;
   generate
      for (b = 0; b < 7; b = b + 1) begin : g_pp
         assign pp[b] = a[b] ? (sext << b) : {W{1'b0}};
      end
      assign pp[7] = a[7] ? ((~sext + {{(W-1){1'b0}}, 1'b1}) << 7)
                          : {W{1'b0}};
   endgenerate

   wire [W-1:0] A0s,A0c,A1s,A1c,A2s,A2c;
   wire [W-1:0] B0s,B0c,B1s,B1c;
   wire [W-1:0] C0s,C0c,C1s,C1c;
   wire [W-1:0] D0s,D0c;
   wire [W-1:0] E0s,E0c;

   // 10 addends -> 7
   csa3 #(.W(W)) u0 (.x(s0),    .y(c0),    .z(pp[0]), .s(A0s), .c(A0c));
   csa3 #(.W(W)) u1 (.x(pp[1]), .y(pp[2]), .z(pp[3]), .s(A1s), .c(A1c));
   csa3 #(.W(W)) u2 (.x(pp[4]), .y(pp[5]), .z(pp[6]), .s(A2s), .c(A2c));
   // 7 -> 5
   csa3 #(.W(W)) u3 (.x(A0s),   .y(A0c),   .z(A1s),   .s(B0s), .c(B0c));
   csa3 #(.W(W)) u4 (.x(A1c),   .y(A2s),   .z(A2c),   .s(B1s), .c(B1c));
   // 5 -> 4
   csa3 #(.W(W)) u5 (.x(B0s),   .y(B0c),   .z(B1s),   .s(C0s), .c(C0c));
   // 4 -> 3
   csa3 #(.W(W)) u6 (.x(C0s),   .y(C0c),   .z(B1c),   .s(D0s), .c(D0c));
   // 3 -> 2
   csa3 #(.W(W)) u7 (.x(D0s),   .y(D0c),   .z(pp[7]), .s(E0s), .c(E0c));

   assign so = E0s;
   assign co = E0c;
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
   // Per-column fused multiply-accumulate cascade + output de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire [AW-1:0] ns [0:N-1];
         wire [AW-1:0] nc [0:N-1];
         reg  [AW-1:0] sr [0:N-1];
         reg  [AW-1:0] cr [0:N-1];

         for (i = 0; i < N; i = i + 1) begin : g_st
            tpu_csa #(.W(AW)) u_s
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

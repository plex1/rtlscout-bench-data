// v11: exact golden staging, fused multiply-accumulate as a shallow Wallace
// tree over explicit partial products (no separate multiplier carry chain).
//
// Golden per-PE staging (a reload freezes every partial sum, so the linear
// 8-stage chain is mandatory):
//   acc(i,j) = acc(i-1,j) + A(T-i-j)[i]*B[i][j]
// Accumulators are kept in redundant (sum,carry) form, so a stage is a single
// compressor tree fed by {s_prev, c_prev, 8 partial products}.

module csa3 #(parameter W = 19)
   (input [W-1:0] x, y, z,
    output [W-1:0] s,
    output [W-1:0] c);
   wire [W-1:0] m = (x & y) | (x & z) | (y & z);
   assign s = x ^ y ^ z;
   assign c = {m[W-2:0], 1'b0};
endmodule

// One fused multiply-accumulate stage: (s0,c0) + a*w  ->  (so,co)
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
      assign pp[7] = a[7] ? (~(sext << 7) + {{(W-1){1'b0}}, 1'b1})
                          : {W{1'b0}};
   endgenerate

   // 10-input Wallace reduction
   wire [W-1:0] a0s, a0c, a1s, a1c, a2s, a2c;
   wire [W-1:0] b0s, b0c, b1s, b1c;
   wire [W-1:0] c0s, c0c;
   wire [W-1:0] d0s, d0c;
   wire [W-1:0] e0s, e0c;

   csa3 #(.W(W)) u0 (.x(s0),   .y(c0),   .z(pp[0]), .s(a0s), .c(a0c));
   csa3 #(.W(W)) u1 (.x(pp[1]),.y(pp[2]), .z(pp[3]), .s(a1s), .c(a1c));
   csa3 #(.W(W)) u2 (.x(pp[4]),.y(pp[5]), .z(pp[6]), .s(a2s), .c(a2c));

   csa3 #(.W(W)) u3 (.x(a0s),  .y(a0c),  .z(a1s),   .s(b0s), .c(b0c));
   csa3 #(.W(W)) u4 (.x(a1c),  .y(a2s),  .z(a2c),   .s(b1s), .c(b1c));

   csa3 #(.W(W)) u5 (.x(b0s),  .y(b0c),  .z(b1s),   .s(c0s), .c(c0c));

   csa3 #(.W(W)) u6 (.x(c0s),  .y(c0c),  .z(b1c),   .s(d0s), .c(d0c));
   csa3 #(.W(W)) u7 (.x(d0s),  .y(d0c),  .z(pp[7]), .s(e0s), .c(e0c));

   assign so = e0s;
   assign co = e0c;
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
   localparam AW = 19;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         wire [AW-1:0] sv [0:N-1];
         wire [AW-1:0] cv [0:N-1];

         tpu_csa #(.W(AW)) u_s0
           (.s0({AW{1'b0}}), .c0({AW{1'b0}}),
            .a(a_tap[0][j]), .w(wreg[0][j]), .so(sv[0]), .co(cv[0]));

         for (i = 1; i < N; i = i + 1) begin : g_st
            reg [AW-1:0] sr, cr;
            wire [AW-1:0] ns, nc;
            tpu_csa #(.W(AW)) u_s
              (.s0(sv[i-1]), .c0(cv[i-1]),
               .a(a_tap[i][j]), .w(wreg[i][j]), .so(ns), .co(nc));
            assign sv[i] = sr;
            assign cv[i] = cr;
            always @(posedge clk) begin
               if (rst) begin
                  sr <= {AW{1'b0}};
                  cr <= {AW{1'b0}};
               end
               else if (en) begin
                  sr <= ns;
                  cr <= nc;
               end
            end
         end
         // stage 0 registers
         reg [AW-1:0] s0r, c0r;
         always @(posedge clk) begin
            if (rst) begin
               s0r <= {AW{1'b0}};
               c0r <= {AW{1'b0}};
            end
            else if (en) begin
               s0r <= sv[0];
               c0r <= cv[0];
            end
         end
         wire [AW-1:0] fin = sv[N-1] + cv[N-1];

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

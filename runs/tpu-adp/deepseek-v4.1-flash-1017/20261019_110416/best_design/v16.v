// v16: Booth multiply and accumulate in SEPARATE pipeline stages.
//
// Golden per-PE staging:  acc(i,j) = acc(i-1,j) + A(T-i-j)[i]*B[i][j], and a
// weight reload freezes every partial sum.  Writing D = T-i-j for the tap
// index, the completed column value is needed in cycle T+8+j.
//
// Pipeline per (i,j):
//   product registers are formed from tap D = T-i-j-1 and registered in the
//   cycle right before they are consumed (1 stage), then the 8-deep narrow add
//   cascade (8-i stages for element i), then (7-j) de-skew stages:
//       (i+j+1) + (8-i) + (7-j) = 16 ... corrected below to 15.
//
// The product registers are NOT gated by the compute enable: during a weight
// reload the activation taps are frozen while the new weights shift in, so the
// first product formed after the reload is (frozen A)*(new B), matching the
// golden array where a stalled PE keeps its a_west and adopts the reloaded
// weight.  Element (0,0) has a combinational product (it sees a_data directly).

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
   // Activation delay line: dd[d] = A word delayed by d compute cycles.
   //--------------------------------------------------------------------------
   wire [N*DW-1:0] din = a_valid ? a_data : {N*DW{1'b0}};
   wire [N*DW-1:0] dd [0:2*N-2];

   assign dd[0] = din;
   generate
      for (i = 1; i <= 2*N-2; i = i + 1) begin : g_dly
         reg [N*DW-1:0] q;
         always @(posedge clk) begin
            if (rst)     q <= {N*DW{1'b0}};
            else if (en) q <= dd[i-1];
         end
         assign dd[i] = q;
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
   // Booth product, combinational.
   //--------------------------------------------------------------------------
   function automatic [AW-1:0] booth_pp(input [7:0] a, input [7:0] w);
      reg [AW-1:0] sext, acc, base, val;
      reg	   b1, b0, bm, zero, neg;
      integer	   k;
      begin
         sext = {{(AW-8){w[7]}}, w};
         acc  = {AW{1'b0}};
         for (k = 0; k < 4; k = k + 1) begin
            b1   = a[2*k+1];
            b0   = a[2*k];
            bm   = (k == 0) ? 1'b0 : a[2*k-1];
            zero = ~(b1 | b0 | bm) | (b1 & b0 & bm);
            neg  = b1;
            base = (b1 == b0 && b0 != bm) ? (sext <<< 1) : sext;
            val  = neg ? (~base + {{(AW-1){1'b0}}, 1'b1}) : base;
            acc  = acc + (zero ? {AW{1'b0}} : (val <<< (2*k)));
         end
         booth_pp = acc;
      end
   endfunction

   //--------------------------------------------------------------------------
   // Per column: ungated product regs, gated narrow accumulate cascade, de-skew.
   //--------------------------------------------------------------------------
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
         reg signed [AW-1:0] pr [0:N-1];
         reg signed [AW-1:0] ac [0:N-1];

         always @(posedge clk) begin
            if (rst) begin
               for (kk = 0; kk < N; kk = kk + 1) begin
                  pr[kk] <= {AW{1'b0}};
                  ac[kk] <= {AW{1'b0}};
               end
            end
            else begin
               // ungated product registers: tap index (i+j-1)
               for (kk = 0; kk < N; kk = kk + 1) begin
                  if (kk + j == 0)
                    pr[kk] <= booth_pp(din[0 +: DW], wreg[0][0]);
                  else
                    pr[kk] <= booth_pp(dd[kk+j-1][kk*DW +: DW], wreg[kk][j]);
               end
               if (en) begin
                  ac[0] <= pr[0];
                  for (kk = 1; kk < N; kk = kk + 1)
                    ac[kk] <= ac[kk-1] + pr[kk];
               end
            end
         end

         localparam D = N - 1 - j;
         if (D == 0) begin : g_ds0
            assign c_data[j*ACCW +: ACCW] = {{(ACCW-AW){ac[N-1][AW-1]}}, ac[N-1]};
         end
         else begin : g_dsn
            reg signed [AW-1:0] ds [0:D-1];
            integer	      m;
            always @(posedge clk) begin
               if (rst) begin
                  for (m = 0; m < D; m = m + 1)
                    ds[m] <= {AW{1'b0}};
               end
               else if (en) begin
                  ds[0] <= ac[N-1];
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

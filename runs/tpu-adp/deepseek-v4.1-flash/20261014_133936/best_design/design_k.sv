// Variant K: activations delayed uniformly by (LAT-2)=13 compute cycles, then
// all 8 products of a column formed simultaneously, summed with a
// combinational 3-level signed adder tree, and captured in ONE sum register per
// column -> available exactly at cycle n+LAT.  No output de-skew chains, no
// systolic skew network: just a single (LAT-2)-deep activation pipe.
module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		clk,
    input		rst,
    input		w_valid,
    input [N*DW-1:0]	w_data,
    input		a_valid,
    input [N*DW-1:0]	a_data,
    output		c_valid,
    output [N*ACCW-1:0]	c_data
    );

   localparam LAT = 2*N - 1;		// 15
   localparam DLY = LAT - 2;		// activation pipe depth (13)
   localparam IW  = 19;
   localparam PW  = 2*DW;

   wire ld = w_valid;
   wire en = ~w_valid;

   integer i, j, k;
   genvar gi, gj;

   // ---- stationary weights -------------------------------------------------
   reg [DW-1:0] wq [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) wq[i][j] <= {DW{1'b0}};
      end
      else if (ld) begin
         for (j = 0; j < N; j = j + 1) wq[0][j] <= w_data[j*DW +: DW];
         for (i = 1; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) wq[i][j] <= wq[i-1][j];
      end
   end

   // ---- activation pipe: zeros while a_valid low, DLY-cycle delay ---------
   wire [DW-1:0] din [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   reg  [DW-1:0] ap [0:N-1][0:DLY-1];
   wire [DW-1:0] aw [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_aw
         assign aw[gi] = ap[gi][DLY-1];
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (k = 0; k < DLY; k = k + 1) ap[i][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1) begin
            ap[i][0] <= din[i];
            for (k = 1; k < DLY; k = k + 1) ap[i][k] <= ap[i][k-1];
         end
      end
   end

   // ---- products, registered ----------------------------------------------
   reg [PW-1:0] P [0:N-1][0:N-1];	// P[j][i]
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1) P[j][i] <= {PW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N; i = i + 1)
             P[j][i] <= $signed(aw[i]) * $signed(wq[i][j]);
      end
   end

   // ---- combinational 3-level signed adder tree ---------------------------
   wire signed [IW-1:0] t0 [0:N-1], t1 [0:N-1], t2 [0:N-1], t3 [0:N-1];
   wire signed [IW-1:0] u0 [0:N-1], u1 [0:N-1];
   wire signed [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_tr
         assign t0[gi] = $signed(P[gi][0]) + $signed(P[gi][1]);
         assign t1[gi] = $signed(P[gi][2]) + $signed(P[gi][3]);
         assign t2[gi] = $signed(P[gi][4]) + $signed(P[gi][5]);
         assign t3[gi] = $signed(P[gi][6]) + $signed(P[gi][7]);
         assign u0[gi] = t0[gi] + t1[gi];
         assign u1[gi] = t2[gi] + t3[gi];
         assign csum[gi] = u0[gi] + u1[gi];
      end
   endgenerate

   // ---- one sum register per column (available at cycle n+LAT) -------------
   reg [IW-1:0] sumreg [0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) sumreg[j] <= csum[j];
      end
   end

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_out
         assign c_data[gi*ACCW +: ACCW] =
                  {{(ACCW-IW){sumreg[gi][IW-1]}}, sumreg[gi]};
      end
   endgenerate

   // ---- valid pipeline -----------------------------------------------------
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

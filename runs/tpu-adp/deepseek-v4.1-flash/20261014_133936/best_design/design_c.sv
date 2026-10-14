// Variant C: systolic activation chain (one register per PE) + split
// multiply/accumulate pipeline stages. Same cycle-exact product/weight timing.
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

   localparam LAT = 2*N - 1;
   localparam IW  = 19;
   localparam PW  = 2*DW;

   wire ld = w_valid;
   wire en = ~w_valid;

   integer i, j, k;

   // stationary weights (shift south while loading)
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

   // activation injection (zeros when a_valid low)
   wire [DW-1:0] din [0:N-1];
   genvar gi, gj;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_din
         assign din[gi] = a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
      end
   endgenerate

   // row skew: arow[i] = din[i] delayed i compute cycles
   reg  [DW-1:0] sk [0:N-1][0:N-1];
   wire [DW-1:0] arow [0:N-1];
   assign arow[0] = din[0];
   generate
      for (gi = 1; gi < N; gi = gi + 1) begin : g_arow
         assign arow[gi] = sk[gi][gi-1];
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 1; i < N; i = i + 1)
           for (k = 0; k < i; k = k + 1) sk[i][k] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 1; i < N; i = i + 1) begin
            sk[i][0] <= din[i];
            for (k = 1; k < i; k = k + 1) sk[i][k] <= sk[i][k-1];
         end
      end
   end

   // east activation chain: a_h[i][0]=arow[i], a_h[i][j+1]=reg(a_h[i][j])
   wire [DW-1:0] a_h [0:N-1][0:N];
   reg  [DW-1:0] areg [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_ah
         assign a_h[gi][0] = arow[gi];
         for (gj = 0; gj < N; gj = gj + 1) begin : g_ahj
            assign a_h[gi][gj+1] = areg[gi][gj];
         end
      end
   endgenerate
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) areg[i][j] <= {DW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) areg[i][j] <= a_h[i][j];
      end
   end

   // products, registered at the PE stage
   wire [PW-1:0] prodc [0:N-1][0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_pi
         for (gj = 0; gj < N; gj = gj + 1) begin : g_pj
            assign prodc[gi][gj] = $signed(a_h[gi][gj]) * $signed(wq[gi][gj]);
         end
      end
   endgenerate

   reg [PW-1:0] P [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) P[i][j] <= {PW{1'b0}};
      end
      else if (en) begin
         for (i = 0; i < N; i = i + 1)
           for (j = 0; j < N; j = j + 1) P[i][j] <= prodc[i][j];
      end
   end

   // per-column accumulation chain (7 registered stages)
   reg [IW-1:0] acc [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (i = 0; i < N-1; i = i + 1) acc[j][i] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            acc[j][0] <= {{(IW-PW){P[0][j][PW-1]}}, P[0][j]};
            for (i = 1; i < N-1; i = i + 1)
              acc[j][i] <= acc[j][i-1] + {{(IW-PW){P[i][j][PW-1]}}, P[i][j]};
         end
      end
   end

   wire [IW-1:0] csum [0:N-1];
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_csum
         assign csum[gi] = acc[gi][N-2] + {{(IW-PW){P[N-1][gi][PW-1]}}, P[N-1][gi]};
      end
   endgenerate

   // output de-skew: column j delayed N-1-j
   reg [IW-1:0] ds [0:N-1][0:N-1];
   always @(posedge clk) begin
      if (rst) begin
         for (j = 0; j < N; j = j + 1)
           for (k = 0; k < N-1-j; k = k + 1) ds[j][k] <= {IW{1'b0}};
      end
      else if (en) begin
         for (j = 0; j < N; j = j + 1) begin
            if (N-1-j > 0) begin
               ds[j][0] <= csum[j];
               for (k = 1; k < N-1-j; k = k + 1) ds[j][k] <= ds[j][k-1];
            end
         end
      end
   end

   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_out
         if (N-1-gi == 0)
           assign c_data[gi*ACCW +: ACCW] = {{(ACCW-IW){csum[gi][IW-1]}}, csum[gi]};
         else
           assign c_data[gi*ACCW +: ACCW] =
                    {{(ACCW-IW){ds[gi][N-2-gi][IW-1]}}, ds[gi][N-2-gi]};
      end
   endgenerate

   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst) vpipe <= {LAT{1'b0}};
      else if (en) vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

// Parallel-pipeline TPU: same cycle-accurate I/O contract as the logikbench
// systolic tpu, but internally scheduled as: A rows delayed 11 compute cycles,
// then per column: 8 multipliers -> 3 registered 2-input add levels.
// c_valid = a_valid delayed exactly 15 compute cycles.

module tpu
  #(parameter N = 8,
    parameter DW = 8,
    parameter ACCW = 32)
   (
    input		     clk,
    input		     rst,
    input		     w_valid,
    input [N*DW-1:0]	     w_data,
    input		     a_valid,
    input [N*DW-1:0]	     a_data,
    output		     c_valid,
    output [N*ACCW-1:0]	     c_data
    );

   localparam LAT = 2*N - 1;   // 15
   localparam DL  = LAT - 4;   // 11 input delay stages
   localparam PW  = 2*DW;      // 16, product width
   localparam S1W = PW+1;      // 17
   localparam S2W = PW+2;      // 18
   localparam S3W = PW+3;      // 19

   wire en = ~w_valid;   // compute enable (array stalls during weight load)
   wire ld = w_valid;

   //######################################################################
   // Weight storage: per-column shift chain, top loads w_data, shifts south.
   // After 8 load cycles wcell[k][j] holds B[7-k][j], so B[i][j] lives in
   // wcell[N-1-i][j].
   //######################################################################
   reg [DW-1:0] wcell [0:N-1][0:N-1];
   integer k2, j2;
   always @(posedge clk) begin
      if (rst) begin
	 for (k2 = 0; k2 < N; k2 = k2 + 1)
	   for (j2 = 0; j2 < N; j2 = j2 + 1)
	     wcell[k2][j2] <= {DW{1'b0}};
      end
      else if (ld) begin
	 for (j2 = 0; j2 < N; j2 = j2 + 1) begin
	    wcell[0][j2] <= w_data[j2*DW +: DW];
	    for (k2 = 1; k2 < N; k2 = k2 + 1)
	      wcell[k2][j2] <= wcell[k2-1][j2];
	 end
      end
   end

   //######################################################################
   // Activation row delay line (zeros injected while a_valid is low).
   //######################################################################
   reg [N*DW-1:0] adel [0:DL-1];
   integer k3;
   always @(posedge clk) begin
      if (rst) begin
	 for (k3 = 0; k3 < DL; k3 = k3 + 1)
	   adel[k3] <= {N*DW{1'b0}};
      end
      else if (en) begin
	 adel[0] <= a_valid ? a_data : {N*DW{1'b0}};
	 for (k3 = 1; k3 < DL; k3 = k3 + 1)
	   adel[k3] <= adel[k3-1];
      end
   end

   wire [N*DW-1:0] arow = adel[DL-1];

   //######################################################################
   // Compute pipeline: products -> 2-input add tree (3 levels), per column.
   //######################################################################
   reg [PW-1:0]  pr [0:N*N-1]; // pr[j*N+i]
   reg [S1W-1:0] r1 [0:N*N-1]; // r1[j*N+k]
   reg [S2W-1:0] r2 [0:N*N-1]; // r2[j*N+k]
   reg [S3W-1:0] r3 [0:N-1];

   genvar i, j;
   generate
      for (j = 0; j < N; j = j + 1) begin : g_col
	 wire [DW-1:0]  wv   [0:N-1];
	 wire [PW-1:0]  prod [0:N-1];
	 for (i = 0; i < N; i = i + 1) begin : g_pe
	    assign wv[i]   = wcell[i][j];
	    assign prod[i] = $signed(arow[i*DW +: DW]) * $signed(wv[i]);
	 end
	 wire [PW-1:0] p1 [0:N-1];
	 for (i = 0; i < N; i = i + 1) begin : g_pr
	    assign p1[i] = pr[j*N+i];
	 end
	 integer ii;
	 always @(posedge clk) begin
	    if (rst) begin
	       for (ii = 0; ii < N*N; ii = ii + 1) begin
		  pr[ii] <= {PW{1'b0}};
		  r1[ii] <= {S1W{1'b0}};
		  r2[ii] <= {S2W{1'b0}};
	       end
	       r3[j]     <= {S3W{1'b0}};
	    end
	    else if (en) begin
	       for (ii = 0; ii < N; ii = ii + 1)
		 pr[j*N+ii] <= prod[ii];
	       // each add level consumes registered values from the previous
	       r1[j*N+0] <= {{1{p1[0][PW-1]}}, p1[0]} + {{1{p1[1][PW-1]}}, p1[1]};
	       r1[j*N+1] <= {{1{p1[2][PW-1]}}, p1[2]} + {{1{p1[3][PW-1]}}, p1[3]};
	       r1[j*N+2] <= {{1{p1[4][PW-1]}}, p1[4]} + {{1{p1[5][PW-1]}}, p1[5]};
	       r1[j*N+3] <= {{1{p1[6][PW-1]}}, p1[6]} + {{1{p1[7][PW-1]}}, p1[7]};
	       r2[j*N+0] <= {{1{r1[j*N+0][S1W-1]}}, r1[j*N+0]} +
			    {{1{r1[j*N+1][S1W-1]}}, r1[j*N+1]};
	       r2[j*N+1] <= {{1{r1[j*N+2][S1W-1]}}, r1[j*N+2]} +
			    {{1{r1[j*N+3][S1W-1]}}, r1[j*N+3]};
	       r3[j] <= {{1{r2[j*N+0][S2W-1]}}, r2[j*N+0]} +
			{{1{r2[j*N+1][S2W-1]}}, r2[j*N+1]};
	    end
	 end

	 assign c_data[j*ACCW +: ACCW] = {{(ACCW-S3W){r3[j][S3W-1]}}, r3[j]};
      end
   endgenerate

   //######################################################################
   // Valid pipeline: c_valid = a_valid delayed LAT compute cycles.
   //######################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
	vpipe <= {LAT{1'b0}};
      else if (en)
	vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

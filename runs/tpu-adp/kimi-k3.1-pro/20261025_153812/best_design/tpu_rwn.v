module tpu_pe #(parameter DW=8,parameter ACCW=19,parameter BOTTOM=0,parameter OW=19) (
  input clk,input rst,input ld,input en,
  input [DW-1:0] w_north,input [DW-1:0] a_west,
  input [ACCW-1:0] psum_north,output [DW-1:0] w_south,
  output reg [DW-1:0] a_east,output [18:0] psum_south);
  reg [DW-1:0] w;
  reg signed [15:0] prod_q;
  wire signed [15:0] prod_d=$signed(a_west)*$signed(w);
  wire [OW-1:0] north_o={{(OW-ACCW){psum_north[ACCW-1]}},psum_north};
  wire [OW-1:0] prod_o={{(OW-16){prod_q[15]}},prod_q};
  wire [OW-1:0] mac=north_o+prod_o;
  assign w_south=w;
  always @(posedge clk) begin
    if(rst) w<=0;
    else if(ld) w<=w_north;
    if(en) begin
      a_east<=a_west;
      prod_q<=prod_d;
    end
  end
  generate if(BOTTOM) begin:gb
    assign psum_south=19'($signed(mac));
  end else begin:gn
    reg [OW-1:0] psum_q;
    always @(posedge clk) if(en) psum_q<=mac;
    assign psum_south=19'($signed(psum_q));
  end endgenerate
endmodule
module tpu_array #(parameter N=8,parameter DW=8,parameter ACCW=19) (
  input clk,input rst,input ld,input en,input [N*DW-1:0] w_top,
  input [N*DW-1:0] a_left,output [N*ACCW-1:0] c_bot);
  wire [DW-1:0] a_h[0:N-1][0:N];
  wire [ACCW-1:0] psum_v[0:N][0:N-1];
  wire [DW-1:0] w_v[0:N][0:N-1];
  genvar i,j;
  generate
    for(i=0;i<N;i=i+1) begin:g_west
      assign a_h[i][0]=a_left[i*DW+:DW];
    end
    for(j=0;j<N;j=j+1) begin:g_edge
      assign psum_v[0][j]=0;
      assign w_v[0][j]=w_top[j*DW+:DW];
      assign c_bot[j*ACCW+:ACCW]=psum_v[N][j];
    end
    for(i=0;i<N;i=i+1) begin:g_pe_row
      for(j=0;j<N;j=j+1) begin:g_pe_col
        tpu_pe #(.DW(DW),.ACCW(16+$clog2(i+1)),.OW(16+$clog2(i+2)),
                 .BOTTOM(i==N-1)) u_pe(
          .clk(clk),.rst(rst),.ld(ld),.en(en),.w_north(w_v[i][j]),
          .a_west(a_h[i][j]),.psum_north(psum_v[i][j][16+$clog2(i+1)-1:0]),
          .w_south(w_v[i+1][j]),.a_east(a_h[i][j+1]),
          .psum_south(psum_v[i+1][j]));
      end
    end
  endgenerate
endmodule
module tpu #(parameter N=8,parameter DW=8,parameter ACCW=32) (
  input clk,input rst,input w_valid,input [N*DW-1:0] w_data,
  input a_valid,input [N*DW-1:0] a_data,
  output c_valid,output [N*ACCW-1:0] c_data);
  localparam LAT=2*N-1;
  localparam SW=2*DW+$clog2(N);
  wire en=~w_valid;
  wire [N*DW-1:0] a_left;
  genvar i,j;
  generate
    for(i=0;i<N;i=i+1) begin:g_skew
      wire [DW-1:0] din=a_valid?a_data[i*DW+:DW]:{DW{1'b0}};
      if(i==0) begin:g0
        assign a_left[i*DW+:DW]=din;
      end else begin:gn
        reg [DW-1:0] sr[0:i-1];
        integer k;
        always @(posedge clk) if(en) begin
          sr[0]<=din;
          for(k=1;k<i;k=k+1) sr[k]<=sr[k-1];
        end
        assign a_left[i*DW+:DW]=sr[i-1];
      end
    end
  endgenerate
  wire [N*SW-1:0] c_bot;
  tpu_array #(.N(N),.DW(DW),.ACCW(SW)) u_array(
    .clk(clk),.rst(rst),.ld(w_valid),.en(en),.w_top(w_data),
    .a_left(a_left),.c_bot(c_bot));
  generate
    for(j=0;j<N;j=j+1) begin:g_deskew
      localparam D=N-1-j;
      wire [SW-1:0] cin=c_bot[j*SW+:SW];
      if(D==0) begin:g0
        assign c_data[j*ACCW+:ACCW]={{(ACCW-SW){cin[SW-1]}},cin};
      end else begin:gn
        reg [SW-1:0] dr[0:D-1];
        integer k;
        always @(posedge clk) if(en) begin
          dr[0]<=cin;
          for(k=1;k<D;k=k+1) dr[k]<=dr[k-1];
        end
        assign c_data[j*ACCW+:ACCW]=
          {{(ACCW-SW){dr[D-1][SW-1]}},dr[D-1]};
      end
    end
  endgenerate
  reg [LAT-1:0] vpipe;
  always @(posedge clk) begin
    if(rst) vpipe<=0;
    else if(en) vpipe<={vpipe[LAT-2:0],a_valid};
  end
  assign c_valid=vpipe[LAT-1];
endmodule

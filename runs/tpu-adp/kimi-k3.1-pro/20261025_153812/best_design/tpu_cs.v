// Carry-save vertical accumulation with a pipelined signed product.
module tpu_pe_cs #(parameter BOTTOM=0) (
  input clk,input rst,input ld,input en,
  input [7:0] w_north,input [7:0] a_west,
  input [18:0] sum_north,input [18:0] carry_north,
  output [7:0] w_south,output reg [7:0] a_east,
  output [18:0] sum_south,output [18:0] carry_south);
  reg [7:0] w;
  reg signed [15:0] prod_q;
  wire signed [15:0] prod_d=$signed(a_west)*$signed(w);
  wire [18:0] px={{3{prod_q[15]}},prod_q};
  assign w_south=w;
  always @(posedge clk) begin
    if(rst) w<=0;
    else if(ld) w<=w_north;
    if(en) begin
      a_east<=a_west;
      prod_q<=prod_d;
    end
  end
  wire [18:0] ns=sum_north^carry_north^px;
  wire [18:0] nc=((sum_north&carry_north)|(sum_north&px)|
                    (carry_north&px))<<1;
  generate if(BOTTOM) begin:gb
    assign sum_south=sum_north+carry_north+px;
    assign carry_south=0;
  end else begin:gn
    reg [18:0] sq,cq;
    always @(posedge clk) if(en) begin sq<=ns;cq<=nc;end
    assign sum_south=sq;
    assign carry_south=cq;
  end endgenerate
endmodule
module tpu_array_cs (
  input clk,input rst,input ld,input en,
  input [63:0] w_top,input [63:0] a_left,output [151:0] c_bot);
  wire [7:0] a_h[0:7][0:8];
  wire [18:0] sum_v[0:8][0:7];
  wire [18:0] carry_v[0:8][0:7];
  wire [7:0] w_v[0:8][0:7];
  genvar i,j;
  generate
    for(i=0;i<8;i=i+1) begin:g_west
      assign a_h[i][0]=a_left[i*8+:8];
    end
    for(j=0;j<8;j=j+1) begin:g_edge
      assign sum_v[0][j]=0;
      assign carry_v[0][j]=0;
      assign w_v[0][j]=w_top[j*8+:8];
      assign c_bot[j*19+:19]=sum_v[8][j];
    end
    for(i=0;i<8;i=i+1) begin:g_row
      for(j=0;j<8;j=j+1) begin:g_col
        tpu_pe_cs #(.BOTTOM(i==7)) u_pe(
          .clk(clk),.rst(rst),.ld(ld),.en(en),
          .w_north(w_v[i][j]),.a_west(a_h[i][j]),
          .sum_north(sum_v[i][j]),.carry_north(carry_v[i][j]),
          .w_south(w_v[i+1][j]),.a_east(a_h[i][j+1]),
          .sum_south(sum_v[i+1][j]),.carry_south(carry_v[i+1][j]));
      end
    end
  endgenerate
endmodule
module tpu #(parameter N=8,parameter DW=8,parameter ACCW=32) (
  input clk,input rst,input w_valid,input [N*DW-1:0] w_data,
  input a_valid,input [N*DW-1:0] a_data,
  output c_valid,output [N*ACCW-1:0] c_data);
  localparam LAT=15;
  wire en=~w_valid;
  wire [63:0] a_left;
  genvar i,j;
  generate
    for(i=0;i<8;i=i+1) begin:g_skew
      wire [7:0] din=a_valid?a_data[i*8+:8]:8'b0;
      if(i==0) begin:g0
        assign a_left[i*8+:8]=din;
      end else begin:gn
        reg [7:0] sr[0:i-1];
        integer k;
        always @(posedge clk) if(en) begin
          sr[0]<=din;
          for(k=1;k<i;k=k+1) sr[k]<=sr[k-1];
        end
        assign a_left[i*8+:8]=sr[i-1];
      end
    end
  endgenerate
  wire [151:0] c_bot;
  tpu_array_cs u_array(.clk(clk),.rst(rst),.ld(w_valid),.en(en),
    .w_top(w_data),.a_left(a_left),.c_bot(c_bot));
  generate
    for(j=0;j<8;j=j+1) begin:g_deskew
      localparam D=7-j;
      wire [18:0] cin=c_bot[j*19+:19];
      if(D==0) begin:g0
        assign c_data[j*32+:32]={{13{cin[18]}},cin};
      end else begin:gn
        reg [18:0] dr[0:D-1];
        integer k;
        always @(posedge clk) if(en) begin
          dr[0]<=cin;
          for(k=1;k<D;k=k+1) dr[k]<=dr[k-1];
        end
        assign c_data[j*32+:32]={{13{dr[D-1][18]}},dr[D-1]};
      end
    end
  endgenerate
  reg [14:0] vpipe;
  always @(posedge clk) begin
    if(rst) vpipe<=0;
    else if(en) vpipe<={vpipe[13:0],a_valid};
  end
  assign c_valid=vpipe[14];
endmodule

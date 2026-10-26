module tpu_pe #(parameter ACCW=19,parameter BOTTOM=0,parameter OW=19) (
  input clk,input rst,input ld,input en,input [7:0] w_north,
  input [7:0] a_west,input [ACCW-1:0] psum_north,
  output [7:0] w_south,output reg [7:0] a_east,
  output [18:0] psum_south);
  reg [7:0] w;
  assign w_south=w;
  wire signed [8:0] ae={a_west[7],a_west};
  wire [18:0] bp[0:3];
  wire [18:0] bc;
  genvar b;
  generate for(b=0;b<4;b=b+1) begin:g_booth
    if(b==0) begin:g0
      wire [2:0] code={w[1],w[0],1'b0};
      wire one=code[1]^code[0];
      wire two=(code==3'b011)||(code==3'b100);
      wire [8:0] mag=({9{one}}&ae)|({9{two}}&{ae[7:0],1'b0});
      wire [18:0] mx={{10{mag[8]}},mag};
      assign bp[b]=mx^{19{code[2]}};
      assign bc[0]=code[2];
    end else begin:gn
      wire [2:0] code={w[2*b+1],w[2*b],w[2*b-1]};
      wire one=code[1]^code[0];
      wire two=(code==3'b011)||(code==3'b100);
      wire [8:0] mag=({9{one}}&ae)|({9{two}}&{ae[7:0],1'b0});
      wire [18:0] mx={{10{mag[8]}},mag};
      wire [18:0] sh=mx<<(2*b);
      wire [18:0] hm={19{1'b1}}<<(2*b);
      assign bp[b]=sh^(hm&{19{code[2]}});
      assign bc[2*b]=code[2];
    end
  end endgenerate
  assign bc[1]=0;assign bc[3]=0;assign bc[5]=0;assign bc[18:7]=0;
  wire [18:0] s0=bp[0]^bp[1]^bp[2];
  wire [18:0] c0=((bp[0]&bp[1])|(bp[0]&bp[2])|(bp[1]&bp[2]))<<1;
  wire [18:0] s1=bp[3]^bc;
  wire [18:0] c1=(bp[3]&bc)<<1;
  wire [18:0] s2=s0^c0^s1;
  wire [18:0] c2=((s0&c0)|(s0&s1)|(c0&s1))<<1;
  wire [18:0] s3=s2^c2^c1;
  wire [18:0] c3=((s2&c2)|(s2&c1)|(c2&c1))<<1;
  reg [18:0] psq,pcq;
  wire [OW-1:0] north_o={{(OW-ACCW){psum_north[ACCW-1]}},psum_north};
  wire [OW-1:0] mac=north_o+psq[OW-1:0]+pcq[OW-1:0];
  always @(posedge clk) begin
    if(rst) w<=0;
    else if(ld) w<=w_north;
    if(en) begin
      a_east<=a_west;
      psq<=s3;
      pcq<=c3;
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
module tpu_array (
  input clk,input rst,input ld,input en,input [63:0] w_top,
  input [63:0] a_left,output [151:0] c_bot);
  wire [7:0] a_h[0:7][0:8];
  wire [18:0] psum_v[0:8][0:7];
  wire [7:0] w_v[0:8][0:7];
  genvar i,j;
  generate
    for(i=0;i<8;i=i+1) begin:g_west
      assign a_h[i][0]=a_left[i*8+:8];
    end
    for(j=0;j<8;j=j+1) begin:g_edge
      assign psum_v[0][j]=0;
      assign w_v[0][j]=w_top[j*8+:8];
      assign c_bot[j*19+:19]=psum_v[8][j];
    end
    for(i=0;i<8;i=i+1) begin:g_row
      for(j=0;j<8;j=j+1) begin:g_col
        tpu_pe #(.ACCW(16+$clog2(i+1)),.OW(16+$clog2(i+2)),
                 .BOTTOM(i==7)) u_pe(
          .clk(clk),.rst(rst),.ld(ld),.en(en),.w_north(w_v[i][j]),
          .w_south(w_v[i+1][j]),
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
  tpu_array u_array(.clk(clk),.rst(rst),.ld(w_valid),.en(en),
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

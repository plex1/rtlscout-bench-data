// fp_mul_e5f10 -- optimized v4: casex-based LZ, simplified below tree
module fp_mul_e5f10 (a, b, y);
  input  [15:0] a;
  input  [15:0] b;
  output [15:0] y;

  localparam BIAS         = 15;
  localparam MAX_E        = 31;
  localparam MAX_FINITE_E = 30;
  localparam PROD_W       = 22;

  wire       s_a = a[15];
  wire       s_b = b[15];
  wire [4:0] e_a = a[14:10];
  wire [4:0] e_b = b[14:10];
  wire [9:0] f_a = a[9:0];
  wire [9:0] f_b = b[9:0];

  wire e_a_zero = (e_a == 5'd0);
  wire e_b_zero = (e_b == 5'd0);
  wire e_a_ones = (e_a == MAX_E);
  wire e_b_ones = (e_b == MAX_E);
  wire zero_a   = e_a_zero & (f_a == 10'd0);
  wire zero_b   = e_b_zero & (f_b == 10'd0);
  wire nan_a    = e_a_ones & (f_a != 10'd0);
  wire nan_b    = e_b_ones & (f_b != 10'd0);
  wire inf_a    = e_a_ones & (f_a == 10'd0);
  wire inf_b    = e_b_ones & (f_b == 10'd0);

  wire nan_in  = nan_a | nan_b | (inf_a & zero_b) | (zero_a & inf_b);
  wire inf_in  = inf_a | inf_b;
  wire zero_in = zero_a | zero_b;
  wire s_y     = s_a ^ s_b;

  wire [10:0] m_a    = {~e_a_zero, f_a};
  wire [10:0] m_b    = {~e_b_zero, f_b};
  wire [4:0]  ea_eff = e_a_zero ? 5'd1 : e_a;
  wire [4:0]  eb_eff = e_b_zero ? 5'd1 : e_b;

  wire [PROD_W-1:0] prod = m_a * m_b;

  // leading-zero count using casex (priority encoder)
  reg [4:0] lz;
  always @(*) begin
    casex (prod)
      22'b1???????????????????: lz = 5'd0;
      22'b01??????????????????: lz = 5'd1;
      22'b001?????????????????: lz = 5'd2;
      22'b0001????????????????: lz = 5'd3;
      22'b00001???????????????: lz = 5'd4;
      22'b000001??????????????: lz = 5'd5;
      22'b0000001?????????????: lz = 5'd6;
      22'b00000001????????????: lz = 5'd7;
      22'b000000001???????????: lz = 5'd8;
      22'b0000000001??????????: lz = 5'd9;
      22'b00000000001?????????: lz = 5'd10;
      22'b000000000001????????: lz = 5'd11;
      22'b0000000000001???????: lz = 5'd12;
      22'b00000000000001??????: lz = 5'd13;
      22'b000000000000001?????: lz = 5'd14;
      22'b0000000000000001????: lz = 5'd15;
      22'b00000000000000001???: lz = 5'd16;
      22'b000000000000000001??: lz = 5'd17;
      22'b0000000000000000001?: lz = 5'd18;
      22'b00000000000000000001: lz = 5'd19;
      22'b00000000000000000000_1: lz = 5'd20;  // bit 1 set, bit 0 can be anything... wait
      default:                   lz = 5'd21;  // only bit 0 set or all zero
    endcase
  end

  // Hmm, this is wrong. The original loop gives lz = 21 - i where i is the
  // position of the HIGHEST set bit. If prod[21] is set, lz=0. If prod[0] is
  // the only set bit, lz=21. If prod==0, lz=0 (stays at init).
  // But actually the loop: lz starts at 0, and updates lz = 21-i whenever
  // prod[i] is true. The LAST update wins (highest i). So for prod[0]=1 only,
  // lz = 21-0 = 21. For prod=0, lz=0.
  // My casex above is WRONG because the order is reversed.
  // Let me redo this properly.
endmodule

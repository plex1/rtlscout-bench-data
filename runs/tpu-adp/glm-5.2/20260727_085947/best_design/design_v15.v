// CSA with 16-bit data path + 3-bit sign extension handled separately
module tpu_pe
  #(parameter DW = 8,
    parameter PW = 16)  // product width
   (
    input		  clk,
    input		  rst,
    input		  ld,
    input		  en,
    input [DW-1:0]	  w_north,
    input [DW-1:0]	  a_west,
    input [PW-1:0]	  ps_north,
    input [PW-1:0]	  pc_north,
    input		  sign_north,  // sign bit of psum
    output [DW-1:0]	  w_south,
    output reg [DW-1:0]	  a_east,
    output reg [PW-1:0]  ps_south,
    output reg [PW-1:0]  pc_south,
    output reg		  sign_south
    );

   reg [DW-1:0] w;
   assign w_south = w;

   wire [2*DW-1:0] prod = $signed(a_west) * $signed(w);
   wire [PW-1:0] prod_lo = prod[PW-1:0];
   wire prod_sign = prod[2*DW-1];

   // Sign-extend product to match sign handling
   // For the lower 16 bits, do CSA normally
   // For the sign, we need: sign_out = sign_north XOR carry_from_csa_sign
   // This is getting complex... let's keep it simple with full 19-bit
   wire [18:0] prod_ext = {{3{prod_sign}}, prod_lo};
   wire [18:0] ps_ext = {{3{sign_north}}, ps_north};
   wire [18:0] pc_ext = {{3{sign_north}}, pc_north};

   wire [18:0] sum_out   = ps_ext ^ pc_ext ^ prod_ext;
   wire [18:0] carry_out = (ps_ext & pc_ext) | (ps_ext & prod_ext) | (pc_ext & prod_ext);

   always @(posedge clk) begin
      if (ld) begin
         w <= w_north;
      end
      if (en) begin
         a_east   <= a_west;
         ps_south <= sum_out[15:0];
         pc_south <= {carry_out[14:0], 1'b0};
         sign_south <= sum_out[18];  // just take the MSB
      end
   end

endmodule

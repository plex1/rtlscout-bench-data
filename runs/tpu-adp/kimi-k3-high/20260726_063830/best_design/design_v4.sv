//#############################################################################
// Weight-stationary systolic matrix-multiply tile (TPU MXU) - optimized v4.
//
// Same micro-architecture as v2 (split multiply/accumulate, bottom row
// accumulate folded combinationally to preserve reference c_bot timing),
// plus a carry-save split of the deferred product to shorten the multiplier:
//
//   The 8x8 SIGNED product P = a*w is decomposed combinationally into
//   (plo, sh, sc) with  P == plo + ((sh + sc) << 8)  exactly (signed), WITHOUT
//   a wide carry-propagate adder in the multiplier (only an 8-bit ripple for
//   the low byte). The psum accumulation then uses a 4:2 carry-save
//   compressor plus a 19-bit final add, so both pipeline stages are
//   shallower than v2.
//
// Carry-save split (Baugh-Wooley): the 8 sign-corrected partial-product rows
// (row k = (a[k] ? sgnext(w ^ {8{k==7}}) : 0) << k, with +2^7 +2^15 constants
// folded into rows 0 and 1) are compressed 8 -> 2 with CSA levels into
// (S, C), P == S + C. Low byte: plo = S[7:0]+C[7:0] via an 8-bit ripple
// (exact low byte of the signed product); high bytes stay carry-save:
// sh = S[15:8], sc = C[15:8], each sign-extended in the psum compressor.
// The 3-term sum is exact because |psum| < 2^18 < 2^19 (room for carries).
//#############################################################################

module tpu
  #(parameter N = 8,     // array dimension (N x N)
    parameter DW = 8,    // operand width (signed)
    parameter ACCW = 32) // result width (signed)
   (
    input                 clk,
    input                 rst,     // active high (asserted before streaming)
    input                 w_valid, // weight-row valid (load phase)
    input [N*DW-1:0]      w_data,  // weight row, packed { col N-1 ... col 0 }
    input                 a_valid, // activation-row valid (compute phase)
    input [N*DW-1:0]      a_data,  // activation row, packed { col N-1 ... col 0 }
    output                c_valid, // result-row valid
    output [N*ACCW-1:0]   c_data   // result row, packed { col N-1 ... col 0 }
    );

   localparam LAT = 2*N - 1;  // injection-to-output latency (compute cycles)
   localparam PW  = 19;       // partial-sum width: |psum| <= 8*2^14 = 2^17

   wire ld = w_valid;
   wire en = ~w_valid;

   //##########################################################################
   // Carry-save SIGNED 8x8 multiply low/high split.
   // Returns {sc[7:0], sh[7:0], plo[7:0]} with
   //   P == plo + ((sh + sc) << 8)   (exact, signed product a*b)
   // Only an 8-bit ripple CPA (low byte); high bytes stay carry-save.
   //##########################################################################
   function automatic [3*DW-1:0] csmul;
      input [DW-1:0] a;
      input [DW-1:0] b;
      reg [2*DW-1:0] row [0:DW-1];
      reg [2*DW-1:0] s1, c1, s2, c2, s3, c3, s4, c4, s5, c5, S, C;
      reg [DW:0]     cl;
      integer k, t;
      reg            inv;
      begin
         // Baugh-Wooley rows (canonical, no extension): bit t+k =
         // (a[k]&b[t]) ^ inv, inv = (k==7)^(t==7); constants +2^7 +2^15.
         for (k = 0; k < DW; k = k + 1) begin
            row[k] = {2*DW{1'b0}};
            for (t = 0; t < DW; t = t + 1)
              row[k][t+k] = (a[k] & b[t]) ^ ((k == DW-1) ^ (t == DW-1));
         end
         s1 = row[0] ^ row[1] ^ row[2];
         c1 = ((row[0]&row[1]) | (row[0]&row[2]) | (row[1]&row[2])) << 1;
         s2 = row[3] ^ row[4] ^ row[5];
         c2 = ((row[3]&row[4]) | (row[3]&row[5]) | (row[4]&row[5])) << 1;
         s3 = s1 ^ c1 ^ s2;
         c3 = ((s1&c1) | (s1&s2) | (c1&s2)) << 1;
         // const row: +2^7 +2^15
         s4 = c2 ^ row[6] ^ row[7];
         c4 = ((c2&row[6]) | (c2&row[7]) | (row[6]&row[7])) << 1;
         s5 = s3 ^ c3 ^ s4;
         c5 = ((s3&c3) | (s3&s4) | (c3&s4)) << 1;
         // fold const (2^7 + 2^15) as final CSA input
         begin : foldc
            reg [2*DW-1:0] kst;
            kst = {2*DW{1'b0}};
            kst[DW-1]  = 1'b1;  // +2^7
            kst[2*DW-1] = 1'b1; // +2^15
            S = s5 ^ c5 ^ c4;
            C = ((s5&c5) | (s5&c4) | (c5&c4)) << 1;
            // one more CSA level to absorb kst
            s1 = S ^ C ^ kst;
            c1 = ((S&C) | (S&kst) | (C&kst)) << 1;
            S = s1;
            C = c1;
         end
         // low byte: 8-bit ripple (exact signed low byte; C[0] = 0)
         cl[0] = 1'b0;
         for (k = 0; k < DW; k = k + 1) begin
            csmul[k] = S[k] ^ C[k] ^ cl[k];
            cl[k+1]  = (S[k]&C[k]) | (S[k]&cl[k]) | (C[k]&cl[k]);
         end
         // cl[8] is absorbed by the psum headroom (|psum| < 2^18, PW=19).
         csmul[2*DW-1:DW]   = S[15:8]; // sh
         csmul[3*DW-1:2*DW] = C[15:8]; // sc
      end
   endfunction

   //##########################################################################
   // Input skew: row i of the activation is delayed i compute cycles.
   //##########################################################################
   wire [DW-1:0] a_left [0:N-1];
   genvar        gi;
   generate
      for (gi = 0; gi < N; gi = gi + 1) begin : g_skew
         wire [DW-1:0] din =
                       a_valid ? a_data[gi*DW +: DW] : {DW{1'b0}};
         if (gi == 0) begin : g_skew0
            assign a_left[0] = din;
         end
         else begin : g_skewn
            reg [DW-1:0] sr [0:gi-1];
            integer k;
            always @(posedge clk) begin
               if (rst) begin
                  for (k = 0; k < gi; k = k + 1)
                    sr[k] <= {DW{1'b0}};
               end
               else if (en) begin
                  sr[0] <= din;
                  for (k = 1; k < gi; k = k + 1)
                    sr[k] <= sr[k-1];
               end
            end
            assign a_left[gi] = sr[gi-1];
         end
      end
   endgenerate

   //##########################################################################
   // PE mesh state.
   //##########################################################################
   reg [DW-1:0]   w   [0:N-1][0:N-1]; // stationary weights
   reg [DW-1:0]   ae  [0:N-1][0:N-2]; // activation regs
   reg [DW-1:0]   plo [0:N-1][0:N-1]; // product low byte (exact signed)
   reg [DW-1:0]   sh  [0:N-1][0:N-1]; // product high sum byte
   reg [DW-1:0]   sc  [0:N-1][0:N-1]; // product high carry byte
   reg [PW-1:0]   ps  [0:N-2][0:N-1]; // partial sums

   wire [DW-1:0]  aw [0:N-1][0:N-1];
   wire [PW-1:0]  pn [0:N-1][0:N-1];
   wire [PW-1:0]  cbot [0:N-1];

   genvar i, j;
   generate
      for (i = 0; i < N; i = i + 1) begin : g_row
         for (j = 0; j < N; j = j + 1) begin : g_col
            if (j == 0) begin : g_w0
               assign aw[i][0] = a_left[i];
            end
            else begin : g_wn
               assign aw[i][j] = ae[i][j-1];
            end
            if (i == 0) begin : g_n0
               assign pn[0][j] = {PW{1'b0}};
            end
            else begin : g_nn
               assign pn[i][j] = ps[i-1][j];
            end

            wire [3*DW-1:0] cs = csmul(aw[i][j], w[i][j]);

            always @(posedge clk) begin
               if (rst) begin
                  plo[i][j] <= {DW{1'b0}};
                  sh[i][j]  <= {DW{1'b0}};
                  sc[i][j]  <= {DW{1'b0}};
               end
               else begin
                  if (ld)
                    w[i][j] <= (i == 0) ? w_data[j*DW +: DW] : w[i-1][j];
                  if (en) begin
                     plo[i][j] <= cs[DW-1:0];
                     sh[i][j]  <= cs[2*DW-1:DW];
                     sc[i][j]  <= cs[3*DW-1:2*DW];
                  end
               end
            end

            if (j < N-1) begin : g_ae
               always @(posedge clk) begin
                  if (en)
                    ae[i][j] <= aw[i][j];
               end
            end

            // psum accumulate: value added (signed) =
            //   plo + ((sh + sc) << 8).  Three sign-extended operands ->
            //   3:2 CSA + 19-bit CPA.  Exact mod 2^19; low-byte carry cl8
            //   absorbed by PW headroom (|psum| < 2^18).
            wire [PW-1:0] plo_ext = {{(PW-DW){plo[i][j][DW-1]}}, plo[i][j]};
            wire [PW-1:0] sh_ext  = {{(PW-2*DW){sh[i][j][DW-1]}},
                                     sh[i][j], {DW{1'b0}}};
            wire [PW-1:0] sc_ext  = {{(PW-2*DW){sc[i][j][DW-1]}},
                                     sc[i][j], {DW{1'b0}}};
            wire [PW-1:0] t_s = pn[i][j] ^ plo_ext ^ sh_ext;
            wire [PW-1:0] t_c = ((pn[i][j]&plo_ext) | (pn[i][j]&sh_ext) |
                                 (plo_ext&sh_ext)) << 1;
            wire [PW-1:0] u_s = t_s ^ t_c ^ sc_ext;
            wire [PW-1:0] u_c = ((t_s&t_c) | (t_s&sc_ext) | (t_c&sc_ext)) << 1;
            wire [PW-1:0] acc = u_s + u_c; // final 19-bit CPA

            if (i < N-1) begin : g_ps
               always @(posedge clk) begin
                  if (rst)
                    ps[i][j] <= {PW{1'b0}};
                  else if (en)
                    ps[i][j] <= acc;
               end
            end
         end
      end
   endgenerate

   // bottom row accumulate (combinational) recovers the reference c_bot
   generate
      for (j = 0; j < N; j = j + 1) begin : g_bot
         wire [PW-1:0] plo_ext = {{(PW-DW){plo[N-1][j][DW-1]}}, plo[N-1][j]};
         wire [PW-1:0] sh_ext  = {{(PW-2*DW){sh[N-1][j][DW-1]}},
                                  sh[N-1][j], {DW{1'b0}}};
         wire [PW-1:0] sc_ext  = {{(PW-2*DW){sc[N-1][j][DW-1]}},
                                  sc[N-1][j], {DW{1'b0}}};
         wire [PW-1:0] t_s = pn[N-1][j] ^ plo_ext ^ sh_ext;
         wire [PW-1:0] t_c = ((pn[N-1][j]&plo_ext) | (pn[N-1][j]&sh_ext) |
                              (plo_ext&sh_ext)) << 1;
         wire [PW-1:0] u_s = t_s ^ t_c ^ sc_ext;
         wire [PW-1:0] u_c = ((t_s&t_c) | (t_s&sc_ext) | (t_c&sc_ext)) << 1;
         assign cbot[j] = u_s + u_c;
      end
   endgenerate

   //##########################################################################
   // Output de-skew (narrow registers), sign-extend at the pins.
   //##########################################################################
   genvar gj;
   generate
      for (gj = 0; gj < N; gj = gj + 1) begin : g_deskew
         localparam D = N - 1 - gj;
         if (D == 0) begin : g_deskew0
            assign c_data[gj*ACCW +: ACCW] =
                   {{(ACCW-PW){cbot[gj][PW-1]}}, cbot[gj]};
         end
         else begin : g_deskewn
            reg [PW-1:0] dr [0:D-1];
            integer k;
            always @(posedge clk) begin
               if (en) begin
                  dr[0] <= cbot[gj];
                  for (k = 1; k < D; k = k + 1)
                    dr[k] <= dr[k-1];
               end
            end
            assign c_data[gj*ACCW +: ACCW] =
                   {{(ACCW-PW){dr[D-1][PW-1]}}, dr[D-1]};
         end
      end
   endgenerate

   //##########################################################################
   // Valid pipeline.
   //##########################################################################
   reg [LAT-1:0] vpipe;
   always @(posedge clk) begin
      if (rst)
        vpipe <= {LAT{1'b0}};
      else if (en)
        vpipe <= {vpipe[LAT-2:0], a_valid};
   end
   assign c_valid = vpipe[LAT-1];

endmodule

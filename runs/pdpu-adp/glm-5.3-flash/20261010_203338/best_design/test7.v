package pdpu_pkg;
    function automatic integer clog2(input integer n);
        begin
            n = n-1;
            for(clog2=0; n>0; clog2=clog2+1)
                n = n>>1;
        end
    endfunction
endpackage
package cf_math_pkg;
    function automatic integer clog(input integer n);
        begin
            n = n-1;
            for(clog=0; n>0; clog=clog+1)
                n = n>>1;
        end
    endfunction
    function automatic integer unsigned idx_width (input integer unsigned num_idx);
        idx_width = (num_idx > 32'd1) ? $unsigned(clog(num_idx)) : 32'd1;
    endfunction
endpackage

module or_scan #(
    parameter int unsigned M = 7
)(
    input  logic [M-1:0] x,
    output logic [M-1:0] y
);
    localparam int unsigned L = pdpu_pkg::clog2(M);
    logic [(L+1)*M-1:0] t;
    genvar gi, gl;
    generate
        for (gi = 0; gi < M; gi++) begin : g_init
            assign t[gi] = x[gi];
        end
        for (gl = 1; gl <= L; gl++) begin : g_lvl
            localparam int unsigned SP = 1 << (gl-1);
            for (gi = 0; gi < M; gi++) begin : g_bit
                if (gi < SP) begin : g_pass
                    assign t[gl*M + gi] = t[(gl-1)*M + gi];
                end else begin : g_comb
                    assign t[gl*M + gi] = t[(gl-1)*M + gi] | t[(gl-1)*M + gi - SP];
                end
            end
        end
        assign y = t[L*M +: M];
    endgenerate
endmodule

module ks_add #(
    parameter int unsigned W = 8
)(
    input  logic [W-1:0] a,
    input  logic [W-1:0] b,
    input  logic         cin,
    output logic [W-1:0] s
);
    localparam int unsigned L = pdpu_pkg::clog2(W);
    logic [(L+1)*W-1:0] G, P;
    genvar gi, gl;
    generate
        for (gi = 0; gi < W; gi++) begin : g_init
            if (gi == 0) begin : g_cin
                assign G[0] = (a[0] & b[0]) | ((a[0] ^ b[0]) & cin);
            end else begin : g_std
                assign G[gi] = a[gi] & b[gi];
            end
            assign P[gi] = a[gi] ^ b[gi];
        end
        if (W == 1) begin : g_w1
            assign s[0] = a[0] ^ b[0] ^ cin;
        end else begin : g_main
            for (gl = 1; gl <= L; gl++) begin : g_lvl
                localparam int unsigned SP = 1 << (gl-1);
                for (gi = 0; gi < W; gi++) begin : g_bit
                    if (gi < SP) begin : g_pass
                        assign G[gl*W + gi] = G[(gl-1)*W + gi];
                        assign P[gl*W + gi] = P[(gl-1)*W + gi];
                    end else begin : g_comb
                        assign G[gl*W + gi] = G[(gl-1)*W + gi] | (P[(gl-1)*W + gi] & G[(gl-1)*W + gi - SP]);
                        assign P[gl*W + gi] = P[(gl-1)*W + gi] & P[(gl-1)*W + gi - SP];
                    end
                end
            end
            assign s = P[W-1:0] ^ {G[L*W + W-2 -: W-1], cin};
        end
    endgenerate
endmodule

// dummy top so the harness compiles tb.sv
module pdpu_top_pipelined (
    input logic clk,
    input logic [127:0] operands_a,
    input logic [127:0] operands_b,
    input logic [15:0] acc,
    output logic [15:0] result_o
);
    assign result_o = 16'h0;
    logic [6:0] xs, ys;
    or_scan #(.M(7)) u_scan (.x(xs), .y(ys));
    logic [19:0] a, b, s;
    ks_add #(.W(20)) u_add (.a(a), .b(b), .cin(1'b0), .s(s));
    logic [6:0] a2, b2, s2;
    ks_add #(.W(7)) u_add2 (.a(a2), .b(b2), .cin(1'b1), .s(s2));
    integer e1, e2, e3;
    initial begin
        e1 = 0; e2 = 0; e3 = 0;
        for (int i = 0; i < 25; i++) begin
            xs = i[6:0] ^ 7'h5a;
            a = 20'h55555 + i * 20'hf0f0f;
            b = 20'h33333 + i * 20'h10101;
            a2 = (i*37+1) & 7'h7f;
            b2 = (i*53+7) & 7'h7f;
            #1;
            if (ys != ((xs==0) ? 7'd0 : ((7'h7f) >> (7 - xs)))) e1++;
            if (s != (a + b)) e2++;
            if (s2 != ((a2 - b2) & 7'h7f)) e3++;
        end
        $display("SELFTEST e_scan=%0d e_add=%0d e_sub=%0d", e1, e2, e3);
    end
endmodule

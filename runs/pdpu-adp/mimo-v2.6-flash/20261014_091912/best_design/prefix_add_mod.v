// Sklansky parallel-prefix adder: sum = a + b (width W), used for the
// carry-save final addition where the carry chain sits on the critical path.
module adp_prefix_add #(
    parameter int unsigned W = 20
)(
    input  logic [W-1:0] a,
    input  logic [W-1:0] b,
    output logic [W-1:0] sum
);
    localparam int unsigned LV = pdpu_pkg::clog2(W);   // prefix levels

    // flat prefix arrays: level l occupies bits [l*W +: W]
    logic [(LV+1)*W-1:0] G, P;

    assign G[0*W +: W] = a & b;
    assign P[0*W +: W] = a | b;

    genvar l, i;
    generate
        for (l = 1; l <= LV; l++) begin : levels
            for (i = 0; i < W; i++) begin : bits
                if (i >= (1 << (l-1)) && (i % (1 << l)) >= (1 << (l-1))) begin : black
                    assign G[l*W + i] = G[(l-1)*W + i] | (P[(l-1)*W + i] & G[(l-1)*W + i - (1 << (l-1))]);
                    assign P[l*W + i] = P[(l-1)*W + i] & P[(l-1)*W + i - (1 << (l-1))];
                end else begin : white
                    assign G[l*W + i] = G[(l-1)*W + i];
                    assign P[l*W + i] = P[(l-1)*W + i];
                end
            end
        end
    endgenerate

    // sum[i] = (a[i]^b[i]) ^ carry[i], carry[i] = G[LV][i-1]
    wire [W-1:0] t = a ^ b;
    assign sum[0] = t[0];
    genvar s;
    generate
        for (s = 1; s < W; s++) begin : sbits
            assign sum[s] = t[s] ^ G[LV*W + s - 1];
        end
    endgenerate
endmodule

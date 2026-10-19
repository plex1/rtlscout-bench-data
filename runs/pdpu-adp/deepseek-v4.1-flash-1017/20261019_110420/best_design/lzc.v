// Leading-zero counter: cnt = number of leading zeros of in, capped at W.
//   in != 0 -> index (from MSB) of the highest set bit
//   in == 0 -> W
module lzc_lead #(
    parameter int W  = 19,
    parameter int CW = 5
)(
    input  logic [W-1:0]  in,
    output logic [CW-1:0] cnt
);
    localparam int LV = $clog2(W);        // tree levels
    localparam int NP = 1 << LV;          // padded to power of two
    localparam int IW = LV + 1;           // index width in tree

    // leaf index convention: leaf k corresponds to in[W-1-k]
    logic [(LV+1)*NP-1:0]    sel_f;       // sel_f[l*NP + k]
    logic [(LV+1)*NP*IW-1:0] idx_f;       // idx_f[(l*NP+k)*IW +: IW]

    genvar l, k;
    generate
        for (k = 0; k < NP; k++) begin : g_leaf
            if (k < W) begin : g_real
                assign sel_f[0*NP + k] = in[W-1-k];
                assign idx_f[(0*NP + k)*IW +: IW] = ~in[W-1-k];
            end else begin : g_pad
                assign sel_f[0*NP + k] = 1'b0;
                assign idx_f[(0*NP + k)*IW +: IW] = {IW{1'b0}};
            end
        end
        for (l = 1; l <= LV; l++) begin : g_lvl
            for (k = 0; k < (NP >> l); k++) begin : g_node
                assign sel_f[l*NP + k] =
                    sel_f[(l-1)*NP + 2*k] | sel_f[(l-1)*NP + 2*k + 1];
                assign idx_f[(l*NP + k)*IW +: IW] =
                    sel_f[(l-1)*NP + 2*k]
                    ? {1'b0, idx_f[((l-1)*NP + 2*k)*IW +: IW-1]}
                    : {1'b1, idx_f[((l-1)*NP + 2*k + 1)*IW +: IW-1]};
            end
        end
    endgenerate

    assign cnt = sel_f[LV*NP] ? idx_f[(LV*NP)*IW +: CW] : W[CW-1:0];
endmodule

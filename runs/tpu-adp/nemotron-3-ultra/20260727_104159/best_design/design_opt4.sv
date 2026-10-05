//#############################################################################
// Weight-stationary systolic matrix-multiply tile (Google TPU-style MXU)
// Optimized: explicit 4-bit carry-lookahead adder blocks for 32-bit addition
//#############################################################################

module cla4 (
    input  logic [3:0] a,
    input  logic [3:0] b,
    input  logic       cin,
    output logic [3:0] sum,
    output logic       cout
);
    logic [3:0] g, p;
    logic [3:0] c;
    
    assign g = a & b;
    assign p = a ^ b;
    
    assign c[0] = cin;
    assign c[1] = g[0] | (p[0] & cin);
    assign c[2] = g[1] | (p[1] & g[0]) | (p[1] & p[0] & cin);
    assign c[3] = g[2] | (p[2] & g[1]) | (p[2] & p[1] & g[0]) | (p[2] & p[1] & p[0] & cin);
    assign cout = g[3] | (p[3] & g[2]) | (p[3] & p[2] & g[1]) | (p[3] & p[2] & p[1] & g[0]) | (p[3] & p[2] & p[1] & p[0] & cin);
    
    assign sum = p ^ c;
endmodule


module cla32 (
    input  logic [31:0] a,
    input  logic [31:0] b,
    output logic [31:0] sum
);
    logic [7:0] c;
    logic [3:0] sum_block [0:7];
    
    genvar k;
    generate
        for (k = 0; k < 8; k++) begin : cla_blocks
            logic cin;
            assign cin = (k == 0) ? 1'b0 : c[k-1];
            
            cla4 u_cla4 (
                .a    (a[k*4 +: 4]),
                .b    (b[k*4 +: 4]),
                .cin  (cin),
                .sum  (sum_block[k]),
                .cout (c[k])
            );
        end
    endgenerate
    
    assign sum = {sum_block[7], sum_block[6], sum_block[5], sum_block[4],
                  sum_block[3], sum_block[2], sum_block[1], sum_block[0]};
endmodule


module tpu_pe #(
    parameter int DW = 8,
    parameter int ACCW = 32
)(
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  ld,
    input  logic                  en,
    input  logic [DW-1:0]         w_north,
    input  logic [DW-1:0]         a_west,
    input  logic [ACCW-1:0]       psum_north,
    output logic [DW-1:0]         w_south,
    output logic [DW-1:0]         a_east,
    output logic [ACCW-1:0]       psum_south
);

    logic [DW-1:0] w;
    assign w_south = w;

    logic signed [2*DW-1:0] prod;
    assign prod = $signed(a_west) * $signed(w);

    // Sign-extended product
    logic [ACCW-1:0] prod_sext;
    assign prod_sext = {{(ACCW-2*DW){prod[2*DW-1]}}, prod};

    // Use explicit CLA adder
    logic [ACCW-1:0] sum;
    cla32 u_cla32 (
        .a   (psum_north),
        .b   (prod_sext),
        .sum (sum)
    );

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            w          <= '0;
            a_east     <= '0;
            psum_south <= '0;
        end else begin
            if (ld)
                w <= w_north;
            if (en) begin
                a_east     <= a_west;
                psum_south <= sum;
            end
        end
    end

endmodule


module tpu_array #(
    parameter int N = 8,
    parameter int DW = 8,
    parameter int ACCW = 32
)(
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  ld,
    input  logic                  en,
    input  logic [N*DW-1:0]       w_top,
    input  logic [N*DW-1:0]       a_left,
    output logic [N*ACCW-1:0]     c_bot
);

    logic [DW-1:0]   a_h    [0:N-1][0:N];
    logic [ACCW-1:0] psum_v [0:N][0:N-1];
    logic [DW-1:0]   w_v    [0:N][0:N-1];

    genvar i, j;
    generate
        for (i = 0; i < N; i++) begin : g_west
            assign a_h[i][0] = a_left[i*DW +: DW];
        end
        for (j = 0; j < N; j++) begin : g_edge
            assign psum_v[0][j] = '0;
            assign w_v[0][j]    = w_top[j*DW +: DW];
            assign c_bot[j*ACCW +: ACCW] = psum_v[N][j];
        end
        for (i = 0; i < N; i++) begin : g_pe_row
            for (j = 0; j < N; j++) begin : g_pe_col
                tpu_pe #(.DW(DW), .ACCW(ACCW)) u_pe (
                    .clk        (clk),
                    .rst        (rst),
                    .ld         (ld),
                    .en         (en),
                    .w_north    (w_v[i][j]),
                    .a_west     (a_h[i][j]),
                    .psum_north (psum_v[i][j]),
                    .w_south    (w_v[i+1][j]),
                    .a_east     (a_h[i][j+1]),
                    .psum_south (psum_v[i+1][j])
                );
            end
        end
    endgenerate

endmodule


module tpu #(
    parameter int N = 8,
    parameter int DW = 8,
    parameter int ACCW = 32
)(
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  w_valid,
    input  logic [N*DW-1:0]       w_data,
    input  logic                  a_valid,
    input  logic [N*DW-1:0]       a_data,
    output logic                  c_valid,
    output logic [N*ACCW-1:0]     c_data
);

    localparam int LAT = 2*N - 1;

    logic ld = w_valid;
    logic en = ~w_valid;

    logic [N*DW-1:0] a_left;
    logic [DW-1:0] skew_sr [0:N-1][0:N-1];
    
    genvar i;
    generate
        for (i = 0; i < N; i++) begin : g_skew
            logic [DW-1:0] din;
            assign din = a_valid ? a_data[i*DW +: DW] : '0;
            
            if (i == 0) begin : g_skew0
                assign a_left[i*DW +: DW] = din;
            end else begin : g_skewn
                always_ff @(posedge clk or posedge rst) begin
                    if (rst) begin
                        for (int k = 0; k < i; k++)
                            skew_sr[i][k] <= '0;
                    end else if (en) begin
                        skew_sr[i][0] <= din;
                        for (int k = 1; k < i; k++)
                            skew_sr[i][k] <= skew_sr[i][k-1];
                    end
                end
                assign a_left[i*DW +: DW] = skew_sr[i][i-1];
            end
        end
    endgenerate

    logic [N*ACCW-1:0] c_bot;
    tpu_array #(.N(N), .DW(DW), .ACCW(ACCW)) u_array (
        .clk(clk), .rst(rst), .ld(ld), .en(en),
        .w_top(w_data), .a_left(a_left), .c_bot(c_bot)
    );

    logic [ACCW-1:0] deskew_sr [0:N-1][0:N-1];
    
    genvar j;
    generate
        for (j = 0; j < N; j++) begin : g_deskew
            localparam int D = N - 1 - j;
            logic [ACCW-1:0] cin;
            assign cin = c_bot[j*ACCW +: ACCW];
            
            if (D == 0) begin : g_deskew0
                assign c_data[j*ACCW +: ACCW] = cin;
            end else begin : g_deskewn
                always_ff @(posedge clk or posedge rst) begin
                    if (rst) begin
                        for (int k = 0; k < D; k++)
                            deskew_sr[j][k] <= '0;
                    end else if (en) begin
                        deskew_sr[j][0] <= cin;
                        for (int k = 1; k < D; k++)
                            deskew_sr[j][k] <= deskew_sr[j][k-1];
                    end
                end
                assign c_data[j*ACCW +: ACCW] = deskew_sr[j][D-1];
            end
        end
    endgenerate

    logic [LAT-1:0] vpipe;
    always_ff @(posedge clk or posedge rst) begin
        if (rst)
            vpipe <= '0;
        else if (en)
            vpipe <= {vpipe[LAT-2:0], a_valid};
    end
    assign c_valid = vpipe[LAT-1];

endmodule

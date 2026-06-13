#include "qarde_ems.hpp"
#include <hls_stream.h>

static void read_intrinsic_llr_ems(
    hls::stream<LLR_TYPE> &intrinsic_LLR,
    LLR_TYPE               L_buf[LDPC_N][GF_Q])
{
#pragma HLS INLINE off

read_L_outer:
    for (int v = 0; v < LDPC_N; ++v) {
    read_L_inner:
        for (int a = 0; a < GF_Q; ++a) {
#pragma HLS PIPELINE II=1
            L_buf[v][a] = intrinsic_LLR.read();
        }
    }
}

static void write_decisions_ems(
    GF_TYPE  x_hat_buf[LDPC_N],
    bool     synd_buf,
    GF_TYPE *decide,
    bool    &synd)
{
// Clang-format off
    #pragma HLS INLINE off
    // Clang-format on

write_decide_loop:
    for (int v = 0; v < LDPC_N; ++v) {
// Clang-format off
        #pragma HLS PIPELINE II=1
        // Clang-format on
        decide[v] = x_hat_buf[v];
    }
    synd = synd_buf;
}

void qarde_accel_ems(
    hls::stream<LLR_TYPE> &intrinsic_LLR,
    LLR_TYPE               alpha,
    LLR_TYPE               offset,
    LLR_TYPE               damp,
    int                    max_iter,
    ems_corr_mode_t        corr_mode,
    bool                  &synd,
    GF_TYPE               *decide)
{
    #pragma HLS INTERFACE axis      port=intrinsic_LLR
    #pragma HLS INTERFACE m_axi     port=decide      offset=slave bundle=gmem depth=LDPC_N
    #pragma HLS INTERFACE s_axilite port=alpha       bundle=control
    #pragma HLS INTERFACE s_axilite port=offset      bundle=control
    #pragma HLS INTERFACE s_axilite port=damp        bundle=control
    #pragma HLS INTERFACE s_axilite port=max_iter    bundle=control
    #pragma HLS INTERFACE s_axilite port=corr_mode   bundle=control
    #pragma HLS INTERFACE s_axilite port=decide      bundle=control
    #pragma HLS INTERFACE s_axilite port=synd        bundle=control
    #pragma HLS INTERFACE s_axilite port=return      bundle=control

    LLR_TYPE L_buf[LDPC_N][GF_Q];
    GF_TYPE  x_hat_buf[LDPC_N];
    bool     synd_buf;

// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=L_buf dim=2 cyclic factor=GF_FACTOR
    #pragma HLS BIND_STORAGE variable=L_buf     type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=x_hat_buf type=ram_t2p impl=bram
    #pragma HLS DATAFLOW
    // Clang-format on

    read_intrinsic_llr_ems(intrinsic_LLR, L_buf);

    qarde_decoder_ems<
        GF_Q, GF_FACTOR,
        LDPC_E, LDPC_N, LDPC_M,
        EMS_NM, NB_FACTOR,
        DEG_C, DEG_V,
        LDR_MIN, LDR_MAX,
        BUBBLE_HALF, BUBBLE, NBOPER>(
        L_buf,
        corr_mode,
        max_iter,
        alpha,
        offset,
        damp,
        x_hat_buf,
        synd_buf);

    write_decisions_ems(x_hat_buf, synd_buf, decide, synd);
}

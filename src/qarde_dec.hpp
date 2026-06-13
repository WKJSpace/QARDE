#ifndef QARDE_DEC_HPP
#define QARDE_DEC_HPP

#include "qarde_tools.hpp"
#include "qarde_gf.hpp"

template<int FP_N, int FP_M, int FP_Q, int Q_FACTOR, int FP_E, int FP_DEG_V, int FP_DEG_C>
struct qarde_dec {
    inline static void fp_post(
        LLR_TYPE L[FP_N][FP_Q],
        const EDGE_TYPE LDPC_adj_v[FP_N][FP_DEG_V],
        const LLR_TYPE LDPC_V_pv[FP_E][FP_Q],
        GF_TYPE x_hat[FP_N],
        GF_TYPE x_hat_synd[FP_DEG_C][FP_N]) {
// Clang-format off
        #pragma HLS INLINE off
        #pragma HLS DEPENDENCE dependent=false variable=LDPC_V_pv type=inter
        #pragma HLS DEPENDENCE dependent=false variable=LDPC_V_pv type=intra
        #pragma HLS ARRAY_PARTITION variable=x_hat_synd dim=1 complete
        // Clang-format on

    VAR_LOOP:
        for (int v = 0; v < FP_N; ++v) {
// Clang-format off
            #pragma HLS PIPELINE
            // Clang-format on

            LLR_TYPE post[FP_Q];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=post cyclic factor=Q_FACTOR
            // Clang-format on

            // Start from intrinsic L[v][a]
            for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                post[a] = L[v][a];
            }

            // Add all incoming C->V messages
            for (int t = 0; t < FP_DEG_V; ++t) {

                int e = (int)LDPC_adj_v[v][t];

                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    post[a] += LDPC_V_pv[e][a];
                }
            }

            // Argmax over symbols
            ValIdx<LLR_TYPE> top1[1]; // Largest first
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=top1 complete
            // Clang-format on

            find_topN<LLR_TYPE, FP_Q, 1>(post, top1, false);
            GF_TYPE hard_decision = top1[0].idx;
            x_hat[v] = hard_decision;

        SYND_BANK_LOOP:
            for (int t = 0; t < FP_DEG_C; ++t) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                x_hat_synd[t][v] = hard_decision;
            }
        }
    }

    static inline bool fp_check_syndrome(
        const GF_TYPE x_hat_synd[FP_DEG_C][FP_N],
        const NODE_TYPE LDPC_chk_v[FP_M][FP_DEG_C],
        const GF_TYPE LDPC_chk_h[FP_M][FP_DEG_C]) {
// Clang-format off
        #pragma HLS INLINE off
        #pragma HLS ARRAY_PARTITION variable=x_hat_synd dim=1 complete
        #pragma HLS ARRAY_PARTITION variable=LDPC_chk_v dim=2 complete
        #pragma HLS ARRAY_PARTITION variable=LDPC_chk_h dim=2 complete
        // Clang-format on

        bool all_ok = true;

    CHK_LOOP:
        for (int j = 0; j < FP_M; ++j) {
// Clang-format off
            #pragma HLS PIPELINE II=1
            // Clang-format on
            GF_TYPE sum = 0;

        EDGE_LOOP:
            for (int t = 0; t < FP_DEG_C; ++t) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                NODE_TYPE v = LDPC_chk_v[j][t];
                GF_TYPE h = LDPC_chk_h[j][t];
                sum = qarde_gf<FP_Q, Q_FACTOR>::gf_add(
                    sum,
                    qarde_gf<FP_Q, Q_FACTOR>::gf_mul(h, x_hat_synd[t][(int)v])
                );
            }

            if (sum != 0)
                all_ok = false;    // parity check failed
        }
        return all_ok;            // all checks passed
    }
};

#endif // QARDE_DEC_HPP

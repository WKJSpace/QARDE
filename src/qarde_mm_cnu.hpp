#ifndef QARDE_MM_CNU_HPP
#define QARDE_MM_CNU_HPP

#include "qarde_gf.hpp"
#include "qarde_tools.hpp"

template <int FP_Q, int Q_FACTOR, int FP_E, int FP_M,
          int FP_DEG_C, int FP_LDR_MIN, int FP_LDR_MAX>
struct qarde_mm_cnu
{
    // -----------------------------
    // Compile-time sanity checks
    // bitonicSort/find_topN require power-of-two sizes
    // -----------------------------
    static inline LLR_TYPE max2(LLR_TYPE a, LLR_TYPE b) { return (a > b) ? a : b; }

    // -----------------------------------------
    // Normalize wrt symbol-0 + clip
    // -----------------------------------------
    static inline void postprocess_ldr(LLR_TYPE U[FP_Q])
    {
        LLR_TYPE z = U[0];
        LLR_TYPE U_res[Q_FACTOR];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=U_res complete
        // Clang-format on
        for (int a = 0; a < (FP_Q / Q_FACTOR); ++a) {
            for (int b = 0; b < Q_FACTOR; b++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                LLR_TYPE v = U[a * Q_FACTOR + b] - z;
                if (v < (LLR_TYPE)FP_LDR_MIN) v = (LLR_TYPE)FP_LDR_MIN;
                if (v > (LLR_TYPE)FP_LDR_MAX) v = (LLR_TYPE)FP_LDR_MAX;
                U_res[b] = v;
            }
            for (int b = 0; b < Q_FACTOR; b++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                U[a * Q_FACTOR + b] = U_res[b];
            }
        }
    }

    // -----------------------------------------
    // Min-Max "convolution"
    //   C[a] = min_x max(A[x], B[a (+) x])
    // -----------------------------------------
    static void combine_minmax(const LLR_TYPE A[FP_Q],
                          const LLR_TYPE B[FP_Q],
                          LLR_TYPE       C[FP_Q])
{
// Clang-format off
    #pragma HLS INLINE off
    // Clang-format on
    static constexpr int NBLOCKS = FP_Q / Q_FACTOR;

row_loop:
    for (int a_out = 0; a_out < FP_Q; ++a_out) {
// Clang-format off
        #pragma HLS PIPELINE II=1
        // Clang-format on

        // Special case: Q_FACTOR == FP_Q  => only one block
        if constexpr (NBLOCKS == 1) {
            LLR_TYPE lane[Q_FACTOR];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=lane complete
            // Clang-format on

            for (int u = 0; u < Q_FACTOR; ++u) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                const int x = u;
                const int y = qarde_gf<FP_Q, Q_FACTOR>::gf_add_idx((GF_TYPE)a_out, (GF_TYPE)x);
                lane[u] = max2(A[x], B[y]);
            }

            ValIdx<LLR_TYPE> lane_sorted[Q_FACTOR];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=lane_sorted complete
            // Clang-format on
            bitonicSort<LLR_TYPE, Q_FACTOR>(lane, lane_sorted);
            C[a_out] = lane_sorted[0].val;
        } else {
            LLR_TYPE block_min_val[NBLOCKS];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=block_min_val complete
            // Clang-format on

        base_loop:
            for (int base = 0; base < FP_Q; base += Q_FACTOR) {
// Clang-format off
                #pragma HLS PIPELINE
                // Clang-format on
                LLR_TYPE lane[Q_FACTOR];
// Clang-format off
                #pragma HLS ARRAY_PARTITION variable=lane complete
                // Clang-format on

                for (int u = 0; u < Q_FACTOR; ++u) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    const int x = base + u;
                    const int y = qarde_gf<FP_Q, Q_FACTOR>::gf_add_idx((GF_TYPE)a_out, (GF_TYPE)x);
                    lane[u] = max2(A[x], B[y]);
                }

                ValIdx<LLR_TYPE> lane_sorted[Q_FACTOR];
// Clang-format off     
                #pragma HLS ARRAY_PARTITION variable=lane_sorted complete
                // Clang-format on
                bitonicSort<LLR_TYPE, Q_FACTOR>(lane, lane_sorted);
                block_min_val[base / Q_FACTOR] = lane_sorted[0].val;
            }

            ValIdx<LLR_TYPE> best_block[1];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=best_block complete
            // Clang-format on
            find_topN<LLR_TYPE, NBLOCKS, 1>(block_min_val, best_block, true);
            C[a_out] = best_block[0].val;
        }
    }
}


    // -----------------------------------------
    // Min-Max CNU for one CN: compute all edge outputs
    // Out[t][a] is the extrinsic msg for excluded edge t
    // -----------------------------------------
    static void cnu_minmax_check(const LLR_TYPE A[FP_DEG_C][FP_Q],
                                LLR_TYPE       Out[FP_DEG_C][FP_Q])
    {
// Clang-format off
        #pragma HLS INLINE off
        // Clang-format on

        LLR_TYPE F[FP_DEG_C][FP_Q];
        LLR_TYPE B[FP_DEG_C][FP_Q];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=F dim=1 complete
        #pragma HLS ARRAY_PARTITION variable=B dim=1 complete
        #pragma HLS ARRAY_PARTITION variable=F dim=2 cyclic factor=Q_FACTOR
        #pragma HLS ARRAY_PARTITION variable=B dim=2 cyclic factor=Q_FACTOR
        // Clang-format on
        for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            F[0][a] = A[0][a];
        }

        // Forward accumulation
    For_acc:
        for (int i = 1; i < FP_DEG_C; ++i) {
            combine_minmax(F[i - 1], A[i], F[i]);
        }

        // B[d-1] = A[d-1]
        for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            B[FP_DEG_C - 1][a] = A[FP_DEG_C - 1][a];
        }

        // Backward accumulation
    Back_acc:
        for (int i = FP_DEG_C - 2; i >= 0; --i) {
// Clang-format off
            // Clang_format on
            combine_minmax(B[i + 1], A[i], B[i]);
        }

        // Emit per excluded edge
        for (int t = 0; t < FP_DEG_C; ++t) {
            if (FP_DEG_C == 1) {
                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    Out[t][a] = (LLR_TYPE)0;
                }
            } else if (t == 0) {
                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    Out[t][a] = B[1][a];
                }
            } else if (t == FP_DEG_C - 1) {
                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    Out[t][a] = F[FP_DEG_C - 2][a];
                }
            } else {
                combine_minmax(F[t - 1], B[t + 1], Out[t]);
            }

            // normalize/clip only at the end (safe for Min-Max)
            postprocess_ldr(Out[t]);
        }
    }

    // -----------------------------------------
    // Run Min-Max CNU for all CNs/edges
    // -----------------------------------------
    static void cnu_run(EDGE_TYPE  LDPC_adj_c[FP_M][FP_DEG_C],
                        LLR_TYPE   LDPC_U_pc[FP_E][FP_Q],
                        LLR_TYPE   LDPC_V_cp[FP_E][FP_Q],
                        LLR_TYPE   damp)
    {
    CN_LOOP:
        for (int c = 0; c < FP_M; ++c) {
            // Load A[t] from U on all edges of this CN
            LLR_TYPE A[FP_DEG_C][FP_Q];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=A dim=1 complete
            #pragma HLS ARRAY_PARTITION variable=A dim=2 cyclic factor=Q_FACTOR
            // Clang-format on

        LOAD_A:
            for (int t = 0; t < FP_DEG_C; ++t) {
// Clang-format off
                #pragma HLS PIPELINE
                // Clang-format on
                const int e = (int)LDPC_adj_c[c][t];
                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    A[t][a] = LDPC_U_pc[e][a];
                }
            }

            // Compute extrinsic Out[t]
            LLR_TYPE Out[FP_DEG_C][FP_Q];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=Out dim=1 complete
            #pragma HLS ARRAY_PARTITION variable=Out dim=2 cyclic factor=Q_FACTOR
            // Clang-format on

            cnu_minmax_check(A, Out);

            // Damping + store
        WRITE_OUT:
            for (int t = 0; t < FP_DEG_C; ++t) {
// Clang-format off
                #pragma HLS PIPELINE II=3
                // Clang-format on
                const int e_out = (int)LDPC_adj_c[c][t];

                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    const LLR_TYPE oldv = LDPC_V_cp[e_out][a];
                    const LLR_TYPE m1   = ((LLR_TYPE)1.0 - damp) * Out[t][a];
                    const LLR_TYPE m2   = damp * oldv;
                    LDPC_V_cp[e_out][a] = m1 + m2;
                }
            }
        }
    }
};

#endif // QARDE_MM_CNU_HPP

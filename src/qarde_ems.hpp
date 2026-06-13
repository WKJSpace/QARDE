#ifndef QARDE_EMS_HPP
#define QARDE_EMS_HPP

#include "qarde_gf.hpp"
#include "qarde_init.hpp"
#include "qarde_ems_cnu.hpp"
#include "qarde_vnu.hpp"
#include "qarde_perm.hpp"
#include "nb_ldpc.hpp"
#include "qarde_dec.hpp"

const CHECK_TYPE H_rows[LDPC_E] = {0};
const NODE_TYPE H_cols[LDPC_E] = {0};
const GF_TYPE H_vals[LDPC_E] = {0};

template <int FP_Q, int Q_FACTOR,
          int FP_E, int FP_N, int FP_M, int FP_NM, int NM_FACTOR,
          int FP_DEG_C, int FP_DEG_V, int FP_LDR_MIN, int FP_LDR_MAX,
          int FP_BUBBLE_HALF, int FP_BUBBLE, int FP_NBOPER>
static void qarde_decoder_ems(
    LLR_TYPE         L[FP_N][FP_Q],
    ems_corr_mode_t  corr_mode,
    int              max_iter,
    LLR_TYPE         alpha,
    LLR_TYPE         offset,
    LLR_TYPE         damp,
    GF_TYPE          x_hat[FP_N],
    bool             &synd)
{
    static EDGE_TYPE LDPC_adj_v[FP_N][FP_DEG_V];
    static EDGE_TYPE LDPC_adj_c[FP_M][FP_DEG_C];
    static CHECK_TYPE LDPC_edge_c[FP_E];
    static NODE_TYPE LDPC_edge_v[FP_E];
    static GF_TYPE  LDPC_edge_h[FP_E];
    static NODE_TYPE LDPC_chk_v[FP_M][FP_DEG_C];
    static GF_TYPE  LDPC_chk_h[FP_M][FP_DEG_C];
    static GF_TYPE  LDPC_x_hat_synd[FP_DEG_C][FP_N];
    static int      LDPC_perm_to_chk[FP_E][FP_Q];
    static int      LDPC_perm_to_var[FP_E][FP_Q];
    static LLR_TYPE LDPC_U_vp[FP_E][FP_Q];
    static LLR_TYPE LDPC_U_pc[FP_E][FP_Q];
    static LLR_TYPE LDPC_V_cp[FP_E][FP_Q];
    static LLR_TYPE LDPC_V_pv[FP_E][FP_Q];

    static bool init_done = false;
    
// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=LDPC_adj_v       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_adj_c       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_chk_v       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_chk_h       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_x_hat_synd  dim=1 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_perm_to_chk dim=2 cyclic factor=Q_FACTOR
    #pragma HLS ARRAY_PARTITION variable=LDPC_perm_to_var dim=2 cyclic factor=Q_FACTOR
    #pragma HLS ARRAY_PARTITION variable=LDPC_U_vp        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_U_pc        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_V_cp        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_V_pv        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_edge_v      cyclic factor=FP_DEG_C
    #pragma HLS ARRAY_PARTITION variable=LDPC_edge_h      cyclic factor=FP_DEG_C
    #pragma HLS DEPENDENCE dependent=false variable=LDPC_U_pc type=inter
    #pragma HLS DEPENDENCE dependent=false variable=LDPC_U_pc type=intra
    #pragma HLS DEPENDENCE dependent=false variable=LDPC_V_pv type=inter
    #pragma HLS DEPENDENCE dependent=false variable=LDPC_V_pv type=intra

    // Clang-format on

    if (!init_done) {
        ldpc_graph_init(
            H_rows, H_cols, H_vals,
            LDPC_adj_v, LDPC_adj_c,
            LDPC_edge_c, LDPC_edge_v, LDPC_edge_h,
            LDPC_perm_to_chk, LDPC_perm_to_var,
            LDPC_U_vp, LDPC_U_pc, LDPC_V_cp, LDPC_V_pv
        );

    INIT_CHK_DIRECT:
        for (int j = 0; j < FP_M; ++j) {
// Clang-format off
            // Clang-format on
        INIT_CHK_DIRECT_INNER:
            for (int t = 0; t < FP_DEG_C; ++t) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                int e = (int)LDPC_adj_c[j][t];
                LDPC_chk_v[j][t] = LDPC_edge_v[e];
                LDPC_chk_h[j][t] = LDPC_edge_h[e];
            }
        }

        init_done = true;
    }

    bool conv = false;

ITER_LOOP:
    for (int it = 0; it < max_iter; ++it) {
// Clang-format off
        #pragma HLS LOOP_TRIPCOUNT min=1 max=10
        // Clang-format on

        if (conv) {
            break;
        }

        qarde_vnu<FP_Q, Q_FACTOR, FP_E,
                  FP_DEG_V, FP_N, FP_LDR_MIN, FP_LDR_MAX>::vnu_update(
            L, corr_mode, alpha, offset, LDPC_adj_v, LDPC_V_pv, LDPC_U_vp
        );

    PERM_V2C:
        for (int e = 0; e < FP_E; ++e) {
            qarde_perm<FP_Q, Q_FACTOR, FP_E>::perm_var_to_check_q(
                e, LDPC_perm_to_var, LDPC_U_vp, LDPC_U_pc
            );
        }

        qarde_ems_cnu<FP_Q, Q_FACTOR, FP_E, FP_M,
                      FP_NM, NM_FACTOR, FP_DEG_C,
                      FP_LDR_MIN, FP_LDR_MAX>::cnu_run(
            LDPC_adj_c, LDPC_U_pc, LDPC_V_cp, damp
        );

    PERM_C2V:
        for (int e = 0; e < FP_E; ++e) {
            qarde_perm<FP_Q, Q_FACTOR, FP_E>::perm_check_to_var_q(
                e, LDPC_perm_to_chk, LDPC_V_cp, LDPC_V_pv
            );
        }

        qarde_dec<FP_N, FP_M, FP_Q, Q_FACTOR, FP_E, FP_DEG_V, FP_DEG_C>::fp_post(
            L, LDPC_adj_v, LDPC_V_pv, x_hat, LDPC_x_hat_synd
        );

        if (qarde_dec<FP_N, FP_M, FP_Q, Q_FACTOR,
                      FP_E, FP_DEG_V, FP_DEG_C>::fp_check_syndrome(
                LDPC_x_hat_synd, LDPC_chk_v, LDPC_chk_h)) {
            conv = true;
        }
    }

    synd = !conv;
}

#endif // QARDE_EMS_HPP

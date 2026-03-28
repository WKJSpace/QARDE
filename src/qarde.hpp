#ifndef QARDE_HPP
#define QARDE_HPP

#include "qarde_gf.hpp"
#include "qarde_init.hpp"
#include "qarde_ems_cnu.hpp"
#include "qarde_mm_cnu.hpp"
#include "qarde_fbems_cnu.hpp"
#include "qarde_vnu.hpp"
#include "qarde_perm.hpp"
#include "nb_ldpc.hpp"
#include "qarde_dec.hpp"

const GF_TYPE H_rows[LDPC_E] = {0};
const GF_TYPE H_cols[LDPC_E] = {0};
const GF_TYPE H_vals[LDPC_E] = {0};

template <int FP_Q, int Q_FACTOR,
          int FP_E, int FP_N, int FP_M, int FP_NM, int NM_FACTOR,
          int FP_DEG_C, int FP_DEG_V, int FP_LDR_MIN, int FP_LDR_MAX,
          int FP_BUBBLE_HALF, int FP_BUBBLE, int FP_NBOPER>
static void qarde_decoder(
    LLR_TYPE         L[FP_N][FP_Q],
    ems_corr_mode_t  corr_mode,
    qarde_cnu_mode_t cnu_mode,
    int              max_iter,
    LLR_TYPE         alpha,
    LLR_TYPE         offset,
    LLR_TYPE         damp,
    GF_TYPE          x_hat[FP_N],
    bool             &synd)
{
    static GF_TYPE  LDPC_adj_v[FP_N][FP_DEG_V];
    static GF_TYPE  LDPC_adj_c[FP_M][FP_DEG_C];
    static GF_TYPE  LDPC_edge_c[FP_E];
    static GF_TYPE  LDPC_edge_v[FP_E];
    static GF_TYPE  LDPC_edge_h[FP_E];
    static GF_TYPE  LDPC_perm_to_chk[FP_E][FP_Q];
    static GF_TYPE  LDPC_perm_to_var[FP_E][FP_Q];
    static LLR_TYPE LDPC_U_vp[FP_E][FP_Q];
    static LLR_TYPE LDPC_U_pc[FP_E][FP_Q];
    static LLR_TYPE LDPC_V_cp[FP_E][FP_Q];
    static LLR_TYPE LDPC_V_pv[FP_E][FP_Q];

    static bool init_done = false;
    
    corr_mode = EMS_CORR_NONE;
    cnu_mode = QARDE_CNU_FBEMS;
    
// Clang-format off
    #pragma HLS ARRAY_PARTITION variable=LDPC_adj_v       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_adj_c       dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_perm_to_chk dim=2 cyclic factor=Q_FACTOR
    #pragma HLS ARRAY_PARTITION variable=LDPC_perm_to_var dim=2 cyclic factor=Q_FACTOR
    #pragma HLS ARRAY_PARTITION variable=LDPC_U_vp        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_U_pc        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_V_cp        dim=2 complete
    #pragma HLS ARRAY_PARTITION variable=LDPC_V_pv        dim=2 complete

    #pragma HLS BIND_STORAGE variable=LDPC_perm_to_chk type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_perm_to_var type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_U_pc        type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_V_pv        type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_edge_c      type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_edge_v      type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=LDPC_edge_h      type=ram_t2p impl=bram
    #pragma HLS RESET variable=init_done
    // Clang-format on

    if (!init_done) {
        qarde_gf<FP_Q, Q_FACTOR>::gf_init(0x43);

        ldpc_graph_init(
            H_rows, H_cols, H_vals,
            LDPC_adj_v, LDPC_adj_c,
            LDPC_edge_c, LDPC_edge_v, LDPC_edge_h,
            LDPC_perm_to_chk, LDPC_perm_to_var,
            LDPC_U_vp, LDPC_U_pc, LDPC_V_cp, LDPC_V_pv
        );

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
// Clang-format off
            #pragma HLS DATAFLOW
            // Clang-format on
            qarde_perm<FP_Q, Q_FACTOR, FP_E>::perm_var_to_check_q(
                e, LDPC_perm_to_chk, LDPC_U_vp, LDPC_U_pc
            );
        }

        switch (cnu_mode) {
        case QARDE_CNU_EMS:
            qarde_ems_cnu<FP_Q, Q_FACTOR, FP_E, FP_M,
                          FP_NM, NM_FACTOR, FP_DEG_C,
                          FP_LDR_MIN, FP_LDR_MAX>::cnu_run(
                LDPC_adj_c, LDPC_U_pc, LDPC_V_cp, damp
            );
            break;

        case QARDE_CNU_MM:
            qarde_mm_cnu<FP_Q, Q_FACTOR, FP_E, FP_M,
                         FP_DEG_C, FP_LDR_MIN, FP_LDR_MAX>::cnu_run(
                LDPC_adj_c, LDPC_U_pc, LDPC_V_cp, damp
            );
            break;

        case QARDE_CNU_FBEMS:
            qarde_fbems_cnu<FP_Q, Q_FACTOR, FP_E, FP_M,
                            FP_NM, NM_FACTOR, FP_DEG_C,
                            FP_LDR_MIN, FP_LDR_MAX,
                            FP_BUBBLE_HALF, FP_BUBBLE, FP_NBOPER>::cnu_run(
                LDPC_adj_c, LDPC_U_pc, LDPC_V_cp, damp
            );
            break;

        default:
            qarde_ems_cnu<FP_Q, Q_FACTOR, FP_E, FP_M,
                          FP_NM, NM_FACTOR, FP_DEG_C,
                          FP_LDR_MIN, FP_LDR_MAX>::cnu_run(
                LDPC_adj_c, LDPC_U_pc, LDPC_V_cp, damp
            );
            break;
        }

    PERM_C2V:
        for (int e = 0; e < FP_E; ++e) {
// Clang-format off
            #pragma HLS DATAFLOW
            // Clang-format on
            qarde_perm<FP_Q, Q_FACTOR, FP_E>::perm_check_to_var_q(
                e, LDPC_perm_to_var, LDPC_V_cp, LDPC_V_pv
            );
        }

        qarde_dec<FP_N, FP_M, FP_Q, Q_FACTOR, FP_E, FP_DEG_V, FP_DEG_C>::fp_post(
            L, LDPC_adj_v, LDPC_V_pv, x_hat
        );

        if (qarde_dec<FP_N, FP_M, FP_Q, Q_FACTOR,
                      FP_E, FP_DEG_V, FP_DEG_C>::fp_check_syndrome(
                x_hat, LDPC_adj_c, LDPC_edge_v, LDPC_edge_h)) {
            conv = true;
        }
    }

    synd = !conv;
}

#endif // QARDE_HPP
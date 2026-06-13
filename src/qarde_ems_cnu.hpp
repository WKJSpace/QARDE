#ifndef QARDE_EMS_CNU_HPP
#define QARDE_EMS_CNU_HPP

#include "qarde_gf.hpp"
#include "qarde_tools.hpp"

template <int FP_Q, int Q_FACTOR, int FP_E, int FP_M, int FP_NM, int NM_FACTOR,
    int FP_DEG_C, int FP_LDR_MIN, int FP_LDR_MAX>
struct qarde_ems_cnu {

    // -----------------------------------------
    //  Top-M truncation for one input message
    // -----------------------------------------
    static void fp_topM_costs(const LLR_TYPE U[FP_Q],
                              GF_TYPE        sym[FP_NM],
                              LLR_TYPE       cost[FP_NM],
                              LLR_TYPE       *Gamma_out)
    {
// Clang-format off
        #pragma HLS INLINE off
    	// Clang-format on

        ValIdx<LLR_TYPE> maxU[FP_NM];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=maxU complete
        // Clang-format on

        // Sort U in descending order and take top FP_NM
        find_topN<LLR_TYPE, FP_Q, FP_NM>(U, maxU, false);

        LLR_TYPE g = maxU[0].val;  // maximum value
        *Gamma_out = g;

        for (int m = 0; m < FP_NM; ++m) {
// Clang-format off
            #pragma HLS UNROLL
        	// Clang-format on
            cost[m] = g - maxU[m].val;
            sym[m]  = maxU[m].idx;
        }

        // Ensure index 0 is included
        int has0 = 0;
        for (int m = 0; m < FP_NM; ++m) {
// Clang-format off
            #pragma HLS UNROLL
        	// Clang-format on
            if (sym[m] == 0) {
                has0 = 1;
                break;
            }
        }

        if (!has0) {
            cost[FP_NM - 1] = g - U[0];
            sym[FP_NM - 1]  = 0;
        }
    }

    // -----------------------------------------
    //  Normalize + clip LDR
    // -----------------------------------------
    static inline void postprocess_ldr(LLR_TYPE U[FP_Q])
    {
        LLR_TYPE z = U[0], U_res[Q_FACTOR];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=U_res type=complete
        // Clang-format on
        for (int a = 0; a < (FP_Q / Q_FACTOR); ++a) {
// Clang-format off
        	// Clang-format on
		rec_U:
        	for (int b = 0 ; b < Q_FACTOR; b++) {
// Clang-format off
        		#pragma HLS UNROLL
        		// Clang-format on
        		LLR_TYPE v = U[a * Q_FACTOR + b] - z;
				if (v < (LLR_TYPE)FP_LDR_MIN) v = (LLR_TYPE)FP_LDR_MIN;
				if (v > (LLR_TYPE)FP_LDR_MAX) v = (LLR_TYPE)FP_LDR_MAX;
				U_res[b] = v;
        	}

        	for (int b = 0 ; b < Q_FACTOR; b++) {
// Clang-format off
				#pragma HLS UNROLL
				// Clang-format on
        		U[a * Q_FACTOR + b] = U_res[b];

        	}
        }
    }

    // -----------------------------------------
    //  EMS DP for one CNU output
    //  U_in[L][FP_Q], where L = FP_DEG_C - 1
    // -----------------------------------------
    static void cnu_ems(const LLR_TYPE U_in[FP_DEG_C - 1][FP_Q],
                        LLR_TYPE       V_out[FP_Q])
    {
        LLR_TYPE Gamma[FP_DEG_C];
        LLR_TYPE Gsum = 0.0;
        GF_TYPE  sym[FP_DEG_C][FP_NM];
        LLR_TYPE cost[FP_DEG_C][FP_NM];
        LLR_TYPE dpList[FP_NM];

// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=Gamma  dim=1 complete
        #pragma HLS ARRAY_PARTITION variable=sym    dim=2 cyclic    factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=cost   dim=2 cyclic    factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=dpList dim=1 complete
        // Clang-format on


        // Top-M truncation for each incoming message (except the current edge)
        for (int k = 0; k < (FP_DEG_C - 1); ++k) {
// Clang-format off
        	// Clang-format on
            fp_topM_costs(U_in[k], sym[k], cost[k], &Gamma[k]);
            Gsum += Gamma[k];
        }

        // DP over GF states: dp_cur -> dp_upd for each neighbor k
        LLR_TYPE dp_cur[FP_Q];
        LLR_TYPE dp_upd[FP_Q];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=dp_cur complete
        #pragma HLS ARRAY_PARTITION variable=dp_upd complete
        // Clang-format on

        // Initial DP state
        for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
            #pragma HLS UNROLL
        	// Clang-format on
            dp_cur[s] = (s == 0) ? (LLR_TYPE)0.0 : (LLR_TYPE)FP_LDR_MAX;
        }

	EDGE_LOOP:
        for (int k = 0; k < (FP_DEG_C - 1); ++k) {

            // Reset dp_upd for this edge
		upd_init:
            for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
                #pragma HLS UNROLL
            	// Clang-format on
                dp_upd[s] = (LLR_TYPE)FP_LDR_MAX;
            }

            // DP update
		upd_update:
            for (int ns = 0; ns < FP_Q; ++ns) {
// Clang-format off
                #pragma HLS PIPELINE II=1
            	// Clang-format on

			upd_reduce_m:
                for (int m = 0; m < FP_NM; ++m) {
// Clang-format off
                    #pragma HLS UNROLL
                	// Clang-format on
                    GF_TYPE u  = sym[k][m];
                    int a = qarde_gf<FP_Q, Q_FACTOR>::gf_add_idx((GF_TYPE)ns, u);
                    LLR_TYPE v = dp_cur[a] + cost[k][m];
                    dpList[m] = v;
                }
                ValIdx<LLR_TYPE> dp_best[1];
// Clang-format off
				#pragma HLS ARRAY_PARTITION variable=dp_best complete
                // Clang-format on

                find_topN<LLR_TYPE, FP_NM, 1>(dpList, dp_best, true);
                dp_upd[ns] = dp_best[0].val;
            }

            // Commit dp_upd -> dp_cur
            for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
                #pragma HLS UNROLL
            	// Clang-format on
                dp_cur[s] = dp_upd[s];
            }
        }

        // Final output
        for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
            #pragma HLS UNROLL
        	// Clang-format on
            V_out[s] = Gsum - dp_cur[s];
        }

        postprocess_ldr(V_out);
    }

    // -----------------------------------------
    //  Run EMS CNU for all checks and edges
    // -----------------------------------------
    static void cnu_run(EDGE_TYPE  LDPC_adj_c[FP_M][FP_DEG_C],
                        LLR_TYPE   LDPC_U_pc[FP_E][FP_Q],
                        LLR_TYPE   LDPC_V_cp[FP_E][FP_Q],
                        LLR_TYPE   damp) {
// Clang-format off
    	// Clang-format on
    CN_LOOP:
        for (int c = 0; c < FP_M; ++c) {
		CN_LOOP_DEGC:
            for (int t = 0; t < FP_DEG_C; ++t) {
                int e_out = (int)LDPC_adj_c[c][t];

                // Collect all U_pc except the one on edge e_out
                LLR_TYPE Uin[FP_DEG_C - 1][FP_Q];
// Clang-format off
                #pragma HLS ARRAY_PARTITION variable=Uin dim=2 complete
                // Clang-format on

                int L = 0;
			U_in_GEN:
                for (int r = 0; r < FP_DEG_C; ++r) {
// Clang-format off
                	#pragma HLS PIPELINE
                	// Clang-format on
                    if (r == t) continue;
                    int L_idx = (r < t) ? r : (r - 1);
                    int e_in = LDPC_adj_c[c][r];
                    for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    	#pragma HLS UNROLL
                    	// Clang-format on
                    	Uin[L_idx][a] = LDPC_U_pc[e_in][a];
                    }
                    L++;
                }

                LLR_TYPE tmp[FP_Q];
// Clang-format off
                #pragma HLS ARRAY_PARTITION variable=tmp cyclic factor=Q_FACTOR
                // Clang-format on

                // EMS for this check edge
                cnu_ems(Uin, tmp);

                // Damping: V_lp = (1 - damp) * tmp + damp * V_lp_old
                LLR_TYPE Vlp[FP_Q], Vlp_out[FP_Q];
// Clang-format off
                #pragma HLS ARRAY_PARTITION variable=Vlp     cyclic factor=Q_FACTOR
                #pragma HLS ARRAY_PARTITION variable=Vlp_out cyclic factor=Q_FACTOR
                // Clang-format on

                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                	// Clang-format on
                    Vlp[a] = LDPC_V_cp[e_out][a];
                }

                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                	// Clang-format on

                	LLR_TYPE m1, m2, s;
// Clang-format off
					#pragma HLS BIND_OP variable=m1 op=dmul impl=fulldsp
					#pragma HLS BIND_OP variable=m2 op=dmul impl=fulldsp
					#pragma HLS BIND_OP variable=s  op=dadd impl=fulldsp
                	// Clang-format on

					m1 = ((LLR_TYPE)1.0 - damp) * tmp[a];
					m2 = damp * Vlp[a];
					s  = m1 + m2;

					Vlp_out[a] = s;
                }

                for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
                    #pragma HLS UNROLL
                	// Clang-format on
                    LDPC_V_cp[e_out][a] = Vlp_out[a];
                }
            }
        }
    }
};

#endif // QARDE_EMS_CNU_HPP

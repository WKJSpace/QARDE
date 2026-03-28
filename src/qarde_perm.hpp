#ifndef QARDE_PERM_HPP
#define QARDE_PERM_HPP
#include "qarde_init.hpp"


template <int FP_Q, int Q_FACTOR, int FP_E>
struct qarde_perm {
	static void perm_var_to_check_q(int e,
			GF_TYPE    LDPC_perm_to_chk[FP_E][FP_Q],
	        LLR_TYPE   LDPC_U_vp[FP_E][FP_Q],
	        LLR_TYPE   LDPC_U_pc[FP_E][FP_Q]) {
// Clang-format off
		#pragma HLS INLINE off
		// CLang-format on

// Clang-format off
		#pragma HLS DEPENDENCE variable=LDPC_U_pc inter false
		#pragma HLS DEPENDENCE variable=LDPC_U_pc intra false
		// Clang-format on

	v2c_perm:
		for (int x = 0; x < FP_Q; x++) {
// Clang-format off
			#pragma HLS UNROLL factor=2
			// Clang-format on
			int y = (int)LDPC_perm_to_chk[e][x];
	        LDPC_U_pc[e][y] = LDPC_U_vp[e][x];
		}
	}

	static void perm_check_to_var_q(int e,
			GF_TYPE  LDPC_perm_to_var[FP_E][FP_Q],
			LLR_TYPE LDPC_V_cp[FP_E][FP_Q],
			LLR_TYPE LDPC_V_pv[FP_E][FP_Q]) {
// Clang-format off
		#pragma HLS INLINE off
		// Clang-format on

// Clang-format off
		#pragma HLS DEPENDENCE variable=LDPC_V_pv inter false
		#pragma HLS DEPENDENCE variable=LDPC_V_pv intra false
        // Clang-format on

	c2v_perm:
        for (int x = 0; x < FP_Q; ++x) {
// Clang-format off
            #pragma HLS UNROLL factor=2
        	// Clang-format on
        	int y   = (int)LDPC_perm_to_var[e][x];
			LDPC_V_pv[e][y] = LDPC_V_cp[e][x];
        }
    }
};

#endif
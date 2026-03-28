#ifndef QARDE_VNU_HPP
#define QARDE_VNU_HPP

#include "nb_ldpc.hpp"

template<int FP_Q, int Q_FACTOR, int FP_E, int FP_DEG_V, int FP_N, int FP_LDR_MIN, int FP_LDR_MAX>
struct qarde_vnu {
	static inline LLR_TYPE vn_correct(LLR_TYPE x,
	                                ems_corr_mode_t mode,
	                                LLR_TYPE alpha,
	                                LLR_TYPE offs)
	{
	    switch (mode) {
	    case EMS_CORR_SCALE:
	        if (alpha != 0.0) {
	            x = (LLR_TYPE)(x / alpha);
	        }
	        break;

	    case EMS_CORR_OFFSET:
	        if (x > 0.0) {
	            LLR_TYPE t = (LLR_TYPE)(x - offs);
	            x = (t > 0.0) ? t : (LLR_TYPE)0.0;
	        } else if (x < 0.0) {
	            LLR_TYPE t = (LLR_TYPE)(x + offs);
	            x = (t < 0.0) ? t : (LLR_TYPE)0.0;
	        }
	        break;

	    case EMS_CORR_NONE:
	    default:
	        break;
	    }
	    return x;
	}

	static inline void postprocess_ldr(LLR_TYPE *U)
	{
// Clang-format off
		#pragma HLS INLINE
		// Clang-format on
		ValIdx<LLR_TYPE> min_tmp[1]; // Find global minimum
// Clang-format off
		#pragma HLS ARRAY_PARTITION variable=min_tmp complete
		// Clang-format on
		find_topN<LLR_TYPE, FP_Q, 1>(U, min_tmp, true);
   		LLR_TYPE z = min_tmp[0].val;

	LDR_LOOP:
		for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
			LLR_TYPE v = U[a] - z;

			if (v < FP_LDR_MIN)      v = (LLR_TYPE)FP_LDR_MIN;
			else if (v > FP_LDR_MAX) v = (LLR_TYPE)FP_LDR_MAX;

			U[a] = v;
		}
	}

	/*
	 * VN update:
	 *  - G : Tanner graph (with embedded messages in edges)
	 *  - L : intrinsic LLRs, size [FP_N][FP_Q]
	 */
	static void vnu_update(LLR_TYPE	L[FP_N][FP_Q],
			ems_corr_mode_t mode,
			LLR_TYPE      	alpha,
			LLR_TYPE      	offs,
			GF_TYPE			LDPC_adj_v[FP_N][FP_DEG_V],
			LLR_TYPE		LDPC_V_pv[FP_E][FP_Q],
			LLR_TYPE		LDPC_U_vp[FP_E][FP_Q])
	{

		// Iterate over every variable node
	VNU_LOOP:
		for (int i = 0; i < FP_N; ++i) {
			LLR_TYPE S[FP_Q];
// Clang-format off
			#pragma HLS ARRAY_PARTITION variable=S type=cyclic factor=Q_FACTOR
			// Clang-format on
		VNU_INIT:
			for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
				#pragma HLS UNROLL factor=Q_FACTOR
				// Clang-format on
				S[a] = L[i][a];
			}

			LLR_TYPE Uvc[FP_Q];
// Clang-format off
			#pragma HLS ARRAY_PARTITION variable=Uvc type=cyclic factor=Q_FACTOR
			// Clang-format on
		VNU_ACC:
			for (int t = 0; t < FP_DEG_V; ++t) {
// Clang-format off
				#pragma HLS PIPELINE
				// Clang-format on
				int e_in = (int)LDPC_adj_v[i][t];
			VNU_CORR:
				for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
					#pragma HLS UNROLL factor=Q_FACTOR
					// Clang-format on
					LLR_TYPE vin = LDPC_V_pv[e_in][a];
					S[a] += vn_correct(vin, mode, alpha, offs);
				}
			}

			for (int t = 0; t < FP_DEG_V; ++t) {
				int e_out = (int)LDPC_adj_v[i][t];

				for (int a = 0; a < FP_Q; a++) {
// Clang-format off
					#pragma HLS UNROLL factor=Q_FACTOR
					// Clang-format on
					LLR_TYPE self = LDPC_V_pv[e_out][a];
					Uvc[a] = S[a] - vn_correct(self, mode, alpha, offs);
				}

				postprocess_ldr(Uvc);

				for (int a = 0; a < FP_Q; ++a) {
// Clang-format off
					#pragma HLS UNROLL factor=Q_FACTOR
					// Clang-format on
					LDPC_U_vp[e_out][a] = Uvc[a];
				}
			}
	    }
	}
};

#endif





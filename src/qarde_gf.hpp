#ifndef FP_GF_HPP
#define FP_GF_HPP

#include "nb_ldpc.hpp"

// One global set of tables sized by FP_Q

template<int FP_Q, int Q_FACTOR>
struct qarde_gf {

	static GF_TYPE ALOG_TBL[2 * FP_Q];
	static GF_TYPE LOG_TBL[FP_Q];


	static inline void gf_init(GF_TYPE prim_poly)
	{
// Clang-format off
		#pragma HLS ARRAY_PARTITION variable=ALOG_TBL type=cyclic factor=Q_FACTOR
		#pragma HLS ARRAY_PARTITION variable=LOG_TBL type=cyclic factor=Q_FACTOR
		#pragma HLS DEPENDENCE  dependent=false variable=ALOG_TBL type=inter
		#pragma HLS DEPENDENCE  dependent=false variable=ALOG_TBL type=intra
		// Clang-format on

	    GF_TYPE x = 1;
	    for (int i = 0; i < FP_Q - 1; ++i) {
// Clang-format off
	    	#pragma HLS PIPELINE
	    	// Clang-format on
	        ALOG_TBL[i] = x;
	        LOG_TBL[x]  = i;
	        x <<= 1;
	        if (x & FP_Q) {
	            x ^= prim_poly;
	        }
	    }

	    for (int i = FP_Q - 1; i < (2 * FP_Q - 2); ++i) {
// Clang-format off
	    	#pragma HLS UNROLL factor=Q_FACTOR
	    	// Clang-format on
	        ALOG_TBL[i] = ALOG_TBL[i - (FP_Q - 1)];
	    }

	    LOG_TBL[0] = 0;
	}

	static inline int gf_mul(GF_TYPE a, GF_TYPE b)
	{
	    if (a == 0 || b == 0) return 0;
	    GF_TYPE la  = LOG_TBL[a];
	    GF_TYPE lb  = LOG_TBL[b];
	    GF_TYPE s   = la + lb;
	    GF_TYPE qm1 = FP_Q - 1;
	    if (s >= qm1) s -= qm1;
	    return ALOG_TBL[s];
	}

	static inline GF_TYPE gf_add(GF_TYPE a, GF_TYPE b) { return a ^ b; }

	static inline GF_TYPE gf_inv(GF_TYPE a)
	{
	    if (a == 0) return 0;
	    GF_TYPE la  = LOG_TBL[a];
	    GF_TYPE qm1 = FP_Q - 1;
	    GF_TYPE e   = qm1 - la;
	    if (e >= qm1) e -= qm1;
	    return ALOG_TBL[e];
	}

	static inline GF_TYPE gf_div(GF_TYPE a, GF_TYPE b)
	{
	    if (a == 0) return 0;
	    if (b == 0) return 0;
	    GF_TYPE la  = LOG_TBL[a];
	    GF_TYPE lb  = LOG_TBL[b];
	    GF_TYPE qm1 = FP_Q - 1;
	    GF_TYPE d   = la - lb;
	    if (d < 0) d += qm1;
	    return ALOG_TBL[d];
	}

	static inline void gf_build_edge_perm(GF_TYPE h, GF_TYPE to_chk[FP_Q], GF_TYPE to_var[FP_Q])
	{
// Clang-format off
		#pragma HLS DEPENDENCE  dependent=false variable=to_var type=inter
		#pragma HLS DEPENDENCE  dependent=false variable=to_var type=intra
		// Clang-format on
	    if (h == 0) {
	        for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
	        	#pragma HLS UNROLL factor=Q_FACTOR
	        	// Clang-format on
	            to_chk[s] = s;
	            to_var[s] = s;
	        }
	        return;
	    }
	    for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
			#pragma HLS UNROLL factor=Q_FACTOR
			// Clang-format on
	        GF_TYPE y  = gf_mul(h, s);
	        to_chk[s]  = y;
	        to_var[y]  = s;
	    }
	}
};

template<int FP_Q, int Q_FACTOR>
GF_TYPE qarde_gf<FP_Q, Q_FACTOR>::ALOG_TBL[2 * FP_Q];

template<int FP_Q, int Q_FACTOR>
GF_TYPE qarde_gf<FP_Q, Q_FACTOR>::LOG_TBL[FP_Q];




#endif // FP_GF_HPP

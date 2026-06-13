#ifndef QARDE_GF_HPP
#define QARDE_GF_HPP

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

        int x = 1;
        for (int i = 0; i < FP_Q - 1; ++i) {
// Clang-format off
            #pragma HLS PIPELINE
            // Clang-format on
            ALOG_TBL[i] = (GF_TYPE)x;
            LOG_TBL[x]  = (GF_TYPE)i;
            x <<= 1;
            if (x & FP_Q) {
                x ^= (int)prim_poly;
            }
            x &= (FP_Q - 1);
        }

        for (int i = FP_Q - 1; i < (2 * FP_Q - 2); ++i) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            ALOG_TBL[i] = ALOG_TBL[i - (FP_Q - 1)];
        }

        LOG_TBL[0] = 0;
    }

    static inline int gf_mul(GF_TYPE a, GF_TYPE b)
    {
        int aa  = a.to_uint();
        int bb  = b.to_uint();
        int res = 0;
        const int prim = (int)GF_PRIM_POLY;
        const int mask = FP_Q - 1;
        const int msb  = FP_Q >> 1;
        const int bits = (FP_Q == 64) ? 6 : 3;

    gf_mul_loop:
        for (int i = 0; i < bits; ++i) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            if (bb & 1) res ^= aa;
            bb >>= 1;
            bool carry = (aa & msb) != 0;
            aa <<= 1;
            if (carry) aa ^= prim;
            aa &= mask;
        }
        return res;
    }

    static inline GF_TYPE gf_add(GF_TYPE a, GF_TYPE b) { return a ^ b; }

    static inline int gf_add_idx(GF_TYPE a, GF_TYPE b) { return ((int)a) ^ ((int)b); }

    static inline GF_TYPE gf_inv(GF_TYPE a)
    {
        if (a == 0) return 0;
        GF_TYPE res = 1;
    gf_inv_loop:
        for (int i = 0; i < FP_Q - 2; ++i) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            res = (GF_TYPE)gf_mul(res, a);
        }
        return res;
    }

    static inline GF_TYPE gf_div(GF_TYPE a, GF_TYPE b)
    {
        if (a == 0) return 0;
        if (b == 0) return 0;
        return (GF_TYPE)gf_mul(a, gf_inv(b));
    }

    static inline void gf_build_edge_perm(GF_TYPE h, int to_chk[FP_Q], int to_var[FP_Q])
    {
        GF_TYPE h_inv = gf_inv(h);

        for (int s = 0; s < FP_Q; ++s) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            to_chk[s] = (h == 0) ? s : gf_mul(h,     (GF_TYPE)s);
            to_var[s] = (h == 0) ? s : gf_mul(h_inv, (GF_TYPE)s);
        }
    }
};

template<int FP_Q, int Q_FACTOR>
GF_TYPE qarde_gf<FP_Q, Q_FACTOR>::ALOG_TBL[2 * FP_Q];

template<int FP_Q, int Q_FACTOR>
GF_TYPE qarde_gf<FP_Q, Q_FACTOR>::LOG_TBL[FP_Q];




#endif // QARDE_GF_HPP

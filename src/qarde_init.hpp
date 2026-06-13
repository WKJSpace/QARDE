#ifndef QARDE_INIT_HPP
#define QARDE_INIT_HPP

#include "nb_ldpc.hpp"
#include "qarde_gf.hpp"

// Here we just define the initializer.

static void ldpc_graph_init(const CHECK_TYPE rows[LDPC_E],
		const NODE_TYPE cols[LDPC_E], const GF_TYPE hs[LDPC_E],
		EDGE_TYPE LDPC_adj_v[LDPC_N][DEG_V], EDGE_TYPE LDPC_adj_c[LDPC_M][DEG_C],
		CHECK_TYPE LDPC_edge_c[LDPC_E], NODE_TYPE LDPC_edge_v[LDPC_E], GF_TYPE LDPC_edge_h[LDPC_E],
		int LDPC_perm_to_chk[LDPC_E][GF_Q], int LDPC_perm_to_var[LDPC_E][GF_Q],
		LLR_TYPE LDPC_U_vp[LDPC_E][GF_Q], LLR_TYPE LDPC_U_pc[LDPC_E][GF_Q],
		LLR_TYPE LDPC_V_lp[LDPC_E][GF_Q], LLR_TYPE LDPC_V_pv[LDPC_E][GF_Q]){

    // local degree counters
    int deg_v[LDPC_N];
    int deg_c[LDPC_M];
// Clang-format off
    #pragma HLS BIND_STORAGE variable=deg_v type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=deg_c type=ram_t2p impl=bram
    // Clang-format on

    // init degrees and adjacency
    for (int i = 0; i < LDPC_N; ++i) {
// Clang-format off
    	#pragma HLS PIPELINE
    	// Clang-format on
        deg_v[i] = 0;
        for (int t = 0; t < DEG_V; ++t)
// Clang-format off
        	#pragma HLS UNROLL
        	// Clang-format on
            LDPC_adj_v[i][t] = -1;
    }

    for (int j = 0; j < LDPC_M; ++j) {
// Clang-format off
		#pragma HLS PIPELINE
		// Clang-format on
        deg_c[j] = 0;
        for (int t = 0; t < DEG_C; ++t)
// Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
            LDPC_adj_c[j][t] = -1;
    }

    // build all edges from COO description
    for (int e = 0; e < LDPC_E; ++e) {
// Clang-format off
		#pragma HLS PIPELINE
		// Clang-format on
        CHECK_TYPE r = rows[e];  // check index
        NODE_TYPE c = cols[e];   // variable index
        GF_TYPE h = hs[e];    // GF coefficient

        LDPC_edge_c[e] = r;
        LDPC_edge_v[e] = c;
        LDPC_edge_h[e] = h;

        // permutations for this edge (h * s mapping)
        qarde_gf<GF_Q, GF_FACTOR>::gf_build_edge_perm(h, LDPC_perm_to_chk[e], LDPC_perm_to_var[e]);

        // variable side
        int dv = deg_v[(int)c];
        if (dv < DEG_V) {
            LDPC_adj_v[(int)c][dv] = e;
            deg_v[(int)c] = dv + 1;
        }

        // check side
        int dc = deg_c[(int)r];
        if (dc < DEG_C) {
            LDPC_adj_c[(int)r][dc] = e;
            deg_c[(int)r] = dc + 1;
        }

        // initialise messages on edge e
        for (int a = 0; a < GF_Q; ++a) {
// Clang-format off
			#pragma HLS UNROLL
			// Clang-format on
            LDPC_U_vp[e][a] = (LLR_TYPE)0.0;
            LDPC_U_pc[e][a] = (LLR_TYPE)0.0;
            LDPC_V_lp[e][a] = (LLR_TYPE)0.0;
            LDPC_V_pv[e][a] = (LLR_TYPE)0.0;
        }
    }

}

#endif // QARDE_INIT_HPP

#ifndef QARDE_FBEMS_CNU_HPP
#define QARDE_FBEMS_CNU_HPP

#include "qarde_gf.hpp"
#include "qarde_tools.hpp"
#include "ap_int.h"

// Mapping:
//   FP_Q            -> GF order (FP_GF in your FB code)
//   FP_DEG_C        -> row degree (FP_ROWDEG)
//   FP_NM           -> NBMAX (list size)
//   Q_FACTOR        -> unroll factor for GF loops

template <
    int FP_Q, int Q_FACTOR, int FP_E, int FP_M, int FP_NM, int NM_FACTOR,
    int FP_DEG_C, int FP_LDR_MIN, int FP_LDR_MAX,
    int FP_BUBBLE_HALF,
    int FP_BUBBLE,
    int FP_NBOPER
>
struct qarde_fbems_cnu {
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
            if (sym[m] == 0) { has0 = 1; break; }
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
        #pragma HLS ARRAY_PARTITION variable=U_res complete
        #pragma HLS ARRAY_PARTITION variable=U cyclic factor=Q_FACTOR
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
    //  FB-EMS ElementaryStep helpers
    // -----------------------------------------
    static inline LLR_TYPE INF_VAL() { return (LLR_TYPE)FP_LDR_MAX; }

    static void update_sorted_frontier(
        const ValIdx<LLR_TYPE> updated_entry,
        const int              pos,
        ValIdx<LLR_TYPE>       sorted[FP_BUBBLE]
    ){
// Clang-format off
        #pragma HLS INLINE
        #pragma HLS ARRAY_PARTITION variable=sorted complete
        // Clang-format on

        ValIdx<LLR_TYPE> original[FP_BUBBLE];
        int i_final = pos;
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=original complete
        // Clang-format on

    update_sorted_copy:
        for (int i = 0; i < FP_BUBBLE; i++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            original[i] = sorted[i];
        }

    update_sorted_scan:
        for (int step = 1; step < FP_BUBBLE; step++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            int left  = pos - step;
            int right = pos + step;

            bool use_left  = (left >= 0) && (updated_entry.val < original[left].val);
            bool use_right = (right < FP_BUBBLE) && (updated_entry.val > original[right].val);

            if (use_left) {
                i_final = left;
            } else if (use_right) {
                i_final = right;
            }
        }

    update_sorted_shift:
        for (int i = 0; i < FP_BUBBLE; i++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            if (i_final < pos && i <= pos && i > i_final) {
                sorted[i] = original[i - 1];
            } else if (i_final > pos && i >= pos && i < i_final) {
                sorted[i] = original[i + 1];
            } else if (i != i_final) {
                sorted[i] = original[i];
            }
        }

        sorted[i_final] = updated_entry;
    }

    static void ElementaryStepInit(
        const LLR_TYPE Input1[FP_NM],
        const LLR_TYPE Input2[FP_NM],
        LLR_TYPE       loc_Output[FP_NM],
        GF_TYPE        loc_IndiceOut[FP_NM],
        TabCompEntry<LLR_TYPE> tab_comp[FP_BUBBLE]
    ){
// Clang-format off
        #pragma HLS INLINE
        #pragma HLS ARRAY_PARTITION variable=Input1 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Input2 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=loc_Output cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=loc_IndiceOut cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=tab_comp complete
        // Clang-format on

        for (int i = 0; i < FP_NM; i++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            loc_Output[i]    = INF_VAL();
            loc_IndiceOut[i] = (GF_TYPE)-1;
        }


        // Seed bubble frontier: first column (rows 0..HALF-1)
        for (int j = 0; j < FP_BUBBLE_HALF; j++) {
            tab_comp[j] = { (LLR_TYPE)(Input1[j] + Input2[0]), j, 0 };
        }

        // Seed bubble frontier: middle row
        for (int j = 0; j < FP_BUBBLE_HALF; j++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            int idx = j + FP_BUBBLE_HALF;
            tab_comp[idx] = { (LLR_TYPE)(Input1[FP_BUBBLE_HALF] + Input2[j]), FP_BUBBLE_HALF, j };
        }
    }



    static void ElementaryStepCore(
        const LLR_TYPE  Input1[FP_NM],
        const LLR_TYPE  Input2[FP_NM],
        const GF_TYPE   IndiceInput1[FP_NM],
        const GF_TYPE   IndiceInput2[FP_NM],
        LLR_TYPE        loc_Output[FP_NM],
        GF_TYPE         loc_IndiceOut[FP_NM],
        TabCompEntry<LLR_TYPE> tab_comp[FP_BUBBLE]
    ){
// Clang-format off
        #pragma HLS INLINE
        #pragma HLS ARRAY_PARTITION variable=Input1 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Input2 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=tab_comp complete
        #pragma HLS ARRAY_PARTITION variable=loc_Output cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=loc_IndiceOut cyclic factor=NM_FACTOR
        // Clang-format on

        ValIdx<LLR_TYPE> sorted_validx[FP_BUBBLE];
        LLR_TYPE vals[FP_BUBBLE];
        bool GF_out_mask[FP_Q];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=sorted_validx complete
        #pragma HLS ARRAY_PARTITION variable=vals complete
        #pragma HLS ARRAY_PARTITION variable=GF_out_mask complete
        // Clang-format on

    init_mask_loop:
        for (int i = 0; i < FP_Q; ++i) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            GF_out_mask[i] = false;
        }

        for (int i = 0; i < FP_BUBBLE; i++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on
            vals[i] = tab_comp[i].val;
        }

        // sorted_validx[0] is the smallest entry (assumed by your existing bitonicSort usage)
        bitonicSort<LLR_TYPE, FP_BUBBLE>(vals, sorted_validx);

        constexpr int HALF = (FP_BUBBLE >> 1);
        int s = 0;

        NBOPER_LOOP:
        for (int ss = 0; ss < FP_NBOPER; ss++) {
// Clang-format off
            #pragma HLS PIPELINE II=7
            // Clang-format on

            int pos = sorted_validx[0].idx;
            LLR_TYPE best_val = sorted_validx[0].val;
            int x = tab_comp[pos].x;
            int y = tab_comp[pos].y;

            if ((IndiceInput1[x] == (GF_TYPE)-1) || (IndiceInput2[y] == (GF_TYPE)-1)) break;

            int idx_sum = qarde_gf<FP_Q, Q_FACTOR>::gf_add_idx(IndiceInput1[x], IndiceInput2[y]);

            if (!GF_out_mask[idx_sum]) {
                loc_Output[s]    = best_val;
                loc_IndiceOut[s] = (GF_TYPE)idx_sum;
                GF_out_mask[idx_sum] = true;
                s++;
            }

            if ((s == FP_NM) || (x >= FP_NM - 1) || (y >= FP_NM - 1)) break;

            bool is_vertical = (pos > (HALF - 1));
            int next_x = x + (is_vertical ? 1 : 0);
            int next_y = y + (is_vertical ? 0 : 1);

            tab_comp[pos].x = next_x;
            tab_comp[pos].y = next_y;

            LLR_TYPE new_val = (LLR_TYPE)(Input1[next_x] + Input2[next_y]);
            tab_comp[pos].val = new_val;

            ValIdx<LLR_TYPE> updated_entry = { new_val, pos };
            int insert_pos = -1;

            for (int i = 0; i < FP_BUBBLE; i++) {
// Clang-format off
                #pragma HLS UNROLL
                // // Clang-format off
                if (sorted_validx[i].idx == pos) insert_pos = i;
            }

            update_sorted_frontier(updated_entry, insert_pos, sorted_validx);
        }
    }

    static void ElementaryStep(
        const LLR_TYPE Input1[FP_NM],
        const LLR_TYPE Input2[FP_NM],
        const GF_TYPE  IndiceInput1[FP_NM],
        const GF_TYPE  IndiceInput2[FP_NM],
        LLR_TYPE       Output[FP_NM],
        GF_TYPE        IndiceOutput[FP_NM]
    ){
// Clang-format off
        #pragma HLS INLINE
        // Clang-format on

        LLR_TYPE loc_Output[FP_NM];
        GF_TYPE  loc_IndiceOut[FP_NM];
        TabCompEntry<LLR_TYPE> tab_comp[FP_BUBBLE];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=loc_Output cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=loc_IndiceOut cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=tab_comp complete
        // Clang-format on

        ElementaryStepInit(Input1, Input2, loc_Output, loc_IndiceOut, tab_comp);
        ElementaryStepCore(Input1, Input2, IndiceInput1, IndiceInput2, loc_Output, loc_IndiceOut, tab_comp);

        for (int i = 0; i < FP_NM; i++) {
// Clang-format off
            #pragma HLS UNROLL
            // Clang-format on 
            Output[i]      = loc_Output[i];
            IndiceOutput[i]= loc_IndiceOut[i];
        }
    }

    // -----------------------------------------
    //  FB recursion per check
    // -----------------------------------------
    static void FB_Recursion(
        const LLR_TYPE M_VtoC_LLR[FP_DEG_C][FP_NM],
        const GF_TYPE  M_VtoC_GF [FP_DEG_C][FP_NM],
        LLR_TYPE       MatriceInter    [2*(FP_DEG_C-2)][FP_NM],
        GF_TYPE        MatriceInter_idx[2*(FP_DEG_C-2)][FP_NM],
        LLR_TYPE       OutForward[FP_NM],
        GF_TYPE        OutForward_idx[FP_NM],
        LLR_TYPE       OutBackward[FP_NM],
        GF_TYPE        OutBackward_idx[FP_NM]
    ){
// Clang-format off
        #pragma HLS INLINE off
        // Clang-format on

        LLR_TYPE OutForward1[FP_NM], OutBackward1[FP_NM];
        GF_TYPE  OutForward1_idx[FP_NM], OutBackward1_idx[FP_NM];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=OutForward1 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=OutForward1_idx cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=OutBackward1 cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=OutBackward1_idx cyclic factor=NM_FACTOR
        // Clang-format on
        FB_recursion_loop:
        for (int kk = 1; kk < FP_DEG_C - 1; kk++) {
// Clang-format off
            // #pragma HLS PIPELINE
            // Clang-format on
            int outback   = FP_DEG_C - kk - 1;
            int for_inter = kk - 1;
            int back_inter= (2*(FP_DEG_C-2)) - kk;

            // load next neighbors
            for (int k = 0; k < FP_NM; k++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                OutForward1[k]     = M_VtoC_LLR[kk][k];
                OutForward1_idx[k] = M_VtoC_GF [kk][k];
                OutBackward1[k]    = M_VtoC_LLR[outback][k];
                OutBackward1_idx[k]= M_VtoC_GF [outback][k];
            }

            // store intermediate
            for (int k = 0; k < FP_NM; k++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                MatriceInter[for_inter][k]      = OutForward[k];
                MatriceInter_idx[for_inter][k]  = OutForward_idx[k];
                MatriceInter[back_inter][k]     = OutBackward[k];
                MatriceInter_idx[back_inter][k] = OutBackward_idx[k];
            }

            // update forward/backward
            ElementaryStep(OutForward,  OutForward1,  OutForward_idx,  OutForward1_idx,  OutForward,  OutForward_idx);
            ElementaryStep(OutBackward, OutBackward1, OutBackward_idx, OutBackward1_idx, OutBackward, OutBackward_idx);
        }
    }

    static void MatriceInter_Process(
        const LLR_TYPE MatriceInter    [2*(FP_DEG_C-2)][FP_NM],
        const GF_TYPE  MatriceInter_idx[2*(FP_DEG_C-2)][FP_NM],
        LLR_TYPE       M_CtoV_LLR[FP_DEG_C][FP_NM],
        GF_TYPE        M_CtoV_GF [FP_DEG_C][FP_NM]
    ){
// Clang-format off
        #pragma HLS INLINE off
        // Clang-format on

        LLR_TYPE For_tmp[FP_NM], Back_tmp[FP_NM], Out_tmp[FP_NM];
        GF_TYPE  For_idx[FP_NM], Back_idx[FP_NM], Out_idx[FP_NM];
// Clang-format off
        #pragma HLS ARRAY_PARTITION variable=For_tmp  cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Back_tmp cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Out_tmp  cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=For_idx  cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Back_idx cyclic factor=NM_FACTOR
        #pragma HLS ARRAY_PARTITION variable=Out_idx  cyclic factor=NM_FACTOR
        // Clang-format on

        MatriceInter_loop:
        for (int k = 0; k < (FP_DEG_C - 2); k++) {
// Clang-format off
            // #pragma HLS PIPELINE
            // Clang-format on

            for (int i = 0; i < FP_NM; i++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                For_tmp[i]  = MatriceInter[k][i];
                For_idx[i]  = MatriceInter_idx[k][i];
                Back_tmp[i] = MatriceInter[(FP_DEG_C - 2) + k][i];
                Back_idx[i] = MatriceInter_idx[(FP_DEG_C - 2) + k][i];
            }

            ElementaryStep(For_tmp, Back_tmp, For_idx, Back_idx, Out_tmp, Out_idx);

            for (int g = 0; g < FP_NM; g++) {
                #pragma HLS UNROLL
                M_CtoV_LLR[k + 1][g] = Out_tmp[g];
                M_CtoV_GF [k + 1][g] = Out_idx[g];
            }
        }
    }

    // -----------------------------------------
    //  Expand list-domain CtoV (cost,index) -> full FP_Q vector
    //  (no DIVDEC / coefficient undo here; assumed handled outside)
    // -----------------------------------------
    static void FinalizeCtoV_NoCoeff(
        const LLR_TYPE M_CtoV_LLR_list[FP_DEG_C][FP_NM],
        const GF_TYPE  M_CtoV_GF_list [FP_DEG_C][FP_NM],
        LLR_TYPE       V_full[FP_DEG_C][FP_Q]
    ){
// Clang-format off
        #pragma HLS INLINE off
        #pragma HLS ARRAY_PARTITION variable=V_full dim=2 cyclic factor=Q_FACTOR
        // Clang-format on

        // offset is not in fp cnu_run signature; keep 0 here
        const LLR_TYPE offset = (LLR_TYPE)0.0;

        for (int t = 0; t < FP_DEG_C; t++) {
// Clang-format off
            #pragma HLS PIPELINE II=10
            // Clang-format on

            int Stp = FP_NM;
            for (int i = 0; i < FP_NM; i++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                if (M_CtoV_GF_list[t][i] == (GF_TYPE)-1) { Stp = i; break; }
            }
            if (Stp <= 0) Stp = 1;

            LLR_TYPE thr = M_CtoV_LLR_list[t][Stp - 1];

            LLR_TYPE fill_val = (thr > (LLR_TYPE)0.0) ? (LLR_TYPE)(thr + offset) : (LLR_TYPE)0.0;

            for (int k = 0; k < FP_Q; k++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                LLR_TYPE v = fill_val;
                for (int i = 0; i < FP_NM; i++) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    int idx = M_CtoV_GF_list[t][i].to_int();
                    if (idx == k) {
                        v = M_CtoV_LLR_list[t][i];
                    }
                }
                if ((k == 0) && (v < (LLR_TYPE)0.0)) {
                    v = (LLR_TYPE)0.0;
                }
                V_full[t][k] = v;
            }

            // Normalize + clip to match your original interface
            postprocess_ldr(V_full[t]);
        }
    }

    // -----------------------------------------
    //  Run FB-EMS CNU for all checks and edges
    // -----------------------------------------
    static void cnu_run(EDGE_TYPE  LDPC_adj_c[FP_M][FP_DEG_C],
                        LLR_TYPE LDPC_U_pc[FP_E][FP_Q],
                        LLR_TYPE LDPC_V_cp[FP_E][FP_Q],
                        LLR_TYPE damp)
    {
    CN_LOOP:
        for (int c = 0; c < FP_M; ++c) {

            // Build list-domain incoming messages for this CN
            LLR_TYPE M_VtoC_LLR[FP_DEG_C][FP_NM];
            GF_TYPE  M_VtoC_GF [FP_DEG_C][FP_NM];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=M_VtoC_LLR dim=2 cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=M_VtoC_GF  dim=2 cyclic factor=NM_FACTOR
            // Clang-format on

            for (int t = 0; t < FP_DEG_C; t++) {
// Clang-format off
                #pragma HLS PIPELINE
                // Clang-format on
                int e_in = (int)LDPC_adj_c[c][t];

                LLR_TYPE Gamma_dummy;
                fp_topM_costs(LDPC_U_pc[e_in], M_VtoC_GF[t], M_VtoC_LLR[t], &Gamma_dummy);
            }

            // Init intermediates
            constexpr int FP_FB_STEP = 2 * (FP_DEG_C - 2);

            LLR_TYPE MatriceInter[FP_FB_STEP][FP_NM];
            GF_TYPE  MatriceInter_idx[FP_FB_STEP][FP_NM];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=MatriceInter dim=2 cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=MatriceInter_idx dim=2 cyclic factor=NM_FACTOR
            // Clang-format on
            for (int k = 0; k < FP_FB_STEP; k++) {
// Clang-format off
                #pragma HLS PIPELINE
                // Clang-format on
                for (int i = 0; i < FP_NM; i++) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    MatriceInter[k][i]     = INF_VAL();
                    MatriceInter_idx[k][i] = (GF_TYPE)-1;
                }
            }

            // Forward/backward accumulators
            LLR_TYPE OutForward[FP_NM], OutBackward[FP_NM];
            GF_TYPE  OutForward_idx[FP_NM], OutBackward_idx[FP_NM];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=OutForward cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=OutForward_idx cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=OutBackward cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=OutBackward_idx cyclic factor=NM_FACTOR
            // Clang-format on

            for (int i = 0; i < FP_NM; i++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                OutForward[i]     = M_VtoC_LLR[0][i];
                OutForward_idx[i] = M_VtoC_GF [0][i];
                OutBackward[i]    = M_VtoC_LLR[FP_DEG_C - 1][i];
                OutBackward_idx[i]= M_VtoC_GF [FP_DEG_C - 1][i];
            }

            // Run FB recursion
            FB_Recursion(M_VtoC_LLR, M_VtoC_GF,
                         MatriceInter, MatriceInter_idx,
                         OutForward, OutForward_idx,
                         OutBackward, OutBackward_idx);

            // Collect list-domain CtoV for all branches
            LLR_TYPE M_CtoV_LLR_list[FP_DEG_C][FP_NM];
            GF_TYPE  M_CtoV_GF_list [FP_DEG_C][FP_NM];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=M_CtoV_LLR_list dim=1 complete
            #pragma HLS ARRAY_PARTITION variable=M_CtoV_LLR_list dim=2 cyclic factor=NM_FACTOR
            #pragma HLS ARRAY_PARTITION variable=M_CtoV_GF_list  dim=1 complete
            #pragma HLS ARRAY_PARTITION variable=M_CtoV_GF_list  dim=2 cyclic factor=NM_FACTOR
            // Clang-format on

            // first (t=0) and last (t=deg-1)
            for (int k = 0; k < FP_NM; k++) {
// Clang-format off
                #pragma HLS UNROLL
                // Clang-format on
                M_CtoV_LLR_list[FP_DEG_C - 1][k] = OutForward[k];
                M_CtoV_GF_list [FP_DEG_C - 1][k] = OutForward_idx[k];
                M_CtoV_LLR_list[0][k]            = OutBackward[k];
                M_CtoV_GF_list [0][k]            = OutBackward_idx[k];
            }

            // middle branches via MatriceInter combine
            MatriceInter_Process(MatriceInter, MatriceInter_idx,
                                 M_CtoV_LLR_list, M_CtoV_GF_list);

            // Expand to full vectors + normalize/clip
            LLR_TYPE V_full[FP_DEG_C][FP_Q];
// Clang-format off
            #pragma HLS ARRAY_PARTITION variable=V_full dim=2 cyclic factor=Q_FACTOR
            // Clang-format on

            FinalizeCtoV_NoCoeff(M_CtoV_LLR_list, M_CtoV_GF_list, V_full);

            // Write back to edges with damping
            for (int t = 0; t < FP_DEG_C; t++) {
// Clang-format off
                #pragma HLS PIPELINE II=4
                // Clang-format on
                int e_out = (int)LDPC_adj_c[c][t];

                for (int a = 0; a < FP_Q; a++) {
// Clang-format off
                    #pragma HLS UNROLL
                    // Clang-format on
                    LLR_TYPE oldv = LDPC_V_cp[e_out][a];
                    LLR_TYPE newv = ((LLR_TYPE)1.0 - damp) * V_full[t][a] + damp * oldv;
                    LDPC_V_cp[e_out][a] = newv;
                }
            }
        }
    }
};

#endif // QARDE_FBEMS_CNU_HPP

#include <iostream>
#include <hls_stream.h>
#include "../src/nb_ldpc.hpp"

void qarde_accel_mm(
    hls::stream<LLR_TYPE> &intrinsic_LLR,
    LLR_TYPE               alpha,
    LLR_TYPE               offset,
    LLR_TYPE               damp,
    int                    max_iter,
    ems_corr_mode_t        corr_mode,
    bool                  &synd,
    GF_TYPE               *decide);

int main()
{
    hls::stream<LLR_TYPE> intrinsic_LLR;
    GF_TYPE decide[LDPC_N];
    bool synd = false;

load_llr_outer:
    for (int v = 0; v < LDPC_N; ++v) {
    load_llr_inner:
        for (int a = 0; a < GF_Q; ++a) {
            intrinsic_LLR.write((a == 0) ? (LLR_TYPE)0 : (LLR_TYPE)4);
        }
    }

    qarde_accel_mm(
        intrinsic_LLR,
        (LLR_TYPE)1,
        (LLR_TYPE)0,
        (LLR_TYPE)0,
        0,
        EMS_CORR_NONE,
        synd,
        decide);

    if (!intrinsic_LLR.empty()) {
        std::cerr << "FAIL: top did not consume intrinsic LLR stream\n";
        return 1;
    }

    std::cout << "qarde_accel_mm smoke test passed\n";
    return 0;
}

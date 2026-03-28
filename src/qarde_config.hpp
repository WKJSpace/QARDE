#ifndef QARDE_CONFIG_HPP
#define QARDE_CONFIG_HPP

#include <hls_stream.h>
#include "qarde.hpp"

void qarde_accel(
    hls::stream<LLR_TYPE> &intrinsic_LLR,
    LLR_TYPE               alpha,
    LLR_TYPE               offset,
    LLR_TYPE               damp,
    int                    max_iter,
    ems_corr_mode_t        corr_mode,
    qarde_cnu_mode_t       cnu_mode,
    bool                  &synd,
    GF_TYPE               *decide);

#endif
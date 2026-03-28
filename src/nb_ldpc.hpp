#ifndef NB_LDPC_HPP
#define NB_LDPC_HPP

#include "qarde_tools.hpp"

// ============================================================
// Global design parameters
// ============================================================

// Galois field order GF(2^M) = GF_Q
#ifndef GF_Q
#define GF_Q        64
#endif

#ifndef GF_FACTOR
#define GF_FACTOR   64
#endif

#ifndef GF_ORDER
#define GF_ORDER    6
#endif

// Tanner graph size
#ifndef LDPC_N
#define LDPC_N      1440    // number of variable nodes
#endif

#ifndef LDPC_M
#define LDPC_M      960     // number of check nodes
#endif

#ifndef DEG_C
#define DEG_C       3       // check node degree
#endif

#ifndef DEG_V
#define DEG_V       2       // variable node degree
#endif

#ifndef LDPC_E
#define LDPC_E      (LDPC_M * DEG_C)
#endif

// EMS truncation parameters
#ifndef EMS_NM
#define EMS_NM      16      // top-M list size
#endif

#ifndef NB_FACTOR
#define NB_FACTOR   16     // parallel factor for NM
#endif

#ifndef EMS_NC
#define EMS_NC      16      // state pruning after convolutions
#endif

#ifndef BUBBLE_HALF
#define BUBBLE_HALF 2
#endif

#ifndef BUBBLE
#define BUBBLE      (2 * BUBBLE_HALF)
#endif

#ifndef NBOPER
#define NBOPER      25
#endif

#ifndef FIX_MIN16
#define FIX_MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

// ============================================================
// Decoder mode / correction mode
// ============================================================

// Correction mode in VNU
typedef enum {
    EMS_CORR_NONE = 0,
    EMS_CORR_SCALE,
    EMS_CORR_OFFSET
} ems_corr_mode_t;

// CNU mode selection
typedef enum {
    QARDE_CNU_EMS = 0,
    QARDE_CNU_MM,
    QARDE_CNU_FBEMS
} qarde_cnu_mode_t;

// ============================================================
// Data types
// ============================================================

typedef SIGNED_DATA_TYPE(GF_Q) GF_TYPE;
typedef ap_fixed<16, 7>        LLR_TYPE;

// ============================================================
// LLR clipping range
// ============================================================

static const int LDR_MIN = -50;
static const int LDR_MAX =  50;

// ============================================================
// H matrix in COO form
// ============================================================

extern const GF_TYPE H_rows[LDPC_E];
extern const GF_TYPE H_cols[LDPC_E];
extern const GF_TYPE H_vals[LDPC_E];

#endif // NB_LDPC_HPP
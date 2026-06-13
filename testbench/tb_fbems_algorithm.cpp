
#include <cmath>
#include <iostream>

static int failures = 0;

static void expect_true(bool cond, const char *msg)
{
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        failures++;
    }
}

static void expect_near(double got, double exp, const char *msg)
{
    double diff = got - exp;
    if (diff < 0) diff = -diff;
    if (diff > 0.01) {
        std::cerr << "FAIL: " << msg << " got=" << got << " expected=" << exp << "\n";
        failures++;
    }
}

#include "../src/qarde_gf.hpp"
#include "../src/qarde_vnu.hpp"
#include "../src/qarde_fbems_cnu.hpp"

static void test_gf64_primitive()
{
    qarde_gf<64, 64>::gf_init(GF_PRIM_POLY);
    int seen[64] = {0};
    for (int i = 0; i < 63; ++i) {
        int v = (int)qarde_gf<64, 64>::ALOG_TBL[i];
        expect_true(v > 0 && v < 64, "GF(64) alog entry in range");
        seen[v]++;
    }
    for (int v = 1; v < 64; ++v) expect_true(seen[v] == 1, "GF(64) alog visits each non-zero symbol once");
    expect_true((int)qarde_gf<64, 64>::gf_mul(2, 32) == 3, "GF(64) multiply uses x^6+x+1 primitive polynomial");
}

static void test_vnu_symbol0_normalization()
{
    LLR_TYPE U[64];
    double orig[64];
    for (int i = 0; i < 64; ++i) {
        orig[i] = (double)(i % 11 - 3);
        U[i] = (LLR_TYPE)orig[i];
    }
    double z = orig[0];
    qarde_vnu<64, 64, 1, 1, 1, -50, 50>::postprocess_ldr(U);
    for (int i = 0; i < 64; ++i) expect_near((double)U[i], orig[i] - z, "VNU keeps symbol-0 as LDR reference");
}

static void brute_elementary(const LLR_TYPE A[4], const LLR_TYPE B[4], const GF_TYPE AI[4], const GF_TYPE BI[4], LLR_TYPE O[4], GF_TYPE OI[4])
{
    bool used[64] = {false};
    for (int k = 0; k < 4; ++k) {
        LLR_TYPE best = (LLR_TYPE)50;
        int best_idx = -1;
        int best_sym = -1;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                int sym = (int)qarde_gf<64, 64>::gf_add(AI[i], BI[j]);
                if (used[sym]) continue;
                LLR_TYPE val = A[i] + B[j];
                if (val < best) {
                    best = val;
                    best_idx = i * 4 + j;
                    best_sym = sym;
                }
            }
        }
        expect_true(best_idx >= 0, "brute elementary found candidate");
        O[k] = best;
        OI[k] = (GF_TYPE)best_sym;
        used[best_sym] = true;
    }
}

static void test_fbems_elementary_step_matches_brute_force()
{
    LLR_TYPE A[4] = {0, 1, 4, 7};
    LLR_TYPE B[4] = {0, 2, 3, 6};
    GF_TYPE AI[4] = {0, 1, 2, 4};
    GF_TYPE BI[4] = {0, 3, 5, 7};
    LLR_TYPE got[4], exp[4];
    GF_TYPE got_i[4], exp_i[4];
    qarde_fbems_cnu<64, 64, 3, 1, 4, 4, 3, -50, 50, 2, 4, 16>::ElementaryStep(A, B, AI, BI, got, got_i);
    brute_elementary(A, B, AI, BI, exp, exp_i);
    for (int i = 0; i < 4; ++i) {
        expect_near((double)got[i], (double)exp[i], "FB-EMS ElementaryStep cost equals brute force");
        expect_true((int)got_i[i] == (int)exp_i[i], "FB-EMS ElementaryStep symbol equals brute force");
    }
}

int main()
{
    test_gf64_primitive();
    test_vnu_symbol0_normalization();
    test_fbems_elementary_step_matches_brute_force();
    if (failures) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "algorithm test passed\n";
    return 0;
}


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
#include "../src/qarde_mm_cnu.hpp"

static void test_gf8_primitive()
{
    qarde_gf<8, 8>::gf_init(GF_PRIM_POLY);
    int seen[8] = {0};
    for (int i = 0; i < 7; ++i) {
        int v = (int)qarde_gf<8, 8>::ALOG_TBL[i];
        expect_true(v > 0 && v < 8, "GF(8) alog entry in range");
        seen[v]++;
    }
    for (int v = 1; v < 8; ++v) expect_true(seen[v] == 1, "GF(8) alog visits each non-zero symbol once");
    expect_true((int)qarde_gf<8, 8>::gf_mul(2, 4) == 3, "GF(8) multiply uses x^3+x+1 primitive polynomial");
    expect_true((int)qarde_gf<8, 8>::gf_mul(3, 7) == 2, "GF(8) multiply wraps correctly");
}

static void test_vnu_symbol0_normalization()
{
    LLR_TYPE U[8] = {5, 2, 7, -1, 9, 3, 5, 8};
    qarde_vnu<8, 8, 1, 1, 1, -50, 50>::postprocess_ldr(U);
    double exp[8] = {0, -3, 2, -6, 4, -2, 0, 3};
    for (int i = 0; i < 8; ++i) expect_near((double)U[i], exp[i], "VNU keeps symbol-0 as LDR reference");
}

static void brute_minmax(const LLR_TYPE A[8], const LLR_TYPE B[8], LLR_TYPE C[8])
{
    for (int s = 0; s < 8; ++s) {
        LLR_TYPE best = (LLR_TYPE)50;
        for (int a = 0; a < 8; ++a) {
            int b = (int)qarde_gf<8, 8>::gf_add((GF_TYPE)s, (GF_TYPE)a);
            LLR_TYPE val = (A[a] > B[b]) ? A[a] : B[b];
            if (val < best) best = val;
        }
        C[s] = best;
    }
}

static void test_minmax_combine_matches_brute_force()
{
    LLR_TYPE A[8] = { 0, 3, 1, 7, 2, 4, 8, 5 };
    LLR_TYPE B[8] = { 0, 6, 2, 5, 1, 9, 3, 4 };
    LLR_TYPE got[8];
    LLR_TYPE exp[8];
    qarde_mm_cnu<8, 8, 3, 1, 3, -50, 50>::combine_minmax(A, B, got);
    brute_minmax(A, B, exp);
    for (int s = 0; s < 8; ++s) expect_near((double)got[s], (double)exp[s], "Minmax combine equals brute-force min-max convolution");
}

int main()
{
    test_gf8_primitive();
    test_vnu_symbol0_normalization();
    test_minmax_combine_matches_brute_force();
    if (failures) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "algorithm test passed\n";
    return 0;
}

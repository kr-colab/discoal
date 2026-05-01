#include "unity.h"
#include "discoal.h"
#include "shapes.h"
#include <math.h>
#include <stdlib.h>

extern long seed1, seed2;
extern void setall(long iseed1, long iseed2);
extern double *currentSize;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    /* Initialize the RNG so ranf() is usable from KS tests later */
    seed1 = 12345;
    seed2 = 67890;
    setall(seed1, seed2);

    /* Allocate currentSize if not already; tests for
     * initializeShapesFromGlobals write into currentSize[i]. */
    if (currentSize == NULL) {
        currentSize = (double *)calloc(MAXPOPS, sizeof(double));
    } else {
        for (int i = 0; i < MAXPOPS; i++) currentSize[i] = 0.0;
    }

    /* Reset shape state for every test */
    for (int i = 0; i < MAXPOPS; i++) {
        popShape[i].type = SHAPE_CONSTANT;
        popShape[i].anchor_value = 1.0;
        popShape[i].rate_param = 0.0;
        popShape[i].anchor_time = 0.0;
        for (int j = 0; j < MAXPOPS; j++) {
            migShape[i][j].type = SHAPE_CONSTANT;
            migShape[i][j].anchor_value = 0.0;
            migShape[i][j].rate_param = 0.0;
            migShape[i][j].anchor_time = 0.0;
        }
    }
}

void tearDown(void) { }
#endif

void test_sizeAt_constant_returns_anchor(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.5;
    popShape[0].anchor_time = 10.0;
    /* For CONSTANT, t should be irrelevant */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 50.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 1e6));
}

void test_sizeAt_constant_per_population(void) {
    popShape[1].type = SHAPE_CONSTANT;
    popShape[1].anchor_value = 0.5;
    popShape[2].type = SHAPE_CONSTANT;
    popShape[2].anchor_value = 4.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, sizeAt(1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 4.0, sizeAt(2, 0.0));
}

void test_sizeAt_exponential_at_anchor_time(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 10.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, sizeAt(0, 10.0));
}

void test_sizeAt_exponential_decays_backward(void) {
    /* alpha > 0 means N decreases as t grows past anchor_time */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    /* t=2: N = 1 * exp(-0.5*2) = exp(-1) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-1.0), sizeAt(0, 2.0));
    /* t=4: N = exp(-2) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-2.0), sizeAt(0, 4.0));
}

void test_sizeAt_exponential_grows_forward(void) {
    /* For t < anchor_time, N is larger */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 5.0;
    /* t=3: N = 1 * exp(-0.5 * (3-5)) = exp(1) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(1.0), sizeAt(0, 3.0));
}

void test_sizeAt_exponential_zero_alpha_equals_constant(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 3.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 1e6));
}

void test_sizeAt_linear_at_anchor(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.5;
    popShape[0].rate_param = 0.25;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.5, sizeAt(0, 0.0));
}

void test_sizeAt_linear_declines_backward_with_positive_gamma(void) {
    /* gamma > 0 = forward growth = backward decline (past smaller) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 5.0;
    popShape[0].rate_param = 2.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, sizeAt(0, 2.0));   /* 5 - 2*2 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 1.0));   /* 5 - 2*1 */
}

void test_sizeAt_linear_grows_backward_with_negative_gamma(void) {
    /* gamma < 0 = forward decline = backward growth (past larger) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = -1.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, sizeAt(0, 1.0));   /* 1 - (-1)*1 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 5.0, sizeAt(0, 4.0));   /* 1 - (-1)*4 */
}

void test_migAt_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, migAt(0, 1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, migAt(0, 1, 100.0));
}

void test_migAt_exponential(void) {
    /* m(t) = 1.0 * exp(-0.5 * (t - 0)). At t=2: exp(-1). */
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-1.0), migAt(0, 1, 2.0));
}

void test_migAt_linear_declines_backward_with_positive_delta(void) {
    /* delta > 0 = forward growth = backward decline.
     * m(t) = 0.5 - 0.05*t. At t=4: 0.5 - 0.2 = 0.3. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migAt(0, 1, 4.0));
}

void test_migAt_linear_grows_backward_with_negative_delta(void) {
    /* delta < 0 = forward decline = backward growth.
     * m(t) = 0.1 - (-0.05)*t = 0.1 + 0.05 t. At t=10: 0.6. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.6, migAt(0, 1, 10.0));
}

void test_migAt_pair_isolation(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.3;
    migShape[1][0].type = SHAPE_CONSTANT;
    migShape[1][0].anchor_value = 0.7;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migAt(0, 1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.7, migAt(1, 0, 0.0));
}

void test_integratedHazardSize_constant(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.0;
    /* k=2, T=4: H = 1*4/2 = 2 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, integratedHazardSize(0, 0.0, 4.0, 2));
    /* k=4, T=10: H = 6*10/2 = 30 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 30.0, integratedHazardSize(0, 0.0, 10.0, 4));
}

void test_integratedHazardSize_exponential_at_anchor(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    /* k=2, T=4, alpha=0.5, N_0 at t=0 = 2.0 */
    /* H = 1/(2*0.5) * (exp(2) - 1) = (exp(2)-1) */
    double expected = (exp(2.0) - 1.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardSize(0, 0.0, 4.0, 2));
}

void test_integratedHazardSize_exponential_offset_t0(void) {
    /* If t0 != anchor_time, the relevant N_0 is sizeAt(t0) */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.2;
    popShape[0].anchor_time = 0.0;
    /* At t0=5, N(t0) = exp(-1). Then H over T=2 with that as anchor:
     * H = 1/(N(t0)*0.2) * (exp(0.2*2) - 1) = 1/(exp(-1)*0.2) * (exp(0.4)-1) */
    double N_t0 = exp(-1.0);
    double expected = 1.0 * (exp(0.4) - 1.0) / (N_t0 * 0.2);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardSize(0, 5.0, 2.0, 2));
}

void test_integratedHazardSize_exponential_alpha_zero(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 3.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    /* alpha=0 should match constant: H = 1*5/3 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 5.0/3.0, integratedHazardSize(0, 0.0, 5.0, 2));
}

void test_integratedHazardSize_exponential_quadrature_match(void) {
    /* High-resolution numerical quadrature should match the closed form */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.5;
    popShape[0].rate_param = 0.7;
    popShape[0].anchor_time = 2.0;
    double t0 = 3.5;
    double T = 4.0;
    int k = 5;
    double pairs = k*(k-1)/2.0;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_lo = i * dt;
        double s_hi = (i+1) * dt;
        double s_mid = 0.5 * (s_lo + s_hi);
        sum += pairs / sizeAt(0, t0 + s_mid) * dt;  /* midpoint rule */
    }
    double closed = integratedHazardSize(0, t0, T, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, sum, closed);
}

void test_integratedHazardSize_linear_basic(void) {
    /* gamma = -0.5 (forward decline => backward growth):
     * N(s) = 1 - (-0.5)*s = 1 + 0.5 s, k=2, T=4
     * H = 1/(-0.5) * log(1/(1+0.5*4)) = -2 * log(1/3) = 2*log(3) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = -0.5;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 2.0 * log(3.0), integratedHazardSize(0, 0.0, 4.0, 2));
}

void test_integratedHazardSize_linear_zero_gamma(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    /* matches constant */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.0 * 5.0 / 2.0, integratedHazardSize(0, 0.0, 5.0, 2));
}

void test_integratedHazardSize_linear_diverges_at_zero_crossing(void) {
    /* gamma = 0.5 (forward growth => backward decline):
     * N(t) = 1 - 0.5 t hits zero at t=2. Integration to T=2.5 should diverge. */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double H = integratedHazardSize(0, 0.0, 2.5, 2);
    TEST_ASSERT_TRUE(isinf(H) || H > 1e15);
}

void test_integratedHazardSize_linear_quadrature_match(void) {
    /* Tolerance is 1e-5, not 1e-6, because the integrand 1/N(s) ranges from
     * 0.54 to 20 over the interval (N at the far end is 0.05), and midpoint-
     * rule O(dt^2) error at N=16384 is ~4e-6 for this parameterization. The
     * closed form is exact within float precision; this test bounds residual
     * algebraic error well below any sign/factor mistake. The EXP analog
     * passes 1e-6 because its integrand is smoother. */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.3;
    popShape[0].anchor_time = 1.0;
    double t0 = 1.5;
    double T = 6.0;
    int k = 4;
    double pairs = k*(k-1)/2.0;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_mid = (i + 0.5) * dt;
        sum += pairs / sizeAt(0, t0 + s_mid) * dt;
    }
    double closed = integratedHazardSize(0, t0, T, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-5, sum, closed);
}

void test_integratedHazardMig_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    /* k=3, T=4: H = 3*0.5*4 = 6 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 6.0, integratedHazardMig(0, 1, 0.0, 4.0, 3));
}

void test_integratedHazardMig_exponential(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    /* k=2, T=4: m(s) = exp(-0.5 s), H = 2 * (1 - exp(-2))/0.5 */
    double expected = 2.0 * (1.0 - exp(-2.0)) / 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardMig(0, 1, 0.0, 4.0, 2));
}

void test_integratedHazardMig_linear(void) {
    /* delta = -0.05 (forward decline => backward growth):
     * m(s) = 0.1 - (-0.05)*s = 0.1 + 0.05 s, k=2, T=4
     * H = 2*(0.1*4 - (-0.05)*16/2) = 2*(0.4 + 0.4) = 1.6 */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.6, integratedHazardMig(0, 1, 0.0, 4.0, 2));
}

void test_integratedHazardMig_quadrature_match(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.3;
    migShape[0][1].rate_param = 0.7;
    migShape[0][1].anchor_time = 1.0;
    double t0 = 1.5;
    double T = 3.0;
    int k = 4;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_mid = (i + 0.5) * dt;
        sum += k * migAt(0, 1, t0 + s_mid) * dt;
    }
    double closed = integratedHazardMig(0, 1, t0, T, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, sum, closed);
}

void test_drawWaitingTimeSize_constant_inverts(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.0;
    /* xi = 3, k=2: T = 3*2/1 = 6 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 6.0, drawWaitingTimeSize(0, 0.0, 3.0, 2));
}

void test_drawWaitingTimeSize_exponential_inverts(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double xi = 2.0;
    int k = 2;
    /* T such that H(T)=xi: T = log(1 + N0*alpha*xi/pairs)/alpha = log(1+1)/0.5 = 2*log(2) */
    double T = drawWaitingTimeSize(0, 0.0, xi, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 2.0 * log(2.0), T);
    /* Verify round trip: H(T) == xi */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardSize(0, 0.0, T, k));
}

void test_drawWaitingTimeSize_linear_inverts(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double xi = 1.0;
    int k = 2;
    double T = drawWaitingTimeSize(0, 0.0, xi, k);
    /* Round trip */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardSize(0, 0.0, T, k));
}

void test_drawWaitingTimeSize_returns_negative_for_k_lt_2(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 1.0;
    TEST_ASSERT_TRUE(drawWaitingTimeSize(0, 0.0, 1.0, 1) < 0.0);
    TEST_ASSERT_TRUE(drawWaitingTimeSize(0, 0.0, 1.0, 0) < 0.0);
}

void test_drawWaitingTimeSize_linear_zero_crossing(void) {
    /* gamma = 0.5 (forward growth => backward decline):
     * N(t) = 1 - 0.5 t hits zero at t=2. H(2) = +inf, so any finite xi maps to T<2. */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    int k = 2;
    /* Try a battery of xi values; T must always be < 2. */
    for (double xi = 0.1; xi < 100.0; xi *= 1.5) {
        double T = drawWaitingTimeSize(0, 0.0, xi, k);
        TEST_ASSERT_TRUE(T > 0.0);
        TEST_ASSERT_TRUE(T < 2.0);
    }
}

void test_drawWaitingTimeMig_constant_inverts(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    /* xi=2, k=4: T = 2/(4*0.5) = 1 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, drawWaitingTimeMig(0, 1, 0.0, 2.0, 4));
}

void test_drawWaitingTimeMig_exponential_round_trip(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.3;
    migShape[0][1].anchor_time = 0.0;
    double xi = 1.0;
    int k = 3;
    double T = drawWaitingTimeMig(0, 1, 0.0, xi, k);
    TEST_ASSERT_TRUE(T > 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardMig(0, 1, 0.0, T, k));
}

void test_drawWaitingTimeMig_exponential_unreachable_xi(void) {
    /* m(s) = 0.5 exp(-0.3 s). Total integrated hazard over [0,inf) for k=3:
     * 3 * 0.5 / 0.3 = 5. Any xi >= 5 should be unreachable. */
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.3;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 5.5, 3) < 0.0);
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 100.0, 3) < 0.0);
}

void test_drawWaitingTimeMig_linear_round_trip(void) {
    /* delta = -0.05 (forward decline => backward growth) so m(s) = 0.1+0.05s
     * grows monotonically and the round-trip is well-defined for any xi. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    double xi = 0.8;
    int k = 4;
    double T = drawWaitingTimeMig(0, 1, 0.0, xi, k);
    TEST_ASSERT_TRUE(T > 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardMig(0, 1, 0.0, T, k));
}

void test_drawWaitingTimeMig_zero_rate_returns_negative(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.0;
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 1.0, 4) < 0.0);
}

extern double ranf(void);  /* from xoshiro256pp_compat */

void test_drawWaitingTimeSize_linear_zero_crossing_stress(void) {
    /* Battery of (gamma, N0) configurations that drive N -> 0 going backward.
     * Under msprime forward-time convention, gamma > 0 (forward growth) means
     * past was smaller, hitting zero at T_cross = N0/gamma. */
    struct { double N0; double gamma; } configs[] = {
        {1.0, 0.1}, {1.0, 1.0}, {1.0, 10.0},
        {0.01, 0.001}, {100.0, 50.0}
    };
    for (int c = 0; c < (int)(sizeof(configs)/sizeof(configs[0])); c++) {
        popShape[0].type = SHAPE_LINEAR;
        popShape[0].anchor_value = configs[c].N0;
        popShape[0].rate_param = configs[c].gamma;
        popShape[0].anchor_time = 0.0;
        double T_cross = configs[c].N0 / configs[c].gamma;
        for (int trial = 0; trial < 1000; trial++) {
            double xi = -log(ranf());
            double T = drawWaitingTimeSize(0, 0.0, xi, 2);
            TEST_ASSERT_TRUE(T > 0.0);
            TEST_ASSERT_TRUE(T < T_cross);
            TEST_ASSERT_FALSE(isnan(T));
            TEST_ASSERT_FALSE(isinf(T));
        }
    }
}

static int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static double ks_statistic_size(int popID, int k, int n) {
    double *T = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        double xi = -log(ranf());
        T[i] = drawWaitingTimeSize(popID, 0.0, xi, k);
    }
    qsort(T, n, sizeof(double), compare_doubles);
    double max_d = 0.0;
    for (int i = 0; i < n; i++) {
        double F_emp = (double)(i + 1) / n;
        double F_theory = 1.0 - exp(-integratedHazardSize(popID, 0.0, T[i], k));
        double d = fabs(F_emp - F_theory);
        if (d > max_d) max_d = d;
    }
    free(T);
    return max_d;
}

void test_drawWaitingTimeSize_distribution_constant(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 1.0;
    int n = 10000;
    double D = ks_statistic_size(0, 2, n);
    /* For n=10000, KS critical at 5% is ~0.0136. Use 0.05 as forgiving threshold. */
    TEST_ASSERT_TRUE(D < 0.05);
}

void test_drawWaitingTimeSize_distribution_exponential(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    int n = 10000;
    double D = ks_statistic_size(0, 4, n);
    TEST_ASSERT_TRUE(D < 0.05);
}

void test_drawWaitingTimeSize_distribution_linear(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.3;
    popShape[0].anchor_time = 0.0;
    int n = 10000;
    double D = ks_statistic_size(0, 3, n);
    TEST_ASSERT_TRUE(D < 0.05);
}

static double ks_statistic_mig(int src, int dst, int k, int n) {
    double *T = malloc(sizeof(double) * n);
    int valid = 0;
    for (int i = 0; i < n; i++) {
        double xi = -log(ranf());
        double t = drawWaitingTimeMig(src, dst, 0.0, xi, k);
        if (t > 0.0) T[valid++] = t;
    }
    if (valid < n / 2) { free(T); return 1.0; }  /* too many unreachable; fail */
    qsort(T, valid, sizeof(double), compare_doubles);
    double max_d = 0.0;
    for (int i = 0; i < valid; i++) {
        double F_emp = (double)(i + 1) / valid;
        double F_theory = 1.0 - exp(-integratedHazardMig(src, dst, 0.0, T[i], k));
        /* Note: F_theory here is conditional-on-finite, which matches valid sample */
        double d = fabs(F_emp - F_theory);
        if (d > max_d) max_d = d;
    }
    free(T);
    return max_d;
}

void test_drawWaitingTimeMig_distribution_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.3;
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 4, 10000) < 0.05);
}

void test_drawWaitingTimeMig_distribution_exponential(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    /* H(infty) = k*m0/beta = 4*1/0.5 = 8. Most xi will resolve. */
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 4, 10000) < 0.05);
}

void test_drawWaitingTimeMig_distribution_linear(void) {
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.1;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 3, 10000) < 0.05);
}

extern double migMatConst[MAXPOPS][MAXPOPS];
extern int npops;

void test_initializeShapesFromGlobals_copies_currentSize(void) {
    /* Set globals as if parameters were just parsed */
    npops = 3;
    currentSize[0] = 1.0;
    currentSize[1] = 0.5;
    currentSize[2] = 2.0;

    /* Pre-condition: popShape entries set to defaults by setUp */
    /* Pre-condition: migMatConst entries are zero by default */

    initializeShapesFromGlobals();

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[0].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, popShape[0].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, popShape[0].rate_param);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, popShape[0].anchor_time);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[1].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, popShape[1].anchor_value);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[2].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, popShape[2].anchor_value);
}

void test_initializeShapesFromGlobals_copies_migMatConst(void) {
    npops = 2;
    currentSize[0] = 1.0;
    currentSize[1] = 1.0;
    migMatConst[0][1] = 0.7;
    migMatConst[1][0] = 0.3;
    migMatConst[0][0] = 0.0;  /* diagonals are zero */
    migMatConst[1][1] = 0.0;

    initializeShapesFromGlobals();

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, migShape[0][1].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.7, migShape[0][1].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][1].rate_param);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][1].anchor_time);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, migShape[1][0].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migShape[1][0].anchor_value);

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][0].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[1][1].anchor_value);
}

void test_allShapesConstant_default_state_is_true(void) {
    /* setUp resets all shapes to SHAPE_CONSTANT */
    npops = 3;
    TEST_ASSERT_EQUAL(1, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pop_is_exponential(void) {
    npops = 3;
    popShape[1].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pop_is_linear(void) {
    npops = 3;
    popShape[2].type = SHAPE_LINEAR;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pair_is_exponential(void) {
    npops = 2;
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_ignores_pops_outside_npops(void) {
    /* Shapes for pops 5..MAXPOPS that are not in use should not affect the answer */
    npops = 2;
    popShape[5].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(1, allShapesConstant());
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_sizeAt_constant_returns_anchor);
    RUN_TEST(test_sizeAt_constant_per_population);
    RUN_TEST(test_sizeAt_exponential_at_anchor_time);
    RUN_TEST(test_sizeAt_exponential_decays_backward);
    RUN_TEST(test_sizeAt_exponential_grows_forward);
    RUN_TEST(test_sizeAt_exponential_zero_alpha_equals_constant);
    RUN_TEST(test_sizeAt_linear_at_anchor);
    RUN_TEST(test_sizeAt_linear_declines_backward_with_positive_gamma);
    RUN_TEST(test_sizeAt_linear_grows_backward_with_negative_gamma);
    RUN_TEST(test_migAt_constant);
    RUN_TEST(test_migAt_exponential);
    RUN_TEST(test_migAt_linear_declines_backward_with_positive_delta);
    RUN_TEST(test_migAt_linear_grows_backward_with_negative_delta);
    RUN_TEST(test_migAt_pair_isolation);
    RUN_TEST(test_integratedHazardSize_constant);
    RUN_TEST(test_integratedHazardSize_exponential_at_anchor);
    RUN_TEST(test_integratedHazardSize_exponential_offset_t0);
    RUN_TEST(test_integratedHazardSize_exponential_alpha_zero);
    RUN_TEST(test_integratedHazardSize_exponential_quadrature_match);
    RUN_TEST(test_integratedHazardSize_linear_basic);
    RUN_TEST(test_integratedHazardSize_linear_zero_gamma);
    RUN_TEST(test_integratedHazardSize_linear_diverges_at_zero_crossing);
    RUN_TEST(test_integratedHazardSize_linear_quadrature_match);
    RUN_TEST(test_integratedHazardMig_constant);
    RUN_TEST(test_integratedHazardMig_exponential);
    RUN_TEST(test_integratedHazardMig_linear);
    RUN_TEST(test_integratedHazardMig_quadrature_match);
    RUN_TEST(test_drawWaitingTimeSize_constant_inverts);
    RUN_TEST(test_drawWaitingTimeSize_exponential_inverts);
    RUN_TEST(test_drawWaitingTimeSize_linear_inverts);
    RUN_TEST(test_drawWaitingTimeSize_returns_negative_for_k_lt_2);
    RUN_TEST(test_drawWaitingTimeSize_linear_zero_crossing);
    RUN_TEST(test_drawWaitingTimeMig_constant_inverts);
    RUN_TEST(test_drawWaitingTimeMig_exponential_round_trip);
    RUN_TEST(test_drawWaitingTimeMig_exponential_unreachable_xi);
    RUN_TEST(test_drawWaitingTimeMig_linear_round_trip);
    RUN_TEST(test_drawWaitingTimeMig_zero_rate_returns_negative);
    RUN_TEST(test_drawWaitingTimeSize_distribution_constant);
    RUN_TEST(test_drawWaitingTimeSize_distribution_exponential);
    RUN_TEST(test_drawWaitingTimeSize_distribution_linear);
    RUN_TEST(test_drawWaitingTimeMig_distribution_constant);
    RUN_TEST(test_drawWaitingTimeMig_distribution_exponential);
    RUN_TEST(test_drawWaitingTimeMig_distribution_linear);
    RUN_TEST(test_drawWaitingTimeSize_linear_zero_crossing_stress);
    RUN_TEST(test_initializeShapesFromGlobals_copies_currentSize);
    RUN_TEST(test_initializeShapesFromGlobals_copies_migMatConst);
    RUN_TEST(test_allShapesConstant_default_state_is_true);
    RUN_TEST(test_allShapesConstant_false_when_a_pop_is_exponential);
    RUN_TEST(test_allShapesConstant_false_when_a_pop_is_linear);
    RUN_TEST(test_allShapesConstant_false_when_a_pair_is_exponential);
    RUN_TEST(test_allShapesConstant_ignores_pops_outside_npops);
    return UNITY_END();
}
#endif

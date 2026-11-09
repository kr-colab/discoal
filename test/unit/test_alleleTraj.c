#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqGeneral_reduces_to_detSweepFreq_constant(void) {
    /* For constant alpha, detSweepFreqGeneral(alpha, alpha*tau) == detSweepFreq(tau, alpha). */
    double alpha = 200.0;
    for (double tau = 0.0; tau < 0.1; tau += 0.005) {
        double x_old = detSweepFreq(tau, alpha);
        double x_new = detSweepFreqGeneral(alpha, alpha * tau);
        TEST_ASSERT_DOUBLE_WITHIN(1e-12, x_old, x_new);
    }
}

void test_detSweepFreqGeneral_handles_zero_A(void) {
    /* At A=0, x should be near 1 (start of sweep). */
    double x = detSweepFreqGeneral(200.0, 0.0);
    TEST_ASSERT_TRUE(x > 0.99);
    TEST_ASSERT_TRUE(x < 1.0);
}

void test_detSweepFreqGeneral_handles_large_A(void) {
    /* At A = -2*log(epsilon) = sweep duration, x ≈ epsilon. */
    double alpha = 200.0;
    double epsilon = 0.05 / alpha;
    double A_ts = -2.0 * log(epsilon);
    double x = detSweepFreqGeneral(alpha, A_ts);
    TEST_ASSERT_DOUBLE_WITHIN(epsilon * 0.01, epsilon, x);
}

void setall(long iseed1, long iseed2);

/* The neutral step's noise has the diffusion's variance x(1-x)dt (with a
 * uniform draw on (-1, 1) it had a third of it). */
void test_neutral_step_has_diffusion_variance(void) {
    const double x = 0.3, dt = 1e-4;
    const int n = 200000;
    double sum_sq = 0.0;
    setall(11, 12);
    for (int i = 0; i < n; i++) {
        double step = neutralStochasticOptimized(dt, x) - (x - x * dt);
        sum_sq += step * step;
    }
    double expected = x * (1.0 - x) * dt;
    TEST_ASSERT_DOUBLE_WITHIN(0.02 * expected, expected, sum_sq / n);
}

/* Paths that start at f0 and step back in time with
 * neutralStochasticBelowCeiling until they reach 0 stay below f0, and
 * their mean duration is that of a neutral allele drifting from 0 up to
 * its first passage at f0: 2[1 + (1 - f0) ln(1 - f0) / f0] units of 2N
 * generations (diffusion limit). discoal's default step, 1/(40N) with
 * N = 1000, gives about 1% more. */
static double mean_time_below_ceiling(double f0, int n_paths) {
    const double dt = 1.0 / (40.0 * 1000.0);
    double total = 0.0;
    for (int i = 0; i < n_paths; i++) {
        double x = f0, t = 0.0;
        while (x > 0.0) {
            x = neutralStochasticBelowCeiling(dt, x, f0);
            TEST_ASSERT_TRUE(x < f0);
            t += dt;
        }
        total += t;
    }
    return total / n_paths;
}

void test_neutral_phase_below_ceiling_has_first_passage_duration(void) {
    const double ceilings[] = {0.05, 0.2};
    setall(21, 22);
    for (int k = 0; k < 2; k++) {
        double f0 = ceilings[k];
        double expected = 2.0 * (1.0 + (1.0 - f0) * log(1.0 - f0) / f0);
        double mean = mean_time_below_ceiling(f0, 10000);
        TEST_ASSERT_DOUBLE_WITHIN(0.03 * expected, expected, mean);
    }
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqGeneral_reduces_to_detSweepFreq_constant);
    RUN_TEST(test_detSweepFreqGeneral_handles_zero_A);
    RUN_TEST(test_detSweepFreqGeneral_handles_large_A);
    RUN_TEST(test_neutral_step_has_diffusion_variance);
    RUN_TEST(test_neutral_phase_below_ceiling_has_first_passage_duration);
    return UNITY_END();
}
#endif

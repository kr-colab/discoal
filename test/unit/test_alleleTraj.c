#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqEuler_against_closed_form(void) {
    /* Existing test from Phase 2 — kept until detSweepFreqEuler is removed in Task 11 */
    double alpha = 200.0;
    double epsilon = 0.05 / alpha;
    double ts = -2.0 * log(epsilon) / alpha;
    int N_steps = 100000;
    double dt = ts / N_steps;
    double x_euler = detSweepFreq(0.0, alpha);
    for (int i = 0; i < N_steps; i++) {
        x_euler = detSweepFreqEuler(x_euler, dt, alpha);
    }
    double x_closed = detSweepFreq(ts, alpha);
    TEST_ASSERT_DOUBLE_WITHIN(1e-3 * x_closed, x_closed, x_euler);
}

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

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqEuler_against_closed_form);
    RUN_TEST(test_detSweepFreqGeneral_reduces_to_detSweepFreq_constant);
    RUN_TEST(test_detSweepFreqGeneral_handles_zero_A);
    RUN_TEST(test_detSweepFreqGeneral_handles_large_A);
    return UNITY_END();
}
#endif

#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqEuler_against_closed_form(void) {
    /* Backward-time integration: start at ttau=0 (x near 1 just after the
     * sweep finishes fixing), walk forward in backward-time toward ttau=ts
     * (x near epsilon, sweep starting). Both detSweepFreq and
     * detSweepFreqEuler should track the same trajectory. */
    double alpha = 200.0;
    double epsilon = 0.05 / alpha;
    double ts = -2.0 * log(epsilon) / alpha;
    int N_steps = 100000;
    double dt = ts / N_steps;
    /* Start at the present (ttau=0); detSweepFreq returns ~1 here. */
    double x_euler = detSweepFreq(0.0, alpha);
    for (int i = 0; i < N_steps; i++) {
        x_euler = detSweepFreqEuler(x_euler, dt, alpha);
    }
    /* End at ttau=ts; detSweepFreq returns ~epsilon. */
    double x_closed = detSweepFreq(ts, alpha);
    /* Euler with 100k steps should match closed form to 1e-3 relative. */
    TEST_ASSERT_DOUBLE_WITHIN(1e-3 * x_closed, x_closed, x_euler);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqEuler_against_closed_form);
    return UNITY_END();
}
#endif

#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqEuler_against_closed_form(void) {
    /* Run both methods over the entire sweep and compare endpoints.
     *
     * detSweepFreq(t, alpha) is the closed-form solution evaluated at
     * absolute time t in [0, ts]. By construction it runs backward in
     * time relative to forward sweep time: it starts near 1 (fixation)
     * at t = 0 and decreases to epsilon at t = ts.
     *
     * detSweepFreqEuler advances dx/dt = +alpha * x * (1-x), i.e. the
     * forward-time logistic growth from low frequency to fixation.
     *
     * To verify agreement, we initialize the Euler trajectory at the
     * value of detSweepFreq at t = ts (i.e., x = epsilon) and integrate
     * forward by N steps of size dt = ts / N. The final value should
     * equal detSweepFreq evaluated at t = 0 (i.e., the value at the
     * other end of the sweep, near 1). */
    double alpha = 200.0;
    /* closed-form sweep duration */
    double epsilon = 0.05 / alpha;
    double ts = -2.0 * log(epsilon) / alpha;
    /* fine Euler grid */
    int N_steps = 100000;
    double dt = ts / N_steps;
    /* start Euler at detSweepFreq(ts, alpha) == epsilon */
    double x_euler = detSweepFreq(ts, alpha);
    for (int i = 0; i < N_steps; i++) {
        x_euler = detSweepFreqEuler(x_euler, dt, alpha);
    }
    /* compare endpoint to detSweepFreq(0, alpha) (near fixation) */
    double x_closed = detSweepFreq(0.0, alpha);
    /* Euler with 100k steps should match closed form to 1e-3 relative */
    TEST_ASSERT_DOUBLE_WITHIN(1e-3 * x_closed, x_closed, x_euler);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqEuler_against_closed_form);
    return UNITY_END();
}
#endif

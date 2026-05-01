#include "unity.h"
#include "discoal.h"
#include "shapes.h"
#include <math.h>
#include <stdlib.h>

extern long seed1, seed2;
extern void setall(long iseed1, long iseed2);

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    /* Initialize the RNG so ranf() is usable from KS tests later */
    seed1 = 12345;
    seed2 = 67890;
    setall(seed1, seed2);

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
    return UNITY_END();
}
#endif

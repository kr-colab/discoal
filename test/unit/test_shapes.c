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

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_sizeAt_constant_returns_anchor);
    RUN_TEST(test_sizeAt_constant_per_population);
    return UNITY_END();
}
#endif

#include "unity.h"
#include "fenwick.h"
#include <stdlib.h>

static Fenwick f;

void setUp(void) { fenwickInit(&f); }
void tearDown(void) { fenwickFree(&f); }

/* fenwickFind against a linear scan of the weights, for every position. */
static void checkAllPositions(const int64_t *w, int n) {
    int64_t total = 0;
    for (int i = 0; i < n; i++) total += w[i];
    TEST_ASSERT_EQUAL_INT64(total, fenwickTotal(&f));
    int64_t before = 0;
    for (int i = 0; i < n; i++) {
        for (int64_t p = before; p < before + w[i]; p++) {
            int64_t offset = -1;
            TEST_ASSERT_EQUAL_INT(i, fenwickFind(&f, p, &offset));
            TEST_ASSERT_EQUAL_INT64(p - before, offset);
        }
        before += w[i];
    }
}

void test_find_matches_linear_scan_with_zero_weights(void) {
    int64_t w[37];
    for (int i = 0; i < 37; i++) {
        w[i] = (i % 3 == 0) ? 0 : (i % 7) + 1;  /* zeros must never be drawn */
        fenwickSet(&f, i, w[i]);
    }
    checkAllPositions(w, 37);
}

void test_updates_and_growth(void) {
    int64_t w[300] = {0};
    srand(7);
    for (int round = 0; round < 2000; round++) {
        int i = rand() % 300;          /* indices beyond the capacity grow it */
        w[i] = rand() % 50;
        fenwickSet(&f, i, w[i]);
    }
    checkAllPositions(w, 300);
    for (int i = 0; i < 300; i++) TEST_ASSERT_EQUAL_INT64(w[i], fenwickGet(&f, i));
}

void test_clear(void) {
    fenwickSet(&f, 3, 10);
    fenwickSet(&f, 40, 5);
    fenwickClear(&f);
    TEST_ASSERT_EQUAL_INT64(0, fenwickTotal(&f));
    fenwickSet(&f, 1, 4);
    int64_t offset;
    TEST_ASSERT_EQUAL_INT(1, fenwickFind(&f, 3, &offset));
    TEST_ASSERT_EQUAL_INT64(3, offset);
}

void test_large_weights(void) {
    /* Sums above 2^31 (e.g. many lineages spanning 10^8 sites). */
    fenwickSet(&f, 0, 3000000000LL);
    fenwickSet(&f, 1, 5000000000LL);
    int64_t offset;
    TEST_ASSERT_EQUAL_INT(1, fenwickFind(&f, 4000000000LL, &offset));
    TEST_ASSERT_EQUAL_INT64(1000000000LL, offset);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_find_matches_linear_scan_with_zero_weights);
    RUN_TEST(test_updates_and_growth);
    RUN_TEST(test_clear);
    RUN_TEST(test_large_weights);
    return UNITY_END();
}

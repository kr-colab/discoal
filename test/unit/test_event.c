#include "unity.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include <string.h>

// Test fixtures
event testEvent;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Initialize test event
    testEvent.time = 0.0;
    testEvent.popnSize = 0.0;
    testEvent.type = 'N';
    testEvent.popID = 0;
    testEvent.popID2 = 0;
    testEvent.popID3 = 0;
    testEvent.lineageNumber = 0;
    testEvent.admixProp = 0.0;
}

void tearDown(void) {
    // Nothing to clean up for this struct
}
#endif

void test_event_initialization(void) {
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, testEvent.time);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, testEvent.popnSize);
    TEST_ASSERT_EQUAL('N', testEvent.type);
    TEST_ASSERT_EQUAL(0, testEvent.popID);
    TEST_ASSERT_EQUAL(0, testEvent.popID2);
    TEST_ASSERT_EQUAL(0, testEvent.popID3);
    TEST_ASSERT_EQUAL(0, testEvent.lineageNumber);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, testEvent.admixProp);
}

void test_event_set_properties(void) {
    testEvent.time = 12.5;
    testEvent.popnSize = 1000.0;
    testEvent.type = 'A';
    testEvent.popID = 1;
    testEvent.popID2 = 2;
    testEvent.popID3 = 3;
    testEvent.lineageNumber = 5;
    testEvent.admixProp = 0.42;

    TEST_ASSERT_FLOAT_WITHIN(0.0001, 12.5, testEvent.time);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 1000.0, testEvent.popnSize);
    TEST_ASSERT_EQUAL('A', testEvent.type);
    TEST_ASSERT_EQUAL(1, testEvent.popID);
    TEST_ASSERT_EQUAL(2, testEvent.popID2);
    TEST_ASSERT_EQUAL(3, testEvent.popID3);
    TEST_ASSERT_EQUAL(5, testEvent.lineageNumber);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.42, testEvent.admixProp);
}

static event makeEvent(double time, char type, int popID) {
    event e;
    memset(&e, 0, sizeof(e));
    e.time = time;
    e.type = type;
    e.popID = popID;
    return e;
}

/* A growth event anchors at the size in force when it runs, so at a tied
 * time it must come after a size change, whatever order it was given in. */
void test_sort_events_puts_growth_after_size_change_at_same_time(void) {
    event e[4] = {
        makeEvent(0.5, 'g', 0),
        makeEvent(0.5, 'n', 0),
        makeEvent(0.2, 'l', 0),
        makeEvent(0.2, 'n', 0),
    };
    sortEventArray(e, 4);
    TEST_ASSERT_EQUAL('n', e[0].type);
    TEST_ASSERT_EQUAL('l', e[1].type);
    TEST_ASSERT_EQUAL('n', e[2].type);
    TEST_ASSERT_EQUAL('g', e[3].type);
    TEST_ASSERT_FLOAT_WITHIN(1e-12, 0.2, e[0].time);
    TEST_ASSERT_FLOAT_WITHIN(1e-12, 0.5, e[3].time);
}

/* Other tied events keep the order they were given in. */
void test_sort_events_is_stable_for_other_ties(void) {
    event e[4] = {
        makeEvent(1.0, 'p', 1),
        makeEvent(0.3, 'n', 0),
        makeEvent(1.0, 'm', 2),
        makeEvent(1.0, 'n', 3),
    };
    sortEventArray(e, 4);
    TEST_ASSERT_EQUAL('n', e[0].type);
    TEST_ASSERT_EQUAL(1, e[1].popID);
    TEST_ASSERT_EQUAL(2, e[2].popID);
    TEST_ASSERT_EQUAL(3, e[3].popID);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_event_initialization);
    RUN_TEST(test_event_set_properties);
    RUN_TEST(test_sort_events_puts_growth_after_size_change_at_same_time);
    RUN_TEST(test_sort_events_is_stable_for_other_ties);
    
    return UNITY_END();
}
#endif 
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

/* sortEventArray orders by time and keeps the given order for ties, as in
 * ms: for tied -eg and -en, the event given last decides whether the
 * population grows. */
void test_sort_events_is_stable(void) {
    event e[6] = {
        makeEvent(0.5, 'g', 0),
        makeEvent(0.5, 'n', 1),
        makeEvent(0.2, 'n', 2),
        makeEvent(0.2, 'l', 3),
        makeEvent(1.0, 'p', 4),
        makeEvent(0.5, 'm', 5),
    };
    sortEventArray(e, 6);
    int expected[6] = {2, 3, 0, 1, 5, 4};
    for (int i = 0; i < 6; i++) {
        TEST_ASSERT_EQUAL(expected[i], e[i].popID);
    }
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_event_initialization);
    RUN_TEST(test_event_set_properties);
    RUN_TEST(test_sort_events_is_stable);
    
    return UNITY_END();
}
#endif 
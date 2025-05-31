#include "unity.h"
#include "discoal.h"
#include <string.h>

// Test fixtures
event testEvent;

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

void test_event_initialization(void) {
    TEST_ASSERT_FLOAT_EQUAL(0.0, testEvent.time, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testEvent.popnSize, 0.0001);
    TEST_ASSERT_EQUAL('N', testEvent.type);
    TEST_ASSERT_EQUAL(0, testEvent.popID);
    TEST_ASSERT_EQUAL(0, testEvent.popID2);
    TEST_ASSERT_EQUAL(0, testEvent.popID3);
    TEST_ASSERT_EQUAL(0, testEvent.lineageNumber);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testEvent.admixProp, 0.0001);
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

    TEST_ASSERT_FLOAT_EQUAL(12.5, testEvent.time, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(1000.0, testEvent.popnSize, 0.0001);
    TEST_ASSERT_EQUAL('A', testEvent.type);
    TEST_ASSERT_EQUAL(1, testEvent.popID);
    TEST_ASSERT_EQUAL(2, testEvent.popID2);
    TEST_ASSERT_EQUAL(3, testEvent.popID3);
    TEST_ASSERT_EQUAL(5, testEvent.lineageNumber);
    TEST_ASSERT_FLOAT_EQUAL(0.42, testEvent.admixProp, 0.0001);
}

int main(void) {
    printf("\nRunning Event Struct Tests\n");
    printf("------------------------\n");
    RUN_TEST(test_event_initialization);
    RUN_TEST(test_event_set_properties);
    printf("\nAll event struct tests passed!\n");
    return 0;
} 
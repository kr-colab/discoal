#include "unity.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include <stdlib.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Minimal setup
    sampleSize = 10;
    nSites = 100;
}

void tearDown(void) {
    // Minimal cleanup
}
#endif

void test_basicNodeCreation(void) {
    // Test basic node allocation and initialization
    rootedNode* node = (rootedNode*)malloc(sizeof(rootedNode));
    TEST_ASSERT_NOT_NULL(node);
    
    // Initialize basic fields that actually exist in the struct
    node->time = 1.0;
    node->population = 1;
    node->sweepPopn = -1;
    node->tskit_node_id = TSK_NULL;
    
    // Test that we can access the basic fields
    TEST_ASSERT_EQUAL_FLOAT(1.0, node->time);
    TEST_ASSERT_EQUAL(1, node->population);
    TEST_ASSERT_EQUAL(-1, node->sweepPopn);
    TEST_ASSERT_EQUAL(TSK_NULL, node->tskit_node_id);
    
    free(node);
}

void test_nodeStateTracking(void) {
    // Test node state tracking fields
    rootedNode* node = (rootedNode*)malloc(sizeof(rootedNode));
    TEST_ASSERT_NOT_NULL(node);
    
    // Initialize state tracking fields
    node->parentsRecorded = 0;
    node->isFullyRecorded = 0;
    node->inActiveSet = 1;
    node->carriesSweepMutation = 0;
    
    // Test state tracking
    TEST_ASSERT_EQUAL(0, node->parentsRecorded);
    TEST_ASSERT_EQUAL(0, node->isFullyRecorded);
    TEST_ASSERT_EQUAL(1, node->inActiveSet);
    TEST_ASSERT_EQUAL(0, node->carriesSweepMutation);
    
    free(node);
}

void test_ancestryFields(void) {
    // Test ancestry-related fields
    rootedNode* node = (rootedNode*)malloc(sizeof(rootedNode));
    TEST_ASSERT_NOT_NULL(node);
    
    // Initialize ancestry fields
    node->nancSites = 10;
    node->lLim = 0;
    node->rLim = 100;
    node->ancestryRoot = NULL;
    node->popListIndex = -1;
    
    // Test ancestry fields
    TEST_ASSERT_EQUAL(10, node->nancSites);
    TEST_ASSERT_EQUAL(0, node->lLim);
    TEST_ASSERT_EQUAL(100, node->rLim);
    TEST_ASSERT_NULL(node->ancestryRoot);
    TEST_ASSERT_EQUAL(-1, node->popListIndex);
    
    free(node);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_basicNodeCreation);
    RUN_TEST(test_nodeStateTracking);
    RUN_TEST(test_ancestryFields);
    
    return UNITY_END();
}
#endif 
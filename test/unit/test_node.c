#include "unity.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include <math.h>
#include <stdlib.h>

// Test fixtures
rootedNode* testNode;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Initialize test node
    testNode = (rootedNode*)malloc(sizeof(rootedNode));
    testNode->leftParent = NULL;
    testNode->rightParent = NULL;
    testNode->leftChild = NULL;
    testNode->rightChild = NULL;
    testNode->time = 0.0;
    testNode->branchLength = 0.0;
    testNode->nancSites = 0;
    testNode->lLim = 0;
    testNode->rLim = 0;
    testNode->id = 0;
    testNode->population = 0;
    testNode->sweepPopn = 0;
    testNode->ancestryRoot = NULL;
    testNode->tskit_node_id = TSK_NULL;
    testNode->parentsRecorded = 0;
    testNode->isFullyRecorded = 0;
    testNode->inActiveSet = 0;
    testNode->carriesSweepMutation = 0;
    testNode->popListIndex = -1;
}

void tearDown(void) {
    // Clean up test node
    free(testNode);
}
#endif

void test_node_initialization(void) {
    TEST_ASSERT_NOT_NULL(testNode);
    TEST_ASSERT_NULL(testNode->leftParent);
    TEST_ASSERT_NULL(testNode->rightParent);
    TEST_ASSERT_NULL(testNode->leftChild);
    TEST_ASSERT_NULL(testNode->rightChild);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, testNode->time);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, testNode->branchLength);
    TEST_ASSERT_EQUAL(0, testNode->nancSites);
    TEST_ASSERT_EQUAL(0, testNode->lLim);
    TEST_ASSERT_EQUAL(0, testNode->rLim);
    TEST_ASSERT_EQUAL(0, testNode->id);
    TEST_ASSERT_EQUAL(0, testNode->population);
    TEST_ASSERT_EQUAL(0, testNode->sweepPopn);
    TEST_ASSERT_NULL(testNode->ancestryRoot);
    TEST_ASSERT_EQUAL(TSK_NULL, testNode->tskit_node_id);
    TEST_ASSERT_EQUAL(0, testNode->parentsRecorded);
    TEST_ASSERT_EQUAL(0, testNode->isFullyRecorded);
    TEST_ASSERT_EQUAL(0, testNode->inActiveSet);
    TEST_ASSERT_EQUAL(0, testNode->carriesSweepMutation);
    TEST_ASSERT_EQUAL(-1, testNode->popListIndex);
}

void test_node_set_properties(void) {
    // Set properties
    testNode->time = 1.5;
    testNode->branchLength = 2.0;
    testNode->id = 1;
    testNode->population = 2;
    
    // Verify properties
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 1.5, testNode->time);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 2.0, testNode->branchLength);
    TEST_ASSERT_EQUAL(1, testNode->id);
    TEST_ASSERT_EQUAL(2, testNode->population);
}

void test_newRootedNode_creation(void) {
    double cTime = 3.14;
    int popn = 7;
    rootedNode* node = newRootedNode(cTime, popn);
    TEST_ASSERT_NOT_NULL(node);
    TEST_ASSERT_NULL(node->leftParent);
    TEST_ASSERT_NULL(node->rightParent);
    TEST_ASSERT_NULL(node->leftChild);
    TEST_ASSERT_NULL(node->rightChild);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, cTime, node->time);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.0, node->branchLength);
    TEST_ASSERT_EQUAL(popn, node->population);
    TEST_ASSERT_EQUAL(-1, node->sweepPopn);
    TEST_ASSERT_NULL(node->ancestryRoot);
    TEST_ASSERT_EQUAL(TSK_NULL, node->tskit_node_id);
    TEST_ASSERT_EQUAL(0, node->parentsRecorded);
    TEST_ASSERT_EQUAL(0, node->isFullyRecorded);
    TEST_ASSERT_EQUAL(0, node->inActiveSet);
    TEST_ASSERT_EQUAL(0, node->carriesSweepMutation);
    TEST_ASSERT_EQUAL(-1, node->popListIndex);
    free(node);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_node_initialization);
    RUN_TEST(test_node_set_properties);
    RUN_TEST(test_newRootedNode_creation);
    
    return UNITY_END();
}
#endif 
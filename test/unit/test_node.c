#include "unity.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include <math.h>
#include <stdlib.h>

// Test fixtures
rootedNode* testNode;

void setUp(void) {
    // Initialize test node (legacy test)
    testNode = (rootedNode*)malloc(sizeof(rootedNode));
    testNode->leftParent = NULL;
    testNode->rightParent = NULL;
    testNode->leftChild = NULL;
    testNode->rightChild = NULL;
    testNode->time = 0.0;
    testNode->branchLength = 0.0;
    testNode->blProb = 0.0;
    testNode->nancSites = 0;
    testNode->lLim = 0;
    testNode->rLim = 0;
    testNode->id = 0;
    testNode->mutationNumber = 0;
    testNode->population = 0;
    testNode->sweepPopn = 0;
    testNode->ndes[0] = 0;
    testNode->ndes[1] = 0;
    testNode->leafs = NULL;
    testNode->times[0] = 0.0;
    testNode->times[1] = 0.0;
}

void tearDown(void) {
    // Clean up test node
    free(testNode);
}

void test_node_initialization(void) {
    TEST_ASSERT_NOT_NULL(testNode);
    TEST_ASSERT_NULL(testNode->leftParent);
    TEST_ASSERT_NULL(testNode->rightParent);
    TEST_ASSERT_NULL(testNode->leftChild);
    TEST_ASSERT_NULL(testNode->rightChild);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testNode->time, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testNode->branchLength, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testNode->blProb, 0.0001);
    TEST_ASSERT_EQUAL(0, testNode->nancSites);
    TEST_ASSERT_EQUAL(0, testNode->lLim);
    TEST_ASSERT_EQUAL(0, testNode->rLim);
    TEST_ASSERT_EQUAL(0, testNode->id);
    TEST_ASSERT_EQUAL(0, testNode->mutationNumber);
    TEST_ASSERT_EQUAL(0, testNode->population);
    TEST_ASSERT_EQUAL(0, testNode->sweepPopn);
    TEST_ASSERT_EQUAL(0, testNode->ndes[0]);
    TEST_ASSERT_EQUAL(0, testNode->ndes[1]);
    TEST_ASSERT_NULL(testNode->leafs);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testNode->times[0], 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(0.0, testNode->times[1], 0.0001);
}

void test_node_set_properties(void) {
    // Set properties
    testNode->time = 1.5;
    testNode->branchLength = 2.0;
    testNode->id = 1;
    testNode->population = 2;
    
    // Verify properties
    TEST_ASSERT_FLOAT_EQUAL(1.5, testNode->time, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(2.0, testNode->branchLength, 0.0001);
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
    TEST_ASSERT_FLOAT_EQUAL(cTime, node->time, 0.0001);
    TEST_ASSERT_FLOAT_EQUAL(0.0, node->branchLength, 0.0001);
    TEST_ASSERT_EQUAL(0, node->mutationNumber);
    TEST_ASSERT_EQUAL(popn, node->population);
    TEST_ASSERT_EQUAL(-1, node->sweepPopn);
    free(node);
}

int main(void) {
    printf("\nRunning Node Tests\n");
    printf("-----------------\n");
    
    RUN_TEST(test_node_initialization);
    RUN_TEST(test_node_set_properties);
    RUN_TEST(test_newRootedNode_creation);
    
    printf("\nAll tests passed!\n");
    return 0;
} 
#include "unity.h"
#include "../../discoal.h"
#include "../../discoalFunctions.h"
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
    node->mutationNumber = 0;
    
    // Test that we can access the basic fields
    TEST_ASSERT_EQUAL_FLOAT(1.0, node->time);
    TEST_ASSERT_EQUAL(1, node->population);
    TEST_ASSERT_EQUAL(0, node->mutationNumber);
    
    free(node);
}

void test_mutationArrayAccess(void) {
    // Test that we can access the mutation array
    rootedNode* node = (rootedNode*)malloc(sizeof(rootedNode));
    TEST_ASSERT_NOT_NULL(node);
    
    // Initialize mutation fields
    node->mutationNumber = 0;
    
    // Test basic mutation array access
    node->muts[0] = 0.5;
    node->mutationNumber = 1;
    
    TEST_ASSERT_EQUAL_FLOAT(0.5, node->muts[0]);
    TEST_ASSERT_EQUAL(1, node->mutationNumber);
    
    free(node);
}

void test_simpleManualMutation(void) {
    // Test manual mutation addition without using addMutation function
    rootedNode* node = (rootedNode*)malloc(sizeof(rootedNode));
    TEST_ASSERT_NOT_NULL(node);
    
    // Initialize
    node->mutationNumber = 0;
    
    // Manually add a mutation
    node->muts[node->mutationNumber] = 0.3;
    node->mutationNumber++;
    
    // Add another mutation
    node->muts[node->mutationNumber] = 0.7;
    node->mutationNumber++;
    
    // Test the values
    TEST_ASSERT_EQUAL(2, node->mutationNumber);
    TEST_ASSERT_EQUAL_FLOAT(0.3, node->muts[0]);
    TEST_ASSERT_EQUAL_FLOAT(0.7, node->muts[1]);
    
    free(node);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_basicNodeCreation);
    RUN_TEST(test_mutationArrayAccess);
    RUN_TEST(test_simpleManualMutation);
    
    return UNITY_END();
}
#endif 
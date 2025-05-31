#include "unity.h"
#include "../../discoal.h"
#include "../../discoalFunctions.h"

void setUp(void) {
    // Initialize any test setup
    initialize();
}

void tearDown(void) {
    // Clean up after tests
    // Free any allocated memory
}

void test_newRootedNode(void) {
    // Test basic node creation
    rootedNode* node = newRootedNode(1.0, 1);
    TEST_ASSERT_NOT_NULL(node);
    TEST_ASSERT_EQUAL(1.0, node->time);
    TEST_ASSERT_EQUAL(1, node->population);
    TEST_ASSERT_NULL(node->leftParent);
    TEST_ASSERT_NULL(node->rightParent);
    TEST_ASSERT_NULL(node->leftChild);
    TEST_ASSERT_NULL(node->rightChild);
    TEST_ASSERT_EQUAL(0, node->mutationNumber);
    TEST_ASSERT_EQUAL(0, node->nancSites);
}

void test_addRemoveNode(void) {
    // Test node list management
    rootedNode* node = newRootedNode(1.0, 1);
    
    // Test adding node
    addNode(node);
    TEST_ASSERT_EQUAL(1, nodePopnSize(1));
    
    // Test removing node
    removeNode(node);
    TEST_ASSERT_EQUAL(0, nodePopnSize(1));
}

void test_pickNodePopn(void) {
    // Test node selection from population
    rootedNode* node1 = newRootedNode(1.0, 1);
    rootedNode* node2 = newRootedNode(2.0, 1);
    
    addNode(node1);
    addNode(node2);
    
    // Test picking node from population
    rootedNode* picked = pickNodePopn(1);
    TEST_ASSERT_NOT_NULL(picked);
    TEST_ASSERT_EQUAL(1, picked->population);
    
    // Cleanup
    removeNode(node1);
    removeNode(node2);
}

void test_nodePopnSize(void) {
    // Test population size calculations
    TEST_ASSERT_EQUAL(0, nodePopnSize(1));  // Empty population
    
    rootedNode* node1 = newRootedNode(1.0, 1);
    rootedNode* node2 = newRootedNode(2.0, 1);
    
    addNode(node1);
    TEST_ASSERT_EQUAL(1, nodePopnSize(1));  // One node
    
    addNode(node2);
    TEST_ASSERT_EQUAL(2, nodePopnSize(1));  // Two nodes
    
    // Cleanup
    removeNode(node1);
    removeNode(node2);
}

int main(void) {
    printf("\nRunning Node Operations Tests\n");
    printf("----------------------------\n");
    
    RUN_TEST(test_newRootedNode);
    RUN_TEST(test_addRemoveNode);
    RUN_TEST(test_pickNodePopn);
    RUN_TEST(test_nodePopnSize);
    
    printf("\nAll tests passed!\n");
    return 0;
} 
#include "unity.h"
#include "../../discoal.h"
#include "../../discoalFunctions.h"
#include "../../ranlib.h"
#include <stdlib.h>
#include <string.h>

// Test fixtures
rootedNode *testNode;
int originalBreakPointsCapacity;
int originalNodesCapacity;
int originalAllNodesCapacity;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Initialize npops to avoid issues
    npops = 1;
    for (int i = 0; i < MAXPOPS; i++) {
        popnSizes[i] = 0;
        sweepPopnSizes[i] = 0;
    }
    
    // Save original values
    originalBreakPointsCapacity = breakPointsCapacity;
    originalNodesCapacity = nodesCapacity;
    originalAllNodesCapacity = allNodesCapacity;
    
    // Reset globals
    breakPoints = NULL;
    breakPointsCapacity = 0;
    breakNumber = 0;
    nodes = NULL;
    nodesCapacity = 0;
    allNodes = NULL;
    allNodesCapacity = 0;
    alleleNumber = 0;
    totNodeNumber = 0;
    
    testNode = NULL;
}

void tearDown(void) {
    // Clean up any allocations
    cleanupBreakPoints();
    cleanupNodeArrays();
    
    if (testNode) {
        cleanupMuts(testNode);
        free(testNode);
        testNode = NULL;
    }
    
    // Restore original values
    breakPointsCapacity = originalBreakPointsCapacity;
    nodesCapacity = originalNodesCapacity;
    allNodesCapacity = originalAllNodesCapacity;
}
#endif

// Test breakPoints initialization
void test_initializeBreakPoints_basic(void) {
    initializeBreakPoints();
    
    TEST_ASSERT_NOT_NULL(breakPoints);
    TEST_ASSERT_EQUAL(1000, breakPointsCapacity);  // INITIAL_BREAKPOINTS_CAPACITY
    TEST_ASSERT_EQUAL(0, breakNumber);
    TEST_ASSERT_EQUAL(666, breakPoints[0]);  // Check marker value
}

// Test breakPoints cleanup
void test_cleanupBreakPoints_basic(void) {
    initializeBreakPoints();
    TEST_ASSERT_NOT_NULL(breakPoints);
    
    cleanupBreakPoints();
    
    TEST_ASSERT_NULL(breakPoints);
    TEST_ASSERT_EQUAL(0, breakPointsCapacity);
    TEST_ASSERT_EQUAL(0, breakNumber);
}

// Test breakPoints capacity growth
void test_ensureBreakPointsCapacity_growth(void) {
    initializeBreakPoints();
    int initialCapacity = breakPointsCapacity;
    
    // Fill to capacity
    for (int i = 0; i < initialCapacity; i++) {
        addBreakPoint(i * 10);
    }
    
    TEST_ASSERT_EQUAL(initialCapacity, breakNumber);
    
    // Add one more to trigger growth
    addBreakPoint(9999);
    
    TEST_ASSERT_EQUAL(initialCapacity * 2, breakPointsCapacity);
    TEST_ASSERT_EQUAL(initialCapacity + 1, breakNumber);
    
    // Verify all values are preserved
    for (int i = 0; i < initialCapacity; i++) {
        TEST_ASSERT_EQUAL(i * 10, breakPoints[i]);
    }
    TEST_ASSERT_EQUAL(9999, breakPoints[initialCapacity]);
}

// Test addBreakPoint
void test_addBreakPoint_multiple(void) {
    initializeBreakPoints();
    
    addBreakPoint(100);
    addBreakPoint(200);
    addBreakPoint(300);
    
    TEST_ASSERT_EQUAL(3, breakNumber);
    TEST_ASSERT_EQUAL(100, breakPoints[0]);
    TEST_ASSERT_EQUAL(200, breakPoints[1]);
    TEST_ASSERT_EQUAL(300, breakPoints[2]);
}

// Test node arrays initialization
void test_initializeNodeArrays_basic(void) {
    initializeNodeArrays();
    
    TEST_ASSERT_NOT_NULL(nodes);
    TEST_ASSERT_NOT_NULL(allNodes);
    TEST_ASSERT_EQUAL(1000, nodesCapacity);  // Initial capacity
    TEST_ASSERT_EQUAL(1000, allNodesCapacity);
}

// Test node arrays cleanup
void test_cleanupNodeArrays_basic(void) {
    initializeNodeArrays();
    TEST_ASSERT_NOT_NULL(nodes);
    TEST_ASSERT_NOT_NULL(allNodes);
    
    cleanupNodeArrays();
    
    TEST_ASSERT_NULL(nodes);
    TEST_ASSERT_NULL(allNodes);
    TEST_ASSERT_EQUAL(0, nodesCapacity);
    TEST_ASSERT_EQUAL(0, allNodesCapacity);
}

// Test nodes array capacity growth
void test_ensureNodesCapacity_growth(void) {
    initializeNodeArrays();
    int initialCapacity = nodesCapacity;
    
    // Request size that requires growth
    ensureNodesCapacity(initialCapacity + 100);
    
    // Should double capacity
    TEST_ASSERT_EQUAL(initialCapacity * 2, nodesCapacity);
    TEST_ASSERT_NOT_NULL(nodes);
}

// Test allNodes array capacity growth
void test_ensureAllNodesCapacity_growth(void) {
    initializeNodeArrays();
    int initialCapacity = allNodesCapacity;
    
    // Request size that requires multiple doublings
    ensureAllNodesCapacity(initialCapacity * 3);
    
    // Should grow to at least 4x initial capacity
    TEST_ASSERT_TRUE(allNodesCapacity >= initialCapacity * 4);
    TEST_ASSERT_NOT_NULL(allNodes);
}

// Test adding nodes
void test_addNode_with_growth(void) {
    initializeNodeArrays();
    
    // Initialize population sizes
    for (int i = 0; i < MAXPOPS; i++) {
        popnSizes[i] = 0;
        sweepPopnSizes[i] = 0;
    }
    
    // Add many nodes to test growth
    for (int i = 0; i < 1500; i++) {
        rootedNode *node = malloc(sizeof(rootedNode));
        node->population = i % 3;  // Distribute across 3 populations
        node->sweepPopn = 0;
        addNode(node);
    }
    
    TEST_ASSERT_EQUAL(1500, alleleNumber);
    TEST_ASSERT_EQUAL(1500, totNodeNumber);
    TEST_ASSERT_TRUE(nodesCapacity >= 1500);
    TEST_ASSERT_TRUE(allNodesCapacity >= 1500);
    
    // Check population counts
    TEST_ASSERT_EQUAL(500, popnSizes[0]);
    TEST_ASSERT_EQUAL(500, popnSizes[1]);
    TEST_ASSERT_EQUAL(500, popnSizes[2]);
    
    // Clean up nodes
    for (int i = 0; i < 1500; i++) {
        free(allNodes[i]);
    }
}

// Test mutations initialization
void test_initializeMuts_basic(void) {
    testNode = malloc(sizeof(rootedNode));
    testNode->muts = NULL;
    testNode->mutsCapacity = 0;
    
    initializeMuts(testNode, 20);
    
    TEST_ASSERT_NOT_NULL(testNode->muts);
    TEST_ASSERT_EQUAL(20, testNode->mutsCapacity);
}

// Test mutations with zero capacity
void test_initializeMuts_zero_capacity(void) {
    testNode = malloc(sizeof(rootedNode));
    testNode->muts = NULL;
    testNode->mutsCapacity = 0;
    
    initializeMuts(testNode, 0);
    
    TEST_ASSERT_NOT_NULL(testNode->muts);
    TEST_ASSERT_EQUAL(10, testNode->mutsCapacity);  // Default capacity
}

// Test mutations cleanup
void test_cleanupMuts_basic(void) {
    testNode = malloc(sizeof(rootedNode));
    initializeMuts(testNode, 20);
    TEST_ASSERT_NOT_NULL(testNode->muts);
    
    cleanupMuts(testNode);
    
    TEST_ASSERT_NULL(testNode->muts);
    TEST_ASSERT_EQUAL(0, testNode->mutsCapacity);
}

// Test mutations capacity growth
void test_ensureMutsCapacity_growth(void) {
    testNode = malloc(sizeof(rootedNode));
    initializeMuts(testNode, 10);
    
    // Add mutations to test capacity
    testNode->muts[0] = 0.1;
    testNode->muts[1] = 0.2;
    
    // Request more capacity
    ensureMutsCapacity(testNode, 25);
    
    // Should grow to at least 32 (power of 2)
    TEST_ASSERT_TRUE(testNode->mutsCapacity >= 32);
    
    // Original values should be preserved
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.1, testNode->muts[0]);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 0.2, testNode->muts[1]);
}

// Test mutation management with many mutations
void test_mutations_stress_test(void) {
    testNode = malloc(sizeof(rootedNode));
    initializeMuts(testNode, 1);  // Start very small
    
    // Add many mutations
    for (int i = 0; i < 1000; i++) {
        ensureMutsCapacity(testNode, i + 1);
        testNode->muts[i] = (double)i / 1000.0;
    }
    
    // Capacity should have grown
    TEST_ASSERT_TRUE(testNode->mutsCapacity >= 1000);
    
    // Verify values
    for (int i = 0; i < 1000; i++) {
        TEST_ASSERT_FLOAT_WITHIN(0.0001, (double)i / 1000.0, testNode->muts[i]);
    }
}

// Test reinitializing already initialized structures
void test_reinitialize_breakPoints(void) {
    // First initialization
    initializeBreakPoints();
    int *oldPtr = breakPoints;
    addBreakPoint(123);
    
    // Reinitialize - should free old and allocate new
    initializeBreakPoints();
    
    TEST_ASSERT_NOT_NULL(breakPoints);
    TEST_ASSERT_EQUAL(0, breakNumber);
    TEST_ASSERT_EQUAL(666, breakPoints[0]);  // Reset to marker
}

// Test cleanup on NULL pointers
void test_cleanup_null_safety(void) {
    // These should not crash
    cleanupBreakPoints();
    cleanupNodeArrays();
    // Note: cleanupMuts doesn't handle NULL node, so we skip it
}

// Test integrated memory management
void test_integrated_memory_usage(void) {
    // Initialize all systems
    initializeBreakPoints();
    initializeNodeArrays();
    
    // Create and add nodes with mutations
    for (int i = 0; i < 100; i++) {
        rootedNode *node = malloc(sizeof(rootedNode));
        node->population = 0;
        node->sweepPopn = 0;
        initializeMuts(node, 5);
        
        // Add some mutations
        for (int j = 0; j < 3; j++) {
            ensureMutsCapacity(node, j + 1);
            node->muts[j] = (double)(i * 10 + j) / 1000.0;
        }
        node->mutationNumber = 3;
        
        addNode(node);
    }
    
    // Add breakpoints
    for (int i = 0; i < 50; i++) {
        addBreakPoint(i * 100);
    }
    
    // Verify everything is allocated
    TEST_ASSERT_EQUAL(100, alleleNumber);
    TEST_ASSERT_EQUAL(50, breakNumber);
    
    // Clean up nodes
    for (int i = 0; i < 100; i++) {
        cleanupMuts(allNodes[i]);
        free(allNodes[i]);
    }
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    printf("Starting memory management tests...\n");
    UNITY_BEGIN();
    
    RUN_TEST(test_initializeBreakPoints_basic);
    RUN_TEST(test_cleanupBreakPoints_basic);
    RUN_TEST(test_ensureBreakPointsCapacity_growth);
    RUN_TEST(test_addBreakPoint_multiple);
    RUN_TEST(test_initializeNodeArrays_basic);
    RUN_TEST(test_cleanupNodeArrays_basic);
    RUN_TEST(test_ensureNodesCapacity_growth);
    RUN_TEST(test_ensureAllNodesCapacity_growth);
    RUN_TEST(test_addNode_with_growth);
    RUN_TEST(test_initializeMuts_basic);
    RUN_TEST(test_initializeMuts_zero_capacity);
    RUN_TEST(test_cleanupMuts_basic);
    RUN_TEST(test_ensureMutsCapacity_growth);
    RUN_TEST(test_mutations_stress_test);
    RUN_TEST(test_reinitialize_breakPoints);
    RUN_TEST(test_cleanup_null_safety);
    RUN_TEST(test_integrated_memory_usage);
    
    return UNITY_END();
}
#endif
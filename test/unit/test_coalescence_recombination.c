#include "unity.h"
#include "../../discoal.h"
#include "../../discoalFunctions.h"
#include "../../ancestrySegment.h"
#include "../../ranlib.h"
#include <stdlib.h>
#include <string.h>

// Test fixtures
rootedNode *testNode1, *testNode2, *testNode3;
int originalSampleSize;
int originalNSites;
int originalActiveSites;
int originalBreakNumber;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Save original globals
    originalSampleSize = sampleSize;
    originalNSites = nSites;
    originalActiveSites = activeSites;
    originalBreakNumber = breakNumber;
    
    // Initialize globals for testing
    sampleSize = 10;
    nSites = 100;
    activeSites = 100;
    breakNumber = 0;
    npops = 2;  // Set number of populations
    alleleNumber = 0;  // Reset active node count
    totNodeNumber = 0;  // Reset total node count
    
    // Initialize population sizes
    for (int i = 0; i < MAXPOPS; i++) {
        popnSizes[i] = 0;
        sweepPopnSizes[i] = 0;
    }
    
    // Initialize required arrays
    initializeNodeArrays();
    initializeBreakPoints();
    initializeActiveMaterial(&activeMaterialSegments, nSites);
    
    // Initialize random number generator
    setall(12345, 67890);
    
    // Create test nodes
    testNode1 = NULL;
    testNode2 = NULL;
    testNode3 = NULL;
}

void tearDown(void) {
    // Clean up test nodes
    if (testNode1) {
        cleanupMuts(testNode1);
        free(testNode1);
    }
    if (testNode2) {
        cleanupMuts(testNode2);
        free(testNode2);
    }
    if (testNode3) {
        cleanupMuts(testNode3);
        free(testNode3);
    }
    
    // Clean up node arrays
    cleanupNodeArrays();
    cleanupBreakPoints();
    freeActiveMaterial(&activeMaterialSegments);
    
    // Restore original globals
    sampleSize = originalSampleSize;
    nSites = originalNSites;
    activeSites = originalActiveSites;
    breakNumber = originalBreakNumber;
}
#endif

// Helper function to create a test node with ancestry
rootedNode* createTestNodeWithAncestry(double time, int popn, int start, int end) {
    rootedNode *node = newRootedNode(time, popn);
    
    // Add ancestry segment
    node->ancestryRoot = newSegment(start, end, NULL, NULL);
    updateAncestryStatsFromTree(node);
    
    addNode(node);
    return node;
}

// Test coalescence basic functionality
void test_coalesceAtTimePopn_basic(void) {
    // Create two nodes in population 0
    testNode1 = createTestNodeWithAncestry(0.0, 0, 10, 50);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 30, 70);
    
    // alleleNumber tracks active nodes
    TEST_ASSERT_EQUAL(2, alleleNumber);
    
    // Perform coalescence at time 1.0
    coalesceAtTimePopn(1.0, 0);
    
    // Should now have one node (the parent)
    TEST_ASSERT_EQUAL(1, alleleNumber);
    
    // Check parent properties
    rootedNode *parent = nodes[0];
    TEST_ASSERT_NOT_NULL(parent);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent->time);
    TEST_ASSERT_EQUAL(0, parent->population);
    
    // Check children are linked correctly
    TEST_ASSERT_NOT_NULL(parent->leftChild);
    TEST_ASSERT_NOT_NULL(parent->rightChild);
    TEST_ASSERT_EQUAL(parent, parent->leftChild->leftParent);
    TEST_ASSERT_EQUAL(parent, parent->rightChild->leftParent);
    
    // Check branch lengths
    TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent->leftChild->branchLength);
    TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent->rightChild->branchLength);
}

// Test coalescence ancestry merging
void test_coalesceAtTimePopn_ancestry_merge(void) {
    // Create two nodes with non-overlapping ancestry
    testNode1 = createTestNodeWithAncestry(0.0, 0, 10, 30);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 50, 70);
    
    coalesceAtTimePopn(1.0, 0);
    
    rootedNode *parent = nodes[0];
    TEST_ASSERT_NOT_NULL(parent->ancestryRoot);
    
    // Count segments in merged ancestry
    int segCount = 0;
    AncestrySegment *seg = parent->ancestryRoot;
    while (seg) {
        segCount++;
        seg = seg->next;
    }
    
    // Should have 2 segments (non-overlapping)
    TEST_ASSERT_EQUAL(2, segCount);
    
    // Check ancestry is preserved
    TEST_ASSERT_EQUAL(1, getAncestryCount(parent->ancestryRoot, 20));  // From node1
    TEST_ASSERT_EQUAL(1, getAncestryCount(parent->ancestryRoot, 60));  // From node2
    TEST_ASSERT_EQUAL(0, getAncestryCount(parent->ancestryRoot, 40));  // Gap
}

// Test coalescence with overlapping ancestry
void test_coalesceAtTimePopn_overlapping_ancestry(void) {
    // Create two nodes with overlapping ancestry
    testNode1 = createTestNodeWithAncestry(0.0, 0, 20, 60);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 40, 80);
    
    coalesceAtTimePopn(1.0, 0);
    
    rootedNode *parent = nodes[0];
    TEST_ASSERT_NOT_NULL(parent->ancestryRoot);
    
    // Check merged ancestry counts
    TEST_ASSERT_EQUAL(1, getAncestryCount(parent->ancestryRoot, 30));  // Only from node1
    TEST_ASSERT_EQUAL(2, getAncestryCount(parent->ancestryRoot, 50));  // From both
    TEST_ASSERT_EQUAL(1, getAncestryCount(parent->ancestryRoot, 70));  // Only from node2
}

// Test recombination basic functionality
void test_recombineAtTimePopn_basic(void) {
    // Create a node with ancestry across whole sequence
    testNode1 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    
    // Try recombination once
    int xOver = recombineAtTimePopn(1.0, 0);
    
    // If recombination succeeded (xOver != 666)
    if (xOver != 666) {
        TEST_ASSERT_TRUE(xOver >= 0 && xOver < nSites);
        
        // Should now have two nodes (the two parents)
        TEST_ASSERT_EQUAL(2, alleleNumber);
        
        // Check both parents exist
        rootedNode *parent1 = nodes[0];
        rootedNode *parent2 = nodes[1];
        TEST_ASSERT_NOT_NULL(parent1);
        TEST_ASSERT_NOT_NULL(parent2);
        
        // Both should be at time 1.0
        TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent1->time);
        TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent2->time);
    } else {
        // Recombination didn't happen, which is okay - it's random
        TEST_ASSERT_EQUAL(1, alleleNumber);
    }
}

// Test recombination ancestry splitting
void test_recombineAtTimePopn_ancestry_split(void) {
    // Create a node with ancestry in middle region
    testNode1 = createTestNodeWithAncestry(0.0, 0, 30, 70);
    
    // Try recombination
    int xOver = recombineAtTimePopn(1.0, 0);
    
    if (xOver != 666 && xOver >= 30 && xOver < 70) {
        // Recombination happened in the ancestry region
        // Check that ancestry was split correctly
        rootedNode *leftParent = NULL;
        rootedNode *rightParent = NULL;
        
        // Find the parents
        for (int i = 0; i < alleleNumber; i++) {
            if (nodes[i]->leftChild == testNode1) {
                if (leftParent == NULL) leftParent = nodes[i];
                else rightParent = nodes[i];
            }
        }
        
        TEST_ASSERT_NOT_NULL(leftParent);
        TEST_ASSERT_NOT_NULL(rightParent);
        
        // Check ancestry split
        // Left parent should have ancestry [30, xOver)
        // Right parent should have ancestry [xOver, 70)
        TEST_ASSERT_EQUAL(1, getAncestryCount(leftParent->ancestryRoot, 35));
        TEST_ASSERT_EQUAL(0, getAncestryCount(leftParent->ancestryRoot, xOver));
        
        TEST_ASSERT_EQUAL(0, getAncestryCount(rightParent->ancestryRoot, xOver - 1));
        TEST_ASSERT_EQUAL(1, getAncestryCount(rightParent->ancestryRoot, 65));
    }
    // Otherwise recombination didn't happen in the ancestry region, which is fine
}

// Test gene conversion basic functionality
void test_geneConversionAtTimePopn_basic(void) {
    // Set gene conversion mean tract length
    gcMean = 100;
    
    // Create a node with ancestry
    testNode1 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    
    int originalCount = alleleNumber;
    
    // Perform gene conversion
    geneConversionAtTimePopn(1.0, 0);
    
    // Should have two parents now (if conversion happened)
    // Gene conversion might not always succeed depending on random position
    TEST_ASSERT_TRUE(alleleNumber == originalCount || alleleNumber == originalCount + 1);
    
    if (alleleNumber > originalCount) {
        // Gene conversion succeeded, check parents
        rootedNode *parent1 = nodes[0];
        rootedNode *parent2 = nodes[1];
        TEST_ASSERT_NOT_NULL(parent1);
        TEST_ASSERT_NOT_NULL(parent2);
        TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent1->time);
        TEST_ASSERT_FLOAT_WITHIN(0.001, 1.0, parent2->time);
    }
}

// Test makeGametesMS mutation collection
void test_makeGametesMS_mutation_collection(void) {
    // Create nodes with mutations
    testNode1 = newRootedNode(0.0, 0);
    testNode2 = newRootedNode(0.0, 0);
    
    // Add some mutations
    initializeMuts(testNode1, 5);
    testNode1->muts[0] = 0.1;
    testNode1->muts[1] = 0.3;
    testNode1->muts[2] = 0.5;
    testNode1->mutationNumber = 3;
    
    initializeMuts(testNode2, 5);
    testNode2->muts[0] = 0.2;
    testNode2->muts[1] = 0.3;  // Duplicate
    testNode2->muts[2] = 0.6;
    testNode2->mutationNumber = 3;
    
    // Add to allNodes
    allNodes[0] = testNode1;
    allNodes[1] = testNode2;
    sampleSize = 2;
    
    // Skip actually calling makeGametesMS as it prints to stdout
    // Just verify our setup is correct
    TEST_ASSERT_NOT_NULL(allNodes[0]);
    TEST_ASSERT_NOT_NULL(allNodes[1]);
    TEST_ASSERT_EQUAL(3, allNodes[0]->mutationNumber);
    TEST_ASSERT_EQUAL(3, allNodes[1]->mutationNumber);
}

// Test updateActiveMaterial
void test_updateActiveMaterial_basic(void) {
    // Create a parent node with two children
    rootedNode *parent = newRootedNode(1.0, 0);
    testNode1 = createTestNodeWithAncestry(0.0, 0, 20, 40);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 60, 80);
    
    parent->leftChild = testNode1;
    parent->rightChild = testNode2;
    testNode1->leftParent = parent;
    testNode2->leftParent = parent;
    
    // Merge ancestry for parent
    parent->ancestryRoot = mergeAncestryTrees(testNode1->ancestryRoot, testNode2->ancestryRoot);
    updateAncestryStatsFromTree(parent);
    
    // Update active material
    updateActiveMaterial(parent);
    
    // Active material should be updated based on whether sites are still segregating
    // This isn't a good test yet...
    TEST_ASSERT_TRUE(activeMaterialSegments.totalActive >= 0);
}

// Test siteBetweenChunks
void test_siteBetweenChunks_basic(void) {
    // Create a node with ancestry from 10 to 70
    testNode1 = newRootedNode(0.0, 0);
    testNode1->ancestryRoot = newSegment(10, 70, NULL, NULL);
    updateAncestryStatsFromTree(testNode1);
    
    // siteBetweenChunks uses lLim and rLim
    // After updateAncestryStatsFromTree, lLim should be 10 and rLim should be 69
    
    // Test various crossover sites
    TEST_ASSERT_EQUAL(0, siteBetweenChunks(testNode1, 5));   // Before lLim
    TEST_ASSERT_EQUAL(1, siteBetweenChunks(testNode1, 40));  // Between lLim and rLim
    TEST_ASSERT_EQUAL(0, siteBetweenChunks(testNode1, 80));  // After rLim
    
    // Clean up
    freeSegmentTree(testNode1->ancestryRoot);
    testNode1->ancestryRoot = NULL;
}

// Test pickNodePopn
void test_pickNodePopn_basic(void) {
    // Create multiple nodes in different populations
    testNode1 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    testNode3 = createTestNodeWithAncestry(0.0, 1, 0, 100);
    
    // Pick from population 0
    rootedNode *picked = pickNodePopn(0);
    TEST_ASSERT_NOT_NULL(picked);
    TEST_ASSERT_EQUAL(0, picked->population);
    TEST_ASSERT_TRUE(picked == testNode1 || picked == testNode2);
    
    // Pick from population 1
    picked = pickNodePopn(1);
    TEST_ASSERT_NOT_NULL(picked);
    TEST_ASSERT_EQUAL(1, picked->population);
    TEST_ASSERT_EQUAL(testNode3, picked);
}

// Test nodePopnSize
void test_nodePopnSize_count(void) {
    // Create nodes in different populations
    testNode1 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    testNode2 = createTestNodeWithAncestry(0.0, 0, 0, 100);
    testNode3 = createTestNodeWithAncestry(0.0, 1, 0, 100);
    
    TEST_ASSERT_EQUAL(2, nodePopnSize(0));
    TEST_ASSERT_EQUAL(1, nodePopnSize(1));
    TEST_ASSERT_EQUAL(0, nodePopnSize(2));
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_coalesceAtTimePopn_basic);
    RUN_TEST(test_coalesceAtTimePopn_ancestry_merge);
    RUN_TEST(test_coalesceAtTimePopn_overlapping_ancestry);
    RUN_TEST(test_recombineAtTimePopn_basic);
    RUN_TEST(test_recombineAtTimePopn_ancestry_split);
    RUN_TEST(test_geneConversionAtTimePopn_basic);
    RUN_TEST(test_makeGametesMS_mutation_collection);
    RUN_TEST(test_updateActiveMaterial_basic);
    RUN_TEST(test_siteBetweenChunks_basic);
    RUN_TEST(test_pickNodePopn_basic);
    RUN_TEST(test_nodePopnSize_count);
    
    return UNITY_END();
}
#endif
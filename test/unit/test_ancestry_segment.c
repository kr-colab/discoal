#include "unity.h"
#include "../../ancestrySegment.h"
#include <stdlib.h>

// Test fixtures
AncestrySegment* testSegment;

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    testSegment = NULL;
}

void tearDown(void) {
    if (testSegment) {
        freeSegmentTree(testSegment);
        testSegment = NULL;
    }
}
#endif

// Test basic segment creation
void test_newSegment_creates_valid_segment(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    TEST_ASSERT_NOT_NULL(testSegment);
    TEST_ASSERT_EQUAL(10, testSegment->start);
    TEST_ASSERT_EQUAL(50, testSegment->end);
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    TEST_ASSERT_NULL(testSegment->left);
    TEST_ASSERT_NULL(testSegment->right);
    TEST_ASSERT_NULL(testSegment->avlTree);
}

// Test segment creation with invalid range
void test_newSegment_handles_invalid_range(void) {
    testSegment = newSegment(50, 10, NULL, NULL);  // start > end
    
    TEST_ASSERT_NOT_NULL(testSegment);
    // The implementation doesn't validate, so check what we get
    TEST_ASSERT_EQUAL(50, testSegment->start);
    TEST_ASSERT_EQUAL(10, testSegment->end);
}

// Test reference counting
void test_reference_counting_retain_release(void) {
    testSegment = newSegment(0, 100, NULL, NULL);
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    
    // Retain should increment
    retainSegment(testSegment);
    TEST_ASSERT_EQUAL(2, testSegment->refCount);
    
    // Release should decrement
    releaseSegment(testSegment);
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    
    // Final release should free (we'll check by setting to NULL after)
    releaseSegment(testSegment);
    testSegment = NULL; // Prevent double-free in tearDown
}

// Test shallow copy with reference counting
void test_shallowCopySegment_shares_references(void) {
    testSegment = newSegment(20, 80, NULL, NULL);
    testSegment->refCount = 1;
    
    AncestrySegment* copy = shallowCopySegment(testSegment);
    
    TEST_ASSERT_NOT_NULL(copy);
    TEST_ASSERT_EQUAL(testSegment, copy);  // Should be same pointer
    TEST_ASSERT_EQUAL(2, testSegment->refCount);  // Ref count should increase
    
    releaseSegment(copy);
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
}

// Test ancestry count queries
void test_getAncestryCount_single_segment(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    // Inside segment [10, 50)
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 20));
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 10));  // Left boundary (inclusive)
    TEST_ASSERT_EQUAL(0, getAncestryCount(testSegment, 50));  // Right boundary (exclusive)
    
    // Outside segment
    TEST_ASSERT_EQUAL(0, getAncestryCount(testSegment, 5));
    TEST_ASSERT_EQUAL(0, getAncestryCount(testSegment, 55));
}

// Test hasAncestry queries
void test_hasAncestry_queries(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    TEST_ASSERT_TRUE(hasAncestry(testSegment, 30));
    TEST_ASSERT_FALSE(hasAncestry(testSegment, 5));
    TEST_ASSERT_FALSE(hasAncestry(testSegment, 60));
    TEST_ASSERT_FALSE(hasAncestry(NULL, 30));  // NULL tree
}

// Test segment tree with multiple segments
void test_getAncestryCount_multiple_segments(void) {
    // Create a linked list of segments using 'next' pointer
    testSegment = newSegment(10, 30, NULL, NULL);
    testSegment->next = newSegment(20, 40, NULL, NULL);
    testSegment->next->next = newSegment(15, 25, NULL, NULL);
    
    // Test regions - getAncestryCount returns count from first matching segment
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 12));  // In first [10,30)
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 22));  // In first [10,30)
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 35));  // In second [20,40)
    TEST_ASSERT_EQUAL(0, getAncestryCount(testSegment, 50));  // In none
}

// Test deep copy functionality
void test_copySegmentTree_creates_independent_copy(void) {
    // Single segment with no next - should be shallow copy
    testSegment = newSegment(10, 30, NULL, NULL);
    
    AncestrySegment* copy = copySegmentTree(testSegment);
    
    TEST_ASSERT_NOT_NULL(copy);
    // Should be shallow copy with ref counting
    TEST_ASSERT_EQUAL(testSegment, copy);
    TEST_ASSERT_EQUAL(2, testSegment->refCount);
    
    freeSegmentTree(copy);  // Should decrement ref count
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    
    // Now test with next pointer - should be deep copy
    testSegment->next = newSegment(40, 60, NULL, NULL);
    copy = copySegmentTree(testSegment);
    
    TEST_ASSERT_NOT_NULL(copy);
    TEST_ASSERT_NOT_EQUAL(testSegment, copy);  // Deep copy
    TEST_ASSERT_EQUAL(testSegment->start, copy->start);
    TEST_ASSERT_EQUAL(testSegment->end, copy->end);
    
    freeSegmentTree(copy);  // Clean up copy
}

// Test NULL handling
void test_null_safety(void) {
    TEST_ASSERT_EQUAL(0, getAncestryCount(NULL, 10));
    TEST_ASSERT_FALSE(hasAncestry(NULL, 10));
    TEST_ASSERT_NULL(copySegmentTree(NULL));
    
    // These should not crash
    freeSegmentTree(NULL);
    retainSegment(NULL);
    releaseSegment(NULL);
}

// Test segment merging (basic case)
void test_mergeAncestryTrees_non_overlapping(void) {
    AncestrySegment* tree1 = newSegment(10, 30, NULL, NULL);
    AncestrySegment* tree2 = newSegment(40, 60, NULL, NULL);
    
    AncestrySegment* merged = mergeAncestryTrees(tree1, tree2);
    
    TEST_ASSERT_NOT_NULL(merged);
    TEST_ASSERT_EQUAL(1, getAncestryCount(merged, 20));  // From tree1
    TEST_ASSERT_EQUAL(1, getAncestryCount(merged, 50));  // From tree2
    TEST_ASSERT_EQUAL(0, getAncestryCount(merged, 35));  // Gap between
    
    freeSegmentTree(merged);
    freeSegmentTree(tree1);
    freeSegmentTree(tree2);
}

// Test splitting operations
void test_splitLeft_basic(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    AncestrySegment* leftTree = splitLeft(testSegment, 30);
    
    TEST_ASSERT_NOT_NULL(leftTree);
    TEST_ASSERT_EQUAL(10, leftTree->start);
    TEST_ASSERT_EQUAL(30, leftTree->end);  // Half-open interval [10, 30)
    
    freeSegmentTree(leftTree);
}

void test_splitRight_basic(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    AncestrySegment* rightTree = splitRight(testSegment, 30);
    
    TEST_ASSERT_NOT_NULL(rightTree);
    TEST_ASSERT_EQUAL(30, rightTree->start);  // Should start at breakpoint
    TEST_ASSERT_EQUAL(50, rightTree->end);
    
    freeSegmentTree(rightTree);
}

// Test edge cases for splitting
void test_split_edge_cases(void) {
    testSegment = newSegment(10, 50, NULL, NULL);
    
    // Split at boundaries
    AncestrySegment* leftEdge = splitLeft(testSegment, 10);
    TEST_ASSERT_NULL(leftEdge);  // Nothing to the left
    
    AncestrySegment* rightEdge = splitRight(testSegment, 51);
    TEST_ASSERT_NULL(rightEdge);  // Nothing to the right
    
    // Split outside range
    AncestrySegment* leftOut = splitLeft(testSegment, 5);
    TEST_ASSERT_NULL(leftOut);
    
    AncestrySegment* rightOut = splitRight(testSegment, 60);
    TEST_ASSERT_NULL(rightOut);  // Nothing to the right of 60 when segment is [10,50)
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_newSegment_creates_valid_segment);
    RUN_TEST(test_newSegment_handles_invalid_range);
    RUN_TEST(test_reference_counting_retain_release);
    RUN_TEST(test_shallowCopySegment_shares_references);
    RUN_TEST(test_getAncestryCount_single_segment);
    RUN_TEST(test_hasAncestry_queries);
    RUN_TEST(test_getAncestryCount_multiple_segments);
    RUN_TEST(test_copySegmentTree_creates_independent_copy);
    RUN_TEST(test_null_safety);
    RUN_TEST(test_mergeAncestryTrees_non_overlapping);
    RUN_TEST(test_splitLeft_basic);
    RUN_TEST(test_splitRight_basic);
    RUN_TEST(test_split_edge_cases);
    
    return UNITY_END();
}
#endif
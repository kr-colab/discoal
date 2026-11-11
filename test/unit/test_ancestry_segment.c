#include "unity.h"
#include "ancestrySegment.h"
#include "segmentPool.h"
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
    testSegment = newSegment(10, 50, TSK_NULL, NULL, NULL);
    
    TEST_ASSERT_NOT_NULL(testSegment);
    TEST_ASSERT_EQUAL(10, testSegment->start);
    TEST_ASSERT_EQUAL(50, testSegment->end);
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    TEST_ASSERT_NULL(testSegment->left);
    TEST_ASSERT_NULL(testSegment->right);
}

// Test segment creation with invalid range
void test_newSegment_handles_invalid_range(void) {
    testSegment = newSegment(50, 10, TSK_NULL, NULL, NULL);  // start > end
    
    TEST_ASSERT_NOT_NULL(testSegment);
    // The implementation doesn't validate, so check what we get
    TEST_ASSERT_EQUAL(50, testSegment->start);
    TEST_ASSERT_EQUAL(10, testSegment->end);
}

// Test reference counting
void test_reference_counting_retain_release(void) {
    testSegment = newSegment(0, 100, TSK_NULL, NULL, NULL);
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
    testSegment = newSegment(20, 80, TSK_NULL, NULL, NULL);
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
    testSegment = newSegment(10, 50, TSK_NULL, NULL, NULL);
    
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
    testSegment = newSegment(10, 50, TSK_NULL, NULL, NULL);
    
    TEST_ASSERT_TRUE(hasAncestry(testSegment, 30));
    TEST_ASSERT_FALSE(hasAncestry(testSegment, 5));
    TEST_ASSERT_FALSE(hasAncestry(testSegment, 60));
    TEST_ASSERT_FALSE(hasAncestry(NULL, 30));  // NULL tree
}

// Test segment tree with multiple segments
void test_getAncestryCount_multiple_segments(void) {
    // Create a linked list of segments using 'next' pointer
    testSegment = newSegment(10, 30, TSK_NULL, NULL, NULL);
    testSegment->next = newSegment(20, 40, TSK_NULL, NULL, NULL);
    testSegment->next->next = newSegment(15, 25, TSK_NULL, NULL, NULL);
    
    // Test regions - getAncestryCount returns count from first matching segment
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 12));  // In first [10,30)
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 22));  // In first [10,30)
    TEST_ASSERT_EQUAL(1, getAncestryCount(testSegment, 35));  // In second [20,40)
    TEST_ASSERT_EQUAL(0, getAncestryCount(testSegment, 50));  // In none
}

// Test deep copy functionality
void test_copySegmentTree_creates_independent_copy(void) {
    // Single segment with no next - should be shallow copy
    testSegment = newSegment(10, 30, TSK_NULL, NULL, NULL);
    
    AncestrySegment* copy = copySegmentTree(testSegment);
    
    TEST_ASSERT_NOT_NULL(copy);
    // Should be shallow copy with ref counting
    TEST_ASSERT_EQUAL(testSegment, copy);
    TEST_ASSERT_EQUAL(2, testSegment->refCount);
    
    freeSegmentTree(copy);  // Should decrement ref count
    TEST_ASSERT_EQUAL(1, testSegment->refCount);
    
    // Now test with next pointer - should be deep copy
    testSegment->next = newSegment(40, 60, TSK_NULL, NULL, NULL);
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
    AncestrySegment* tree1 = newSegment(10, 30, TSK_NULL, NULL, NULL);
    AncestrySegment* tree2 = newSegment(40, 60, TSK_NULL, NULL, NULL);
    
    AncestrySegment* merged = mergeAncestryTrees(tree1, tree2, TSK_NULL);
    
    TEST_ASSERT_NOT_NULL(merged);
    TEST_ASSERT_EQUAL(1, getAncestryCount(merged, 20));  // From tree1
    TEST_ASSERT_EQUAL(1, getAncestryCount(merged, 50));  // From tree2
    TEST_ASSERT_EQUAL(0, getAncestryCount(merged, 35));  // Gap between
    
    freeSegmentTree(merged);
    freeSegmentTree(tree1);
    freeSegmentTree(tree2);
}


/* A list [0,10) [20,30) [40,50) with counts 1, 2, 3 and tskit nodes 5, 6, 7. */
static AncestrySegment* threeSegments(void) {
    AncestrySegment *a = newSegment(0, 10, 5, NULL, NULL);
    AncestrySegment *b = newSegment(20, 30, 6, NULL, NULL);
    AncestrySegment *c = newSegment(40, 50, 7, NULL, NULL);
    a->count = 1; b->count = 2; c->count = 3;
    a->next = b; b->next = c;
    return a;
}

/* The pieces of a list inside [lo, hi), as new segments: the halves that
 * a split at a breakpoint must give. */
static AncestrySegment* copyPart(AncestrySegment *list, int lo, int hi) {
    AncestrySegment *head = NULL, *tail = NULL;
    for (AncestrySegment *s = list; s; s = s->next) {
        int start = s->start > lo ? s->start : lo;
        int end = s->end < hi ? s->end : hi;
        if (start >= end) continue;
        AncestrySegment *c = newSegment(start, end, s->tskit_node_id, NULL, NULL);
        c->count = s->count;
        if (tail) tail->next = c; else head = c;
        tail = c;
    }
    return head;
}

/* splitInPlace must give the pieces of the list left and right of the
 * breakpoint. */
static void assertSameList(AncestrySegment *expected, AncestrySegment *actual) {
    while (expected && actual) {
        TEST_ASSERT_EQUAL(expected->start, actual->start);
        TEST_ASSERT_EQUAL(expected->end, actual->end);
        TEST_ASSERT_EQUAL(expected->count, actual->count);
        TEST_ASSERT_EQUAL(expected->tskit_node_id, actual->tskit_node_id);
        expected = expected->next;
        actual = actual->next;
    }
    TEST_ASSERT_NULL(expected);
    TEST_ASSERT_NULL(actual);
}

/* The span that splitInPlace() reports equals a walk of the half. With a
 * sample size of 3 the segment [40,50) (count 3) has found its MRCA and
 * does not count. */
static void assertSpan(AncestrySegment *half, SegmentSpan span) {
    SegmentSpan walk = {0, 100, 0};
    for (AncestrySegment *s = half; s; s = s->next)
        spanAdd(&walk, s->start, s->end, s->count, 3);
    TEST_ASSERT_EQUAL(walk.nancSites, span.nancSites);
    TEST_ASSERT_EQUAL(walk.lLim, span.lLim);
    TEST_ASSERT_EQUAL(walk.rLim, span.rLim);
}

void test_splitInPlace_matches_copying_split(void) {
    /* Inside a segment, at a segment's start, in a gap, at a segment's end,
     * near the ends, and before or after all the material. */
    int breakpoints[] = {25, 20, 35, 30, 1, 49, 0, 50};
    for (int k = 0; k < 8; k++) {
        AncestrySegment *ref = threeSegments();
        AncestrySegment *wantLeft = copyPart(ref, 0, breakpoints[k]);
        AncestrySegment *wantRight = copyPart(ref, breakpoints[k], 100);
        AncestrySegment *list = threeSegments();
        AncestrySegment *left = NULL, *right = NULL;
        SegmentSpan lSpan = {0, 100, 0}, rSpan = {0, 100, 0};
        SegmentSpan whole = {0, 100, 0};
        for (AncestrySegment *s = list; s; s = s->next)
            spanAdd(&whole, s->start, s->end, s->count, 3);
        splitInPlace(list, breakpoints[k], 3, &whole, &left, &right, &lSpan, &rSpan);
        assertSameList(wantLeft, left);
        assertSameList(wantRight, right);
        assertSpan(left, lSpan);
        assertSpan(right, rSpan);
        freeSegmentTree(left);
        freeSegmentTree(right);
        freeSegmentTree(wantLeft);
        freeSegmentTree(wantRight);
        freeSegmentTree(ref);
    }
}


/* Free every segment and reset the pool after it has expanded: the free
 * list must hold every segment of the pool, not only those free when it
 * last grew. A pool started at 1000 segments (the size discoal uses for
 * small models) holds 20000, then 60000 after one more expansion; 50000
 * segments reach that expansion. Run under valgrind or AddressSanitizer,
 * a write past the end of the free list is reported. */
void test_segment_pool_free_list_holds_every_segment(void) {
    enum { COUNT = 50000 };
    AncestrySegment **segs = malloc(COUNT * sizeof(AncestrySegment *));
    TEST_ASSERT_NOT_NULL(segs);
    destroySegmentPool();
    initSegmentPool(1000);
    PoolStats before = getPoolStatistics();
    for (int i = 0; i < COUNT; i++) {
        segs[i] = allocSegmentFromPool();
        TEST_ASSERT_NOT_NULL(segs[i]);
    }
    PoolStats grown = getPoolStatistics();
    TEST_ASSERT_EQUAL(before.expansions + 1, grown.expansions);
    TEST_ASSERT_TRUE(grown.pool_capacity > COUNT);
    for (int i = 0; i < COUNT; i++)
        freeSegmentToPool(segs[i]);
    resetSegmentPool();
    PoolStats after = getPoolStatistics();
    TEST_ASSERT_EQUAL(0, after.current_in_use);
    /* Every segment can be allocated again. */
    for (int i = 0; i < COUNT; i++)
        segs[i] = allocSegmentFromPool();
    for (int i = 0; i < COUNT; i++)
        freeSegmentToPool(segs[i]);
    free(segs);
    destroySegmentPool();
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
    RUN_TEST(test_splitInPlace_matches_copying_split);
    RUN_TEST(test_segment_pool_free_list_holds_every_segment);
    
    return UNITY_END();
}
#endif
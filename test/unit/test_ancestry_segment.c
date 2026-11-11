#include "unity.h"
#include "ancestrySegment.h"
#include "segmentPool.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

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

/* ---- The merge against a reference ---- */

/* The merge as it was before it reused segments: it walks every interval
 * between the boundaries of both lists, copies each piece into a new
 * segment, and leaves the children's lists alone. The new merge must give
 * the same list, edges, fixed runs and span. */
static void refAppend(AncestrySegment **head, AncestrySegment **tail, int start,
                      int end, uint16_t count, tsk_id_t node) {
    if (*tail && (*tail)->end == start && (*tail)->count == count &&
        (*tail)->tskit_node_id == node) {
        (*tail)->end = end;
        return;
    }
    AncestrySegment *seg = newSegment(start, end, node, NULL, NULL);
    seg->count = count;
    if (*tail) (*tail)->next = seg; else *head = seg;
    *tail = seg;
}

static AncestrySegment* referenceMerge(AncestrySegment *l, AncestrySegment *r,
        tsk_id_t (*newParentNode)(void *ctx),
        void (*addEdge)(tsk_id_t parent, tsk_id_t child, double left, double right),
        int fullCount, void (*onFixed)(int start, int end), SegmentSpan *span) {
    AncestrySegment *result = NULL, *tail = NULL;
    tsk_id_t parent = TSK_NULL;
    int fixedStart = -1, fixedEnd = -1;
    int pos = INT_MAX;
    if (l) pos = l->start;
    if (r && r->start < pos) pos = r->start;
    while (l || r) {
        int nextPos = INT_MAX;
        if (l) {
            if (l->start > pos && l->start < nextPos) nextPos = l->start;
            if (l->end > pos && l->end < nextPos) nextPos = l->end;
        }
        if (r) {
            if (r->start > pos && r->start < nextPos) nextPos = r->start;
            if (r->end > pos && r->end < nextPos) nextPos = r->end;
        }
        if (nextPos == INT_MAX) break;
        AncestrySegment *lc = (l && pos >= l->start && pos < l->end) ? l : NULL;
        AncestrySegment *rc = (r && pos >= r->start && pos < r->end) ? r : NULL;
        if (lc && rc) {
            if (parent == TSK_NULL) parent = newParentNode(NULL);
            refAppend(&result, &tail, pos, nextPos, lc->count + rc->count, parent);
            spanAdd(span, pos, nextPos, lc->count + rc->count, fullCount);
            addEdge(parent, lc->tskit_node_id, pos, nextPos);
            addEdge(parent, rc->tskit_node_id, pos, nextPos);
            if (lc->count + rc->count == fullCount) {
                if (fixedEnd == pos) {
                    fixedEnd = nextPos;
                } else {
                    if (fixedStart >= 0) onFixed(fixedStart, fixedEnd);
                    fixedStart = pos;
                    fixedEnd = nextPos;
                }
            }
        } else if (lc || rc) {
            AncestrySegment *c = lc ? lc : rc;
            refAppend(&result, &tail, pos, nextPos, c->count, c->tskit_node_id);
            spanAdd(span, pos, nextPos, c->count, fullCount);
        }
        pos = nextPos;
        while (l && l->end <= pos) l = l->next;
        while (r && r->end <= pos) r = r->next;
    }
    if (fixedStart >= 0) onFixed(fixedStart, fixedEnd);
    return result;
}

/* A record of the callbacks: edges and fixed runs, in order. */
static int eventLog[128][4];
static int eventCount;
static tsk_id_t logParent(void *ctx) { (void) ctx; return 99; }
static void logEdge(tsk_id_t parent, tsk_id_t child, double left, double right) {
    int *e = eventLog[eventCount++];
    e[0] = (int) parent; e[1] = (int) child; e[2] = (int) left; e[3] = (int) right;
}
static void logFixed(int start, int end) {
    int *e = eventLog[eventCount++];
    e[0] = -1; e[1] = -1; e[2] = start; e[3] = end;
}

/* A random canonical list (adjacent segments differ in count or node),
 * with gaps, counts 1 to 3 and tskit nodes 1 to 3. */
static AncestrySegment* randomList(void) {
    AncestrySegment *head = NULL, *tail = NULL;
    int pos = rand() % 4;
    int n = rand() % 8;
    for (int k = 0; k < n; k++) {
        int len = 1 + rand() % 6;
        refAppend(&head, &tail, pos, pos + len, 1 + rand() % 3, 1 + rand() % 3);
        pos += len + (rand() % 3 == 0 ? 0 : rand() % 5);
    }
    return head;
}

/* The list without its pieces whose sites have found their MRCA
 * (count == fullCount), which the merge does not keep. A dropped piece
 * leaves a gap, so no neighbours can be joined across it. */
static AncestrySegment* withoutFixed(AncestrySegment *list, int fullCount) {
    AncestrySegment *head = NULL, *tail = NULL;
    while (list) {
        AncestrySegment *next = list->next;
        if (list->count == fullCount) {
            list->next = NULL;
            freeSegmentTree(list);
        } else {
            list->next = NULL;
            if (tail) tail->next = list; else head = list;
            tail = list;
        }
        list = next;
    }
    return head;
}

static SegmentSpan spanOf(AncestrySegment *list, int fullCount) {
    SegmentSpan span = {0, 1000, 0};
    for (; list; list = list->next)
        spanAdd(&span, list->start, list->end, list->count, fullCount);
    return span;
}

void test_merge_matches_reference(void) {
    int fullCount = 6;  /* two counts of 3 make a newly fixed piece */
    srand(11);
    for (int iter = 0; iter < 5000; iter++) {
        /* The merge uses up its lists, so it gets copies. */
        AncestrySegment *l1 = randomList(), *r1 = randomList();
        AncestrySegment *l2 = copyPart(l1, 0, INT_MAX);
        AncestrySegment *r2 = copyPart(r1, 0, INT_MAX);

        SegmentSpan wantSpan = {0, 1000, 0}, gotSpan = {0, 1000, 0};
        eventCount = 0;
        AncestrySegment *want = withoutFixed(referenceMerge(l1, r1, logParent,
            logEdge, fullCount, logFixed, &wantSpan), fullCount);
        int wantLog[128][4], wantCount = eventCount;
        memcpy(wantLog, eventLog, sizeof wantLog);

        SegmentSpan lSpan = spanOf(l2, fullCount), rSpan = spanOf(r2, fullCount);
        eventCount = 0;
        AncestrySegment *got = mergeAncestryTreesCoalescing(l2, r2, logParent,
            logEdge, NULL, fullCount, logFixed, &lSpan, &rSpan, &gotSpan);

        assertSameList(want, got);
        TEST_ASSERT_EQUAL(wantCount, eventCount);
        for (int k = 0; k < eventCount; k++)
            TEST_ASSERT_EQUAL_INT_ARRAY(wantLog[k], eventLog[k], 4);
        TEST_ASSERT_EQUAL(wantSpan.nancSites, gotSpan.nancSites);
        TEST_ASSERT_EQUAL(wantSpan.lLim, gotSpan.lLim);
        TEST_ASSERT_EQUAL(wantSpan.rLim, gotSpan.rLim);

        freeSegmentTree(want);
        freeSegmentTree(got);
        freeSegmentTree(l1);
        freeSegmentTree(r1);
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
    RUN_TEST(test_merge_matches_reference);
    RUN_TEST(test_segment_pool_free_list_holds_every_segment);
    
    return UNITY_END();
}
#endif
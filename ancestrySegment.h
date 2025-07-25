#ifndef __ANCESTRY_SEGMENT_H__
#define __ANCESTRY_SEGMENT_H__

#include <stdint.h>

typedef struct AncestrySegment {
    int start, end;  // genomic interval [start, end)
    struct AncestrySegment *left, *right;  // child segments (for tree structure)
    struct AncestrySegment *next;  // linked list for segments at same level
    uint16_t count;  // number of lineages
    int isLeaf;  // 1 if this is a leaf segment, 0 otherwise
    int refCount;  // Reference count for sharing
    void *avlTree;  // Optional AVL tree for fast lookups (only on root)
} AncestrySegment;

// Basic operations
AncestrySegment* newSegment(int start, int end, AncestrySegment *left, AncestrySegment *right);
void freeSegmentTree(AncestrySegment *root);
AncestrySegment* copySegmentTree(AncestrySegment *root);

// Reference counting operations
AncestrySegment* retainSegment(AncestrySegment *seg);
void releaseSegment(AncestrySegment *seg);
AncestrySegment* shallowCopySegment(AncestrySegment *seg);

// Query operations
uint16_t getAncestryCount(AncestrySegment *root, int site);
int hasAncestry(AncestrySegment *root, int site);

// Tree operations for coalescence and recombination
AncestrySegment* mergeAncestryTrees(AncestrySegment *left, AncestrySegment *right);
AncestrySegment* splitLeft(AncestrySegment *root, int breakpoint);
AncestrySegment* splitRight(AncestrySegment *root, int breakpoint);

// Structures for split operations
typedef struct {
    AncestrySegment *left;
    AncestrySegment *right;
} splitResult;

typedef struct {
    AncestrySegment *converted;
    AncestrySegment *unconverted;
} gcSplitResult;

// Gene conversion specific split
gcSplitResult splitSegmentTreeForGeneConversion(AncestrySegment *root, int startPos, int endPos);

// Debug operations
void printSegmentTree(AncestrySegment *root, int depth);
int verifySegmentTree(AncestrySegment *root, int nSites);

#endif
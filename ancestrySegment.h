#ifndef __ANCESTRY_SEGMENT_H__
#define __ANCESTRY_SEGMENT_H__

#include <stdint.h>
#include <tskit.h>

typedef struct AncestrySegment {
    int start, end;  // genomic interval [start, end)
    struct AncestrySegment *left, *right;  // child segments (for tree structure)
    struct AncestrySegment *next;  // linked list for segments at same level
    void *avlTree;  // Optional AVL tree for fast lookups (only on root)
    tsk_id_t tskit_node_id;  // tskit node ID this segment represents (like msprime's segment_t.value)
    uint16_t count;  // number of lineages
    uint16_t refCount;  // Reference count for sharing (max 65535)
    uint8_t isLeaf;  // 1 if this is a leaf segment, 0 otherwise
    uint8_t isRecorded;  // 1 if this segment has been recorded as edge in tskit, 0 otherwise
} AncestrySegment;

// Basic operations
AncestrySegment* newSegment(int start, int end, tsk_id_t tskit_node_id, AncestrySegment *left, AncestrySegment *right);
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
AncestrySegment* mergeAncestryTrees(AncestrySegment *left, AncestrySegment *right, tsk_id_t parent_tskit_node_id);
AncestrySegment* splitLeft(AncestrySegment *root, int breakpoint);
AncestrySegment* splitRight(AncestrySegment *root, int breakpoint);

// Structure for gene conversion split operation
typedef struct {
    AncestrySegment *converted;
    AncestrySegment *unconverted;
} gcSplitResult;

// Gene conversion specific split
gcSplitResult splitSegmentTreeForGeneConversion(AncestrySegment *root, int startPos, int endPos);

// Debug operations
void printSegmentTree(AncestrySegment *root, int depth);
int verifySegmentTree(AncestrySegment *root, int nSites);

// Recording status operations
int areAllSegmentsRecorded(AncestrySegment *root);
void markSegmentRecorded(AncestrySegment *seg);

#endif
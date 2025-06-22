#ifndef __ANCESTRY_SEGMENT_AVL_H__
#define __ANCESTRY_SEGMENT_AVL_H__

#include "ancestrySegment.h"

// AVL tree node for fast segment lookup
typedef struct AVLNode {
    AncestrySegment *segment;
    struct AVLNode *left;
    struct AVLNode *right;
    int height;
} AVLNode;

// AVL tree structure
typedef struct {
    AVLNode *root;
    int size;
} AVLTree;

// AVL tree operations
AVLTree* createAVLTree(void);
void freeAVLTree(AVLTree *tree);
void insertSegment(AVLTree *tree, AncestrySegment *segment);
AncestrySegment* findSegmentContaining(AVLTree *tree, int site);
AVLTree* buildAVLFromList(AncestrySegment *listHead);

// Convert between representations
AncestrySegment* convertAVLToList(AVLTree *tree);

#endif
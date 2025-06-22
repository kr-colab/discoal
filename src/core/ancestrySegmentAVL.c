#include <stdio.h>
#include <stdlib.h>
#include "ancestrySegmentAVL.h"

// Helper functions for AVL tree
static int max(int a, int b) {
    return (a > b) ? a : b;
}

static int height(AVLNode *node) {
    return node ? node->height : 0;
}

static int getBalance(AVLNode *node) {
    return node ? height(node->left) - height(node->right) : 0;
}

static AVLNode* createAVLNode(AncestrySegment *segment) {
    AVLNode *node = (AVLNode*)malloc(sizeof(AVLNode));
    if (!node) return NULL;
    
    node->segment = segment;
    node->left = NULL;
    node->right = NULL;
    node->height = 1;
    return node;
}

// Right rotate
static AVLNode* rotateRight(AVLNode *y) {
    AVLNode *x = y->left;
    AVLNode *T2 = x->right;
    
    x->right = y;
    y->left = T2;
    
    y->height = max(height(y->left), height(y->right)) + 1;
    x->height = max(height(x->left), height(x->right)) + 1;
    
    return x;
}

// Left rotate
static AVLNode* rotateLeft(AVLNode *x) {
    AVLNode *y = x->right;
    AVLNode *T2 = y->left;
    
    y->left = x;
    x->right = T2;
    
    x->height = max(height(x->left), height(x->right)) + 1;
    y->height = max(height(y->left), height(y->right)) + 1;
    
    return y;
}

// Insert a segment into AVL tree
static AVLNode* insertAVLNode(AVLNode *node, AncestrySegment *segment) {
    // Standard BST insertion
    if (!node) {
        return createAVLNode(segment);
    }
    
    // Use segment start as key
    if (segment->start < node->segment->start) {
        node->left = insertAVLNode(node->left, segment);
    } else {
        node->right = insertAVLNode(node->right, segment);
    }
    
    // Update height
    node->height = 1 + max(height(node->left), height(node->right));
    
    // Get balance factor
    int balance = getBalance(node);
    
    // Left Left Case
    if (balance > 1 && segment->start < node->left->segment->start) {
        return rotateRight(node);
    }
    
    // Right Right Case
    if (balance < -1 && segment->start > node->right->segment->start) {
        return rotateLeft(node);
    }
    
    // Left Right Case
    if (balance > 1 && segment->start > node->left->segment->start) {
        node->left = rotateLeft(node->left);
        return rotateRight(node);
    }
    
    // Right Left Case
    if (balance < -1 && segment->start < node->right->segment->start) {
        node->right = rotateRight(node->right);
        return rotateLeft(node);
    }
    
    return node;
}

// Find segment containing a site (O(log n))
static AncestrySegment* findInAVL(AVLNode *node, int site) {
    if (!node) return NULL;
    
    AncestrySegment *seg = node->segment;
    
    // Check if site is in current segment
    if (site >= seg->start && site < seg->end) {
        return seg;
    }
    
    // Search left or right
    if (site < seg->start) {
        return findInAVL(node->left, site);
    } else {
        return findInAVL(node->right, site);
    }
}

// Free AVL nodes
static void freeAVLNodes(AVLNode *node) {
    if (!node) return;
    freeAVLNodes(node->left);
    freeAVLNodes(node->right);
    free(node);
}

// Public functions
AVLTree* createAVLTree(void) {
    AVLTree *tree = (AVLTree*)malloc(sizeof(AVLTree));
    if (!tree) return NULL;
    
    tree->root = NULL;
    tree->size = 0;
    return tree;
}

void freeAVLTree(AVLTree *tree) {
    if (!tree) return;
    freeAVLNodes(tree->root);
    free(tree);
}

void insertSegment(AVLTree *tree, AncestrySegment *segment) {
    if (!tree) return;
    tree->root = insertAVLNode(tree->root, segment);
    tree->size++;
}

AncestrySegment* findSegmentContaining(AVLTree *tree, int site) {
    if (!tree) return NULL;
    return findInAVL(tree->root, site);
}

// Build AVL tree from linked list of segments
AVLTree* buildAVLFromList(AncestrySegment *listHead) {
    AVLTree *tree = createAVLTree();
    if (!tree) return NULL;
    
    AncestrySegment *current = listHead;
    while (current) {
        insertSegment(tree, current);
        current = current->next;
    }
    
    return tree;
}

// In-order traversal to convert AVL back to sorted list
static void inOrderToList(AVLNode *node, AncestrySegment **tail) {
    if (!node) return;
    
    inOrderToList(node->left, tail);
    
    // Add current segment to list
    if (*tail) {
        (*tail)->next = node->segment;
    }
    *tail = node->segment;
    node->segment->next = NULL;
    
    inOrderToList(node->right, tail);
}

AncestrySegment* convertAVLToList(AVLTree *tree) {
    if (!tree || !tree->root) return NULL;
    
    AncestrySegment *head = NULL;
    AncestrySegment *tail = NULL;
    
    // Find leftmost node as head
    AVLNode *current = tree->root;
    while (current->left) {
        current = current->left;
    }
    head = current->segment;
    
    // Convert tree to list
    inOrderToList(tree->root, &tail);
    
    return head;
}
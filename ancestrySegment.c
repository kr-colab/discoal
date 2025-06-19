#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "ancestrySegment.h"
#include "ancestrySegmentAVL.h"

AncestrySegment* newSegment(int start, int end, AncestrySegment *left, AncestrySegment *right) {
    AncestrySegment *seg = (AncestrySegment*)calloc(1, sizeof(AncestrySegment));
    if (!seg) {
        fprintf(stderr, "Memory allocation failed for AncestrySegment\n");
        exit(1);
    }
    seg->start = start;
    seg->end = end;
    seg->left = left;
    seg->right = right;
    seg->next = NULL;
    seg->isLeaf = (left == NULL && right == NULL) ? 1 : 0;
    seg->refCount = 1;  // Initial reference count
    seg->avlTree = NULL;  // AVL tree only used on root segments
    
    // For leaf segments, count is 1; for internal nodes, sum of children
    if (seg->isLeaf) {
        seg->count = 1;
    } else if (left && right) {
        seg->count = left->count + right->count;
    }
    
    return seg;
}

void freeSegmentTree(AncestrySegment *root) {
    if (!root) return;
    releaseSegment(root);
}

AncestrySegment* copySegmentTree(AncestrySegment *root) {
    if (!root) return NULL;
    
    // If the tree is immutable (no next pointers), we can share it
    if (!root->next && root->refCount > 0) {
        return retainSegment(root);
    }
    
    // Otherwise, create a deep copy
    AncestrySegment *newRoot = NULL;
    AncestrySegment *lastNew = NULL;
    AncestrySegment *current = root;
    
    while (current) {
        // For child segments, retain them instead of copying
        AncestrySegment *newSeg = newSegment(current->start, current->end, 
                                            retainSegment(current->left), 
                                            retainSegment(current->right));
        newSeg->count = current->count;
        newSeg->isLeaf = current->isLeaf;
        
        if (!newRoot) {
            newRoot = newSeg;
        } else {
            lastNew->next = newSeg;
        }
        lastNew = newSeg;
        current = current->next;
    }
    
    return newRoot;
}

uint16_t getAncestryCount(AncestrySegment *root, int site) {
    // Use AVL tree for O(log n) lookup if available
    if (root && root->avlTree) {
        AncestrySegment *found = findSegmentContaining((AVLTree*)root->avlTree, site);
        return found ? found->count : 0;
    }
    
    // Fall back to linear search
    AncestrySegment *current = root;
    while (current) {
        if (site >= current->start && site < current->end) {
            return current->count;
        }
        current = current->next;
    }
    
    return 0;  // No ancestry at this site
}

int hasAncestry(AncestrySegment *root, int site) {
    return getAncestryCount(root, site) > 0;
}

// Reference counting operations
AncestrySegment* retainSegment(AncestrySegment *seg) {
    if (seg) {
        seg->refCount++;
    }
    return seg;
}

void releaseSegment(AncestrySegment *seg) {
    if (!seg) return;
    
    seg->refCount--;
    if (seg->refCount <= 0) {
        // Release next segment in the list
        if (seg->next) {
            releaseSegment(seg->next);
        }
        // Release child segments
        if (seg->left) {
            releaseSegment(seg->left);
        }
        if (seg->right) {
            releaseSegment(seg->right);
        }
        // Free AVL tree if present
        if (seg->avlTree) {
            freeAVLTree((AVLTree*)seg->avlTree);
            seg->avlTree = NULL;
        }
        free(seg);
    }
}

// Create a shallow copy that shares the segment data
AncestrySegment* shallowCopySegment(AncestrySegment *seg) {
    if (!seg) return NULL;
    
    // For immutable segments, we can just retain and return
    return retainSegment(seg);
}

// Helper function to create a segment list covering an interval
static AncestrySegment* createSegmentList(int start, int end, uint16_t count) {
    if (start >= end) return NULL;
    return newSegment(start, end, NULL, NULL);
}

// Helper function to coalesce adjacent segments with same count
static AncestrySegment* coalesceSegments(AncestrySegment *head) {
    if (!head || !head->next) return head;
    
    AncestrySegment *current = head;
    while (current && current->next) {
        if (current->end == current->next->start && 
            current->count == current->next->count) {
            // Merge with next segment
            AncestrySegment *toRemove = current->next;
            current->end = toRemove->end;
            current->next = toRemove->next;
            free(toRemove);
            // Don't advance current, check if we can merge with the new next
        } else {
            current = current->next;
        }
    }
    
    return head;
}

// Helper function to add a segment to the result list
static void addSegmentToResult(AncestrySegment **result, AncestrySegment **tail, 
                               int start, int end, uint16_t count) {
    if (start >= end || count == 0) return;  // Skip empty or zero-count segments
    
    // Check if we can merge with the previous segment
    if (*tail && (*tail)->end == start && (*tail)->count == count) {
        (*tail)->end = end;  // Extend previous segment
    } else {
        // Create new segment
        AncestrySegment *newSeg = newSegment(start, end, NULL, NULL);
        newSeg->count = count;
        
        if (!*result) {
            *result = newSeg;
            *tail = newSeg;
        } else {
            (*tail)->next = newSeg;
            *tail = newSeg;
        }
    }
}

// Merge two ancestry trees (for coalescence)
AncestrySegment* mergeAncestryTrees(AncestrySegment *leftTree, AncestrySegment *rightTree) {
    if (!leftTree) return copySegmentTree(rightTree);
    if (!rightTree) return copySegmentTree(leftTree);
    
    // Optimized case: both trees have single segments covering all sites
    if (leftTree && rightTree && leftTree->start == 0 && rightTree->start == 0 && 
        leftTree->end == rightTree->end && !leftTree->next && !rightTree->next) {
        // Common case: both cover all sites with no recombination
        AncestrySegment *merged = newSegment(0, leftTree->end, 
                                           retainSegment(leftTree), 
                                           retainSegment(rightTree));
        merged->count = leftTree->count + rightTree->count;
        return merged;
    }
    
    // General interval merging using a sweep-line approach
    AncestrySegment *result = NULL;
    AncestrySegment *tail = NULL;
    AncestrySegment *l = leftTree;
    AncestrySegment *r = rightTree;
    
    // Start position should be the minimum of both trees
    int pos = INT_MAX;
    if (l && l->start < pos) pos = l->start;
    if (r && r->start < pos) pos = r->start;
    if (pos == INT_MAX) return NULL;  // Both trees empty
    
    while (l || r) {
        // Find the next position where something changes
        int nextPos = INT_MAX;
        
        // Consider starts and ends of segments
        if (l) {
            if (l->start > pos && l->start < nextPos) nextPos = l->start;
            if (l->end > pos && l->end < nextPos) nextPos = l->end;
        }
        if (r) {
            if (r->start > pos && r->start < nextPos) nextPos = r->start;
            if (r->end > pos && r->end < nextPos) nextPos = r->end;
        }
        
        if (nextPos == INT_MAX) break;  // No more events
        
        // Calculate total count for interval [pos, nextPos)
        uint16_t count = 0;
        if (l && pos >= l->start && pos < l->end) count += l->count;
        if (r && pos >= r->start && pos < r->end) count += r->count;
        
        // Add segment if count > 0
        if (count > 0) {
            addSegmentToResult(&result, &tail, pos, nextPos, count);
        }
        
        // Move to next position
        pos = nextPos;
        
        // Advance pointers past segments that have ended
        while (l && l->end <= pos) l = l->next;
        while (r && r->end <= pos) r = r->next;
    }
    
    // Build AVL tree if result has enough segments for it to be worthwhile
    if (result) {
        int segmentCount = 0;
        AncestrySegment *temp = result;
        while (temp && segmentCount < 3) {  // Count up to threshold
            segmentCount++;
            temp = temp->next;
        }
        
        // Build AVL tree for fast lookups if we have 3 or more segments
        if (segmentCount >= 30) {
            result->avlTree = buildAVLFromList(result);
        }
    }
    
    return result;
}

// Split ancestry tree at breakpoint (for recombination) - left side
AncestrySegment* splitLeft(AncestrySegment *root, int breakpoint) {
    if (!root || breakpoint <= 0) return NULL;
    
    // Special case: if breakpoint is beyond all segments, share the entire tree
    AncestrySegment *last = root;
    while (last->next) last = last->next;
    if (breakpoint >= last->end) {
        // Can share the entire tree
        if (!root->next) {
            return retainSegment(root);
        } else {
            return copySegmentTree(root);
        }
    }
    
    // Check if we can share a prefix of the tree
    if (!root->next && root->start == 0 && breakpoint >= root->end) {
        return retainSegment(root);
    }
    
    AncestrySegment *result = NULL;
    AncestrySegment *tail = NULL;
    AncestrySegment *current = root;
    
    while (current && current->start < breakpoint) {
        if (current->end <= breakpoint) {
            // Entire segment is to the left
            addSegmentToResult(&result, &tail, current->start, current->end, current->count);
        } else {
            // Segment spans the breakpoint
            addSegmentToResult(&result, &tail, current->start, breakpoint, current->count);
        }
        current = current->next;
    }
    
    return result;
}

// Split ancestry tree at breakpoint (for recombination) - right side
AncestrySegment* splitRight(AncestrySegment *root, int breakpoint) {
    if (!root) return NULL;
    
    // Special case: if breakpoint is at or before the start, share the entire tree
    if (breakpoint <= root->start) {
        if (!root->next) {
            return retainSegment(root);
        } else {
            return copySegmentTree(root);
        }
    }
    
    // Check if we can share a suffix of the tree
    AncestrySegment *firstRight = root;
    while (firstRight && firstRight->end <= breakpoint) {
        firstRight = firstRight->next;
    }
    
    if (firstRight && !firstRight->next && firstRight->start >= breakpoint) {
        // Single segment entirely to the right
        return retainSegment(firstRight);
    }
    
    AncestrySegment *result = NULL;
    AncestrySegment *tail = NULL;
    AncestrySegment *current = root;
    
    while (current) {
        if (current->end > breakpoint) {
            if (current->start >= breakpoint) {
                // Entire segment is to the right
                addSegmentToResult(&result, &tail, current->start, current->end, current->count);
            } else {
                // Segment spans the breakpoint
                addSegmentToResult(&result, &tail, breakpoint, current->end, current->count);
            }
        }
        current = current->next;
    }
    
    return result;
}

void printSegmentTree(AncestrySegment *root, int depth) {
    AncestrySegment *current = root;
    while (current) {
        for (int i = 0; i < depth; i++) printf("  ");
        printf("[%d,%d) count=%d%s\n", current->start, current->end, 
               current->count, current->isLeaf ? " (leaf)" : "");
        current = current->next;
    }
}

int verifySegmentTree(AncestrySegment *root, int nSites) {
    if (!root) return 1;  // Empty tree is valid
    
    AncestrySegment *current = root;
    int lastEnd = 0;
    
    while (current) {
        // Check for gaps or overlaps
        if (current->start != lastEnd) {
            fprintf(stderr, "Gap or overlap in segments: %d to %d\n", lastEnd, current->start);
            return 0;
        }
        
        // Check bounds
        if (current->start < 0 || current->end > nSites) {
            fprintf(stderr, "Segment out of bounds: [%d,%d)\n", current->start, current->end);
            return 0;
        }
        
        lastEnd = current->end;
        current = current->next;
    }
    
    return 1;
}

// Check if all segments in the tree have been recorded
int areAllSegmentsRecorded(AncestrySegment *root) {
    if (!root) return 1;  // Empty tree is considered fully recorded
    
    AncestrySegment *current = root;
    while (current) {
        if (!current->isRecorded) {
            return 0;  // Found an unrecorded segment
        }
        current = current->next;
    }
    
    return 1;  // All segments are recorded
}

// Mark a segment as recorded
void markSegmentRecorded(AncestrySegment *seg) {
    if (seg && seg->refCount == 1) {
        // Only mark as recorded if this segment is not shared
        seg->isRecorded = 1;
    }
}

// Split ancestry tree for gene conversion
// Returns converted tract and everything else
gcSplitResult splitSegmentTreeForGeneConversion(AncestrySegment *root, int startPos, int endPos) {
    gcSplitResult result;
    result.converted = NULL;
    result.unconverted = NULL;
    
    if (!root || startPos >= endPos) {
        result.unconverted = copySegmentTree(root);
        return result;
    }
    
    AncestrySegment *convTail = NULL;
    AncestrySegment *unconvTail = NULL;
    AncestrySegment *current = root;
    
    while (current) {
        if (current->end <= startPos || current->start >= endPos) {
            // Segment is completely outside conversion tract
            addSegmentToResult(&result.unconverted, &unconvTail, 
                             current->start, current->end, current->count);
        } else if (current->start >= startPos && current->end <= endPos) {
            // Segment is completely inside conversion tract
            addSegmentToResult(&result.converted, &convTail,
                             current->start, current->end, current->count);
        } else {
            // Segment overlaps with conversion tract boundary
            if (current->start < startPos) {
                // Part before conversion tract
                addSegmentToResult(&result.unconverted, &unconvTail,
                                 current->start, startPos, current->count);
            }
            
            // Part inside conversion tract
            int convStart = (current->start > startPos) ? current->start : startPos;
            int convEnd = (current->end < endPos) ? current->end : endPos;
            if (convStart < convEnd) {
                addSegmentToResult(&result.converted, &convTail,
                                 convStart, convEnd, current->count);
            }
            
            if (current->end > endPos) {
                // Part after conversion tract
                addSegmentToResult(&result.unconverted, &unconvTail,
                                 endPos, current->end, current->count);
            }
        }
        current = current->next;
    }
    
    return result;
}
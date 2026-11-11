#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "ancestrySegment.h"
#include "segmentPool.h"

AncestrySegment* newSegment(int start, int end, tsk_id_t tskit_node_id, AncestrySegment *left, AncestrySegment *right) {
    AncestrySegment *seg = allocSegmentFromPool();
    if (!seg) {
        fprintf(stderr, "Memory allocation failed for AncestrySegment\n");
        exit(1);
    }
    seg->start = start;
    seg->end = end;
    seg->tskit_node_id = tskit_node_id;  // Store the tskit node ID
    seg->left = left;
    seg->right = right;
    seg->next = NULL;
    seg->isLeaf = (left == NULL && right == NULL) ? 1 : 0;
    seg->refCount = 1;  // Initial reference count
    
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
        AncestrySegment *newSeg = newSegment(current->start, current->end, current->tskit_node_id,
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

/* The segment of a lineage that covers `site`, or NULL: a walk of the
 * sorted list. Only rare paths look up a site (migration that depends on
 * the sweep site, partial-sweep carriers), so no search structure is kept
 * up to date at every coalescence. */
AncestrySegment* segmentAt(AncestrySegment *root, int site) {
    for (AncestrySegment *current = root; current; current = current->next) {
        if (site >= current->start && site < current->end) {
            return current;
        }
    }
    return NULL;
}

uint16_t getAncestryCount(AncestrySegment *root, int site) {
    AncestrySegment *found = segmentAt(root, site);
    return found ? found->count : 0;
}

int hasAncestry(AncestrySegment *root, int site) {
    return getAncestryCount(root, site) > 0;
}

// Reference counting operations
AncestrySegment* retainSegment(AncestrySegment *seg) {
    if (seg) {
        if (seg->refCount >= UINT16_MAX) {
            fprintf(stderr, "Error: AncestrySegment reference count overflow\n");
            exit(1);
        }
        seg->refCount++;
    }
    return seg;
}

// Stack node for iterative traversal
typedef struct SegmentStackNode {
    AncestrySegment *segment;
    struct SegmentStackNode *next;
} SegmentStackNode;

void releaseSegment(AncestrySegment *seg) {
    if (!seg) return;
    
    // Initialize stack with root segment
    SegmentStackNode *stack = NULL;
    SegmentStackNode *node = (SegmentStackNode*)malloc(sizeof(SegmentStackNode));
    node->segment = seg;
    node->next = NULL;
    stack = node;
    
    while (stack) {
        // Pop from stack
        SegmentStackNode *top = stack;
        AncestrySegment *current = top->segment;
        stack = top->next;
        free(top);
        
        if (!current) continue;
        
        // Decrement reference count
        current->refCount--;
        
        if (current->refCount <= 0) {
            // Process the entire 'next' chain for this segment
            AncestrySegment *chain_current = current;
            int chain_length = 0;
            
            while (chain_current) {
                chain_length++;
                AncestrySegment *next = chain_current->next;
                
                // Push children onto stack for later processing
                if (chain_current->left) {
                    node = (SegmentStackNode*)malloc(sizeof(SegmentStackNode));
                    node->segment = chain_current->left;
                    node->next = stack;
                    stack = node;
                }
                if (chain_current->right) {
                    node = (SegmentStackNode*)malloc(sizeof(SegmentStackNode));
                    node->segment = chain_current->right;
                    node->next = stack;
                    stack = node;
                }
                
                // Free current segment
                freeSegmentToPool(chain_current);
                
                // Handle next segment in chain
                if (next && next->refCount > 0) {
                    next->refCount--;
                    if (next->refCount > 0) {
                        break;  // Next segment still has references
                    }
                }
                chain_current = next;
            }
        }
    }
}

// Create a shallow copy that shares the segment data
AncestrySegment* shallowCopySegment(AncestrySegment *seg) {
    if (!seg) return NULL;
    
    // For immutable segments, we can just retain and return
    return retainSegment(seg);
}

// Helper function to add a segment to the result list
static void addSegmentToResult(AncestrySegment **result, AncestrySegment **tail, 
                               int start, int end, uint16_t count, tsk_id_t parent_tskit_node_id) {
    if (start >= end || count == 0) return;  // Skip empty or zero-count segments
    
    // Extend the previous segment if this one continues it: adjacent, with
    // the same count and the same tskit node
    if (*tail && (*tail)->end == start && (*tail)->count == count &&
        (*tail)->tskit_node_id == parent_tskit_node_id) {
        (*tail)->end = end;  // Extend previous segment
    } else {
        // Create new segment
        AncestrySegment *newSeg = newSegment(start, end, parent_tskit_node_id, NULL, NULL);
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
AncestrySegment* mergeAncestryTrees(AncestrySegment *leftTree, AncestrySegment *rightTree, tsk_id_t parent_tskit_node_id) {
    if (!leftTree && !rightTree) return NULL;
    
    // If only one tree exists, create segments with parent node ID
    if (!leftTree) {
        AncestrySegment *result = NULL;
        AncestrySegment *tail = NULL;
        AncestrySegment *current = rightTree;
        while (current) {
            addSegmentToResult(&result, &tail, current->start, current->end, current->count, parent_tskit_node_id);
            current = current->next;
        }
        return result;
    }
    
    if (!rightTree) {
        AncestrySegment *result = NULL;
        AncestrySegment *tail = NULL;
        AncestrySegment *current = leftTree;
        while (current) {
            addSegmentToResult(&result, &tail, current->start, current->end, current->count, parent_tskit_node_id);
            current = current->next;
        }
        return result;
    }
    
    // Optimized case: both trees have single segments covering all sites
    if (leftTree && rightTree && leftTree->start == 0 && rightTree->start == 0 && 
        leftTree->end == rightTree->end && !leftTree->next && !rightTree->next) {
        // Common case: both cover all sites with no recombination
        AncestrySegment *merged = newSegment(0, leftTree->end, parent_tskit_node_id,
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
            addSegmentToResult(&result, &tail, pos, nextPos, count, parent_tskit_node_id);
        }
        
        // Move to next position
        pos = nextPos;
        
        // Advance pointers past segments that have ended
        while (l && l->end <= pos) l = l->next;
        while (r && r->end <= pos) r = r->next;
    }
    
    return result;
}

/* Merge the ancestry of two lineages at a common-ancestor event, recording
 * only what coalesces, as msprime does (msp_merge_two_ancestors):
 *  - where only one lineage has ancestral material, the parent's segment
 *    keeps that lineage's tskit node: nothing coalesces there, so no node
 *    and no edge are needed;
 *  - where both have material, they coalesce: the parent needs a tskit
 *    node, which newParentNode(ctx) makes at the first such interval, and
 *    an edge from it to each child's node over the interval (addEdge).
 * So a common-ancestor event of lineages whose material does not overlap
 * adds nothing to the tree sequence.
 * An overlap whose count reaches fullCount (all samples) has found its
 * MRCA here, and only here: onFixed(start, end) reports each such run of
 * sites once, so the caller can update the sites still active without
 * walking the parent's other segments. If span is not NULL, the span of
 * the result is added to it, so the caller need not walk the result. */
AncestrySegment* mergeAncestryTreesCoalescing(AncestrySegment *leftTree,
        AncestrySegment *rightTree, tsk_id_t (*newParentNode)(void *ctx),
        void (*addEdge)(tsk_id_t parent, tsk_id_t child, double left, double right),
        void *ctx, int fullCount, void (*onFixed)(int start, int end),
        SegmentSpan *span) {
    AncestrySegment *result = NULL, *tail = NULL;
    AncestrySegment *l = leftTree, *r = rightTree;
    tsk_id_t parent = TSK_NULL;
    int fixedStart = -1, fixedEnd = -1;  // a run of newly fixed sites
    int pos = INT_MAX;
    if (l) pos = l->start;
    if (r && r->start < pos) pos = r->start;

    while (l || r) {
        /* The next position where a segment starts or ends. */
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
            if (parent == TSK_NULL) parent = newParentNode(ctx);
            addSegmentToResult(&result, &tail, pos, nextPos,
                               lc->count + rc->count, parent);
            if (span) spanAdd(span, pos, nextPos, lc->count + rc->count, fullCount);
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
            addSegmentToResult(&result, &tail, pos, nextPos, c->count,
                               c->tskit_node_id);
            if (span) spanAdd(span, pos, nextPos, c->count, fullCount);
        }

        pos = nextPos;
        while (l && l->end <= pos) l = l->next;
        while (r && r->end <= pos) r = r->next;
    }

    if (fixedStart >= 0) onFixed(fixedStart, fixedEnd);

    return result;
}

// Split ancestry tree at breakpoint (for recombination) - left side
/* Split a segment list at `breakpoint` in place, as msprime does: the
 * segments left of it stay in *left, the rest go to *right, and a segment
 * that spans it is cut, with one new segment for its right part. This
 * moves the list into the two halves, so the caller must drop its own
 * pointer to it. Only for a list that nothing else shares (every segment
 * has a reference count of 1 and no child segments, as in minimal tree
 * sequence mode); for any other list it changes nothing and returns 0,
 * and the caller copies the halves with splitLeft() and splitRight().
 * The spans of the two halves are added to *leftSpan and *rightSpan
 * (valid only when it returns 1). */
int splitInPlace(AncestrySegment *root, int breakpoint, int fullCount,
                 AncestrySegment **left, AncestrySegment **right,
                 SegmentSpan *leftSpan, SegmentSpan *rightSpan) {
    AncestrySegment *lastLeft = NULL;   /* last segment that starts left of it */
    for (AncestrySegment *s = root; s; s = s->next) {
        if (s->refCount != 1 || s->left || s->right) return 0;
        if (s->start < breakpoint) lastLeft = s;
        if (s->end <= breakpoint) {
            spanAdd(leftSpan, s->start, s->end, s->count, fullCount);
        } else if (s->start >= breakpoint) {
            spanAdd(rightSpan, s->start, s->end, s->count, fullCount);
        } else {
            spanAdd(leftSpan, s->start, breakpoint, s->count, fullCount);
            spanAdd(rightSpan, breakpoint, s->end, s->count, fullCount);
        }
    }
    if (lastLeft == NULL) {             /* everything is right of it */
        *left = NULL;
        *right = root;
        return 1;
    }
    AncestrySegment *rightHead;
    if (lastLeft->end > breakpoint) {   /* cut the spanning segment */
        rightHead = newSegment(breakpoint, lastLeft->end,
                               lastLeft->tskit_node_id, NULL, NULL);
        rightHead->count = lastLeft->count;
        rightHead->next = lastLeft->next;
        lastLeft->end = breakpoint;
    } else {
        rightHead = lastLeft->next;
    }
    lastLeft->next = NULL;
    *left = root;
    *right = rightHead;
    return 1;
}

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
            addSegmentToResult(&result, &tail, current->start, current->end, current->count, current->tskit_node_id);
        } else {
            // Segment spans the breakpoint
            addSegmentToResult(&result, &tail, current->start, breakpoint, current->count, current->tskit_node_id);
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
                addSegmentToResult(&result, &tail, current->start, current->end, current->count, current->tskit_node_id);
            } else {
                // Segment spans the breakpoint
                addSegmentToResult(&result, &tail, breakpoint, current->end, current->count, current->tskit_node_id);
            }
        }
        current = current->next;
    }
    
    return result;
}

// Full ARG mode split - left side with parent node ID
AncestrySegment* splitLeftWithParent(AncestrySegment *root, int breakpoint, tsk_id_t parent_tskit_node_id) {
    if (!root || breakpoint <= root->start) return NULL;
    
    AncestrySegment *result = NULL;
    AncestrySegment *tail = NULL;
    AncestrySegment *current = root;
    
    while (current && current->start < breakpoint) {
        if (current->end <= breakpoint) {
            // Entire segment is to the left - use parent's tskit_node_id
            addSegmentToResult(&result, &tail, current->start, current->end, current->count, parent_tskit_node_id);
        } else {
            // Segment spans the breakpoint - use parent's tskit_node_id
            addSegmentToResult(&result, &tail, current->start, breakpoint, current->count, parent_tskit_node_id);
        }
        current = current->next;
    }
    
    return result;
}

// Full ARG mode split - right side with parent node ID
AncestrySegment* splitRightWithParent(AncestrySegment *root, int breakpoint, tsk_id_t parent_tskit_node_id) {
    if (!root) return NULL;
    
    // Special case: if breakpoint is at or before the start, share the entire tree but update IDs
    if (breakpoint <= root->start) {
        AncestrySegment *result = NULL;
        AncestrySegment *tail = NULL;
        AncestrySegment *current = root;
        
        while (current) {
            addSegmentToResult(&result, &tail, current->start, current->end, current->count, parent_tskit_node_id);
            current = current->next;
        }
        return result;
    }
    
    AncestrySegment *result = NULL;
    AncestrySegment *tail = NULL;
    AncestrySegment *current = root;
    
    while (current) {
        if (current->end > breakpoint) {
            if (current->start >= breakpoint) {
                // Entire segment is to the right - use parent's tskit_node_id
                addSegmentToResult(&result, &tail, current->start, current->end, current->count, parent_tskit_node_id);
            } else {
                // Segment spans the breakpoint - use parent's tskit_node_id
                addSegmentToResult(&result, &tail, breakpoint, current->end, current->count, parent_tskit_node_id);
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
                             current->start, current->end, current->count, current->tskit_node_id);
        } else if (current->start >= startPos && current->end <= endPos) {
            // Segment is completely inside conversion tract
            addSegmentToResult(&result.converted, &convTail,
                             current->start, current->end, current->count, current->tskit_node_id);
        } else {
            // Segment overlaps with conversion tract boundary
            if (current->start < startPos) {
                // Part before conversion tract
                addSegmentToResult(&result.unconverted, &unconvTail,
                                 current->start, startPos, current->count, current->tskit_node_id);
            }
            
            // Part inside conversion tract
            int convStart = (current->start > startPos) ? current->start : startPos;
            int convEnd = (current->end < endPos) ? current->end : endPos;
            if (convStart < convEnd) {
                addSegmentToResult(&result.converted, &convTail,
                                 convStart, convEnd, current->count, current->tskit_node_id);
            }
            
            if (current->end > endPos) {
                // Part after conversion tract
                addSegmentToResult(&result.unconverted, &unconvTail,
                                 endPos, current->end, current->count, current->tskit_node_id);
            }
        }
        current = current->next;
    }
    
    return result;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "activeSegment.h"

// Create a new active segment
ActiveSegment* newActiveSegment(int start, int end) {
    ActiveSegment *seg = (ActiveSegment*)malloc(sizeof(ActiveSegment));
    if (!seg) {
        fprintf(stderr, "Memory allocation failed for ActiveSegment\n");
        exit(1);
    }
    seg->start = start;
    seg->end = end;
    seg->next = NULL;
    return seg;
}

// Free a single segment
void freeActiveSegment(ActiveSegment *seg) {
    if (seg) free(seg);
}

// Free entire segment list
static void freeSegmentList(ActiveSegment *head) {
    while (head) {
        ActiveSegment *next = head->next;
        freeActiveSegment(head);
        head = next;
    }
}

// Initialize active material - all sites start as active
void initializeActiveMaterial(ActiveMaterial *am, int nSites) {
    if (!am) return;
    
    // Start with single segment covering all sites
    am->segments = newActiveSegment(0, nSites);
    am->avlTree = NULL;
    am->totalActive = nSites;
}

// Free all memory associated with active material
void freeActiveMaterial(ActiveMaterial *am) {
    if (!am) return;
    
    freeSegmentList(am->segments);
    am->segments = NULL;
    
    if (am->avlTree) {
        freeAVLTree((AVLTree*)am->avlTree);
        am->avlTree = NULL;
    }
    
    am->totalActive = 0;
}

// Check if a site is active - uses AVL tree if available
int isActiveSite(ActiveMaterial *am, int site) {
    if (!am || !am->segments) return 0;
    
    // Use AVL tree for fast lookup if available
    if (am->avlTree) {
        ActiveSegment *found = (ActiveSegment*)findSegmentContaining((AVLTree*)am->avlTree, site);
        return found != NULL;
    }
    
    // Linear search through segments
    ActiveSegment *seg = am->segments;
    while (seg) {
        if (site >= seg->start && site < seg->end) {
            return 1;
        }
        seg = seg->next;
    }
    
    return 0;
}

// Get total count of active sites
int getActiveSiteCount(ActiveMaterial *am) {
    return am ? am->totalActive : 0;
}

// Coalesce adjacent segments with no gap between them
ActiveSegment* coalesceActiveSegments(ActiveSegment *head) {
    if (!head || !head->next) return head;
    
    ActiveSegment *current = head;
    while (current && current->next) {
        if (current->end == current->next->start) {
            // Merge with next segment
            ActiveSegment *toRemove = current->next;
            current->end = toRemove->end;
            current->next = toRemove->next;
            freeActiveSegment(toRemove);
            // Don't advance - check if we can merge with new next
        } else {
            current = current->next;
        }
    }
    
    return head;
}

// Remove a fixed region from the active segments
void removeFixedRegion(ActiveMaterial *am, int start, int end) {
    if (!am || !am->segments || start >= end) return;
    
    ActiveSegment *prev = NULL;
    ActiveSegment *current = am->segments;
    
    while (current) {
        // Case 1: Fixed region completely before this segment
        if (end <= current->start) {
            break;  // No more segments affected
        }
        
        // Case 2: Fixed region completely after this segment
        if (start >= current->end) {
            prev = current;
            current = current->next;
            continue;
        }
        
        // Case 3: Fixed region overlaps with this segment
        if (start <= current->start && end >= current->end) {
            // Entire segment is fixed - remove it
            am->totalActive -= (current->end - current->start);
            
            if (prev) {
                prev->next = current->next;
            } else {
                am->segments = current->next;
            }
            
            ActiveSegment *toRemove = current;
            current = current->next;
            freeActiveSegment(toRemove);
        } else if (start > current->start && end < current->end) {
            // Fixed region splits this segment
            int oldEnd = current->end;
            current->end = start;  // Truncate current segment
            
            // Create new segment for the part after fixed region
            ActiveSegment *newSeg = newActiveSegment(end, oldEnd);
            newSeg->next = current->next;
            current->next = newSeg;
            
            am->totalActive -= (end - start);
            
            prev = newSeg;  // Skip the new segment
            current = newSeg->next;
        } else if (start <= current->start) {
            // Fixed region covers start of segment
            am->totalActive -= (end - current->start);
            current->start = end;
            prev = current;
            current = current->next;
        } else {
            // Fixed region covers end of segment
            am->totalActive -= (current->end - start);
            current->end = start;
            prev = current;
            current = current->next;
        }
    }
    
    // Invalidate AVL tree after modifications
    if (am->avlTree) {
        freeAVLTree((AVLTree*)am->avlTree);
        am->avlTree = NULL;
    }
}

// Update active material based on ancestry segments
void updateActiveMaterialFromAncestry(ActiveMaterial *am, AncestrySegment *ancestry, 
                                     int sampleSize, int nSites) {
    if (!am || !ancestry) return;
    
    // Walk through ancestry segments and remove fixed regions
    AncestrySegment *seg = ancestry;
    while (seg) {
        if (seg->count == sampleSize) {
            // This region is fixed - remove from active segments
            removeFixedRegion(am, seg->start, seg->end);
        }
        seg = seg->next;
    }
    
    // Coalesce adjacent segments
    am->segments = coalesceActiveSegments(am->segments);
    
    // Rebuild AVL tree if we have many segments
    if (am->avlTree) {
        freeAVLTree((AVLTree*)am->avlTree);
        am->avlTree = NULL;
    }
    
    // Count segments and potentially build AVL tree
    int segmentCount = 0;
    ActiveSegment *temp = am->segments;
    while (temp && segmentCount < 10) {
        segmentCount++;
        temp = temp->next;
    }
    
    if (segmentCount >= 10) {
        // Build AVL tree for fast lookups
        AVLTree *tree = createAVLTree();
        if (tree) {
            temp = am->segments;
            while (temp) {
                // Insert using the segment pointer as data
                insertSegment(tree, (AncestrySegment*)temp);
                temp = temp->next;
            }
            am->avlTree = tree;
        }
    }
}

// Print active segments for debugging
void printActiveSegments(ActiveMaterial *am) {
    if (!am) {
        printf("ActiveMaterial: NULL\n");
        return;
    }
    
    printf("ActiveMaterial: %d total active sites\n", am->totalActive);
    ActiveSegment *seg = am->segments;
    int count = 0;
    while (seg) {
        printf("  Segment %d: [%d, %d)\n", count++, seg->start, seg->end);
        seg = seg->next;
    }
    if (am->avlTree) {
        printf("  AVL tree: present\n");
    }
}

// Verify integrity of active material structure
int verifyActiveMaterial(ActiveMaterial *am, int nSites) {
    if (!am) return 0;
    
    int totalCount = 0;
    ActiveSegment *seg = am->segments;
    int lastEnd = -1;
    
    while (seg) {
        // Check bounds
        if (seg->start < 0 || seg->end > nSites || seg->start >= seg->end) {
            fprintf(stderr, "Invalid segment bounds: [%d, %d)\n", seg->start, seg->end);
            return 0;
        }
        
        // Check ordering
        if (seg->start <= lastEnd) {
            fprintf(stderr, "Segments overlap or out of order\n");
            return 0;
        }
        
        totalCount += (seg->end - seg->start);
        lastEnd = seg->end;
        seg = seg->next;
    }
    
    // Check total count
    if (totalCount != am->totalActive) {
        fprintf(stderr, "Total count mismatch: counted %d, stored %d\n", 
                totalCount, am->totalActive);
        return 0;
    }
    
    return 1;
}
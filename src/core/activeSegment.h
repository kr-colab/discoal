#ifndef __ACTIVE_SEGMENT_H__
#define __ACTIVE_SEGMENT_H__

#include "ancestrySegment.h"
#include "ancestrySegmentAVL.h"

// Segment representing a contiguous region of active material
typedef struct ActiveSegment {
    int start;
    int end;
    struct ActiveSegment *next;
} ActiveSegment;

// Main structure for tracking active material
typedef struct {
    ActiveSegment *segments;  // Linked list of active segments
    void *avlTree;           // Unused; NULL
    int totalActive;         // Total count of active sites
    /* The list as sorted arrays of starts and ends, for binary search in
     * isActiveSite(). Rebuilt from the list when a lookup finds it stale
     * (lookupStale), as every change of the list sets it. */
    int *lookupStart, *lookupEnd;
    int lookupCount, lookupCapacity, lookupStale;
} ActiveMaterial;

// Core operations
void initializeActiveMaterial(ActiveMaterial *am, int nSites);
void freeActiveMaterial(ActiveMaterial *am);
int isActiveSite(ActiveMaterial *am, int site);
int getActiveSiteCount(ActiveMaterial *am);

// Update active material based on ancestry information
void updateActiveMaterialFromAncestry(ActiveMaterial *am, AncestrySegment *ancestry, 
                                     int sampleSize, int nSites);

// Utility functions
void printActiveSegments(ActiveMaterial *am);
int verifyActiveMaterial(ActiveMaterial *am, int nSites);

// Internal helpers (exposed for testing)
ActiveSegment* newActiveSegment(int start, int end);
void freeActiveSegment(ActiveSegment *seg);
ActiveSegment* coalesceActiveSegments(ActiveSegment *head);
void removeFixedRegion(ActiveMaterial *am, int start, int end);

#endif
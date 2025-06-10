#ifndef __ANCESTRY_WRAPPER_H__
#define __ANCESTRY_WRAPPER_H__

#include "discoal.h"

// Wrapper functions that use both ancSites array and ancestry tree
// This allows gradual migration from array to tree

// Get ancestry count at a site
static inline uint16_t getAncestryAt(rootedNode *node, int site) {
    #ifdef USE_ANCESTRY_TREE_ONLY
        return node->ancestryRoot ? getAncestryCount(node->ancestryRoot, site) : 0;
    #else
        // During transition, verify consistency
        #ifdef DEBUG_ANCESTRY
        if (node->ancestryRoot) {
            uint16_t treeCount = getAncestryCount(node->ancestryRoot, site);
            uint16_t arrayCount = node->ancSites[site];
            if (treeCount != arrayCount) {
                fprintf(stderr, "ERROR: Ancestry mismatch at node %d site %d: tree=%d, array=%d\n",
                        node->id, site, treeCount, arrayCount);
            }
        }
        #endif
        return node->ancSites[site];
    #endif
}

// Set ancestry count at a site
static inline void setAncestryAt(rootedNode *node, int site, uint16_t count) {
    #ifndef USE_ANCESTRY_TREE_ONLY
        node->ancSites[site] = count;
    #endif
    // Tree update will be handled by merge/split operations
}

// Check if site has any ancestry
static inline int hasAncestryAt(rootedNode *node, int site) {
    return getAncestryAt(node, site) > 0;
}

// Check if site is polymorphic (has ancestry but not fixed)
static inline int isPolymorphicAt(rootedNode *node, int site, int sampleSize) {
    uint16_t count = getAncestryAt(node, site);
    return count > 0 && count < sampleSize;
}

// Check if site is fixed (MRCA reached)
static inline int isFixedAt(rootedNode *node, int site, int sampleSize) {
    return getAncestryAt(node, site) == sampleSize;
}

#endif
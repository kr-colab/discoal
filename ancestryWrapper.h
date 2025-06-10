#ifndef __ANCESTRY_WRAPPER_H__
#define __ANCESTRY_WRAPPER_H__

#include "discoal.h"

// Wrapper functions for ancestry tree operations

// Get ancestry count at a site
static inline uint16_t getAncestryAt(rootedNode *node, int site) {
    return node->ancestryRoot ? getAncestryCount(node->ancestryRoot, site) : 0;
}

// Note: Setting ancestry at individual sites is no longer supported
// Ancestry is managed through tree merge/split operations

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
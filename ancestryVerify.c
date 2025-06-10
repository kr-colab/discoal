#include <stdio.h>
#include <stdlib.h>
#include "ancestryVerify.h"

extern rootedNode **allNodes;
extern int totNodeNumber;

int verifyAncestryConsistency(rootedNode *node, int nSites) {
    if (!node) return 1;
    if (!node->ancestryRoot && !node->ancSites) return 1; // Both NULL is OK
    
    if (!node->ancestryRoot || !node->ancSites) {
        fprintf(stderr, "ERROR: Node %d has inconsistent ancestry representations\n", node->id);
        return 0;
    }
    
    // Check each site
    for (int i = 0; i < nSites; i++) {
        uint16_t treeCount = getAncestryCount(node->ancestryRoot, i);
        uint16_t arrayCount = node->ancSites[i];
        
        if (treeCount != arrayCount) {
            fprintf(stderr, "ERROR: Node %d site %d mismatch: tree=%d, array=%d\n", 
                    node->id, i, treeCount, arrayCount);
            return 0;
        }
    }
    
    return 1;
}

void verifyAllNodes(int nSites) {
    int errors = 0;
    
    for (int i = 0; i < totNodeNumber; i++) {
        if (allNodes[i] && !verifyAncestryConsistency(allNodes[i], nSites)) {
            errors++;
        }
    }
    
    if (errors > 0) {
        fprintf(stderr, "VERIFICATION FAILED: %d nodes have inconsistent ancestry\n", errors);
        exit(1);
    }
}
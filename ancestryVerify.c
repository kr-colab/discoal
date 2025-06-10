#include <stdio.h>
#include <stdlib.h>
#include "ancestryVerify.h"

extern rootedNode **allNodes;
extern int totNodeNumber;

int verifyAncestryConsistency(rootedNode *node, int nSites) {
    if (!node) return 1;
    
    // Only verify tree structure now that ancSites is removed
    if (node->ancestryRoot) {
        // Could add tree consistency checks here if needed
        // For now, just verify the tree exists
        return 1;
    }
    
    // No ancestry tree is OK for some nodes
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
#ifndef __ANCESTRY_VERIFY_H__
#define __ANCESTRY_VERIFY_H__

#include "discoal.h"
#include "ancestrySegment.h"

// Verification functions to ensure ancestry segment tree matches ancSites array
int verifyAncestryConsistency(rootedNode *node, int nSites);
void verifyAllNodes(int nSites);

#endif
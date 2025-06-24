#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "ancestryWrapper.h"
#include "activeSegment.h"
#include <time.h>
#include "discoal.h"
#include "discoalFunctions.h"
#include "ranlib.h"
#include "alleleTraj.h"
#include "tskitInterface.h"


// Initial capacity for breakPoints array
#define INITIAL_BREAKPOINTS_CAPACITY 1000

void initializeBreakPoints() {
	// Clean up any existing allocation first
	if (breakPoints != NULL) {
		free(breakPoints);
		breakPoints = NULL;
	}
	
	breakPointsCapacity = INITIAL_BREAKPOINTS_CAPACITY;
	breakPoints = malloc(sizeof(int) * breakPointsCapacity);
	if (breakPoints == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for breakPoints array\n");
		exit(1);
	}
	breakPoints[0] = 666; // Initialize with original marker value
	breakNumber = 0;
}

void ensureBreakPointsCapacity() {
	if (breakNumber >= breakPointsCapacity) {
		int newCapacity = breakPointsCapacity * 2;
		int *newBreakPoints = realloc(breakPoints, sizeof(int) * newCapacity);
		if (newBreakPoints == NULL) {
			fprintf(stderr, "Error: Failed to reallocate memory for breakPoints array (capacity: %d -> %d)\n", 
					breakPointsCapacity, newCapacity);
			exit(1);
		}
		breakPoints = newBreakPoints;
		breakPointsCapacity = newCapacity;
	}
}

void cleanupBreakPoints() {
	if (breakPoints != NULL) {
		free(breakPoints);
		breakPoints = NULL;
		breakPointsCapacity = 0;
		breakNumber = 0;
	}
}

void addBreakPoint(int bp) {
	ensureBreakPointsCapacity();
	breakPoints[breakNumber] = bp;
	breakNumber += 1;
}

// Population list management functions for O(1) node selection

void initializePopLists() {
	for (int p = 0; p < MAXPOPS; p++) {
		popLists[p].nodes = NULL;
		popLists[p].count = 0;
		popLists[p].capacity = 0;
	}
}

void cleanupPopLists() {
	for (int p = 0; p < MAXPOPS; p++) {
		if (popLists[p].nodes != NULL) {
			free(popLists[p].nodes);
			popLists[p].nodes = NULL;
			popLists[p].count = 0;
			popLists[p].capacity = 0;
		}
	}
}

void addNodeToPopList(rootedNode *node, int popn) {
	if (popn < 0 || popn >= MAXPOPS) return;  // Handle negative/invalid populations
	
	PopulationNodeList *list = &popLists[popn];
	
	// Grow array if needed
	if (list->count >= list->capacity) {
		int newCap = list->capacity ? list->capacity * 2 : 16;
		rootedNode **newNodes = realloc(list->nodes, newCap * sizeof(rootedNode*));
		if (newNodes == NULL) {
			fprintf(stderr, "Error: Failed to allocate memory for population list\n");
			exit(1);
		}
		list->nodes = newNodes;
		list->capacity = newCap;
	}
	
	node->popListIndex = list->count;
	list->nodes[list->count++] = node;
}

void removeNodeFromPopList(rootedNode *node) {
	int popn = node->population;
	if (popn < 0 || popn >= MAXPOPS) return;  // Handle negative/invalid populations
	
	PopulationNodeList *list = &popLists[popn];
	int index = node->popListIndex;
	
	if (index < 0 || index >= list->count) return;
	
	// Swap with last element
	if (index < list->count - 1) {
		list->nodes[index] = list->nodes[list->count - 1];
		list->nodes[index]->popListIndex = index;
	}
	
	list->count--;
	node->popListIndex = -1;
}

// Fast O(1) node selection from population
rootedNode *pickNodePopnFast(int popn) {
	if (popn < 0 || popn >= MAXPOPS) {
		fprintf(stderr, "Error: Invalid population %d in pickNodePopnFast\n", popn);
		exit(1);
	}
	
	PopulationNodeList *list = &popLists[popn];
	if (list->count == 0) {
		fprintf(stderr, "Error: No nodes in population %d\n", popn);
		exit(1);
	}
	
	int index = ignuin(0, list->count - 1);
	return list->nodes[index];
}

void initialize(){
	int i,j,p, count=0;
	int leafID=0;
	int tmpCount = 0;
	
	/* Initialize the arrays */
	totChunkNumber = 0;
	initializeBreakPoints();
	initializeNodeArrays();
	initializePopLists();  // Initialize population lists for O(1) node selection
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			
			// Always track sample node ID in tskit (create node first)
			tsk_id_t tsk_id = tskit_add_node(0.0, p, 1);  // time=0, is_sample=1
			set_tskit_node_id(nodes[count], tsk_id);
			addSampleNodeId(tsk_id);
			
			// Initialize ancestry segment tree for leaf node with tskit node ID
			nodes[count]->ancestryRoot = newSegment(0, nSites, tsk_id, NULL, NULL);

			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;
			nodes[count]->inActiveSet = 1;  // Sample nodes start in active set
			
			// Add to population list
			addNodeToPopList(nodes[count], p);
			
			count += 1;
		}
	}


	breakNumber = 0;
	// Initialize segment-based active material (all sites start as active)
	initializeActiveMaterial(&activeMaterialSegments, nSites);
	
	//set initial counts
	alleleNumber = sampleSize;
	totNodeNumber = sampleSize;
	freedNodeCount = 0;  // Initialize freed node counter
	//ancient pop samples?
	if(ancSampleFlag == 1){
		//go through ancient sample events
		for(i = 0; i < eventNumber; i++){
			if (events[i].type == 'A'){
				j=0;
				while(nodes[j]->population != events[i].popID){
					j++;
				}
				for(tmpCount=0;tmpCount<events[i].lineageNumber;tmpCount++){
					// Remove from population list before changing population
					removeNodeFromPopList(nodes[tmpCount+j]);
					nodes[tmpCount+j]->population = (events[i].popID + 1) * -1;
					popnSizes[events[i].popID]--;
				}
			}
		}

	}
	
	activeSites = nSites;
	if (npops>1){
		if(tDiv==666 && migFlag == 0){
			fprintf(stderr,"tDiv or migration not set in population split model\n");
			exit(1);
		}
		//initialize migration matrix
		for(i=0;i<npops;i++){
			for(j=0;j<npops;j++){
				migMat[i][j]=migMatConst[i][j];
			}
		}
		eventFlag = 0;
	}
	
	//priors on parameters?
        if(priorUA==1)
          uA = genunf(pUALow,pUAUp);
	if(priorTheta==1)
	  theta = genunf(pThetaLow,pThetaUp);
	if(priorC==1)
	  partialSweepFinalFreq = genunf(pCLow,pCUp);
	if(priorRho==1){
          rho = genunf(pRhoLow,pRhoUp);
        }
        else if(priorRho==2){
	  rho = genexp(pRhoMean);
          if (rho > pRhoUp)
              rho = pRhoUp;
        }
        if(gammaCoRatioMode==1)
            my_gamma = rho*gammaCoRatio;
	if(priorAlpha==1)
	  alpha = genunf(pAlphaLow,pAlphaUp);
	if(priorX==1)
	  sweepSite = genunf(pXLow,pXUp);
	if(priorF0==1)
	  f0 = genunf(pF0Low,pF0Up);
	if(priorTau==1){
	  tau = genunf(pTauLow,pTauUp);
	  for(i=0;i<eventNumber;i++){
	    if(events[i].type=='s')
	      events[i].time=tau;
	  }
	}
	if(priorE1==1){
		events[1].time = genunf(pE1TLow,pE1THigh);
		events[1].popnSize = genunf(pE1SLow,pE1SHigh);
	}
	if(priorE2==1){
		events[2].time = genunf(pE2TLow,pE2THigh);
		events[2].popnSize = genunf(pE2SLow,pE2SHigh);
	}
	sortEventArray(events,eventNumber);
	
}

void initializeTwoSite(){
	int i,j,p, count=0;
	int leafID=0;
	/* initialize the arrays */
	totChunkNumber = 0;
	initializeBreakPoints();
	initializeNodeArrays();
	initializePopLists();  // Initialize population lists for O(1) node selection
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			
			// Track sample node ID in tskit (create node first)
			tsk_id_t tsk_id = tskit_add_node(0.0, p, 1);  // time=0, is_sample=1
			set_tskit_node_id(nodes[count], tsk_id);
			addSampleNodeId(tsk_id);
			
			// Initialize ancestry segment tree for leaf node with tskit node ID
			nodes[count]->ancestryRoot = newSegment(0, nSites, tsk_id, NULL, NULL);
			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;
			nodes[count]->inActiveSet = 1;  // Sample nodes start in active set
			
			// Add to population list
			addNodeToPopList(nodes[count], p);
			
			count += 1;
		}
	}
	
	breakNumber = 0;
	// Initialize segment-based active material (all sites start as active)
	initializeActiveMaterial(&activeMaterialSegments, nSites);
	
	//set initial counts
	alleleNumber = sampleSize;
	totNodeNumber = sampleSize;
	freedNodeCount = 0;  // Initialize freed node counter
	activeSites = nSites;
	if (npops>1){
		if(tDiv==666 && mig[0] == 0 && mig[1] == 0){
			fprintf(stderr,"tDiv not set in population split model\n");
			exit(1);
		}
		
		eventFlag = 0;
	}
	
	//priors on parameters?
        if(priorUA==1)
          uA = genunf(pUALow,pUAUp);
	if(priorTheta==1)
	  theta = genunf(pThetaLow,pThetaUp);
	if(priorRho==1)
          rho = genunf(pRhoLow,pRhoUp);
        else if(priorRho==2){
	  rho = genexp(pRhoMean);
	  if (rho > pRhoUp)
	    rho = pRhoUp;
	}
        if(gammaCoRatioMode==1)
            my_gamma = rho*gammaCoRatio;
	if(priorAlpha==1)
	  alpha = genunf(pAlphaLow,pAlphaUp);
	if(priorX==1)
	  sweepSite = genunf(pXLow,pXUp);
	if(priorF0==1)
	  f0 = genunf(pF0Low,pF0Up);
	if(priorTau==1){
	  tau = genunf(pTauLow,pTauUp);
	  for(i=0;i<eventNumber;i++){
	    if(events[i].type=='s')
	      events[i].time=tau;
	  }
	}
	if(priorE1==1){
		events[1].time = genunf(pE1TLow,pE1THigh);
		events[1].popnSize = genunf(pE1SLow,pE1SHigh);
	}
	if(priorE2==1){
		events[2].time = genunf(pE2TLow,pE2THigh);
		events[2].popnSize = genunf(pE2SLow,pE2SHigh);
	}
	sortEventArray(events,eventNumber);
}


static int total_nodes_allocated = 0;

rootedNode *newRootedNode(double cTime, int popn) {
	rootedNode *temp;

	temp =  malloc(sizeof(rootedNode));
	if(temp==NULL){
		printf("fucked: out o' memory from rootedNode alloc\n");
	}
	total_nodes_allocated++;
	// fprintf(stderr, "DEBUG: Allocated node %p, total allocated: %d\n", temp, total_nodes_allocated);
	temp->rightParent = NULL;
	temp->leftParent = NULL;
	temp->rightChild = NULL;
	temp->leftChild = NULL; 
	temp->time = cTime;
	temp->branchLength=0.0;
	temp->population = popn;
	temp->sweepPopn = -1;
	
	// Initialize node state tracking
	temp->parentsRecorded = 0;
	temp->isFullyRecorded = 0;
	temp->inActiveSet = 0;  // Will be set to 1 for sample nodes
	temp->carriesSweepMutation = 0;  // No sweep mutation by default
	
	// Initialize ancestry segment tree to NULL (will be set during initialization)
	temp->ancestryRoot = NULL;

	// Initialize tskit node ID
	temp->tskit_node_id = TSK_NULL;
	
	// Initialize population list index
	temp->popListIndex = -1;

	return temp;
}

// Forward declaration
void addToRemovedList(rootedNode *node);

// Mark that a parent has recorded its connection to this child
void markParentRecorded(rootedNode *child, rootedNode *parent) {
	if (child == NULL) return;
	
	if (parent == child->leftParent) {
		child->parentsRecorded |= 0x01;  // Set bit 0
	} else if (parent == child->rightParent) {
		child->parentsRecorded |= 0x02;  // Set bit 1
	}
	
	// Check if both parents recorded (or if single parent for coalescence)
	int expectedParents = (child->leftParent != NULL) + (child->rightParent != NULL);
	int recordedParents = ((child->parentsRecorded & 0x01) ? 1 : 0) + 
	                     ((child->parentsRecorded & 0x02) ? 1 : 0);
	
	// Only mark as fully recorded if:
	// 1. All expected parents are recorded
	// 2. Node is not in active set
	// 3. All segments on the node are recorded
	if (recordedParents == expectedParents && !child->inActiveSet) {
		// Check if all segments are recorded
		if (areAllSegmentsRecorded(child->ancestryRoot)) {
			child->isFullyRecorded = 1;
			// If node is fully recorded and not in active set, add to removed list
			addToRemovedList(child);
		}
	}
}

// Try to free a node if it's safe to do so
void tryFreeNode(rootedNode *node) {
	if (node == NULL || !node->isFullyRecorded) return;
	
	// Critical safety check: ensure not in active set
	if (node->inActiveSet) return;
	
	// Additional safety check: ensure not in nodes array
	for (int i = 0; i < alleleNumber; i++) {
		if (nodes[i] == node) return;  // Still active!
	}
	
	// Free ancestry segments using reference counting
	if (node->ancestryRoot) {
		freeSegmentTree(node->ancestryRoot);
		node->ancestryRoot = NULL;
	}
	
	// Free the node
	free(node);
	
	// Update statistics
	freedNodeCount++;
}

// Global list to track nodes that have been removed from active set but not yet freed
// This allows us to safely free them later when we're sure they're not being used
#define MAX_REMOVED_NODES 10000
static rootedNode *removedNodes[MAX_REMOVED_NODES];
static int removedNodeCount = 0;

// Add a node to the removed list for later freeing
void addToRemovedList(rootedNode *node) {
	if (node == NULL || removedNodeCount >= MAX_REMOVED_NODES) return;
	
	// Check if already in removed list
	for (int i = 0; i < removedNodeCount; i++) {
		if (removedNodes[i] == node) return;  // Already tracked
	}
	
	// Mark as not in active set
	node->inActiveSet = 0;
	
	// Add to removed list
	removedNodes[removedNodeCount++] = node;
}

// Sweep through removed nodes and free any that are fully recorded
// This implements a mark-and-sweep approach for safe memory management
int sweepAndFreeRemovedNodes() {
	int freedCount = 0;
	int newRemovedCount = 0;
	
	// Go through removed nodes and free those that are marked as fully recorded
	for (int i = 0; i < removedNodeCount; i++) {
		rootedNode *node = removedNodes[i];
		if (node == NULL) continue;
		
		// Only process nodes that are marked as fully recorded
		// This marking happens in markParentRecorded when all edges are recorded
		if (node->isFullyRecorded) {
			// Extra safety check: ensure it's really not in the active array
			int inActiveArray = 0;
			for (int j = 0; j < alleleNumber; j++) {
				if (nodes[j] == node) {
					inActiveArray = 1;
					fprintf(stderr, "WARNING: Fully recorded node still in active array!\n");
					break;
				}
			}
			
			if (!inActiveArray) {
				// Safe to free - node is fully recorded and not active
				if (node->ancestryRoot) {
					freeSegmentTree(node->ancestryRoot);
					node->ancestryRoot = NULL;
				}
				
				// NULL out any parent pointers to this node before freeing
				if (node->leftParent) {
					if (node->leftParent->leftChild == node) node->leftParent->leftChild = NULL;
					if (node->leftParent->rightChild == node) node->leftParent->rightChild = NULL;
				}
				if (node->rightParent) {
					if (node->rightParent->leftChild == node) node->rightParent->leftChild = NULL;
					if (node->rightParent->rightChild == node) node->rightParent->rightChild = NULL;
				}
				
				// NULL out any child pointers from this node (defensive)
				if (node->leftChild) {
					if (node->leftChild->leftParent == node) node->leftChild->leftParent = NULL;
					if (node->leftChild->rightParent == node) node->leftChild->rightParent = NULL;
				}
				if (node->rightChild) {
					if (node->rightChild->leftParent == node) node->rightChild->leftParent = NULL;
					if (node->rightChild->rightParent == node) node->rightChild->rightParent = NULL;
				}
				
				// Debug: print info about node being freed
				// fprintf(stderr, "Freeing node %p (time=%f, pop=%d)\n", 
				//         node, node->time, node->population);
				free(node);
				freedCount++;
				freedNodeCount++;
			} else {
				// This shouldn't happen - keep for debugging
				removedNodes[newRemovedCount++] = node;
			}
		} else {
			// Not fully recorded yet - this shouldn't happen either
			// since we only add fully recorded nodes to the list
			removedNodes[newRemovedCount++] = node;
		}
	}
	
	// Update the removed count
	removedNodeCount = newRemovedCount;
	
	return freedCount;
}

////////  Migration Stuff ///////////

//migrateAtTime -- this only switches between two populations
void migrateAtTime(double cTime,int srcPopn, int destPopn){
	// Check if source population has any nodes
	if (popLists[srcPopn].count == 0) {
		// No nodes to migrate, this can happen in complex demographic scenarios
		return;
	}
	
	rootedNode *temp;
	temp = pickNodePopn(srcPopn);
	
	// Remove from source population list
	removeNodeFromPopList(temp);
	
	temp->population = destPopn;
	
	// Add to destination population list
	addNodeToPopList(temp, destPopn);

	popnSizes[srcPopn]-=1;
	popnSizes[destPopn]+=1;
	if (srcPopn==0)
		sweepPopnSizes[temp->sweepPopn]-=1;

}

//recurrentMutAtTime -- this adds a non-sweeping node to the sweep or vice-versa
void recurrentMutAtTime(double cTime,int srcPopn, int sp){
        rootedNode *temp;

        temp = pickNodePopnSweep(srcPopn,sp);
        temp->sweepPopn = (sp+1)%2;
        sweepPopnSizes[sp]-=1;
        sweepPopnSizes[temp->sweepPopn]+=1;
}



// migrateExceptSite -- does migration unless specific site is present in allele picked
void migrateExceptSite(double site, double scalar, int srcPopn, int destPopn){
	// Check if source population has any nodes
	if (popLists[srcPopn].count == 0) {
		return;
	}
	
	rootedNode *temp;
	int doMig=0;
	temp = pickNodePopn(srcPopn);
	
	if(isAncestralHere(temp,site) == 0 || ranf() < scalar) doMig = 1; //migration of site is 10% rate of other sites
	if(doMig == 1)
	{
		// Remove from source population list
		removeNodeFromPopList(temp);
		
		temp->population = destPopn;
		
		// Add to destination population list
		addNodeToPopList(temp, destPopn);
		
		popnSizes[srcPopn]-=1;
		popnSizes[destPopn]+=1;
		if (srcPopn==0)
			sweepPopnSizes[temp->sweepPopn]-=1;
	}
	
	
}


////////

// Structure to hold a node and its ancestry segments for edge recording
typedef struct {
    rootedNode *node;
    AncestrySegment *segments;
} NodeWithSegments;

// Helper function to find descendant nodes by tracing through unrecorded nodes
static int findTrueDescendants(rootedNode *node, rootedNode **results, int maxResults) {
    int count = 0;
    
    // If this node has a tskit ID, it's a true descendant
    if (get_tskit_node_id(node) != TSK_NULL) {
        if (count < maxResults) {
            results[count] = node;
            count++;
        }
        return count;
    }
    
    // This node doesn't have a tskit ID - trace through its children
    if (node->leftChild != NULL) {
        int childCount = findTrueDescendants(node->leftChild, results + count, maxResults - count);
        count += childCount;
    }
    
    if (node->rightChild != NULL && node->rightChild != node->leftChild && count < maxResults) {
        int childCount = findTrueDescendants(node->rightChild, results + count, maxResults - count);
        count += childCount;
    }
    
    return count;
}

// Helper function to find descendants with tskit nodes by tracing through unrecorded nodes
// This properly handles genomic intervals through recombination events
static int findDescendantsWithSegments(rootedNode *node, NodeWithSegments *results, int maxResults) {
    int count = 0;
    
    // If this node has no ancestry, there's nothing to trace
    if (node->ancestryRoot == NULL) {
        return 0;
    }
    
    // If this node has a tskit ID, we can record edges to it
    if (get_tskit_node_id(node) != TSK_NULL) {
        if (count < maxResults) {
            results[count].node = node;
            results[count].segments = node->ancestryRoot;
            count++;
        }
        return count;
    }
    
    // This node doesn't have a tskit ID (unrecorded recombination node)
    // Trace through its children
    if (node->leftChild != NULL) {
        // For each child, we need to track which segments flow through
        NodeWithSegments childResults[10];
        int childCount = findDescendantsWithSegments(node->leftChild, childResults, 10);
        
        // Pass through the results from the recursive call
        for (int i = 0; i < childCount && count < maxResults; i++) {
            results[count] = childResults[i];
            count++;
        }
    }
    
    if (node->rightChild != NULL && node->rightChild != node->leftChild) {
        // Handle second child (for recombination, leftChild == rightChild)
        NodeWithSegments childResults[10];
        int childCount = findDescendantsWithSegments(node->rightChild, childResults, 10);
        
        for (int i = 0; i < childCount && count < maxResults; i++) {
            results[count] = childResults[i];
            count++;
        }
    }
    
    return count;
}

/**
 * Records edges from a parent node to its children based on ancestry segments.
 * Handles both full ARG mode and minimal tree sequence mode.
 * 
 * @param parent_tsk_id The tskit ID of the parent node
 * @param child The child node whose segments should be recorded
 * @param child_tsk_id The tskit ID of the child node (used in full ARG mode)
 * @param isMinimalMode Whether we're in minimal tree sequence mode
 */
static void recordEdgesForChild(tsk_id_t parent_tsk_id, rootedNode *child, tsk_id_t child_tsk_id, int isMinimalMode) {
    if (child->ancestryRoot == NULL) {
        return;
    }
    
    AncestrySegment *seg = child->ancestryRoot;
    
    if (!isMinimalMode) {
        // Full ARG mode: record edges from parent to immediate child
        if (child_tsk_id != TSK_NULL) {
            while (seg != NULL) {
                tskit_add_edges(parent_tsk_id, child_tsk_id, 
                    (double)seg->start, (double)seg->end);
                markSegmentRecorded(seg);
                seg = seg->next;
            }
        }
    } else {
        // Minimal mode: use segment IDs to skip recombination nodes
        while (seg != NULL) {
            // Each segment knows which tskit node it represents
            if (seg->tskit_node_id != TSK_NULL && seg->tskit_node_id != parent_tsk_id) {
                tskit_add_edges(parent_tsk_id, seg->tskit_node_id, 
                    (double)seg->start, (double)seg->end);
            }
            markSegmentRecorded(seg);
            seg = seg->next;
        }
    }
}

/**
 * Records all edges from a parent to its children during a coalescent event.
 * Handles child ordering requirements and marks all segments as recorded.
 * 
 * @param parent The parent node
 * @param lChild The left child node
 * @param rChild The right child node
 * @param parent_tsk_id The tskit ID of the parent
 * @param lchild_tsk_id The tskit ID of the left child
 * @param rchild_tsk_id The tskit ID of the right child
 */
static void recordCoalescentEdges(rootedNode *parent, rootedNode *lChild, rootedNode *rChild,
                                  tsk_id_t parent_tsk_id, tsk_id_t lchild_tsk_id, tsk_id_t rchild_tsk_id) {
    if (parent_tsk_id == TSK_NULL) {
        return;
    }
    
    extern int minimalTreeSeq;
    
    // Record edges for both children
    recordEdgesForChild(parent_tsk_id, lChild, lchild_tsk_id, minimalTreeSeq);
    recordEdgesForChild(parent_tsk_id, rChild, rchild_tsk_id, minimalTreeSeq);
    
    // CRITICAL: Mark the parent's segments as recorded too!
    // The parent has new segments from mergeAncestryTrees that need to be marked
    // otherwise the parent node can never be freed
    if (parent->ancestryRoot != NULL) {
        AncestrySegment *seg = parent->ancestryRoot;
        while (seg != NULL) {
            markSegmentRecorded(seg);
            seg = seg->next;
        }
    }
}

/**
 * Records edges for recombination and gene conversion events.
 * In these events, a single child has two parents, each getting different segments.
 * Only used in full ARG mode (not minimal tree sequence mode).
 * 
 * @param child The child node
 * @param lParent The left parent node (gets segments left of crossover/converted tract)
 * @param rParent The right parent node (gets segments right of crossover/unconverted segments)
 * @param child_tsk_id The tskit ID of the child
 * @param lparent_tsk_id The tskit ID of the left parent
 * @param rparent_tsk_id The tskit ID of the right parent
 */
static void recordRecombinationEdges(rootedNode *child, rootedNode *lParent, rootedNode *rParent,
                                     tsk_id_t child_tsk_id, tsk_id_t lparent_tsk_id, tsk_id_t rparent_tsk_id) {
    extern int minimalTreeSeq;
    if (minimalTreeSeq) {
        // In minimal mode, don't create tskit nodes for recombination parents
        // They will remain with tskit_node_id = TSK_NULL
        // During coalescence, we'll trace through to find the actual samples
        return;
    }
    
    // Record edges based on actual ancestry segments
    // Left parent gets segments (left of crossover for recombination, converted tract for gene conversion)
    if (lparent_tsk_id != TSK_NULL && child_tsk_id != TSK_NULL && lParent->ancestryRoot != NULL) {
        AncestrySegment *seg = lParent->ancestryRoot;
        while (seg != NULL) {
            tskit_add_edges(lparent_tsk_id, child_tsk_id,
                (double)seg->start, (double)seg->end);
            markSegmentRecorded(seg);
            seg = seg->next;
        }
    }
    
    // Right parent gets segments (right of crossover for recombination, unconverted segments for gene conversion)
    if (rparent_tsk_id != TSK_NULL && child_tsk_id != TSK_NULL && rParent->ancestryRoot != NULL) {
        AncestrySegment *seg = rParent->ancestryRoot;
        while (seg != NULL) {
            tskit_add_edges(rparent_tsk_id, child_tsk_id,
                (double)seg->start, (double)seg->end);
            markSegmentRecorded(seg);
            seg = seg->next;
        }
    }
}

void coalesceAtTimePopn(double cTime, int popn){
	// Check if population has at least 2 nodes for coalescence
	if (popLists[popn].count < 2) {
		// Cannot coalesce with fewer than 2 nodes
		return;
	}
	
	rootedNode *temp, *lChild, *rChild;
	int i;
	
	temp = newRootedNode(cTime,popn);

	lChild = pickNodePopn(popn);
	temp->leftChild = lChild;
	lChild->leftParent = temp;
	lChild->branchLength = cTime - lChild->time;
	removeNode(lChild);

	rChild =  pickNodePopn(popn);
	temp->rightChild = rChild;
	rChild->leftParent = temp;
	rChild->branchLength = cTime - rChild->time;
	removeNode(rChild);
	
	//deal with ancMaterial
	temp->nancSites = 0;
	temp->lLim = nSites;
	temp->rLim = 0;
	
	// Propagate sweep mutation flag if either child carries it
	if (lChild->carriesSweepMutation || rChild->carriesSweepMutation) {
		temp->carriesSweepMutation = 1;
	}
	
	// Add the parent node to discoal's arrays and tskit BEFORE merging/recording
	addNode(temp);
	
	// Get tskit IDs for all nodes involved in coalescence  
	tsk_id_t parent_tsk_id = get_tskit_node_id(temp);
	tsk_id_t lchild_tsk_id = get_tskit_node_id(lChild);
	tsk_id_t rchild_tsk_id = get_tskit_node_id(rChild);
	
	// Merge ancestry segment trees with parent's tskit node ID
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot, parent_tsk_id);
	
	// Update stats from tree when in tree-only mode
	updateAncestryStatsFromTree(temp);
	
	extern int minimalTreeSeq;
	
	// Record edges using the new convenience function
	recordCoalescentEdges(temp, lChild, rChild, parent_tsk_id, lchild_tsk_id, rchild_tsk_id);
	
	
	// Mark parent recording and try to free nodes if using tskit
	if (tskitOutputMode || minimalTreeSeq) {
		// Use the same first/second child logic as edge recording
		rootedNode *first_child, *second_child;
		if (lchild_tsk_id < rchild_tsk_id) {
			first_child = lChild;
			second_child = rChild;
		} else {
			first_child = rChild;
			second_child = lChild;
		}
		
		// For first child - mark parent recorded
		if (minimalTreeSeq && first_child->ancestryRoot != NULL) {
			rootedNode *trueDescendants[10];
			int numDesc = findTrueDescendants(first_child, trueDescendants, 10);
			
			if (numDesc > 0 && trueDescendants[0] != first_child) {
				// Mark all nodes in the bypass chain
				rootedNode *current = first_child;
				while (current != NULL && current != trueDescendants[0]) {
					markParentRecorded(current, temp);
					rootedNode *next = (current->leftChild != NULL) ? 
					                  current->leftChild : current->rightChild;
					// Temporarily disabled node freeing
					// if (current != first_child) {
					//     tryFreeNode(current);
					// }
					current = next;
				}
				// Mark the true descendants
				for (int d = 0; d < numDesc; d++) {
					markParentRecorded(trueDescendants[d], temp);
				}
			} else {
				markParentRecorded(first_child, temp);
			}
		} else {
			markParentRecorded(first_child, temp);
		}
		
		// For second child - mark parent recorded
		if (minimalTreeSeq && second_child->ancestryRoot != NULL) {
			rootedNode *trueDescendants[10];
			int numDesc = findTrueDescendants(second_child, trueDescendants, 10);
			
			if (numDesc > 0 && trueDescendants[0] != second_child) {
				// Mark all nodes in the bypass chain
				rootedNode *current = second_child;
				while (current != NULL && current != trueDescendants[0]) {
					markParentRecorded(current, temp);
					rootedNode *next = (current->leftChild != NULL) ? 
					                  current->leftChild : current->rightChild;
					// Temporarily disabled node freeing
					// if (current != second_child) {
					//     tryFreeNode(current);
					// }
					current = next;
				}
				// Mark the true descendants
				for (int d = 0; d < numDesc; d++) {
					markParentRecorded(trueDescendants[d], temp);
				}
			} else {
				markParentRecorded(second_child, temp);
			}
		} else {
			markParentRecorded(second_child, temp);
		}
		
		// Try to free children based on their node type
		// A coalescent node has both left and right children
		int lChildIsCoal = isCoalNode(lChild);
		int rChildIsCoal = isCoalNode(rChild);
		
		// Check if nodes are leaves (no children at all)
		int lChildIsLeaf = (lChild->leftChild == NULL && lChild->rightChild == NULL);
		int rChildIsLeaf = (rChild->leftChild == NULL && rChild->rightChild == NULL);
		
		// Leaf nodes need special handling since they have no children to record
		if (lChildIsLeaf) {
			// Mark all segments as recorded since leaf nodes have simple ancestry
			if (lChild->ancestryRoot) {
				markSegmentRecorded(lChild->ancestryRoot);
			}
			lChild->isFullyRecorded = 1;
			addToRemovedList(lChild);
		}
		
		if (rChildIsLeaf) {
			// Mark all segments as recorded since leaf nodes have simple ancestry
			if (rChild->ancestryRoot) {
				markSegmentRecorded(rChild->ancestryRoot);
			}
			rChild->isFullyRecorded = 1;
			addToRemovedList(rChild);
		}
		
		// With periodic sweep approach, we don't need immediate freeing
		// Just mark nodes as removed and let the sweep handle it
	}
	
	//update active anc. material (addNode was moved earlier)
	updateActiveMaterial(temp);
	//popnSizes[popn]--; //decrese popnSize	
}

/*
void coalesceAtTimePopnTrackLeafs(double cTime, int popn){
	rootedNode *temp, *lChild, *rChild;
	int i;

	temp = newRootedNode(cTime,popn);
	temp->leafs = calloc(sizeof(int) * sampleSize);
	lChild = pickNodePopn(popn);
	temp->leftChild = lChild;
	lChild->leftParent = temp;
	lChild->branchLength = cTime - lChild->time;
	removeNode(lChild);

	rChild =  pickNodePopn(popn);
	temp->rightChild = rChild;
	rChild->leftParent = temp;
	rChild->branchLength = cTime - rChild->time;
	removeNode(rChild);
	
	//deal with ancMaterial
	temp->nancSites = 0;
	temp->lLim = nSites;
	temp->rLim = 0;
	// Merge ancestry segment trees
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot);
	//deal with leafs
	for(i=0;i<sampleSize;i++){
		temp->leafs[i] =  lChild->leafs[i] + rChild->leafs[i];
	}
	
	addNode(temp);
	//update active anc. material
	updateActiveMaterial(temp);
	//popnSizes[popn]--; //decrese popnSize	
}*/


// Helper function to update nancSites, lLim, rLim from ancestry tree
void updateAncestryStatsFromTree(rootedNode *node) {
	if (!node || !node->ancestryRoot) return;
	
	node->nancSites = 0;
	node->lLim = nSites;
	node->rLim = 0;
	
	// Walk the segment list instead of querying every site
	AncestrySegment *seg = node->ancestryRoot;
	while (seg) {
		// Check if this segment is polymorphic (has ancestry but not fixed)
		if (seg->count > 0 && seg->count < sampleSize) {
			// Add all sites in this polymorphic segment
			node->nancSites += (seg->end - seg->start);
			
			// Update left boundary
			if (seg->start < node->lLim) {
				node->lLim = seg->start;
			}
			// Update right boundary (end is exclusive, so last site is end-1)
			if (seg->end - 1 > node->rLim) {
				node->rLim = seg->end - 1;
			}
		}
		seg = seg->next;
	}
}

/*updateActiveMaterial- does bookkeeping on the material that has found it's MRCA-- added for efficiency */
void updateActiveMaterial(rootedNode *aNode){
	// Use segment-based approach to update active material
	// This removes regions where the node has ancestry from all samples
	updateActiveMaterialFromAncestry(&activeMaterialSegments, aNode->ancestryRoot, 
	                                sampleSize, nSites);
	
	// Update global active site count
	activeSites = getActiveSiteCount(&activeMaterialSegments);
	
	#ifdef DEBUG_ACTIVE_MATERIAL
	fprintf(stderr, "DEBUG updateActiveMaterial: node at time %f, %d sites still active\n", 
	        aNode->time, activeSites);
	printActiveSegments(&activeMaterialSegments);
	#endif
	
	#ifdef DEBUG_TREE_ONLY
	if (activeSites == 0) {
		fprintf(stderr, "DEBUG: All sites have found MRCA\n");
	}
	#endif
}

/* asks whether a potential breakpoint is in material that hasn't found MRCA yet*/
int isActive(int site){
	return isActiveSite(&activeMaterialSegments, site);
}

int siteBetweenChunks(rootedNode *aNode, int xOverSite){

	if(aNode->rLim < aNode->lLim)return(0);
	if (aNode->lLim < xOverSite && aNode->rLim >= xOverSite){
		return 1;
	}
	else{
		return 0;
	}
}

int isAncestralHere(rootedNode *aNode, float site){
	int bp;
	bp = floor(site * nSites);
	return isPolymorphicAt(aNode, bp, sampleSize);
}

int hasMaterialHere(rootedNode *aNode, float site){
	int bp;
	bp = floor(site * nSites);
	return hasAncestryAt(aNode, bp);
}

int nAncestorsHere(rootedNode *aNode, float site){
	int bp;
	bp = floor(site * nSites);
	return getAncestryAt(aNode, bp);
}

int isLeaf(rootedNode *aNode){
	if(aNode->leftChild == NULL && aNode->rightChild == NULL)
		return(1);
	else
		return(0);
}

int isCoalNode(rootedNode *aNode){
	if(aNode->leftChild != NULL && aNode->rightChild != NULL)
		return(1);
	else
		return(0);
}

int recombineAtTimePopn(double cTime, int popn){
	rootedNode *aNode, *lParent, *rParent;
	int i;
	int xOver;
	
	aNode = pickNodePopn(popn);
//	printf("recombine\n");
	xOver = ignuin(0,nSites-1);
//	printf("xo: %d test1: %d test2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver) );
//	printNode(aNode);
//	for(i=0;i<nSites;i++)printf("%d",activeMaterial[i]);
//	printf("\n");
//	printf("there xover: %d t1: %d t2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver));
		
	if (siteBetweenChunks(aNode, xOver) == 1 &&  isActive(xOver) == 1  ){
		// Early exit if no ancestral material
		if(aNode->nancSites == 0 || aNode->lLim > aNode->rLim) {
			return xOver;
		}
		
		removeNode(aNode); 
		lParent = newRootedNode(cTime,popn);
		rParent = newRootedNode(cTime,popn);
		aNode->leftParent = lParent;
		aNode->rightParent = rParent;
		aNode->branchLength = cTime - aNode->time;
		lParent->leftChild = aNode;
		rParent->leftChild = aNode;
		lParent->population = aNode->population;
		rParent->population = aNode->population;
		
		
		// Split ancestry segment tree at crossover point
		extern int minimalTreeSeq;
		if (!minimalTreeSeq) {
			// Full ARG mode: use parent's tskit_node_id in segments
			lParent->ancestryRoot = splitLeftWithParent(aNode->ancestryRoot, xOver, lParent->tskit_node_id);
			rParent->ancestryRoot = splitRightWithParent(aNode->ancestryRoot, xOver, rParent->tskit_node_id);
		} else {
			// Minimal mode: preserve original tskit_node_id
			lParent->ancestryRoot = splitLeft(aNode->ancestryRoot, xOver);
			rParent->ancestryRoot = splitRight(aNode->ancestryRoot, xOver);
		}
		
		// Always free the original tree since split functions create new segments
		// With our segment-based approach, we don't need to preserve the child's tree
		freeSegmentTree(aNode->ancestryRoot);
		aNode->ancestryRoot = NULL;
		
		// Update stats from tree when in tree-only mode
		updateAncestryStatsFromTree(lParent);
		updateAncestryStatsFromTree(rParent);
		
		// Propagate sweep mutation flag based on which parent gets the sweep site
		if (aNode->carriesSweepMutation) {
			// Determine which parent inherits the sweep mutation based on sweep site position
			extern double sweepSite;
			if (sweepSite >= 0.0) {
				if (sweepSite < (float)xOver / nSites) {
					// Sweep site is on the left, so left parent inherits
					lParent->carriesSweepMutation = 1;
				} else {
					// Sweep site is on the right, so right parent inherits
					rParent->carriesSweepMutation = 1;
				}
			}
		}
		
		addNode(lParent);
		addNode(rParent);
		
		// Record edges in tskit for recombination (only in full mode)
		tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
		tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
		tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
		
		recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
		
	//	printf("reco-> lParent: %u rParent:%u child: %u time:%f xover:%d\n",lParent,rParent,aNode,cTime,xOver);
		
		// Mark that the child node is no longer in active set
		// It needs both parents to record before it can be freed
		aNode->inActiveSet = 0;
		
		// Mark parent recording if using tskit
		if (tskitOutputMode || minimalTreeSeq) {
			markParentRecorded(aNode, lParent);
			markParentRecorded(aNode, rParent);
		}
		
		return xOver;
	}
	return 666;
}

void geneConversionAtTimePopn(double cTime, int popn){
	rootedNode *aNode, *lParent, *rParent;
	int i;
	int xOver;
	int tractL;
	
	aNode = pickNodePopn(popn);
//	printf("GC\n",aNode);
	xOver = ignuin(0,nSites);
	tractL = (int) ceil(log(genunf(0,1))/log(1.0-(1.0/gcMean)));
//	printf("xo: %d test1: %d tractL: %d gcMean: %d\n",xOver,siteBetweenChunks(aNode, xOver),tractL,gcMean );
//	printNode(aNode);
//	for(i=0;i<nSites;i++)printf("%d",activeMaterial[i]);
//	printf("\n");
		
	if (siteBetweenChunks(aNode, xOver) == 1 &&  isActive(xOver) == 1  ){
		removeNode(aNode); 
		lParent = newRootedNode(cTime,popn);
		rParent = newRootedNode(cTime,popn);
		aNode->leftParent = lParent;
		aNode->rightParent = rParent;
		aNode->branchLength = cTime - aNode->time;
		lParent->leftChild = aNode;
		rParent->leftChild = aNode;
		lParent->population = aNode->population;
		rParent->population = aNode->population;
		lParent->nancSites=0;
		rParent->nancSites=0;
		lParent->lLim= nSites;
		lParent->rLim=0;
		rParent->lLim = nSites;
		rParent->rLim=0;
		
		
		// For gene conversion, handle the ancestry segments
		// Gene conversion creates a tract where one parent gets a segment and the other gets the rest
		if (aNode->ancestryRoot) {
			// Split the tree for gene conversion tract
			gcSplitResult gcSplit = splitSegmentTreeForGeneConversion(aNode->ancestryRoot, xOver, xOver + tractL);
			lParent->ancestryRoot = gcSplit.converted;     // Gets the converted tract
			rParent->ancestryRoot = gcSplit.unconverted;   // Gets everything else
			
			// Always free the original tree since split functions create new segments
			// With our segment-based approach, we don't need to preserve the child's tree
			freeSegmentTree(aNode->ancestryRoot);
			aNode->ancestryRoot = NULL;
		} else {
			lParent->ancestryRoot = NULL;
			rParent->ancestryRoot = NULL;
		}
		
		// Update stats from the ancestry trees
		updateAncestryStatsFromTree(lParent);
		updateAncestryStatsFromTree(rParent);
		
		// Propagate sweep mutation flag based on whether sweep site is in converted tract
		if (aNode->carriesSweepMutation) {
			extern double sweepSite;
			if (sweepSite >= 0.0) {
				int sweepSitePos = (int)(sweepSite * nSites);
				if (sweepSitePos >= xOver && sweepSitePos < xOver + tractL) {
					// Sweep site is in the converted tract
					lParent->carriesSweepMutation = 1;
				} else {
					// Sweep site is outside the converted tract
					rParent->carriesSweepMutation = 1;
				}
			}
		}
		
		addNode(lParent);
		addNode(rParent);
		
		// Record edges in tskit for gene conversion (only in full mode)
		tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
		tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
		tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
		
		recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
		
		// Mark that the child node is no longer in active set
		aNode->inActiveSet = 0;
		
		// Mark parent recording if using tskit
		if (tskitOutputMode || minimalTreeSeq) {
			markParentRecorded(aNode, lParent);
			markParentRecorded(aNode, rParent);
		}
	}
}


/*neutralPhase-- do coalescent and recombination events until
specified time. returns endTime. */
double neutralPhase(int *bpArray,double startTime, double endTime, double sizeRatio){
	double cTime, cRate[npops], rRate[npops], gcRate[npops], totRate, waitTime, bp,r;
	int  i;

	if(startTime == endTime){
		return(endTime);
	}
	cTime = 0.0;
	cTime += startTime;
	waitTime = 0.0;

	while (activeSites > 0){
		totRate = 0.0;
//		printf("currPopSize[0]: %d\n",popnSizes[0]);
		for(i=0;i<npops;i++){
			cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeRatio;
			rRate[i] = rho * popnSizes[i] * 0.5 ;
			gcRate[i] = my_gamma * popnSizes[i] * 0.5;
			totRate += cRate[i] + rRate[i] + gcRate[i];
		}
		//printf("%f\n",rRate);

		//find time of next event
		waitTime = genexp(1.0)  * (1.0/ totRate);
		cTime += waitTime;
		if (cTime >= endTime){
			return(endTime);
		}
		//find event type
		else{ 
			r =ranf();
			if (r < (rRate[0] / totRate)){
			//	printf("R %f active: %f\n",cTime, activeChromosomePercent());
				bp = recombineAtTimePopn(cTime,0);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((rRate[0]+gcRate[0])/totRate)){
					geneConversionAtTimePopn(cTime,0);
				}
				else{
					if(r < ((rRate[0]+gcRate[0]+cRate[0]) / totRate)){
						coalesceAtTimePopn(cTime,0);
					}
					else{
						if(r < ((rRate[0]+gcRate[0]+cRate[0]+rRate[1]) / totRate)){
							bp = recombineAtTimePopn(cTime,1);
							if (bp != 666){
								addBreakPoint(bp);
							}
						}
						else{
							if(r< ((rRate[0]+gcRate[0]+cRate[0]+rRate[1]+gcRate[1]) / totRate)){
								geneConversionAtTimePopn(cTime,1);
							}
							else{
								coalesceAtTimePopn(cTime,1);
							}
						}

					}
				}
			}
		}
	}
	return(cTime);
}


/*neutralPhase<ig-- does coalescent,recombination, and migration events until
specified time. returns endTime. */
double neutralPhaseMig(int *bpArray,double startTime, double endTime, double sizeRatio){
	double cTime, cRate[npops], rRate[npops], mRate[npops],totRate, waitTime, bp,r;
	int  i;

	if(startTime == endTime){
		return(endTime);
	}
	cTime = 0.0;
	cTime += startTime;
	waitTime = 0.0;

	while (activeSites > 0){
		totRate = 0.0;
		for(i=0;i<npops;i++){
			cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 1.0/ sizeRatio;
			rRate[i] = rho * popnSizes[i] * 1.0 ;/// sizeRatio;
			mRate[i] = mig[i] * popnSizes[i] * 1.0;/// sizeRatio;
			totRate += cRate[i] + rRate[i] + mRate[i];
		}
		//printf("%f\n",rRate);

		//find time of next event
		waitTime = genexp(1.0)  * (1.0/ totRate);
		cTime += waitTime;
		if (cTime >= endTime){
			return(endTime);
		}
		//find event type
		else{ 
			r =ranf();
			if (r < (rRate[0] / totRate)){
			//	printf("R %f active: %f\n",cTime, activeChromosomePercent());
				bp = recombineAtTimePopn(cTime,0);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((rRate[0]+cRate[0]) / totRate)){
					coalesceAtTimePopn(cTime,0);
				}
				else{
					if(r < ((rRate[0]+cRate[0]+rRate[1]) / totRate)){
						bp = recombineAtTimePopn(cTime,1);
						if (bp != 666){
							addBreakPoint(bp);
						}
					}
					else{
						if(r < ((rRate[0]+cRate[0]+rRate[1]+cRate[1]) / totRate)){
							coalesceAtTimePopn(cTime,1);
						}
						else{
							if(r < ((rRate[0]+cRate[0]+rRate[1]+cRate[1]+mRate[0]) / totRate)){
								migrateAtTime(cTime,0,1);
							}
							else
								migrateAtTime(cTime,1,0);	
						}
					}
				}

			}
		}
	}
	return(cTime);
}

/*neutralPhaseMigExclude-- does coalescent,recombination, and migration events until
specified time. returns endTime. excludes migrations from designated site to mimic local selection*/
double neutralPhaseMigExclude(int *bpArray,double startTime, double endTime, double sizeRatio, double selSite, double migScale){
	double cTime, cRate[npops], rRate[npops], mRate[npops],totRate, waitTime, bp,r;
	int  i;

	if(startTime == endTime){
		return(endTime);
	}
	cTime = 0.0;
	cTime += startTime;
	waitTime = 0.0;

	while (activeSites > 0){
		totRate = 0.0;
		for(i=0;i<npops;i++){
			cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5* sizeRatio;
			rRate[i] = rho * popnSizes[i] * 0.5 ;
			mRate[i] = mig[i] * popnSizes[i] * 0.5;
			totRate += cRate[i] + rRate[i] + mRate[i];
		}
		//printf("%f\n",rRate);

		//find time of next event
		waitTime = genexp(1.0)  * (1.0/ totRate);
		cTime += waitTime;
		if (cTime >= endTime && endTime > 0){
			return(endTime);
		}
		//find event type
		else{ 
			r =ranf();
			if (r < (rRate[0] / totRate)){
			//	printf("R %f active: %f\n",cTime, activeChromosomePercent());
				bp = recombineAtTimePopn(cTime,0);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((rRate[0]+cRate[0]) / totRate)){
					coalesceAtTimePopn(cTime,0);
				}
				else{
					if(r < ((rRate[0]+cRate[0]+rRate[1]) / totRate)){
						bp = recombineAtTimePopn(cTime,1);
						if (bp != 666){
							addBreakPoint(bp);
						}
					}
					else{
						if(r < ((rRate[0]+cRate[0]+rRate[1]+cRate[1]) / totRate)){
							coalesceAtTimePopn(cTime,1);
						}
						else{
							if(r < ((rRate[0]+cRate[0]+rRate[1]+cRate[1]+mRate[0]) / totRate)){
								migrateExceptSite(selSite,migScale,0,1);
							}
							else
								migrateExceptSite(selSite,migScale,1,0);	
						}
					}
				}

			}
		}
	}
	return(cTime);
}


/*neutralPhaseGeneralPopNumber--coalescent, recombination, gc events until
specified time. returns endTime. can handle multiple popns*/
double neutralPhaseGeneralPopNumber(int *bpArray,double startTime, double endTime, double *sizeRatio){
	double cTime, cRate[npops], rRate[npops], gcRate[npops], mRate[npops],totRate, waitTime, bp,r, r2;
	double totCRate, totRRate, totGCRate,totMRate, eSum;
	int  i,j;
	int coalescenceCounter = 0;  // Counter for coalescent events only
	const int SWEEP_INTERVAL = 10;  // Sweep every 10 coalescent events

	if(startTime == endTime){
		return(endTime);
	}
	cTime = 0.0;
	cTime += startTime;
	waitTime = 0.0;

	while (activeSites > 0){
		totRate = 0.0;
		totCRate = 0.0;
		totRRate = 0.0;
		totGCRate = 0.0;
		totMRate = 0.0;
		for(i=0;i<npops;i++){
			mRate[i]=0.0;
		}
		//printf("currPopSize[0]: %d currPopSize[1]: %d alleleNumber: %d totNodeNumber: %d activeSites: %d cTime: %f\n",popnSizes[0],popnSizes[1],alleleNumber,totNodeNumber, activeSites, cTime);
		for(i=0;i<npops;i++){
			cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeRatio[i];
			rRate[i] = rho * popnSizes[i] * 0.5;// * ((float)activeSites/nSites);
			gcRate[i] = my_gamma * popnSizes[i] * 0.5 ;
			
			for(j=0;j<npops;j++) mRate[i]+=migMat[i][j];
			mRate[i] *= popnSizes[i] * 0.5;
			totCRate += cRate[i];
			totRRate += rRate[i];
			totMRate += mRate[i];
			totGCRate += gcRate[i];
			totRate += cRate[i] + rRate[i] + mRate[i] + gcRate[i];
		}
		//printf("totRate: %f totCRate: %f totRRate: %f \n",totRate,totCRate,totRRate);

		//find time of next event
		waitTime = genexp(1.0)  * (1.0/ totRate);
		cTime += waitTime;
		
		if (cTime >= endTime){
			return(endTime);
		}
		//find event type
		else{ 
			r =ranf();
			if (r < (totRRate/ totRate)){
				//pick popn
				eSum = rRate[0];
				i = 0;
				r2 = ranf();
				while(eSum/totRRate < r2) eSum += rRate[++i];
				bp = recombineAtTimePopn(cTime,i);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((totRRate + totGCRate)/totRate)){
					//pick popn
					eSum = gcRate[0];
					i = 0;
					r2 = ranf();
					while(eSum/totGCRate < r2) eSum += gcRate[++i];
					geneConversionAtTimePopn(cTime,i);
				}
				else{
					if(r < ((totMRate+totRRate + totGCRate)/totRate)){
						//pick source popn
						eSum = mRate[0];
						i = 0;
						j = 0;
						r2 = ranf();
						while(eSum/totMRate < r2) eSum += mRate[++i];
						//printf("outer totMRate: %f eSum: %f i:%d\n",totMRate,eSum,i);
						//pick dest popn
						eSum = migMat[i][0]* popnSizes[i] * 0.5;
						//printf("outer eSumNew:%f mRate[%d]:%f\n",eSum,i,mRate[i]);
						r2 = ranf();
						while(eSum/mRate[i] < r2){
						//	printf("eSum: %f mRate[%d]:%f migMat[%d][0]:%f migMat[%d][1]:%f popnSize: %d\n",eSum,i,mRate[i],i,migMat[i][0],i,migMat[i][1], popnSizes[i]);
							eSum += migMat[i][++j] * popnSizes[i] * 0.5;
						//	printf("eSum: %f mRate[%d]:%f migMat[%d][0]:%f migMat[%d][1]:%f popnSize: %d\n",eSum,i,mRate[i],i,migMat[i][0],i,migMat[i][1], popnSizes[i]);
							
						} 
						migrateAtTime(cTime,i,j);
					}
				
					else{
						//coalesce 
						//pick popn
						eSum = cRate[0];
						i = 0;
						r2 = ranf();
						while(eSum/totCRate < r2){
							 eSum += cRate[++i];
							}
						coalesceAtTimePopn(cTime,i);
						
						// Increment coalescence counter and sweep
						coalescenceCounter++;
						if (coalescenceCounter % SWEEP_INTERVAL == 0) {
							// Flush buffered edges before freeing nodes to avoid referencing freed nodes
							extern int tskit_flush_edges_periodic(void);
							tskit_flush_edges_periodic();
							sweepAndFreeRemovedNodes();
						}
					
					
					}
				}
			}
		}
	}
	
	// Final sweep at the end of the phase
	// Flush buffered edges before freeing nodes
	extern int tskit_flush_edges_periodic(void);
	tskit_flush_edges_periodic();
	sweepAndFreeRemovedNodes();
	
	return(cTime);
}

void ensureTrajectoryCapacity(long int requiredSize) {
	// This function is now deprecated - we use file-based trajectories
	// Keeping it for compatibility but it just checks size limits
	if (requiredSize >= 500000000) {  // Match legacy limit
		fprintf(stderr, "trajectory too bigly. step= %ld. killing myself gently\n", requiredSize);
		exit(1);
	}
}


// Function to mmap an accepted trajectory file
void mmapAcceptedTrajectory(const char *filename, long int numSteps) {
	// Clean up any previous trajectory
	if (trajectoryFd != -1) {
		if (currentTrajectory && currentTrajectory != MAP_FAILED) {
			munmap(currentTrajectory, trajectoryFileSize);
		}
		close(trajectoryFd);
		trajectoryFd = -1;
	}
	
	// Open the trajectory file
	trajectoryFd = open(filename, O_RDONLY);
	if (trajectoryFd == -1) {
		perror("Failed to open trajectory file for mmap");
		exit(1);
	}
	
	// Calculate file size
	trajectoryFileSize = numSteps * sizeof(float);
	
	// Memory map the file
	currentTrajectory = (float *)mmap(NULL, trajectoryFileSize, PROT_READ, 
	                                  MAP_PRIVATE, trajectoryFd, 0);
	if (currentTrajectory == MAP_FAILED) {
		perror("Failed to mmap trajectory file");
		close(trajectoryFd);
		exit(1);
	}
}

// Function to clean up rejected trajectory files
void cleanupRejectedTrajectory(const char *filename) {
	if (filename && filename[0] != '\0') {
		unlink(filename);
	}
}

/*proposeTrajectory-- this function creates a sweep trajectory and deals with
complications like changing population size, or soft sweeps, etc 
returns the acceptance probability of the trajectory */
double proposeTrajectory(int currentEventNumber, float *currentTrajectory, double *sizeRatio, char sweepMode, \
double initialFreq, double *finalFreq, double alpha, double f0, double currentTime)
{	
	double tInc, tIncOrig, minF,ttau, N;
	double N_0 = (double) EFFECTIVE_POPN_SIZE;
	double Nmax, localNextTime,localCurrentTime, currentSizeRatio;
	int i, insweepphase;
	long int j;
	float x;
	
	// For sweep simulations, write directly to a temporary file
	char tempFilename[256];
	snprintf(tempFilename, sizeof(tempFilename), "/tmp/discoal_traj_%d_%ld_%d.tmp", 
	         getpid(), time(NULL), rand());
	
	FILE *trajFile = fopen(tempFilename, "wb");
	if (!trajFile) {
		perror("Failed to create trajectory file");
		exit(1);
	}
	
	// Use buffered writes for efficiency
	float writeBuffer[1024];
	int bufferPos = 0;
	
	tIncOrig = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	j=0;
	N = (double) floor(N_0 * sizeRatio[0]);
	Nmax = currentSizeRatio = sizeRatio[0];
	x = initialFreq;
	minF = f0;
	insweepphase = 1;
	ttau=0.0;
	//go through each event until last one
	for(i=currentEventNumber;i<eventNumber;i++){
		localCurrentTime = events[i].time;
		if(i == eventNumber - 1){
			localNextTime = MAXTIME;
			}
		else{
			localNextTime = events[i+1].time;
		}
		if(events[i].type == 'n'){
			currentSizeRatio = events[i].popnSize;
			N = floor(N_0 *events[i].popnSize);
			if(currentSizeRatio > Nmax) Nmax = currentSizeRatio;
		}
		if(minF < 1.0/(2.*N)) minF = 1.0/(2.*N);
		tInc = 1.0 / (deltaTMod * N);
		//iterate until epoch time or sweep freq
		while( x > 1.0/(2.*N) && (currentTime+ttau) < localNextTime){
			ttau += tIncOrig;
			if(x > minF && insweepphase){
				//get next sweep allele freq
				switch(sweepMode){
					case 'd':
					x = detSweepFreq(ttau, alpha * currentSizeRatio);
					break;
					case 's':
					x = 1.0 - genicSelectionStochasticForwardsOptimized(tInc, (1.0 - x), alpha * currentSizeRatio);
					break;
					case 'N':
					x = neutralStochasticOptimized(tInc, x);
					break;
				}
			}
			else{
				insweepphase = 0;
				tInc = 1.0 / (deltaTMod * N );
				x = neutralStochasticOptimized(tInc, x);
			}
			//printf("j: %ld x: %f\n",j,x);
			
			// Check trajectory size to prevent runaway
			if (j >= 500000000) {  // Match legacy limit
				fprintf(stderr, "trajectory too bigly. step= %ld. killing myself gently\n", j);
				fclose(trajFile);
				unlink(tempFilename);
				exit(1);
			}
			
			// Write to buffer
			writeBuffer[bufferPos++] = x;
			if (bufferPos >= 1024) {
				// Flush buffer to file
				fwrite(writeBuffer, sizeof(float), bufferPos, trajFile);
				bufferPos = 0;
			}
			j++;
		}
	}
	
	// Flush any remaining data in buffer
	if (bufferPos > 0) {
		fwrite(writeBuffer, sizeof(float), bufferPos, trajFile);
	}
	fclose(trajFile);
	
	// Store the filename and trajectory length globally
	strncpy(trajectoryFilename, tempFilename, sizeof(trajectoryFilename) - 1);
	currentTrajectoryStep = 0;
	totalTrajectorySteps = j;
	
	// Note: We don't mmap here because this function may be called multiple times
	// during rejection sampling. The accepted trajectory will be mmap'd later.
	
	return(currentSizeRatio/Nmax);
	
}


/*sweepPhaseEventsGeneralPopNumber-- this does the compressed time sweep thing. 
  generalized to account for popnSize changes returns the time after the sweep.
  Same as sweepPhaseEvents but includes recurrent adaptive mutation */
double sweepPhaseEventsGeneralPopNumber(int *bpArray, double startTime, double endTime, double sweepSite,\
double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
double *sizeRatio, char sweepMode,double f0, double uA)
{

	double totRate,bp;
	double  ttau, x, tInc, tIncOrig;
	double pCoalB, pCoalb, pRecB, pRecb, r, sum,eventRand,eventProb, pLeftRecB, pLeftRecb;
	double pRecurMut, pGCB, pGCb;
	double sweepPopTotRate,cRate[npops], rRate[npops], gcRate[npops];
	double totCRate, totRRate, totGCRate, eSum, r2;
	double N = (double) EFFECTIVE_POPN_SIZE;
	double minF;
	double cTime = startTime;
	int insweepphase, i;

	//initialize stuff
	pCoalB = pCoalb = pRecB = pRecb = totRate = pRecurMut = pLeftRecB = pLeftRecb = totGCRate = totCRate = totRRate = 0;
	sweepPopTotRate = pGCB = pGCb = 0;

	N = (double) floor(N * sizeRatio[0]);
	ttau = 0.0;
	x = initialFreq;
	minF = f0;
	if(minF < 1.0/(2.*N))
		minF = 1.0/(2.*N);

	//if new sweep then reset sweepPopnIDs; only popn 0 sweeps, sweepPopn==1 is the beneficial class
	if(!*stillSweeping){
		for(i=0;i<alleleNumber;i++){
			if(nodes[i]->population==0){
				if(partialSweepMode == 1){
					//for partial sweeps choose randomly acccording to final sweep freq
					if(ranf()>partialSweepFinalFreq){
						nodes[i]->sweepPopn = 0;						
					}
					else{
						nodes[i]->sweepPopn = 1;
						if(isAncestralHere(nodes[i],sweepSite) && hidePartialSNP == 0) {
							// Mark that this node carries the sweep mutation
							nodes[i]->carriesSweepMutation = 1;
						}
					}
				}
				else{
			 		nodes[i]->sweepPopn=1;
				}
			}
		}
	*stillSweeping = 1;
	}
	
	//assume that sweep always happens in popn 0!!!
	//using popnSize global to manage bookkeeping
	sweepPopnSizes[1] = nodePopnSweepSize(0,1);
	sweepPopnSizes[0] = nodePopnSweepSize(0,0);
	//  printf("sweepPopnSizes: %d %d \n",sweepPopnSizes[1],sweepPopnSizes[0]);

	//set time increment
	tInc = 1.0 / (deltaTMod * N);
	tIncOrig = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	insweepphase = 1;
	
	//go for epoch time, sweep freq, or root
	while( x > 1.0/(2.*N) && (cTime+ttau) < endTime && popnSizes[0] > 1){ 
		//rejection algorithm of Braverman et al. 1995
		eventRand = ranf();
		eventProb = 1.0;
		//wait for something
		while(eventProb > eventRand && x > (1.0 / (2*N)) && (cTime+ttau) < endTime){
			ttau += tIncOrig;

			if(x > minF && insweepphase)
			{
				//get next sweep allele freq
				switch(sweepMode){
					case 'd':
					x = detSweepFreq(ttau, alpha * sizeRatio[0]);
				//	printf("x here:%f ttau: %f alpha*sizeRatio: %f\n",x,ttau,alpha*sizeRatio);

					break;
					case 's':
					x = 1.0 - genicSelectionStochasticForwardsOptimized(tInc, (1.0 - x), alpha * sizeRatio[0]);
				//	printf("x here:%f ttau: %f alpha*sizeRatio: %f\n",x,ttau,alpha*sizeRatio[0]);
					break;
					case 'N':
				//	printf("here\n");
					x = neutralStochasticOptimized(tInc, x);
					break;
				}
				//	printf("x:%g ttau: %f\n",x,ttau );
			}
			else{
				insweepphase = 0;
				tInc = 1.0 / (deltaTMod * N );
				x = neutralStochastic(tInc, x);
			}

			//calculate event probs
			//first 4 events are probs of events in population 0
			pCoalB = ((sweepPopnSizes[1] * (sweepPopnSizes[1] - 1) ) * 0.5)/x*tIncOrig / sizeRatio[0];
			pCoalb = ((sweepPopnSizes[0] * (sweepPopnSizes[0] - 1) ) * 0.5)/(1-x)*tIncOrig / sizeRatio[0];
			pRecB = rho * sweepPopnSizes[1]*0.5 *tIncOrig; // / sizeRatio[0];
			pRecb = rho * sweepPopnSizes[0]*0.5 *tIncOrig;// / sizeRatio[0];
			pGCB = my_gamma * sweepPopnSizes[1]*0.5 *tIncOrig;// / sizeRatio[0];
			pGCb = my_gamma * sweepPopnSizes[0]*0.5 *tIncOrig;/// sizeRatio[0];
			pRecurMut = (uA * sweepPopnSizes[1]*0.5 *tIncOrig)/x;///sizeRatio[0];
			if (sweepSite < 0.0){
				pLeftRecB = leftRho * sweepPopnSizes[1]*0.5 * tIncOrig * (1-x);
				pLeftRecb = leftRho * sweepPopnSizes[0]*0.5 * tIncOrig * x;
                        }
			sweepPopTotRate = pCoalB + pCoalb + pRecB + pRecb + pGCB + pGCb + pRecurMut + pLeftRecB + pLeftRecb;
			//printf("nB: %d; x: %g; pCoalB: %g; tInc: %g cTime+ttau: %g popnSizes[0]:%d\n",sweepPopnSizes[1],x,pCoalB,tInc,cTime+ttau,  popnSizes[0]);
			//now two events in population 1
			totRate = 0.0;
			totCRate = 0.0;
			totRRate = 0.0;
			totGCRate = 0.0;
			totRate = sweepPopTotRate;
			
			//printf("currPopSize[0]: %d currPopSize[1]: %d\n",popnSizes[0],popnSizes[1]);
			for(i=1;i<npops;i++){
				cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 * tIncOrig / sizeRatio[i];
				rRate[i] = rho * popnSizes[i] * 0.5 * tIncOrig;// / sizeRatio[i];
				gcRate[i] = my_gamma * popnSizes[i] * 0.5 * tIncOrig;// / sizeRatio[i];
				totCRate += cRate[i];
				totRRate += rRate[i];
				totGCRate += gcRate[i];
				totRate += cRate[i] + rRate[i] + gcRate[i];
			}
			

			eventProb *= 1-totRate;
			//printf("x(t): %g t: %g tinc: %g eventProb: %g eventRand: %g totRate: %g swPopnSize1: %d swPopnSize2: %d\n",x,ttau,tInc, eventProb,
			//	eventRand,totRate,sweepPopnSizes[0],sweepPopnSizes[1]);	


		}

		//Decide which event took place
		//first was it in population 0 (sweep) or not
		double eventTime = cTime + ttau;

		if(ranf() < sweepPopTotRate / totRate){
			//event is in sweep pop; choose event
			r = ranf();
			sum = pCoalB;
			//coalescent in B?
			if(r < sum / sweepPopTotRate){
				coalesceAtTimePopnSweep(eventTime, 0,1);
			}
			else{
				sum += pCoalb;
				//coalescent in b?
				if(r < sum / sweepPopTotRate){
					coalesceAtTimePopnSweep(eventTime,0, 0);
				}
				else{
					sum += pRecb;
					//recombination in b, also need bookkeeping
					if(r < sum / sweepPopTotRate){
						bp = recombineAtTimePopnSweep(eventTime,0, 0, sweepSite, (1.0-x));
						if(bp != 666){
							addBreakPoint(bp);
						}
					}
					else{
						sum += pRecB;
						if( r < sum / sweepPopTotRate){
							//recombination in B and bookkeeping
							bp = recombineAtTimePopnSweep(eventTime,0, 1, sweepSite, x);
							if(bp != 666){
								addBreakPoint(bp);
							}
						}
						else{
							sum+= pGCB;
							if( r < sum / sweepPopTotRate){
								geneConversionAtTimePopnSweep(eventTime,0, 1, sweepSite, x);
							}
							else{
								sum+= pGCb;
								if( r < sum / sweepPopTotRate){
									geneConversionAtTimePopnSweep(eventTime,0, 0, sweepSite, x);
								}
								else{
									sum += pLeftRecb;
									if ( r < sum / totRate){
									        recombineToLeftPopnSweep(0, 0, x);
									}
									else{
										sum += pLeftRecB;
										if ( r < sum / totRate){
										recombineToLeftPopnSweep(0, 1, x);
										}
										else{
											//recurrent adaptive mutation:
											//node in pop zero's sweep group exits sweep
											recurrentMutAtTime(eventTime,0, 1);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else{	//event is in another population
			totRate -= sweepPopTotRate; //discount sweepPopnEvent probs from tot
			r = ranf();
			if (r < (totRRate/ totRate)){
				//pick popn
				eSum = rRate[1];
				i = 1;
				r2 = ranf();
				while(eSum/totRRate < r2) eSum += rRate[++i];
				bp = recombineAtTimePopn(eventTime,i);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((totRRate + totGCRate)/totRate)){
					//pick popn
					eSum = gcRate[1];
					i = 1;
					r2 = ranf();
					while(eSum/totGCRate < r2) eSum += gcRate[++i];
					geneConversionAtTimePopn(eventTime,i);
				}
				else{
					//coalesce 
					//pick popn
					eSum = cRate[1];
				//	printf("esum= %f\n",eSum);
					i = 1;
					r2 = ranf();
					while(eSum/totCRate < r2){
						 eSum += cRate[++i];
						}
					coalesceAtTimePopn(eventTime,i);
					
					
				}
			}
		}

	}
	if((cTime+ttau) >= endTime){
		*stillSweeping = 1; 
	}
	else{
		*stillSweeping = 0; 
	}
	if(sweepPopnSizes[1]==0) *stillSweeping = 0; 

	*finalFreq = x;
//	printf("sp: %d %d\n",	sweepPopnSizes[0], 	sweepPopnSizes[1]);
	return(cTime+(ttau));
}

/*sweepPhaseEventsConditionalTrajectory-- does sweep phase with trajectory created externally */
double sweepPhaseEventsConditionalTrajectory(int *bpArray, double startTime, double endTime, double sweepSite,\
double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
double *sizeRatio, char sweepMode,double f0, double uA)
{

	double totRate,bp;
	double  ttau, x, tInc, tIncOrig;
	double pCoalB, pCoalb, pRecB, pRecb, r, sum,eventRand,eventProb, pLeftRecB, pLeftRecb;
	double pRecurMut, pGCB, pGCb;
	double sweepPopTotRate,cRate[npops], rRate[npops], gcRate[npops];
	double totCRate, totRRate, totGCRate, eSum, r2;
	double N = (double) EFFECTIVE_POPN_SIZE;
	double minF;
	double cTime = startTime;
	int insweepphase, i;

	//initialize stuff
	pCoalB = pCoalb = pRecB = pRecb = totRate = pRecurMut = pLeftRecB = pLeftRecb = totGCRate = totCRate = totRRate = 0;
	sweepPopTotRate = pGCB = pGCb = 0;

	N = (double) floor(N * sizeRatio[0]);
	ttau = 0.0;
	x = initialFreq;
	minF = f0;
	if(minF < 1.0/(2.*N))
		minF = 1.0/(2.*N);

	//if new sweep then reset sweepPopnIDs; only popn 0 sweeps, sweepPopn==1 is the beneficial class
	if(!*stillSweeping){
		for(i=0;i<alleleNumber;i++){
			if(nodes[i]->population==0){
				if(partialSweepMode == 1){
					//for partial sweeps choose randomly acccording to final sweep freq
					if(ranf()>partialSweepFinalFreq){
						nodes[i]->sweepPopn = 0;						
					}
					else{
						nodes[i]->sweepPopn = 1;
						if(isAncestralHere(nodes[i],sweepSite) && hidePartialSNP == 0) {
							// Mark that this node carries the sweep mutation
							nodes[i]->carriesSweepMutation = 1;
						}
					}
				}
				else{
			 		nodes[i]->sweepPopn=1;
				}
			}
		}
	*stillSweeping = 1;
	}
	
	//assume that sweep always happens in popn 0!!!
	//using popnSize global to manage bookkeeping
	sweepPopnSizes[1] = nodePopnSweepSize(0,1);
	sweepPopnSizes[0] = nodePopnSweepSize(0,0);
	//  printf("sweepPopnSizes: %d %d \n",sweepPopnSizes[1],sweepPopnSizes[0]);

	//set time increment
	tInc = 1.0 / (deltaTMod * N);
	tIncOrig = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	insweepphase = 1;
	
	//go for epoch time, sweep freq, or root
	while( x > 1.0/(2.*N) && (cTime+ttau) < endTime && popnSizes[0] > 1){ 
		//rejection algorithm of Braverman et al. 1995
		eventRand = ranf();
		eventProb = 1.0;
		//wait for something
		while(eventProb > eventRand && x > (1.0 / (2*N)) && (cTime+ttau) < endTime ){
			ttau += tIncOrig;

			if (currentTrajectoryStep >= totalTrajectorySteps) {
			fprintf(stderr, "Error: trajectory step %ld exceeds total steps %ld\n", 
					currentTrajectoryStep, totalTrajectorySteps);
			exit(1);
		}
		x = currentTrajectory[currentTrajectoryStep++];

			//calculate event probs
			//first 4 events are probs of events in population 0
			pCoalB = ((sweepPopnSizes[1] * (sweepPopnSizes[1] - 1) ) * 0.5)/x*tIncOrig / sizeRatio[0];
			pCoalb = ((sweepPopnSizes[0] * (sweepPopnSizes[0] - 1) ) * 0.5)/(1-x)*tIncOrig / sizeRatio[0];
			pRecB = rho * sweepPopnSizes[1]*0.5 *tIncOrig; // / sizeRatio[0];
			pRecb = rho * sweepPopnSizes[0]*0.5 *tIncOrig;// / sizeRatio[0];
			pGCB = my_gamma * sweepPopnSizes[1]*0.5 *tIncOrig;// / sizeRatio[0];
			pGCb = my_gamma * sweepPopnSizes[0]*0.5 *tIncOrig;/// sizeRatio[0];
			pRecurMut = (uA * sweepPopnSizes[1]*0.5 *tIncOrig)/x;///sizeRatio[0];
			if (sweepSite < 0.0){
				pLeftRecB = leftRho * sweepPopnSizes[1]*0.5 * tIncOrig * (1-x);
				pLeftRecb = leftRho * sweepPopnSizes[0]*0.5 * tIncOrig * x;
                        }
			sweepPopTotRate = pCoalB + pCoalb + pRecB + pRecb + pGCB + pGCb + pRecurMut + pLeftRecB + pLeftRecb;
			//printf("nB: %d; x: %g; pCoalB: %g; tInc: %g cTime+ttau: %g popnSizes[0]:%d\n",sweepPopnSizes[1],x,pCoalB,tInc,cTime+ttau,  popnSizes[0]);
			//now two events in population 1
			totRate = 0.0;
			totCRate = 0.0;
			totRRate = 0.0;
			totGCRate = 0.0;
			totRate = sweepPopTotRate;
			
			//printf("currPopSize[0]: %d currPopSize[1]: %d\n",popnSizes[0],popnSizes[1]);
			for(i=1;i<npops;i++){
				cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 * tIncOrig / sizeRatio[i];
				rRate[i] = rho * popnSizes[i] * 0.5 * tIncOrig;// / sizeRatio[i];
				gcRate[i] = my_gamma * popnSizes[i] * 0.5 * tIncOrig;// / sizeRatio[i];
				totCRate += cRate[i];
				totRRate += rRate[i];
				totGCRate += gcRate[i];
				totRate += cRate[i] + rRate[i] + gcRate[i];
			}
			

			eventProb *= 1-totRate;
			//printf("x(t): %g t: %g tinc: %g eventProb: %g eventRand: %g totRate: %g swPopnSize1: %d swPopnSize2: %d\n",x,ttau,tInc, eventProb,
			//	eventRand,totRate,sweepPopnSizes[0],sweepPopnSizes[1]);	


		}
		if(cTime+ttau >= endTime) return(cTime+ttau);

		//Decide which event took place
		//first was it in population 0 (sweep) or not
		double eventTime = cTime + ttau;
		if(ranf() < sweepPopTotRate / totRate){
			//event is in sweep pop; choose event
			r = ranf();
			sum = pCoalB;
			//coalescent in B?
			if(r < sum / sweepPopTotRate){
				coalesceAtTimePopnSweep(eventTime, 0,1);
			}
			else{
				sum += pCoalb;
				//coalescent in b?
				if(r < sum / sweepPopTotRate){
					coalesceAtTimePopnSweep(eventTime,0, 0);
				}
				else{
					sum += pRecb;
					//recombination in b, also need bookkeeping
					if(r < sum / sweepPopTotRate){
						bp = recombineAtTimePopnSweep(eventTime,0, 0, sweepSite, (1.0-x));
						if(bp != 666){
							addBreakPoint(bp);
						}
					}
					else{
						sum += pRecB;
						if( r < sum / sweepPopTotRate){
							//recombination in B and bookkeeping
							bp = recombineAtTimePopnSweep(eventTime,0, 1, sweepSite, x);
							if(bp != 666){
								addBreakPoint(bp);
							}
						}
						else{
							sum+= pGCB;
							if( r < sum / sweepPopTotRate){
								geneConversionAtTimePopnSweep(eventTime,0, 1, sweepSite, x);
							}
							else{
								sum+= pGCb;
								if( r < sum / sweepPopTotRate){
									geneConversionAtTimePopnSweep(eventTime,0, 0, sweepSite, x);
								}
								else{
									sum += pLeftRecb;
									if ( r < sum / totRate){
									        recombineToLeftPopnSweep(0, 0, x);
									}
									else{
										sum += pLeftRecB;
										if ( r < sum / totRate){
										recombineToLeftPopnSweep(0, 1, x);
										}
										else{
											//recurrent adaptive mutation:
											//node in pop zero's sweep group exits sweep
											recurrentMutAtTime(eventTime,0, 1);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else{	
			if(npops>1){
				//event is in another population
				totRate -= sweepPopTotRate; //discount sweepPopnEvent probs from tot
				if(totRate > 0){ //need to catch cases without recombination here
					r = ranf();
					if (r < (totRRate/ totRate)){
						//pick popn
						eSum = rRate[1];
						i = 1;
						r2 = ranf();
						while(eSum/totRRate < r2) eSum += rRate[++i];
						bp = recombineAtTimePopn(eventTime,i);
						if (bp != 666){
							addBreakPoint(bp);
						}
					}
					else{
						if(r < ((totRRate + totGCRate)/totRate)){
							//pick popn
							eSum = gcRate[1];
							i = 1;
							r2 = ranf();
							while(eSum/totGCRate < r2) eSum += gcRate[++i];
							geneConversionAtTimePopn(eventTime,i);
						}
						else{
							//coalesce 
							//pick popn
							eSum = cRate[1];
						//	printf("esum= %f\n",eSum);
							i = 1;
							r2 = ranf();
							while(eSum/totCRate < r2){
								 eSum += cRate[++i];
								}
							coalesceAtTimePopn(eventTime,i);
					
					
						}
					}
				}
			}
		}

	}
	if((cTime+ttau) >= endTime){
		*stillSweeping = 1; 
	}
	else{
		*stillSweeping = 0; 
	}
	if(sweepPopnSizes[1]==0) *stillSweeping = 0; 

	*finalFreq = x;
//	printf("sp: %d %d\n",	sweepPopnSizes[0], 	sweepPopnSizes[1]);
	return(cTime+(ttau));
}


/*recurrentSweepPhaseGeneralPopNumber--coalescent, recombination, gc events, and sweeps! until
specified time. returns endTime. can handle multiple popns*/
double recurrentSweepPhaseGeneralPopNumber(int *bpArray,double startTime, double endTime, double *finalFreq, double alpha, char sweepMode, double *sizeRatio){
	double cTime, cRate[npops], rRate[npops], gcRate[npops], mRate[npops],totRate, waitTime, bp,r, r2;
	double totCRate, totRRate, totGCRate,totMRate, eSum,curSweepSite, initFreq, probAccept;
	int  i,j;

	if(startTime == endTime){
		return(endTime);
	}
	cTime = 0.0;
	cTime += startTime;
	waitTime = 0.0;

	while (activeSites > 0){
		totRate = 0.0;
		totCRate = 0.0;
		totRRate = 0.0;
		totGCRate = 0.0;
		totMRate = 0.0;
		for(i=0;i<npops;i++){
			mRate[i]=0.0;
		}
	//	printf("currPopSize[0]: %d currPopSize[1]: %d alleleNumber: %d totNodeNumber: %d activeSites: %d cTime: %f\n",popnSizes[0],popnSizes[1],alleleNumber,totNodeNumber, activeSites, cTime);
		for(i=0;i<npops;i++){
			cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeRatio[i];
			rRate[i] = rho * popnSizes[i] * 0.5;// * ((float)activeSites/nSites);
			gcRate[i] = my_gamma * popnSizes[i] * 0.5 ;
			
			for(j=0;j<npops;j++) mRate[i]+=migMat[i][j];
			mRate[i] *= popnSizes[i] * 0.5;
			totCRate += cRate[i];
			totRRate += rRate[i];
			totMRate += mRate[i];
			totGCRate += gcRate[i];
			totRate += cRate[i] + rRate[i] + mRate[i] + gcRate[i];
		}
		//add the recurrent sweep probability
		totRate += recurSweepRate;

		//find time of next event
		waitTime = genexp(1.0)  * (1.0/ totRate);
		cTime += waitTime;
		if (cTime >= endTime){
			return(endTime);
		}
		//find event type
		else{ 
			r =ranf();
			if (r < (totRRate/ totRate)){
				//pick popn
				eSum = rRate[0];
				i = 0;
				r2 = ranf();
				while(eSum/totRRate < r2) eSum += rRate[++i];
				bp = recombineAtTimePopn(cTime,i);
				if (bp != 666){
					addBreakPoint(bp);
				}
			}
			else{
				if(r < ((totRRate + totGCRate)/totRate)){
					//pick popn
					eSum = gcRate[0];
					i = 0;
					r2 = ranf();
					while(eSum/totGCRate < r2) eSum += gcRate[++i];
					geneConversionAtTimePopn(cTime,i);
				}
				else{
					if(r < ((totMRate+totRRate + totGCRate)/totRate)){
						//pick source popn
						eSum = mRate[0];
						i = 0;
						j = 0;
						r2 = ranf();
						while(eSum/totMRate < r2) eSum += mRate[++i];
						//printf("outer totMRate: %f eSum: %f i:%d\n",totMRate,eSum,i);
						//pick dest popn
						eSum = migMat[i][0]* popnSizes[i] * 0.5;
						//printf("outer eSumNew:%f mRate[%d]:%f\n",eSum,i,mRate[i]);
						r2 = ranf();
						while(eSum/mRate[i] < r2){
						//	printf("eSum: %f mRate[%d]:%f migMat[%d][0]:%f migMat[%d][1]:%f popnSize: %d\n",eSum,i,mRate[i],i,migMat[i][0],i,migMat[i][1], popnSizes[i]);
							eSum += migMat[i][++j] * popnSizes[i] * 0.5;
						//	printf("eSum: %f mRate[%d]:%f migMat[%d][0]:%f migMat[%d][1]:%f popnSize: %d\n",eSum,i,mRate[i],i,migMat[i][0],i,migMat[i][1], popnSizes[i]);
							
						} 
						migrateAtTime(cTime,i,j);
					}
				
					else{
						if(r < ((totCRate + totMRate + totRRate + totGCRate)/totRate)){
							//coalesce 
							//pick popn
							eSum = cRate[0];
							i = 0;
							r2 = ranf();
							while(eSum/totCRate < r2){
							 	eSum += cRate[++i];
							}
							coalesceAtTimePopn(cTime,i);
						}
						else{
							if(sweepSite<0){
								curSweepSite = -1.0;
								leftRho = genunf(0.0,2.0 * alpha);
							}
							else{
								curSweepSite = ranf();
							}
							if(partialSweepMode==1){
								initFreq=MIN(partialSweepFinalFreq,1.0-(1.0/(2*sizeRatio[0]*EFFECTIVE_POPN_SIZE)));
							}
							else{
								initFreq=1.0-(1.0/(2*sizeRatio[0]*EFFECTIVE_POPN_SIZE));
							}
							//generate a proposed trajectory
							probAccept = proposeTrajectory(currentEventNumber, currentTrajectory, sizeRatio, sweepMode, initFreq, finalFreq, alpha, f0, cTime);
							while(ranf()>probAccept){
								probAccept = proposeTrajectory(currentEventNumber, currentTrajectory, sizeRatio, sweepMode, initFreq, finalFreq, alpha, f0, cTime);
								//printf("probAccept: %lf\n",probAccept);
							}
							cTime= sweepPhaseEventsConditionalTrajectory(bpArray, cTime, endTime, curSweepSite,\
								initFreq, finalFreq, &activeSweepFlag, alpha,\
								sizeRatio, sweepMode,0, 0);
							
						}
					
					}
				}
			}
		}
	}
	return(cTime);
}



/*recombineAtTimePopnSweep-- preforms recombination on
an individual drawn from a popn and assigns parental
	popn based on sweep site and frequency of popn / 
	beneficial mutation */
int recombineAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq){
	rootedNode *aNode, *lParent, *rParent;
	int i;
	int xOver;
	double r;


		aNode = pickNodePopnSweep(popn, sp);
	//	printf("picked: %p\n",aNode);
		xOver = ignuin(0,nSites-1);
	//	printf("xo: %d test1: %d test2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver) );
	//	printNode(aNode);
	//	for(i=0;i<nSites;i++)printf("%d",activeMaterial[i]);
	//	printf("\n");

		if (siteBetweenChunks(aNode, xOver) == 1 &&  isActive(xOver) == 1  ){
			removeNode(aNode); 
			lParent = newRootedNode(cTime,popn);
			rParent = newRootedNode(cTime,popn);
			aNode->leftParent = lParent;
			aNode->rightParent = rParent;
			aNode->branchLength = cTime - aNode->time;
			lParent->leftChild = aNode;
			rParent->leftChild = aNode;
			lParent->population = aNode->population;
			rParent->population = aNode->population;
			lParent->nancSites=0;
			rParent->nancSites=0;
			lParent->lLim= nSites;
			lParent->rLim=0;
			rParent->lLim = nSites;
			rParent->rLim=0;
			
			//determine sweep popn affinity
			//left side?
		//	printf("here-- popnFreq:%f sweepSite:%f xOver: %d test: %d \n",popnFreq,sweepSite,xOver,sweepSite < (float) xOver / nSites);
			if(sweepSite < (float) xOver / nSites){ //lParent has sweep site set to same as child node,rParent is then random
				lParent->sweepPopn = sp;
				r = ranf();
				if (r < popnFreq)rParent->sweepPopn = sp;
				else rParent->sweepPopn = (sp == 0) ? 1:0;
			}
			else{ //rParent has sweep site, lParent is then random
				rParent->sweepPopn = sp;
				r = ranf();
				if (r < popnFreq) lParent->sweepPopn = sp;
				else lParent->sweepPopn = (sp == 0) ? 1:0;
			}
		//	printf("lParentsp:%d rParentsp:%d\n",lParent->sweepPopn,rParent->sweepPopn);
			
			// Propagate sweep mutation flag to the parent that contains the sweep site
			if (aNode->carriesSweepMutation) {
				if (sweepSite < (float)xOver / nSites) {
					// Sweep site is on the left
					lParent->carriesSweepMutation = 1;
				} else {
					// Sweep site is on the right
					rParent->carriesSweepMutation = 1;
				}
			}
			
			// Split ancestry segment tree at crossover point
			extern int minimalTreeSeq;
			if (!minimalTreeSeq) {
				// Full ARG mode: use parent's tskit_node_id in segments
				lParent->ancestryRoot = splitLeftWithParent(aNode->ancestryRoot, xOver, lParent->tskit_node_id);
				rParent->ancestryRoot = splitRightWithParent(aNode->ancestryRoot, xOver, rParent->tskit_node_id);
			} else {
				// Minimal mode: preserve original tskit_node_id
				lParent->ancestryRoot = splitLeft(aNode->ancestryRoot, xOver);
				rParent->ancestryRoot = splitRight(aNode->ancestryRoot, xOver);
			}
			
			// Update stats from ancestry trees
			updateAncestryStatsFromTree(lParent);
			updateAncestryStatsFromTree(rParent);
			
			//add in the nodes
			addNode(lParent);
			addNode(rParent);
			
			// Record edges in tskit for recombination during sweep (only in full mode)
			tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
			tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
			tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
			
			recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
			
			// Mark that the child node is no longer in active set
			aNode->inActiveSet = 0;
			
			// Mark parent recording if using tskit
			if (tskitOutputMode || minimalTreeSeq) {
				markParentRecorded(aNode, lParent);
				markParentRecorded(aNode, rParent);
			}
			
			//sweepPopnSizes[sp]++;
			return xOver;
		}
		return 666; //debuging
}

int recombineToLeftPopnSweep(int popn, int sp, double popnFreq){
	rootedNode *aNode;

	aNode = pickNodePopnSweep(popn, sp);
	sweepPopnSizes[aNode->sweepPopn]--;

	aNode->sweepPopn = (sp == 0) ? 1:0;
	sweepPopnSizes[aNode->sweepPopn]++;
	return 0;
}

/*geneConversionAtTimePopnSweep-- preforms gene conversion on
an individual drawn from a popn and assigns parental
	popn based on sweep site and frequency of popn / 
	beneficial mutation */
void geneConversionAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq){
	rootedNode *aNode, *lParent, *rParent;
	int i;
	int xOver, tractL;
	double r;


	aNode = pickNodePopnSweep(popn, sp);
	//	printf("picked: %p\n",aNode);
	xOver = ignuin(0,nSites);
	tractL = (int) ceil(log(genunf(0,1))/log(1.0-(1.0/gcMean)));

	//	printf("xo: %d test1: %d test2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver) );
	//	printNode(aNode);
	//	for(i=0;i<nSites;i++)printf("%d",activeMaterial[i]);
	//	printf("\n");

	if (siteBetweenChunks(aNode, xOver) == 1 &&  isActive(xOver) == 1  ){
		removeNode(aNode); 
		lParent = newRootedNode(cTime,popn);
		rParent = newRootedNode(cTime,popn);
		aNode->leftParent = lParent;
		aNode->rightParent = rParent;
		aNode->branchLength = cTime - aNode->time;
		lParent->leftChild = aNode;
		rParent->leftChild = aNode;
		lParent->population = aNode->population;
		rParent->population = aNode->population;
		lParent->nancSites=0;
		rParent->nancSites=0;
		lParent->lLim= nSites;
		lParent->rLim=0;
		rParent->lLim = nSites;
		rParent->rLim=0;
		
		
		// For gene conversion during sweep, handle the ancestry segments
		if (aNode->ancestryRoot) {
			// Split the tree for gene conversion tract
			gcSplitResult gcSplit = splitSegmentTreeForGeneConversion(aNode->ancestryRoot, xOver, xOver + tractL);
			lParent->ancestryRoot = gcSplit.converted;     // Gets the converted tract
			rParent->ancestryRoot = gcSplit.unconverted;   // Gets everything else
			
			// Always free the original tree since split functions create new segments
			// With our segment-based approach, we don't need to preserve the child's tree
			freeSegmentTree(aNode->ancestryRoot);
			aNode->ancestryRoot = NULL;
		} else {
			lParent->ancestryRoot = NULL;
			rParent->ancestryRoot = NULL;
		}
		
		// Update stats from the ancestry trees
		updateAncestryStatsFromTree(lParent);
		updateAncestryStatsFromTree(rParent);
			//determine sweep popn affinity
			//left side?
		//	printf("here-- popnFreq:%f sweepSite:%f xOver: %d test: %d \n",popnFreq,sweepSite,xOver,sweepSite < (float) xOver / nSites);
		if(sweepSite >= (float) xOver / nSites  && sweepSite < (float) (xOver+tractL) / nSites ){ //lParent has sweep site set to same as child node,rParent is then random
			lParent->sweepPopn = sp;
			r = ranf();
			if (r < popnFreq)rParent->sweepPopn = sp;
			else rParent->sweepPopn = (sp == 0) ? 1:0;
			
			// Propagate sweep mutation flag to lParent (has the converted tract with sweep site)
			if (aNode->carriesSweepMutation) {
				lParent->carriesSweepMutation = 1;
			}
		}
		else{ //rParent has sweep site, lParent is then random
			rParent->sweepPopn = sp;
			r = ranf();
			if (r < popnFreq) lParent->sweepPopn = sp;
			else lParent->sweepPopn = (sp == 0) ? 1:0;
			
			// Propagate sweep mutation flag to rParent (has the sweep site)
			if (aNode->carriesSweepMutation) {
				rParent->carriesSweepMutation = 1;
			}
		}
		//	printf("lParentsp:%d rParentsp:%d\n",lParent->sweepPopn,rParent->sweepPopn);
			//add in the nodes
		addNode(lParent);
		addNode(rParent);
		
		// Record edges in tskit for gene conversion during sweep (only in full mode)
		tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
		tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
		tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
		
		recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
		
		// Mark that the child node is no longer in active set
		aNode->inActiveSet = 0;
		
		// Mark parent recording if using tskit
		if (tskitOutputMode || minimalTreeSeq) {
			markParentRecorded(aNode, lParent);
			markParentRecorded(aNode, rParent);
		}
	}
}

void coalesceAtTimePopnSweep(double cTime, int popn, int sp){
	rootedNode *temp, *lChild, *rChild;
	int i;

	temp = newRootedNode(cTime,popn);

	lChild = pickNodePopnSweep(popn,sp);
	temp->leftChild = lChild;
	temp->sweepPopn = sp;
	lChild->leftParent = temp;
	lChild->branchLength = cTime - lChild->time;
	removeNode(lChild);

	rChild =  pickNodePopnSweep(popn,sp);
	temp->rightChild = rChild;
	rChild->leftParent = temp;
	rChild->branchLength = cTime - rChild->time;
	removeNode(rChild);
	
	//deal with ancMaterial
	temp->nancSites = 0;
	temp->lLim = nSites;
	temp->rLim = 0;
	
	// Propagate sweep mutation flag if either child carries it
	if (lChild->carriesSweepMutation || rChild->carriesSweepMutation) {
		temp->carriesSweepMutation = 1;
	}
	
	// Add the parent node to discoal's arrays and tskit BEFORE merging/recording
	addNode(temp);
	
	// Get tskit IDs for all nodes involved in coalescence  
	tsk_id_t parent_tsk_id = get_tskit_node_id(temp);
	
	// Merge ancestry segment trees with parent's tskit node ID
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot, parent_tsk_id);
	
	// Update stats from tree when in tree-only mode
	updateAncestryStatsFromTree(temp);
	
	// Record edges in tskit - one edge per child spanning their full range
	// Get tskit IDs for children (parent already retrieved above)
	tsk_id_t lchild_tsk_id = get_tskit_node_id(lChild);
	tsk_id_t rchild_tsk_id = get_tskit_node_id(rChild);
	
	extern int minimalTreeSeq;
	
	// Record edges using the new convenience function
	recordCoalescentEdges(temp, lChild, rChild, parent_tsk_id, lchild_tsk_id, rchild_tsk_id);
	
	// Mark parent recording and try to free nodes if using tskit
	if (tskitOutputMode || minimalTreeSeq) {
		// Use the same first/second child logic as edge recording
		rootedNode *first_child, *second_child;
		if (lchild_tsk_id < rchild_tsk_id) {
			first_child = lChild;
			second_child = rChild;
		} else {
			first_child = rChild;
			second_child = lChild;
		}
		
		// For first child - mark parent recorded
		if (minimalTreeSeq && first_child->ancestryRoot != NULL) {
			rootedNode *trueDescendants[10];
			int numDesc = findTrueDescendants(first_child, trueDescendants, 10);
			
			if (numDesc > 0 && trueDescendants[0] != first_child) {
				// Mark all nodes in the bypass chain
				rootedNode *current = first_child;
				while (current != NULL && current != trueDescendants[0]) {
					markParentRecorded(current, temp);
					rootedNode *next = (current->leftChild != NULL) ? 
					                  current->leftChild : current->rightChild;
					// Temporarily disabled node freeing
					// if (current != first_child) {
					//     tryFreeNode(current);
					// }
					current = next;
				}
				// Mark the true descendants
				for (int d = 0; d < numDesc; d++) {
					markParentRecorded(trueDescendants[d], temp);
				}
			} else {
				markParentRecorded(first_child, temp);
			}
		} else {
			markParentRecorded(first_child, temp);
		}
		
		// For second child - mark parent recorded
		if (minimalTreeSeq && second_child->ancestryRoot != NULL) {
			rootedNode *trueDescendants[10];
			int numDesc = findTrueDescendants(second_child, trueDescendants, 10);
			
			if (numDesc > 0 && trueDescendants[0] != second_child) {
				// Mark all nodes in the bypass chain
				rootedNode *current = second_child;
				while (current != NULL && current != trueDescendants[0]) {
					markParentRecorded(current, temp);
					rootedNode *next = (current->leftChild != NULL) ? 
					                  current->leftChild : current->rightChild;
					// Temporarily disabled node freeing
					// if (current != second_child) {
					//     tryFreeNode(current);
					// }
					current = next;
				}
				// Mark the true descendants
				for (int d = 0; d < numDesc; d++) {
					markParentRecorded(trueDescendants[d], temp);
				}
			} else {
				markParentRecorded(second_child, temp);
			}
		} else {
			markParentRecorded(second_child, temp);
		}
		
		// Try to free children based on their node type
		// A coalescent node has both left and right children
		int lChildIsCoal = isCoalNode(lChild);
		int rChildIsCoal = isCoalNode(rChild);
		
		// Check if nodes are leaves (no children at all)
		int lChildIsLeaf = (lChild->leftChild == NULL && lChild->rightChild == NULL);
		int rChildIsLeaf = (rChild->leftChild == NULL && rChild->rightChild == NULL);
		
		// Leaf nodes need special handling since they have no children to record
		if (lChildIsLeaf) {
			// Mark all segments as recorded since leaf nodes have simple ancestry
			if (lChild->ancestryRoot) {
				markSegmentRecorded(lChild->ancestryRoot);
			}
			lChild->isFullyRecorded = 1;
			addToRemovedList(lChild);
		}
		
		if (rChildIsLeaf) {
			// Mark all segments as recorded since leaf nodes have simple ancestry
			if (rChild->ancestryRoot) {
				markSegmentRecorded(rChild->ancestryRoot);
			}
			rChild->isFullyRecorded = 1;
			addToRemovedList(rChild);
		}
		
		// With periodic sweep approach, we don't need immediate freeing
		// Just mark nodes as removed and let the sweep handle it
	}
	
	//update active anc. material
	updateActiveMaterial(temp);
	//popnSizes[popn]--; //decrese popnSize
	//sweepPopnSizes[sp]--;	
}



/*******************************************************/
void dropMutations(){
	// Always use edge-based mutation placement with tskit
	int ret = tskit_place_mutations_edge_based(theta);
	if (ret < 0) {
		fprintf(stderr, "Error: Failed to place mutations using edge-based algorithm\n");
	}
}





/* Hash table entry for mutation duplicate detection */
typedef struct MutHashEntry {
	double mutation;
	struct MutHashEntry *next;
} MutHashEntry;

#define MUTATION_HASH_SIZE 40009  /* Prime number larger than MAXMUTS to ensure good performance */

/* Hash function for mutations - uses bit representation of double */
static inline unsigned int hashMutation(double mut) {
	union { double d; uint64_t i; } u;
	u.d = mut;
	/* Mix the bits to get better distribution */
	u.i ^= (u.i >> 32);
	u.i *= 0x9e3779b97f4a7c15ULL;  /* Golden ratio constant */
	return (unsigned int)(u.i % MUTATION_HASH_SIZE);
}




//////////////////////////////////////////////////////////
///////////////
///
///Node Maintenance Utilities
//
//
/*pickNodePopn-- picks an allele at random from active nodes from a specific popn 
*/
rootedNode *pickNodePopn(int popn){
	// Use optimized O(1) population list-based selection
	return pickNodePopnFast(popn);
}

void mergePopns(int popnSrc, int popnDest){
	// Use population lists for efficient merging
	PopulationNodeList *srcList = &popLists[popnSrc];
	
	// Create a temporary copy since we'll be modifying the list
	int nodeCount = srcList->count;
	rootedNode **tempNodes = malloc(nodeCount * sizeof(rootedNode*));
	if (!tempNodes && nodeCount > 0) {
		fprintf(stderr, "Error: Failed to allocate memory for population merge\n");
		exit(1);
	}
	
	// Copy node pointers
	for (int i = 0; i < nodeCount; i++) {
		tempNodes[i] = srcList->nodes[i];
	}
	
	// Move all nodes from source to destination
	for (int i = 0; i < nodeCount; i++) {
		rootedNode *node = tempNodes[i];
		removeNodeFromPopList(node);  // Remove from source list FIRST
		node->population = popnDest;
		addNodeToPopList(node, popnDest);
	}
	
	free(tempNodes);
	
	// Update population sizes
	popnSizes[popnDest] += popnSizes[popnSrc];
	popnSizes[popnSrc] = 0;
	
	//set migration rates to zero
	migMat[popnSrc][popnDest] = 0.0;
	migMat[popnDest][popnSrc] = 0.0;
	
}

void admixPopns(int popnSrc, int popnDest1, int popnDest2, double admixProp){
	// Use population lists for efficient admixture
	PopulationNodeList *srcList = &popLists[popnSrc];
	double rn;
	
	// If source population is empty, nothing to do
	if (srcList->count == 0) {
		// This can happen in backward-time simulations when admixture events
		// occur before any coalescent events in the source population
		return;
	}
	
	// Create a temporary copy of the nodes array since we'll be modifying the list
	int nodeCount = srcList->count;
	rootedNode **tempNodes = malloc(nodeCount * sizeof(rootedNode*));
	if (!tempNodes) {
		fprintf(stderr, "Error: Failed to allocate memory for admixture\n");
		exit(1);
	}
	
	// Copy node pointers
	for (int i = 0; i < nodeCount; i++) {
		tempNodes[i] = srcList->nodes[i];
	}
	
	// Process each node
	for (int i = 0; i < nodeCount; i++) {
		rootedNode *node = tempNodes[i];
		rn = ranf();
		if(rn < admixProp){
			removeNodeFromPopList(node);  // Remove BEFORE changing population
			node->population = popnDest1;
			addNodeToPopList(node, popnDest1);
			popnSizes[popnDest1]++;
			popnSizes[popnSrc]--;
		}
		else{
			removeNodeFromPopList(node);  // Remove BEFORE changing population
			node->population = popnDest2;
			addNodeToPopList(node, popnDest2);
			popnSizes[popnDest2]++;
			popnSizes[popnSrc]--;
		}
	}
	
	free(tempNodes);
}

//addAncientSample -- adds ancient samples by basically flipping population id of already alloced alleles and sets there time
void addAncientSample(int lineageNumber, int popnDest, double addTime, int stillSweeping, double currentFreq){
	int i;
	int count = 0;
	double rn;
	for(i=0; i < alleleNumber && count < lineageNumber; i++){
		if(nodes[i]->population == (popnDest+1) * -1){
			nodes[i]->population = popnDest;
			nodes[i]->time = addTime;
			// Add to destination population list
			addNodeToPopList(nodes[i], popnDest);
			if(stillSweeping == 1){
				rn = ranf();
				if(rn<currentFreq){
					nodes[i]->sweepPopn=1;
				}
			}
			//printf("time for %d: %f\n", i, nodes[i]->time);
			popnSizes[popnDest]++;
			count++;
		}
	}
	
}

void initializeNodeArrays() {
	// Always initialize sample node tracking (needed for each replicate)
	initializeSampleNodeIds();
	
	// Only allocate if not already allocated (prevents leaks in multiple replicates)
	if (nodes != NULL) {
		return;  // Already initialized
	}
	
	int initialCapacity = 1000;  // Start small
	
	nodesCapacity = initialCapacity;
	nodes = calloc(nodesCapacity, sizeof(rootedNode*));
	
	if (nodes == NULL) {
		fprintf(stderr, "Error: Failed to allocate initial node arrays\n");
		exit(1);
	}
}

void ensureNodesCapacity(int requiredSize) {
	if (requiredSize >= nodesCapacity) {
		int newCapacity = nodesCapacity;
		while (newCapacity <= requiredSize) {
			newCapacity *= 2;
		}
		
		rootedNode **newNodes = realloc(nodes, sizeof(rootedNode*) * newCapacity);
		if (newNodes == NULL) {
			fprintf(stderr, "Error: Failed to reallocate nodes array (requested: %d pointers)\n", newCapacity);
			exit(1);
		}
		
		// Initialize the new portion of the array to NULL
		for (int i = nodesCapacity; i < newCapacity; i++) {
			newNodes[i] = NULL;
		}
		
		nodes = newNodes;
		nodesCapacity = newCapacity;
	}
}


// Clean up nodes in the removed list (called between replicates)
void cleanupRemovedNodes() {
	for (int i = 0; i < removedNodeCount; i++) {
		if (removedNodes[i] != NULL) {
			if (removedNodes[i]->ancestryRoot != NULL) {
				freeSegmentTree(removedNodes[i]->ancestryRoot);
				removedNodes[i]->ancestryRoot = NULL;
			}
			free(removedNodes[i]);
			removedNodes[i] = NULL;
		}
	}
	removedNodeCount = 0;
}

void cleanupNodeArrays() {
	// First free any nodes still in the removed list
	cleanupRemovedNodes();
	
	if (nodes != NULL) {
		// Free any remaining nodes and their ancestry segments
		for (int i = 0; i < alleleNumber; i++) {
			if (nodes[i] != NULL) {
				// Free ancestry segments using reference counting
				if (nodes[i]->ancestryRoot != NULL) {
					freeSegmentTree(nodes[i]->ancestryRoot);
					nodes[i]->ancestryRoot = NULL;
				}
				free(nodes[i]);
				nodes[i] = NULL;
			}
		}
		free(nodes);
		nodes = NULL;
		nodesCapacity = 0;
	}
	// Cleanup sample node tracking
	cleanupSampleNodeIds();
	
	// Cleanup population lists
	cleanupPopLists();
}

// Sample node tracking functions (always available)
void initializeSampleNodeIds() {
	// Only allocate if not already allocated (prevents leaks in multiple replicates)
	if (sample_node_ids != NULL) {
		sample_node_count = 0;  // Just reset the count
		return;
	}
	
	sample_node_capacity = 1000;  // Start with reasonable size
	sample_node_count = 0;
	sample_node_ids = calloc(sample_node_capacity, sizeof(tsk_id_t));
	
	if (sample_node_ids == NULL) {
		fprintf(stderr, "Error: Failed to allocate sample node IDs array\n");
		exit(1);
	}
	
	// Initialize all entries to invalid ID
	for (int i = 0; i < sample_node_capacity; i++) {
		sample_node_ids[i] = TSK_NULL;
	}
}

void ensureSampleNodeCapacity(int required_size) {
	if (required_size >= sample_node_capacity) {
		int new_capacity = sample_node_capacity;
		while (new_capacity <= required_size) {
			new_capacity *= 2;
		}
		
		tsk_id_t *new_array = realloc(sample_node_ids, sizeof(tsk_id_t) * new_capacity);
		if (new_array == NULL) {
			fprintf(stderr, "Error: Failed to reallocate sample node IDs array\n");
			exit(1);
		}
		
		// Initialize new entries
		for (int i = sample_node_capacity; i < new_capacity; i++) {
			new_array[i] = TSK_NULL;
		}
		
		sample_node_ids = new_array;
		sample_node_capacity = new_capacity;
	}
}

void addSampleNodeId(tsk_id_t tskit_node_id) {
	ensureSampleNodeCapacity(sample_node_count + 1);
	sample_node_ids[sample_node_count] = tskit_node_id;
	sample_node_count++;
}

void cleanupSampleNodeIds() {
	if (sample_node_ids != NULL) {
		free(sample_node_ids);
		sample_node_ids = NULL;
		sample_node_capacity = 0;
		sample_node_count = 0;
	}
}

void addNode(rootedNode *aNode){
	ensureNodesCapacity(alleleNumber + 1);
	nodes[alleleNumber] = aNode;
	aNode->inActiveSet = 1;  // Mark as being in active set
	
	// Record node to tskit (unless in minimal mode and this is a recombination parent)
	extern int minimalTreeSeq;
	int isRecombParent = (aNode->leftChild != NULL && aNode->rightChild != NULL && 
	                      aNode->leftChild == aNode->rightChild);
	
	if (!minimalTreeSeq || !isRecombParent) {
		// Create tskit node for non-recombination nodes or in full mode
		tsk_id_t tsk_id = tskit_add_node(aNode->time, aNode->population, 0);  // is_sample=0
		set_tskit_node_id(aNode, tsk_id);
	} else {
		// In minimal mode, recombination parents don't get tskit nodes
		set_tskit_node_id(aNode, TSK_NULL);
	}
	
	// Add to population list
	addNodeToPopList(aNode, aNode->population);
	
	alleleNumber += 1;
	totNodeNumber += 1;
	popnSizes[aNode->population]+=1;
	if(aNode->population==0)
		sweepPopnSizes[aNode->sweepPopn]+=1;
}

//removeNodeAt-- removes a node given an index   
void removeNodeAt(int index){
	int i;
	rootedNode *temp;

	for (i = index ; i < alleleNumber - 1; i++){
		temp =  nodes[i + 1];
		nodes[i] = temp;
	}
	// NULL out the last position to avoid duplicate pointers
	nodes[alleleNumber - 1] = NULL;
	alleleNumber -= 1;
	
	
}

//removeNode -- removes a given node, uses above routine
void removeNode(rootedNode *aNode){
	int i = 0;

//find node index
	while(nodes[i] != aNode){
		i++;
	}
	popnSizes[aNode->population]-=1;
	if (aNode->population==0)
		sweepPopnSizes[aNode->sweepPopn]-=1;
	
	// Remove from population list
	removeNodeFromPopList(aNode);
	
	removeNodeAt(i);
	
	// Mark as no longer in active set
	aNode->inActiveSet = 0;
	
	// Note: Node will be added to removed list automatically when it becomes
	// fully recorded (in markParentRecorded), not here
}

/* addNodeAtIndex - variation on the theme here. adds
a rootedNode at a specific index. useful for adding pops
	or outgroups to current pops */
void addNodeAtIndex(rootedNode *aNode, int anIndex){
	nodes[anIndex] = aNode;
}

/*shiftNodes- this moves the nodes over an offset and
similiarly increases the indices alleleNumber, totNodeNumber
	by the offset. It doesn't add nodes though so fill it! */
void shiftNodes(int offset){
	int i;
	rootedNode *temp;

	for(i = alleleNumber - 1; i >= 0; i--){
		temp =  nodes[i];
		nodes[i+offset] = temp;
	}
	alleleNumber += offset;
	totNodeNumber += offset;
}

/*pickNodePopnSweep-- picks an allele at random from active nodes from specific popn and sweepPopn
*/
rootedNode *pickNodePopnSweep(int popn,int sp){
	int i, popnSize, indexArray[alleleNumber], index;

	//get popnSize -- may later decide to keep this is some sort of variable that is passed
	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
                //fprintf(stderr,"node %d in pop %d and sweepPopn %d, looking for %d and %d\n", i,nodes[i]->population,nodes[i]->sweepPopn,popn,sp);
		if(nodes[i]->sweepPopn == sp && nodes[i]->population == popn){
			indexArray[popnSize] = i;  //store indexes as we go
			popnSize++;
		}
	}
	if(popnSize == 0){
		fprintf(stderr,"error encountered in pickNodePopnSweep\n");
		fprintf(stderr,"tried to pick allele from popn %d, sweepPopn %d but popnSize is %d! Rho=%f\n",popn,sp,popnSize,rho);
		printf("popnSizes[0]:%d popnSizes[1]:%d sweepPopnSizes[0]:%d sweepPopnSizes[1]:%d\n", popnSizes[0],popnSizes[1],\
			sweepPopnSizes[0],sweepPopnSizes[1]);
		for(i = 0 ; i < alleleNumber; i++){
			printNode(nodes[i]);
		}
		exit(1);
	}
	//choose random index
	index = indexArray[ignuin(0, popnSize - 1)];

	return(nodes[index]);
}


void printNode(rootedNode *aNode){
	printf("node: %p time: %f lLim: %d rLim: %d nancSites: %d popn: %d sweepPopn: %d\n",aNode, aNode->time,aNode->lLim,\
		aNode->rLim, aNode->nancSites, aNode->population, aNode->sweepPopn);
	// ancSites array removed - ancestry now tracked via tree
}

void freeTree(rootedNode *aNode){
	int i;
	//printf("final nodeNumber = %d\n",totNodeNumber);
	//cleanup nodes
	for (i = 0; i < totNodeNumber; i++){
		if (nodes[i] == NULL) {
			continue;
		}
		// Free ancestry segment tree
		if (nodes[i]->ancestryRoot) {
			freeSegmentTree(nodes[i]->ancestryRoot);
			nodes[i]->ancestryRoot = NULL;
		}
		free(nodes[i]);
		nodes[i] = NULL;
	}
}

/*nodePopnSize-- returns popnSize of popn from
	nodes array */
int nodePopnSize(int popn){
	int i, popnSize;

	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popn){
			popnSize++;
		}
	}
	return(popnSize);
}

/*nodePopnSweepSize-- returns popnSize of popn/sweepPopn
	nodes array */
int nodePopnSweepSize(int popn, int sp){
	int i, popnSize;

	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popn && nodes[i]->sweepPopn == sp){
			popnSize++;
		}
	}
	return(popnSize);
}


void printAllActiveNodes(){
	int i;
	for(i = 0 ; i < alleleNumber; i++){
		printNode(nodes[i]);
	}
}

/********** event stuff */
int compare_events(const void *a,const void *b){
	struct event *pa = (struct event *) a;
	struct event *pb = (struct event *) b;
	if ((pa->time - pb->time) > 0 ){
		return 1;
	}
	else{
		return - 1;
	}
}

void sortEventArray(struct event *eArray, int eNumber){
	qsort(eArray, eNumber, sizeof(eArray[0]), compare_events);
}

//////// MISC
int compare_doubles(const void *a,const void *b){
	double *pa = (double *) a;
	double *pb = (double *) b;
	if ((*pa - *pb) > 0 ){
		return 1;
	}
	else{
		return -1;
	}
}
int compare_floats(const void *a,const void *b){
	float *pa = (float *) a;
	float *pb = (float *) b;
	if ((*pa - *pb) > 0 ){
		return 1;
	}
	else{
		return -1;
	}
}


/* used for getting random number seeds */
unsigned int devrand(void) {
	int fn; 
	unsigned int r; 

	fn = open("/dev/urandom", O_RDONLY); 
	if (fn == -1) 
		exit(-1); /* Failed! */ 
	if (read(fn, &r, 4) != 4) 
		exit(-1); /* Failed! */ 
	close(fn); 
	return r;
}


// Muts dynamic memory management functions

/*makeGametesMS-- MS style sample output using tskit genotype matrix */
void makeGametesMS(int argc, const char *argv[]) {
	if (tsk_tables == NULL) {
		fprintf(stderr, "Error: No tskit tables available for output generation\n");
		return;
	}
	
	// Sort tables before building tree sequence (required when nodes are freed during simulation)
	int ret = tsk_table_collection_sort(tsk_tables, NULL, 0);
	if (ret != 0) {
		fprintf(stderr, "Error: Failed to sort tables: %s\n", tsk_strerror(ret));
		return;
	}
	
	// Build trees from the table collection for genotype calculation
	tsk_treeseq_t ts;
	ret = tsk_treeseq_init(&ts, tsk_tables, TSK_TS_INIT_BUILD_INDEXES);
	if (ret != 0) {
		fprintf(stderr, "Error: Failed to build tree sequence: %s\n", tsk_strerror(ret));
		return;
	}
	
	// Get number of samples and sites
	tsk_size_t num_samples = tsk_treeseq_get_num_samples(&ts);
	tsk_size_t num_sites = tsk_treeseq_get_num_sites(&ts);
	
	
	// Print segsites header
	printf("\n//\nsegsites: %zu", (size_t)num_sites);
	
	if (num_sites > 0) {
		printf("\npositions: ");
		
		// Print positions
		for (tsk_id_t site_id = 0; site_id < (tsk_id_t)num_sites; site_id++) {
			double position = tsk_tables->sites.position[site_id];
			// Convert from absolute position back to [0,1] scale
			double relative_position = position / tsk_tables->sequence_length;
			printf("%6.6f ", relative_position);
		}
		printf("\n");
		
		// Generate and print genotype matrix
		// Use the sample_node_ids array to get the correct sample order
		extern tsk_id_t *sample_node_ids;
		extern int sample_node_count;
		
		if (sample_node_ids == NULL || sample_node_count == 0) {
			fprintf(stderr, "Error: No sample node IDs available\n");
			tsk_treeseq_free(&ts);
			return;
		}
		
		// Use tskit's efficient variant API for genotype generation
		tsk_variant_t variant;
		ret = tsk_variant_init(&variant, &ts, NULL, 0, NULL, 0);
		if (ret != 0) {
			fprintf(stderr, "Error: Failed to initialize variant: %s\n", tsk_strerror(ret));
			tsk_treeseq_free(&ts);
			return;
		}
		
		// Create a matrix to store all genotypes (samples x sites)
		int32_t **genotype_matrix = malloc(num_samples * sizeof(int32_t*));
		if (genotype_matrix == NULL) {
			fprintf(stderr, "Error: Failed to allocate genotype matrix\n");
			tsk_variant_free(&variant);
			tsk_treeseq_free(&ts);
			return;
		}
		for (tsk_size_t i = 0; i < num_samples; i++) {
			genotype_matrix[i] = malloc(num_sites * sizeof(int32_t));
			if (genotype_matrix[i] == NULL) {
				fprintf(stderr, "Error: Failed to allocate genotype row\n");
				// Cleanup already allocated rows
				for (tsk_size_t j = 0; j < i; j++) {
					free(genotype_matrix[j]);
				}
				free(genotype_matrix);
				tsk_variant_free(&variant);
				tsk_treeseq_free(&ts);
				return;
			}
		}
		
		// Decode genotypes for each site
		for (tsk_id_t site_id = 0; site_id < (tsk_id_t)num_sites; site_id++) {
			ret = tsk_variant_decode(&variant, site_id, 0);
			if (ret != 0) {
				fprintf(stderr, "Error: Failed to decode variant at site %d: %s\n", site_id, tsk_strerror(ret));
				// Cleanup
				for (tsk_size_t i = 0; i < num_samples; i++) {
					free(genotype_matrix[i]);
				}
				free(genotype_matrix);
				tsk_variant_free(&variant);
				tsk_treeseq_free(&ts);
				return;
			}
			
			// Copy genotypes for this site
			for (tsk_size_t sample_idx = 0; sample_idx < num_samples; sample_idx++) {
				genotype_matrix[sample_idx][site_id] = variant.genotypes[sample_idx];
			}
		}
		
		// Print genotypes in sample order
		for (tsk_size_t sample_idx = 0; sample_idx < num_samples; sample_idx++) {
			for (tsk_size_t site_idx = 0; site_idx < num_sites; site_idx++) {
				printf("%d", genotype_matrix[sample_idx][site_idx]);
			}
			printf("\n");
		}
		
		// Cleanup
		for (tsk_size_t i = 0; i < num_samples; i++) {
			free(genotype_matrix[i]);
		}
		free(genotype_matrix);
		tsk_variant_free(&variant);
	} else {
		printf("\n");
		// Still need to print sample lines even with no sites
		extern int sampleSize;
		for (int i = 0; i < sampleSize; i++) {
			printf("\n");
		}
	}
	
	tsk_treeseq_free(&ts);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "activeSegment.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "fenwick.h"
#include "ranlib.h"
#include "alleleTraj.h"
#include "tskitInterface.h"
#include "shapes.h"

static int64_t sweepClassWeight(int popn, int sp);

/* The lineages of each population, for O(1) uniform picks, with their
 * recombination weights. Only this file reads the weights, and only
 * through currentWeights() for totals and positions. */
typedef struct {
    rootedNode **nodes;      // Array of node pointers
    int count;              // Number of active nodes in this population
    int capacity;           // Current array capacity
    /* Recombination weight of each lineage, at the lineage's position in
     * nodes[]: the number of breakpoint positions in its span of material
     * without an MRCA (lineageWeight()). Recombination and gene conversion
     * draw a lineage, and a breakpoint in its span, from it. */
    Fenwick weights;
    /* Slots count .. staleEnd-1 may still hold the weight of a lineage
     * that left the list: removal does not clear the last slot, as the
     * next lineage added there usually overwrites it. currentWeights()
     * clears them before the weights are read as a whole. */
    int staleEnd;
} PopulationNodeList;

static PopulationNodeList popLists[MAXPOPS];

/* Slots of nodes[] in use: the active lineages and the empty slots that
 * removeNode() leaves until compactNodes(). */
static int nodesUsed = 0;
static void compactNodes(void);

// Population list management functions for O(1) node selection

/* Empty every population list (the lists start zeroed). */
void initializePopLists() {
	cleanupPopLists();
}

void cleanupPopLists() {
	for (int p = 0; p < MAXPOPS; p++) {
		free(popLists[p].nodes);
		fenwickFree(&popLists[p].weights);
		popLists[p] = (PopulationNodeList){0};
	}
}

/* The recombination weight of a lineage: the number of breakpoint
 * positions x with lLim < x <= rLim, inside its span of material without
 * an MRCA. As in Hudson's model (ms "links", msprime's recombination
 * mass), a breakpoint anywhere in the span splits the lineage, also in a
 * gap between its chunks whose sites have found their MRCA. */
static int64_t lineageWeight(const rootedNode *node) {
	if (node->nancSites <= 0 || node->rLim <= node->lLim) return 0;
	return (int64_t) node->rLim - node->lLim;
}

/* Set a lineage's recombination weight from its span, if it is in a
 * population list. */
static void setLineageWeight(const rootedNode *node) {
	if (node->popListIndex >= 0 && node->population >= 0 && node->population < MAXPOPS)
		fenwickSet(&popLists[node->population].weights, node->popListIndex,
			lineageWeight(node));
}

/* The weights of a population, with the slots of lineages that have left
 * it cleared: every reader of the total or of cumulative positions goes
 * through it. */
static Fenwick *currentWeights(int popn) {
	PopulationNodeList *list = &popLists[popn];
	for (int k = list->count; k < list->staleEnd; k++)
		fenwickSet(&list->weights, k, 0);
	list->staleEnd = list->count;
	return &list->weights;
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
	setLineageWeight(node);
}

void removeNodeFromPopList(rootedNode *node) {
	int popn = node->population;
	if (popn < 0 || popn >= MAXPOPS) return;  // Handle negative/invalid populations
	
	PopulationNodeList *list = &popLists[popn];
	int index = node->popListIndex;
	
	if (index < 0 || index >= list->count) return;
	
	// Swap with last element (and its weight)
	if (index < list->count - 1) {
		list->nodes[index] = list->nodes[list->count - 1];
		list->nodes[index]->popListIndex = index;
		fenwickSet(&list->weights, index, fenwickGet(&list->weights, list->count - 1));
	}
	/* The last slot keeps its weight until a lineage is added there or the
	 * weights are read as a whole (currentWeights()). */
	if (list->staleEnd < list->count) list->staleEnd = list->count;
	
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

/* Append a population size change with uniform priors on its time (in
 * internal units) and its size (relative to N0). Its 'n' event holds the
 * lower bounds until drawPriors() draws the values for a replicate, so the
 * preflight check of the size trajectories sees the lower bounds. `what`
 * names the source in error messages. Returns 0, or -1 after printing an
 * error. */
int addPriorSizeChange(const char *what, double timeLow, double timeHigh,
                       double sizeLow, double sizeHigh, int popID){
	if (numPriorSizeChanges >= MAX_PRIOR_SIZE_CHANGES) {
		fprintf(stderr, "Error: %s: at most %d population size changes "
			"with priors\n", what, MAX_PRIOR_SIZE_CHANGES);
		return -1;
	}
	if (!(timeLow >= 0 && timeLow <= timeHigh && sizeLow > 0 &&
	      sizeLow <= sizeHigh)) {
		fprintf(stderr, "Error: %s: the bounds must satisfy "
			"0 <= time low <= time high and 0 < size low <= size high "
			"(got time %g to %g, size %g to %g)\n", what,
			internalTimeToUser(timeLow), internalTimeToUser(timeHigh),
			sizeLow, sizeHigh);
		return -1;
	}
	int slot = numPriorSizeChanges++;
	priorSizeChanges[slot] = (struct priorSizeChange){timeLow, timeHigh,
		sizeLow, sizeHigh};
	ensureEventsCapacity();
	events[eventNumber].type = 'n';
	events[eventNumber].time = timeLow;
	events[eventNumber].popnSize = sizeLow;
	events[eventNumber].popID = popID;
	events[eventNumber].priorSlot = slot + 1;
	eventNumber++;
	return 0;
}

/* Draw this replicate's values for every parameter that has a prior, then
 * sort the events, whose times may have changed. Prior size changes are
 * drawn in the order they were given, not in event order, so the random
 * stream does not depend on how the previous replicate sorted the events. */
void drawPriors(){
	int i, k;
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
		/* Exponential with mean pRhoMean, capped (not truncated) at pRhoUp. */
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
	if(numPriorSizeChanges > 0){
		double drawnTime[MAX_PRIOR_SIZE_CHANGES], drawnSize[MAX_PRIOR_SIZE_CHANGES];
		for(k = 0; k < numPriorSizeChanges; k++){
			drawnTime[k] = genunf(priorSizeChanges[k].timeLow, priorSizeChanges[k].timeHigh);
			drawnSize[k] = genunf(priorSizeChanges[k].sizeLow, priorSizeChanges[k].sizeHigh);
		}
		for(i = 0; i < eventNumber; i++){
			k = events[i].priorSlot - 1;
			if(k >= 0){
				events[i].time = drawnTime[k];
				events[i].popnSize = drawnSize[k];
			}
		}
	}
	sortEventArray(events,eventNumber);
	/* The check before the run saw the lower bounds of the prior size
	 * changes; check the trajectories this replicate actually uses. */
	if(numPriorSizeChanges > 0 && validateShapeTrajectories(events, eventNumber) != 0)
		exit(1);
}

/* Ancient samples are taken from the start of their population's
 * samples, in event order (the events are sorted by time). Return the
 * index of the 'A' event that takes sample i of population p, or -1 if
 * the sample is taken at time 0. */
static int ancientEventOfSample(int p, int i){
	int k, first = 0;
	for(k = 0; k < eventNumber; k++){
		if(events[k].type != 'A' || events[k].popID != p) continue;
		first += events[k].lineageNumber;
		if(i < first) return k;
	}
	return -1;
}

void initialize(){
	int i,p, count=0;
	int leafID=0;

	/* Initialize the arrays */
	initializeNodeArrays();
	/* The samples go straight into nodes[], which starts with room for
	 * 1000 lineages; addNode() grows it for the lineages that follow. */
	ensureNodesCapacity(sampleSize);
	initializePopLists();  // Initialize population lists for O(1) node selection
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			int ancientEvent = ancientEventOfSample(p, i);
			double sampleTime = (ancientEvent < 0) ? 0.0 : events[ancientEvent].time;
			nodes[count] = newRootedNode(sampleTime,p);
			nodes[count]->nodesIndex = count;
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			
			// Always track sample node ID in tskit (create node first)
			tsk_id_t tsk_id = tskit_add_node(sampleTime, p, 1);
			set_tskit_node_id(nodes[count], tsk_id);
			addSampleNodeId(tsk_id);
			
			// Initialize ancestry segment tree for leaf node with tskit node ID
			nodes[count]->ancestryRoot = newSegment(0, nSites, tsk_id, NULL, NULL);

			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;

			if (ancientEvent < 0) {
				addNodeToPopList(nodes[count], p);
			} else {
				/* Dormant until its 'A' event: a negative population
				 * marks it, and addAncientSample() adds it to p then. */
				nodes[count]->population = (p + 1) * -1;
				popnSizes[p]--;
			}
			
			count += 1;
		}
	}

	// Initialize segment-based active material (all sites start as active);
	// free the previous replicate's first (the global starts zeroed)
	freeActiveMaterial(&activeMaterialSegments);
	initializeActiveMaterial(&activeMaterialSegments, nSites);
	
	//set initial counts
	alleleNumber = sampleSize;
	nodesUsed = sampleSize;
	/* Restore initial sizes; 'n'/'g'/'l' events mutate currentSize per replicate. */
	for (p = 0; p < MAXPOPS; p++) currentSize[p] = currentSizeConst[p];

	//initialize shape state for time-varying parameter framework
	initializeShapesFromGlobals();
	activeSites = nSites;
	if (npops>1){
		if(tDiv==666 && migFlag == 0){
			fprintf(stderr,"tDiv or migration not set in population split model\n");
			exit(1);
		}
	}
	
	drawPriors();
}

/* Lineages freed in minimal mode, kept for reuse: each event frees and
 * makes lineages, and a list of them saves a malloc and a free per
 * lineage. Linked through leftParent; released by releaseFreeLineages(). */
static rootedNode *freeLineages = NULL;

void releaseFreeLineages(void) {
	while (freeLineages) {
		rootedNode *next = freeLineages->leftParent;
		free(freeLineages);
		freeLineages = next;
	}
}

rootedNode *newRootedNode(double cTime, int popn) {
	rootedNode *temp;

	if (freeLineages) {
		temp = freeLineages;
		freeLineages = temp->leftParent;
	} else {
		temp = malloc(sizeof(rootedNode));
		if (temp == NULL) {
			fprintf(stderr, "Error: out of memory for a lineage\n");
			exit(1);
		}
	}
	temp->rightParent = NULL;
	temp->leftParent = NULL;
	temp->rightChild = NULL;
	temp->leftChild = NULL; 
	temp->time = cTime;
	temp->branchLength=0.0;
	temp->population = popn;
	/* -1 means "not yet classified into a sweep background". Sample-init,
	 * sweepPhaseEvents* (via nodePopnSweepSize), and the sweep coalescence/
	 * recombination paths overwrite this with 0 or 1. Consumers that update
	 * sweepPopnSizes[] must guard with `sweepPopn >= 0` to avoid an OOB
	 * write for nodes that haven't been classified yet (everything in a
	 * no-sweep run, plus pre-sweep internal nodes). */
	temp->sweepPopn = -1;
	temp->nodesIndex = -1;
	/* No material yet: an empty span, as emptySpan(). Events that give the
	 * node material set its span (setLineageSpan()). */
	temp->nancSites = 0;
	temp->lLim = nSites;
	temp->rLim = 0;
	
	// Initialize node state tracking
	
	// Initialize ancestry segment tree to NULL (will be set during initialization)
	temp->ancestryRoot = NULL;

	// Initialize tskit node ID
	temp->tskit_node_id = TSK_NULL;
	
	// Initialize population list index
	temp->popListIndex = -1;

	return temp;
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
	if (srcPopn==0 && temp->sweepPopn >= 0)
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

////////

// Structure to hold a node and its ancestry segments for edge recording

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

    // Record edges for both children
    recordEdgesForChild(parent_tsk_id, lChild, lchild_tsk_id, minimalTreeSeq);
    recordEdgesForChild(parent_tsk_id, rChild, rchild_tsk_id, minimalTreeSeq);
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
            seg = seg->next;
        }
    }
    
    // Right parent gets segments (right of crossover for recombination, unconverted segments for gene conversion)
    if (rparent_tsk_id != TSK_NULL && child_tsk_id != TSK_NULL && rParent->ancestryRoot != NULL) {
        AncestrySegment *seg = rParent->ancestryRoot;
        while (seg != NULL) {
            tskit_add_edges(rparent_tsk_id, child_tsk_id,
                (double)seg->start, (double)seg->end);
            seg = seg->next;
        }
    }
}

/* Minimal tree sequence mode records what msprime records: a node and
 * edges only where two lineages' ancestral material coalesces. Each
 * segment carries the tskit node of the lineage it belongs to, so a
 * lineage's node (rootedNode) is not needed once it leaves the active set,
 * and is freed then. */

/* The tskit node of a coalescence parent, made when its children's
 * ancestral material first overlaps (mergeAncestryTreesCoalescing). */
static tsk_id_t coalescenceParentNode(void *ctx){
	rootedNode *parent = ctx;
	tsk_id_t id = tskit_add_node(parent->time, parent->population, 0);
	if (id < 0) {
		fprintf(stderr, "Error: failed to add a tree sequence node\n");
		exit(1);
	}
	set_tskit_node_id(parent, id);
	return id;
}

static void coalescenceEdge(tsk_id_t parent, tsk_id_t child, double left, double right){
	if (tskit_add_edges(parent, child, left, right) != 0) {
		fprintf(stderr, "Error: failed to add a tree sequence edge\n");
		exit(1);
	}
}

/* Set a lineage's span of material without an MRCA (nancSites, lLim,
 * rLim), and its recombination weight if it is already in a population
 * list. */
static void setLineageSpan(rootedNode *node, const SegmentSpan *span) {
	node->nancSites = span->nancSites;
	node->lLim = span->lLim;
	node->rLim = span->rLim;
	setLineageWeight(node);
}

static SegmentSpan emptySpan(void) {
	SegmentSpan span = {0, nSites, 0};
	return span;
}

/* Free a lineage that has left the active set, after the event that
 * removed it has recorded its edges: its segments are freed and its struct
 * goes to the list for reuse by newRootedNode(). */
static void freeLineage(rootedNode *node){
	if (node->ancestryRoot) freeSegmentTree(node->ancestryRoot);
	node->leftParent = freeLineages;
	freeLineages = node;
}

/* Sites whose MRCA a coalescence has just found are no longer active. */
static void sitesFixed(int start, int end){
	removeFixedRegion(&activeMaterialSegments, start, end);
}

/* The parent of a common-ancestor event takes the merged ancestry of its
 * two children (minimal mode); the children are freed. The merge reports
 * the sites that reach their MRCA here, so the active sites are updated
 * without walking the parent's segments (updateActiveMaterial()), whose
 * fixed ones were removed when they were fixed. */
static void coalesceLineagesMinimal(rootedNode *parent, rootedNode *lChild, rootedNode *rChild){
	SegmentSpan span = emptySpan();
	SegmentSpan lSpan = {lChild->nancSites, lChild->lLim, lChild->rLim};
	SegmentSpan rSpan = {rChild->nancSites, rChild->lLim, rChild->rLim};
	parent->ancestryRoot = mergeAncestryTreesCoalescing(lChild->ancestryRoot,
		rChild->ancestryRoot, coalescenceParentNode, coalescenceEdge, parent,
		sampleSize, sitesFixed, &lSpan, &rSpan, &span);
	/* The merge took over the children's segments. */
	lChild->ancestryRoot = rChild->ancestryRoot = NULL;
	activeSites = getActiveSiteCount(&activeMaterialSegments);
	setLineageSpan(parent, &span);
	parent->leftChild = parent->rightChild = NULL;
	freeLineage(lChild);
	freeLineage(rChild);
}

void coalesceAtTimePopn(double cTime, int popn){
	// Check if population has at least 2 nodes for coalescence
	if (popLists[popn].count < 2) {
		// Cannot coalesce with fewer than 2 nodes
		return;
	}
	
	rootedNode *temp, *lChild, *rChild;
	
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

	/* Minimal mode merges first, so that the parent enters its population
	 * list with its span and its weight is set once. The full ARG adds the
	 * parent first: the merge records it as the edges' parent. */
	if (minimalTreeSeq) {
		coalesceLineagesMinimal(temp, lChild, rChild);
		addNode(temp);
		return;
	}
	addNode(temp);
	
	// Get tskit IDs for all nodes involved in coalescence  
	tsk_id_t parent_tsk_id = get_tskit_node_id(temp);
	tsk_id_t lchild_tsk_id = get_tskit_node_id(lChild);
	tsk_id_t rchild_tsk_id = get_tskit_node_id(rChild);
	
	// Merge ancestry segment trees with parent's tskit node ID
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot, parent_tsk_id);
	
	// Update stats from tree when in tree-only mode
	updateAncestryStatsFromTree(temp);

	// Record edges using the new convenience function
	recordCoalescentEdges(temp, lChild, rChild, parent_tsk_id, lchild_tsk_id, rchild_tsk_id);

	/* The children's edges are recorded and the parent has its own copy of
	 * their segments, so nothing reads the children again. */
	temp->leftChild = temp->rightChild = NULL;
	freeLineage(lChild);
	freeLineage(rChild);
	
	//update active anc. material (addNode was moved earlier)
	updateActiveMaterial(temp);
}

/* The span of a lineage from a walk of its segment list. Splits and
 * merges in minimal mode get it during their own walk instead. */
void updateAncestryStatsFromTree(rootedNode *node) {
	if (!node || !node->ancestryRoot) return;
	SegmentSpan span = emptySpan();
	for (AncestrySegment *seg = node->ancestryRoot; seg; seg = seg->next)
		spanAdd(&span, seg->start, seg->end, seg->count, sampleSize);
	setLineageSpan(node, &span);
}

/* Split a recombining lineage's segments at xOver between its two
 * parents, and set the parents' spans. The full ARG gives the parents'
 * segments their own tskit nodes, so it copies the halves. Minimal mode
 * cuts the child's list in place, as msprime does, and the parents then
 * own its segments. In both modes the child's ancestryRoot is NULL after
 * the split. */
static void splitLineage(rootedNode *aNode, rootedNode *lParent, rootedNode *rParent, int xOver){
	if (minimalTreeSeq) {
		SegmentSpan lSpan = emptySpan(), rSpan = emptySpan();
		SegmentSpan whole = {aNode->nancSites, aNode->lLim, aNode->rLim};
		splitInPlace(aNode->ancestryRoot, xOver, sampleSize, &whole,
			&lParent->ancestryRoot, &rParent->ancestryRoot, &lSpan, &rSpan);
		aNode->ancestryRoot = NULL;
		setLineageSpan(lParent, &lSpan);
		setLineageSpan(rParent, &rSpan);
		return;
	}
	lParent->ancestryRoot = splitLeftWithParent(aNode->ancestryRoot, xOver, lParent->tskit_node_id);
	rParent->ancestryRoot = splitRightWithParent(aNode->ancestryRoot, xOver, rParent->tskit_node_id);
	updateAncestryStatsFromTree(lParent);
	updateAncestryStatsFromTree(rParent);
	/* The parents have copies; the child's segments are not needed. */
	freeSegmentTree(aNode->ancestryRoot);
	aNode->ancestryRoot = NULL;
}

/*updateActiveMaterial- does bookkeeping on the material that has found it's MRCA-- added for efficiency */
void updateActiveMaterial(rootedNode *aNode){
	// Use segment-based approach to update active material
	// This removes regions where the node has ancestry from all samples
	updateActiveMaterialFromAncestry(&activeMaterialSegments, aNode->ancestryRoot, 
	                                sampleSize, nSites);
	
	// Update global active site count
	activeSites = getActiveSiteCount(&activeMaterialSegments);
	
}

/* 1 if aNode carries material at site (a fraction of the locus) that has
 * not found its MRCA: some, but not all, samples descend from it there. */
int isAncestralHere(rootedNode *aNode, float site){
	int bp = floor(site * nSites);
	uint16_t count = aNode->ancestryRoot ? getAncestryCount(aNode->ancestryRoot, bp) : 0;
	return count > 0 && count < sampleSize;
}

/* Rates of recombination and gene conversion in a population, per unit of
 * time (2N0 generations). Crossovers and gene conversion tracts start at
 * rho / 2 and gamma / 2 over the locus, so at rho / (2 nSites) and
 * gamma / (2 nSites) per site. A lineage has a crossover or an internal
 * tract at each position of its span (its weight w).
 *
 * It also has tracts that start at or left of its first site and run into
 * its material (left tracts, drawLeftTract()), at gcMean sites' worth of
 * rate. Only those that end in its span (with probability 1 - pc^w, pc =
 * 1 - 1/gcMean) change anything, so their rate is gcMean (1 - pc^w) sites,
 * which is at most min(gcMean, w). Left tracts are drawn at the bound
 * min(gcMean x lineages, total weight) and thinned (leftTractLineage()):
 * with gcMean much larger than the spans, most left tracts would cover
 * all the material and do nothing. */
static double leftTractBound(int64_t spanWeight, int lineages){
	const double perLineage = (double) gcMean * lineages;
	return perLineage < (double) spanWeight ? perLineage : (double) spanWeight;
}

static double gcWeight(int64_t spanWeight, int lineages){
	return (double) spanWeight + leftTractBound(spanWeight, lineages);
}

/* 1 - pc^w: the share of left tracts that end in a span of w positions. */
static double leftTractEffective(int64_t w){
	if (w <= 0) return 0.0;
	if (gcMean <= 1) return 1.0;
	return -expm1((double) w * log1p(-1.0 / gcMean));
}

void lineageEventRates(int popn, double *recombination, double *geneConversion){
	int64_t total = fenwickTotal(currentWeights(popn));
	*recombination = rho * 0.5 * (double) total / nSites;
	*geneConversion = my_gamma * 0.5 * gcWeight(total, popLists[popn].count) / nSites;
}

/* A lineage of population popn and a position in its span, drawn with
 * probability proportional to the lineage's weight: one uniform position
 * over all the population's weights gives both. Returns NULL if no lineage
 * has a span. */
static rootedNode *drawLineageAndPosition(int popn, int *position){
	Fenwick *w = currentWeights(popn);
	int64_t total = fenwickTotal(w);
	if (total <= 0) return NULL;
	int64_t at = (int64_t) (ranf() * (double) total);
	if (at >= total) at = total - 1;
	int64_t offset;
	rootedNode *node = popLists[popn].nodes[fenwickFind(w, at, &offset)];
	*position = node->lLim + 1 + (int) offset;
	return node;
}

int recombineAtTimePopn(double cTime, int popn){
	rootedNode *aNode, *lParent, *rParent;
	int xOver;
	
	/* The breakpoint is in the lineage's span, so both parents get
	 * material without an MRCA. */
	aNode = drawLineageAndPosition(popn, &xOver);
	if (aNode != NULL){
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
		
		splitLineage(aNode, lParent, rParent, xOver);
		
		addNode(lParent);
		addNode(rParent);
		
		// Record edges in tskit for recombination (only in full mode)
		tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
		tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
		tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
		
		recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
		/* The parents hold the child's segments (or, in the full ARG, copies
		 * of them) and its edges are recorded: the child is not needed. */
		lParent->leftChild = rParent->leftChild = NULL;
		freeLineage(aNode);
		return xOver;
	}
	return -1;  /* no lineage has a span */
}

/* A gene conversion tract length: geometric on 1, 2, 3, ... with mean
 * gcMean. With a mean of 1 every tract is one site; the inversion formula
 * gives 0 there (log(1 - 1) is -inf), so the length is at least 1. The
 * uniform is drawn in every case, so that the random number stream does
 * not depend on the mean. A tract longer than the locus has the same
 * effect at any length, so the length is at most nSites + 1; this keeps
 * the tract ends far from integer overflow for a very large mean. */
static int drawTractLength(void){
	double u = genunf(0,1);
	if (gcMean <= 1) return 1;
	double length = ceil(log(u)/log(1.0-(1.0/gcMean)));
	if (length < 1) return 1;
	if (length > nSites + 1) return nSites + 1;
	return (int) length;
}

/* The span of the material without an MRCA in a list of segments. */
static SegmentSpan spanOfSegments(const AncestrySegment *seg){
	SegmentSpan span = emptySpan();
	for (; seg; seg = seg->next)
		spanAdd(&span, seg->start, seg->end, seg->count, sampleSize);
	return span;
}

/* A tract that starts at or left of a lineage's first site without an
 * MRCA (lLim) and ends in its span of w positions, as msprime's left gene
 * conversion events (and ms's): the chromosome goes on to the left of the
 * locus, and of the lineage's material. A tract that starts j sites left
 * of lLim reaches it if its length L > j. Given that, it covers lLim ..
 * lLim + K - 1 and starts at lLim - j, with K >= 1 and j >= 0 independent
 * and geometric with the parameter of L, as the geometric has no memory;
 * K is drawn given K <= w, so that the tract ends in the span. The start
 * matters only in a sweep, for the parent that gets the sweep site. */
static void drawLeftTract(const rootedNode *node, int64_t w, int *start, int *end){
	int k = 1;
	if (gcMean > 1) {
		const double u = genunf(0,1);
		k = 1 + (int) floor(log1p(-leftTractEffective(w) * u) / log1p(-1.0 / gcMean));
		if (k < 1) k = 1;
		if (k > w) k = (int) w;
	}
	*end = node->lLim + k;
	*start = node->lLim - (drawTractLength() - 1);
}

/* A lineage for a left tract, drawn at the bound of leftTractBound() and
 * thinned to the rate gcMean (1 - pc^w): uniformly with uniform() when
 * gcMean x lineages is the bound, else by weight with byWeight(). Returns
 * NULL if the draw is thinned out; *w is the lineage's weight. */
static rootedNode *leftTractLineage(int64_t totalWeight, int lineages,
		rootedNode *(*uniform)(int, int), rootedNode *(*byWeight)(int, int, int *),
		int popn, int sp, int64_t *w){
	rootedNode *aNode;
	double accept;
	int unused;
	if ((double) gcMean * lineages <= (double) totalWeight) {
		aNode = uniform(popn, sp);
		if (aNode == NULL) return NULL;
		*w = lineageWeight(aNode);
		accept = leftTractEffective(*w);
	} else {
		aNode = byWeight(popn, sp, &unused);
		if (aNode == NULL) return NULL;
		*w = lineageWeight(aNode);
		accept = (double) gcMean * leftTractEffective(*w) / (double) *w;
	}
	return (ranf() < accept) ? aNode : NULL;
}

/* Gene conversion of aNode in population popn with the tract [start, end):
 * the converted parent (lParent) gets the material in the tract, the other
 * parent (rParent) the rest. If one of the two would get no material
 * without an MRCA, as when the tract falls in a gap of the material or
 * covers all of it, the event changes nothing, as in msprime, and aNode
 * stays. In a sweep (sp >= 0, population 0) the parent with the sweep site
 * keeps aNode's class sp, and the other is a random draw from the
 * population: it joins class sp with probability popnFreq, the frequency
 * of that class. A sweep outside the locus (sweepSite < 0) is in no tract. */
static void convertTract(double cTime, int popn, rootedNode *aNode, int start, int end,
                         int sp, double sweepSite, double popnFreq){
	rootedNode *lParent, *rParent;

	if (aNode == NULL || aNode->ancestryRoot == NULL) return;
	gcSplitResult split = splitSegmentTreeForGeneConversion(aNode->ancestryRoot, start, end);
	SegmentSpan convertedSpan = spanOfSegments(split.converted);
	SegmentSpan restSpan = spanOfSegments(split.unconverted);
	if (convertedSpan.nancSites == 0 || restSpan.nancSites == 0) {
		freeSegmentTree(split.converted);
		freeSegmentTree(split.unconverted);
		return;
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

	/* The split made new segments; the child's are not needed. */
	lParent->ancestryRoot = split.converted;
	rParent->ancestryRoot = split.unconverted;
	freeSegmentTree(aNode->ancestryRoot);
	aNode->ancestryRoot = NULL;
	setLineageSpan(lParent, &convertedSpan);
	setLineageSpan(rParent, &restSpan);

	if (sp >= 0) {
		const int other = (sp == 0) ? 1 : 0;
		const int inTract = sweepSite >= 0.0 && sweepSite >= (float) start / nSites
			&& sweepSite < (float) end / nSites;
		rootedNode *withSite = inTract ? lParent : rParent;
		rootedNode *drawn = inTract ? rParent : lParent;
		withSite->sweepPopn = sp;
		drawn->sweepPopn = (ranf() < popnFreq) ? sp : other;
	}

	addNode(lParent);
	addNode(rParent);

	tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
	tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
	tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
	recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
	/* As in recombineAtTimePopn(): the child is not needed. */
	lParent->leftChild = rParent->leftChild = NULL;
	freeLineage(aNode);
}

/* Lineage draws of population popn with the signatures of the sweep ones
 * (the class argument is not used), for leftTractLineage(). */
static rootedNode *uniformInPopulation(int popn, int sp){
	(void) sp;
	return pickNodePopn(popn);
}

static rootedNode *byWeightInPopulation(int popn, int sp, int *position){
	(void) sp;
	return drawLineageAndPosition(popn, position);
}

/* A gene conversion event in population popn: a tract that starts in a
 * lineage's span (a position drawn by weight), or one that comes in from
 * the left of its first site (a lineage drawn uniformly), in proportion
 * to their rates (gcWeight()). */
void geneConversionAtTimePopn(double cTime, int popn){
	const double internal = (double) fenwickTotal(currentWeights(popn));
	const double total = gcWeight((int64_t) internal, popLists[popn].count);
	rootedNode *aNode;
	int start, end;

	if (total <= 0) return;
	if (ranf() * total < internal) {
		aNode = drawLineageAndPosition(popn, &start);
		if (aNode == NULL) return;
		end = start + drawTractLength();
	} else {
		int64_t w;
		aNode = leftTractLineage((int64_t) internal, popLists[popn].count,
			uniformInPopulation, byWeightInPopulation, popn, -1, &w);
		if (aNode == NULL) return;
		drawLeftTract(aNode, w, &start, &end);
	}
	convertTract(cTime, popn, aNode, start, end, -1, 0.0, 0.0);
}

/* Draw the next event time across all rate components using NHPP per-component
 * draws. Returns the minimum waiting time and writes the firing event class
 * + index into *winnerKind / *winnerArg.
 *
 * winnerKind values:
 *   0 = recomb       (winnerArg unused; caller samples pop proportional to rRate[i])
 *   1 = gene conv    (winnerArg unused; caller samples pop proportional to gcRate[i])
 *   2 = migration    (winnerArg = src_pop * MAXPOPS + dst_pop)
 *   3 = coalescence  (winnerArg = popID)
 *
 * Recombination and gene-conversion rates are constant in time within an
 * interval (depend only on lineage counts and rho/gamma), so an Exp(total)
 * draw is exact under any shape mix. Coalescence and migration rates can
 * vary in time when shapes are non-constant; we use closed-form
 * drawWaitingTime accessors per component.
 *
 * cRate and mRate are accepted in the signature for symmetry with the
 * caller's rate vector and possible future use; they are not consumed
 * directly by this routine. */
static double drawNHPPWaitingTime(double *cRate, double *rRate, double *gcRate,
                                  double *mRate, double totRRate, double totGCRate,
                                  double t0, int *winnerKind, int *winnerArg) {
    double T_min = 1e300;  /* effectively +infinity */
    int i, j;
    *winnerKind = -1;
    *winnerArg = -1;

    /* Recombination: aggregate exponential draw across all populations */
    if (totRRate > 0.0) {
        double T = -log(ranf()) / totRRate;
        if (T < T_min) {
            T_min = T;
            *winnerKind = 0;
            *winnerArg = -1;
        }
    }

    /* Gene conversion: aggregate exponential draw across all populations */
    if (totGCRate > 0.0) {
        double T = -log(ranf()) / totGCRate;
        if (T < T_min) {
            T_min = T;
            *winnerKind = 1;
            *winnerArg = -1;
        }
    }

    /* Coalescence per population */
    for (i = 0; i < npops; i++) {
        if (popnSizes[i] < 2) continue;
        double xi = -log(ranf());
        double T = drawWaitingTimeSize(i, t0, xi, popnSizes[i]);
        if (T > 0.0 && T < T_min) {
            T_min = T;
            *winnerKind = 3;
            *winnerArg = i;
        }
    }

    /* Migration per ordered pair */
    for (i = 0; i < npops; i++) {
        if (popnSizes[i] < 1) continue;
        for (j = 0; j < npops; j++) {
            if (i == j) continue;
            double m = migAt(i, j, t0);
            if (m <= 0.0) continue;
            double xi = -log(ranf());
            double T = drawWaitingTimeMig(i, j, t0, xi, popnSizes[i]);
            if (T > 0.0 && T < T_min) {
                T_min = T;
                *winnerKind = 2;
                *winnerArg = i * MAXPOPS + j;
            }
        }
    }

    (void)cRate; (void)mRate;
    return T_min;
}

/* Coalescence rate for k lineages in a population of relative size `size`:
 * k choose 2 over size. With fewer than two lineages the rate is 0 even when
 * the size is 0, as it is for an ancestral deme before its lineages arrive;
 * the plain formula would give 0/0 = NaN, and a NaN event time never reaches
 * the end of the interval, so the event loop would never finish. */
static inline double pairCoalRate(int k, double size){
	return (k < 2) ? 0.0 : k * (k - 1) * 0.5 / size;
}

/* Migration reachability in the final interval, after the last event:
 * migReach[i][j] is 1 if a lineage in population i can reach population j
 * through migration with positive rate. No event changes the model after
 * the last one, so it is computed once when the final interval starts. (A
 * linear migration shape that reaches zero is still counted as positive.) */
static unsigned char migReach[MAXPOPS][MAXPOPS];

static void computeMigReach(double t){
	int i, j, k;
	for (i = 0; i < npops; i++)
		for (j = 0; j < npops; j++)
			migReach[i][j] = (i == j) || migAt(i, j, t) > 0.0;
	for (k = 0; k < npops; k++)
		for (i = 0; i < npops; i++)
			if (migReach[i][k])
				for (j = 0; j < npops; j++)
					if (migReach[k][j]) migReach[i][j] = 1;
}

/* 1 if population i holds a lineage with ancestral material (a site that
 * has no common ancestor yet). Other lineages may stay apart forever. */
static int holdsAncestralMaterial(int i){
	for (int k = 0; k < popLists[i].count; k++)
		if (popLists[i].nodes[k]->nancSites > 0) return 1;
	return 0;
}

/* In the final interval, lineages in different populations can coalesce
 * only if some population can be reached from every population that holds
 * lineages with ancestral material. If none can, they never coalesce: the
 * run would never end. Exit with an error, as msprime does ("infinite
 * waiting time"). */
static void requireLineagesCanMeet(double t){
	int c, i;
	int holds[MAXPOPS];
	for (i = 0; i < npops; i++) holds[i] = holdsAncestralMaterial(i);
	for (c = 0; c < npops; c++) {
		for (i = 0; i < npops; i++)
			if (holds[i] && !migReach[i][c]) break;
		if (i == npops) return;
	}
	/* t is internal time (2N generations); report it in the 4N units of the
	 * command line. */
	fprintf(stderr, "Error: after the last demographic event (at time %g in "
	        "units of 4N generations), lineages with ancestral material remain "
	        "in populations", t / USER_TO_INTERNAL_TIME);
	for (i = 0; i < npops; i++)
		if (holds[i]) fprintf(stderr, " %d", i);
	fprintf(stderr, ", and no population can be reached from all of them by "
	        "migration, so they can never coalesce. Check the migration, merge "
	        "and pulse events of the model.\n");
	exit(1);
}

/*neutralPhaseGeneralPopNumber--coalescent, recombination, gc events until
specified time. returns endTime. can handle multiple popns*/
double neutralPhaseGeneralPopNumber(double startTime, double endTime, double *sizeRatio){
	double cTime, cRate[npops], rRate[npops], gcRate[npops], mRate[npops],totRate, waitTime, r, r2;
	double totCRate, totRRate, totGCRate,totMRate, eSum;
	int  i,j;
	int coalescenceCounter = 0;  // Counter for coalescent events only
	const int EDGE_FLUSH_INTERVAL = 10;  // flush the edge buffer every 10 coalescences
	/* Sizes and migration rates change only at events, which come between
	 * calls of this sampler. */
	const int constantShapes = allShapesConstant();

	if(startTime == endTime){
		return(endTime);
	}
	/* The last interval (no later event) is the one with no end. Only sites
	 * without a common ancestor matter: lineages that carry no ancestral
	 * material may stay apart forever. */
	int finalInterval = (isinf(endTime) && npops > 1 && activeSites > 0);
	if (finalInterval) {
		computeMigReach(startTime);
		requireLineagesCanMeet(startTime);
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
		for(i=0;i<npops;i++){
			cRate[i] = pairCoalRate(popnSizes[i], sizeAt(i, currentTime));
			lineageEventRates(i, &rRate[i], &gcRate[i]);

			for(j=0;j<npops;j++) mRate[i]+=migAt(i, j, currentTime);
			mRate[i] *= popnSizes[i] * MIG_LINEAGE_WEIGHT;
			totCRate += cRate[i];
			totRRate += rRate[i];
			totMRate += mRate[i];
			totGCRate += gcRate[i];
			totRate += cRate[i] + rRate[i] + mRate[i] + gcRate[i];
		}

		if (constantShapes) {
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
					recombineAtTimePopn(cTime,i);
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
							//pick dest popn
							eSum = migAt(i, 0, currentTime) * popnSizes[i] * MIG_LINEAGE_WEIGHT;
							r2 = ranf();
							while(eSum/mRate[i] < r2){
								eSum += migAt(i, ++j, currentTime) * popnSizes[i] * MIG_LINEAGE_WEIGHT;

							}
							migrateAtTime(cTime,i,j);
							if (finalInterval) requireLineagesCanMeet(cTime);
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

							/* Flush the buffered edges (squashed per flush) every
							 * EDGE_FLUSH_INTERVAL coalescences. */
							coalescenceCounter++;
							if (coalescenceCounter % EDGE_FLUSH_INTERVAL == 0) {
								extern int tskit_flush_edges_periodic(void);
								tskit_flush_edges_periodic();
							}

						}
					}
				}
			}
		}
		else {
			/* Non-constant shapes: NHPP per-component sampler. Each rate
			 * component draws its own waiting time using the closed-form
			 * inverse hazard for size/migration shapes; recomb / gc are
			 * time-constant within an interval so an Exp(total) draw is
			 * exact for those classes. The minimum across components is
			 * the next event; its identity is returned via winnerKind /
			 * winnerArg, so we dispatch directly without a second
			 * within-class lottery for migration / coalescence. */
			int winnerKind = -1, winnerArg = -1;
			waitTime = drawNHPPWaitingTime(cRate, rRate, gcRate, mRate,
			                               totRRate, totGCRate, cTime,
			                               &winnerKind, &winnerArg);
			cTime += waitTime;

			/* winnerKind < 0: no event can ever occur from here. */
			if (winnerKind < 0 || cTime >= endTime){
				return(endTime);
			}

			if (winnerKind == 0) {
				/* recombination: pick population proportional to rRate */
				eSum = rRate[0];
				i = 0;
				r2 = ranf();
				while(eSum/totRRate < r2) eSum += rRate[++i];
				recombineAtTimePopn(cTime,i);
			}
			else if (winnerKind == 1) {
				/* gene conversion: pick population proportional to gcRate */
				eSum = gcRate[0];
				i = 0;
				r2 = ranf();
				while(eSum/totGCRate < r2) eSum += gcRate[++i];
				geneConversionAtTimePopn(cTime,i);
			}
			else if (winnerKind == 2) {
				/* migration: src/dst already chosen by the per-pair NHPP draw */
				int src = winnerArg / MAXPOPS;
				int dst = winnerArg % MAXPOPS;
				migrateAtTime(cTime, src, dst);
				if (finalInterval) requireLineagesCanMeet(cTime);
			}
			else if (winnerKind == 3) {
				/* coalescence: population already chosen by per-pop NHPP draw */
				int popID = winnerArg;
				coalesceAtTimePopn(cTime, popID);

				/* Flush the buffered edges every EDGE_FLUSH_INTERVAL coalescences. */
				coalescenceCounter++;
				if (coalescenceCounter % EDGE_FLUSH_INTERVAL == 0) {
					extern int tskit_flush_edges_periodic(void);
					tskit_flush_edges_periodic();
				}
			}
			else {
				/* No event drawable (all rates zero / unreachable). Advance
				 * to endTime to terminate the phase cleanly. */
				return(endTime);
			}
		}
	}
	
	/* Flush the buffered edges at the end of the phase. */
	extern int tskit_flush_edges_periodic(void);
	tskit_flush_edges_periodic();
	
	return(cTime);
}

/* Unmap and close the mapped sweep trajectory, if any. Its file is
 * anonymous (tmpfile), so closing it, or the end of the process however it
 * ends, removes it. */
void releaseTrajectory(void) {
	if (trajectoryFd != -1) {
		if (currentTrajectory && currentTrajectory != MAP_FAILED) {
			munmap(currentTrajectory, trajectoryFileSize);
		}
		close(trajectoryFd);
		trajectoryFd = -1;
		currentTrajectory = NULL;
	}
}

/* Map an accepted trajectory for the sweep sampler. fd is an open file of
 * numSteps floats (an anonymous file from proposeTrajectory()); this takes
 * it over, and releases the trajectory mapped before. */
void mmapAcceptedTrajectory(int fd, long int numSteps) {
	releaseTrajectory();
	trajectoryFd = fd;
	trajectoryFileSize = numSteps * sizeof(float);
	/* A trajectory with no steps (the sweep ends at once) has nothing to
	 * read, and mmap cannot map zero bytes. */
	if (numSteps == 0) {
		currentTrajectory = NULL;
		return;
	}
	currentTrajectory = (float *)mmap(NULL, trajectoryFileSize, PROT_READ,
	                                  MAP_PRIVATE, trajectoryFd, 0);
	if (currentTrajectory == MAP_FAILED) {
		perror("Failed to mmap trajectory file");
		exit(1);
	}
}

/* The frequency at which a sweep trajectory, followed back in time, ends:
 * the selected allele's origin. A sweep from a new mutation ends at one
 * copy, 1/(2N). A sweep from standing variation (f0 above 1/(2N) when the
 * sweep starts) ends with the neutral history below f0, which is followed
 * back to frequency 0: that history, reversed, starts at the allele's
 * origin, and stopping at the first return to 1/(2N) would cut it short. */
static double trajectoryEnd(int fromStandingVariation, double N){
	return fromStandingVariation ? 0.0 : 1.0/(2.*N);
}

/* Whether the current sweep is from standing variation (f0 above one copy
 * at the sweep's start). proposeTrajectory() decides it once per sweep, and
 * the sweep sampler reads it, so a size change between the sampler's calls
 * cannot turn a soft sweep into a hard one, which would end the sweep
 * before its trajectory does. */
static int sweepFromStandingVariation = 0;

/* A partial sweep starts (going back in time) with its lineages in
 * population 0 split into the two sweep classes. The selected allele is
 * then at partialSweepFinalFreq, carried by the beneficial class. Add a
 * site at the selected position with a mutation on each beneficial
 * lineage that has ancestral material there; the samples below these
 * lineages at the site carry the allele. Each sweep, recurrent ones
 * included, gets its own site. -h (hidePartialSNP) leaves it out, as does
 * a selected site outside the locus. */
/* The selected site of the partial sweep in progress (as a fraction of
 * the locus), and its tskit site once a mutation has been added there. */
static double partialSweepSite = -1.0;
static tsk_id_t partialSweepSiteId = TSK_NULL;

/* Mark a lineage as a carrier of the partial sweep's selected allele: a
 * mutation at the selected site. In the full ARG every lineage has its own
 * tskit node; in minimal mode a lineage may have none, and the material at
 * the site belongs to the node of the lineage's segment there (the node in
 * which it last coalesced). The site is created with its first mutation. */
static void addPartialSweepCarrier(rootedNode *lineage){
	tsk_id_t node;
	if (minimalTreeSeq) {
		AncestrySegment *seg = segmentAt(lineage->ancestryRoot,
			(int) floor(partialSweepSite * nSites));
		node = seg ? seg->tskit_node_id : TSK_NULL;
	} else {
		node = get_tskit_node_id(lineage);
	}
	if (node == TSK_NULL)
		return;
	if (partialSweepSiteId == TSK_NULL) {
		partialSweepSiteId = tskit_add_site(partialSweepSite * nSites, "0");
		if (partialSweepSiteId < 0) {
			fprintf(stderr, "Error: failed to add the partial sweep site\n");
			exit(1);
		}
	}
	if (tskit_add_mutation(partialSweepSiteId, node, "1") < 0) {
		fprintf(stderr, "Error: failed to add the partial sweep mutation\n");
		exit(1);
	}
}

static void recordPartialSweepSite(double site){
	compactNodes();
	int i;
	partialSweepSite = (hidePartialSNP || site < 0.0) ? -1.0 : site;
	partialSweepSiteId = TSK_NULL;
	if (partialSweepSite < 0.0)
		return;
	for(i = 0; i < alleleNumber; i++){
		if(nodes[i]->population == 0 && nodes[i]->sweepPopn == 1 &&
		   isAncestralHere(nodes[i], site))
			addPartialSweepCarrier(nodes[i]);
	}
}

/* Start a sweep (going back in time): put the lineages of population 0
 * into the sweep classes (sweepPopn 1 is the beneficial class). In a
 * partial sweep each lineage is in the beneficial class with probability
 * partialSweepFinalFreq, and the selected site is recorded; otherwise
 * every lineage is. */
static void startSweep(double sweepSite){
	compactNodes();
	int i;
	for(i=0;i<alleleNumber;i++){
		if(nodes[i]->population!=0)
			continue;
		if(partialSweepMode == 1)
			nodes[i]->sweepPopn = (ranf() > partialSweepFinalFreq) ? 0 : 1;
		else
			nodes[i]->sweepPopn = 1;
	}
	if(partialSweepMode == 1)
		recordPartialSweepSite(sweepSite);
}

/* The selected allele's frequency at the start of a sweep (going back in
 * time): one copy short of fixation in population 0 at time t, or the
 * final frequency of a partial sweep if that is lower. The size is
 * sizeAt(0, t), since currentSize is not updated for growth. */
double sweepStartFreq(double t){
	if (floor(EFFECTIVE_POPN_SIZE * sizeAt(0, t)) < 1.0) {
		fprintf(stderr, "Error: population 0 has fewer than one individual "
			"(N x size = %g) at the start of a sweep (time %g in 4N0 "
			"generations)\n", EFFECTIVE_POPN_SIZE * sizeAt(0, t),
			t / USER_TO_INTERNAL_TIME);
		exit(1);
	}
	const double nearlyFixed = 1.0 - (1.0 / (2.0 * EFFECTIVE_POPN_SIZE * sizeAt(0, t)));
	return partialSweepMode ? MIN(partialSweepFinalFreq, nearlyFixed) : nearlyFixed;
}

/* The -S trajectory, built by buildSweepTrajectoryFromFile(); NULL
 * without -S. */
static float *fileTrajectory = NULL;
static long int fileTrajectorySteps = 0;

/* The file of the last proposed trajectory, not yet accepted; -1 if none. */
static int proposedFd = -1;

/* Proposals for one sweep before discoal stops: at this point the
 * acceptance probability is so small that the run would take hours. */
#define MAX_TRAJECTORY_PROPOSALS 1000000L

/* The bound for the acceptance of the trajectories of a sweep that starts
 * at currentTime (see sampleSweepTrajectory()): the largest size of
 * population 0, relative to N, at any time an origin can fall in. A
 * trajectory has at most TRAJSTEPSTART steps of 1/(deltaTMod N) (the
 * proposal stops with an error after that), so origins are within that
 * span of the start, and the bound is the larger of the size at the start
 * and the size changes of population 0 in the span. Population 0 does
 * not grow continuously at or after a sweep (getParameters() rejects
 * it), so its size changes only at these events. */
static double trajectoryAcceptanceBound(double currentTime){
	double horizon = currentTime +
		(double) TRAJSTEPSTART / (deltaTMod * EFFECTIVE_POPN_SIZE);
	double bound = sizeAt(0, currentTime);
	int i;
	for (i = 0; i < eventNumber; i++) {
		if (events[i].type == 'n' && events[i].popID == 0 &&
		    events[i].time >= currentTime && events[i].time < horizon &&
		    events[i].popnSize > bound)
			bound = events[i].popnSize;
	}
	return bound;
}

/* Draw a sweep trajectory by rejection sampling: propose trajectories
 * until one is accepted, then map the accepted one for the sweep samplers
 * (which also releases the trajectory of the previous sweep, with
 * recurrent sweeps). The arguments are those of proposeTrajectory().
 *
 * New beneficial mutations arise at a rate proportional to the size of
 * population 0, so a trajectory is accepted with probability
 * N(origin) / bound, with N(origin) the size at the allele's origin. For
 * exact rejection sampling the bound must be the same for every proposal
 * and at least N(origin) for every trajectory that can be drawn
 * (trajectoryAcceptanceBound()). A bound over each trajectory alone would
 * be too small for a trajectory whose origin is recent, if population 0
 * was larger further in the past, and would accept such trajectories too
 * often. */
void sampleSweepTrajectory(int currentEventNumber, double *sizeRatio,
	char sweepMode, double initialFreq, double *finalFreq, double alpha,
	double f0, double currentTime){
	/* -S: the trajectory from the file, built once. It is not a mapped
	 * file, so trajectoryFd stays -1 and nothing unmaps it. */
	if (fileTrajectory != NULL) {
		sweepFromStandingVariation = 0;  /* the file ends at the origin */
		currentTrajectory = fileTrajectory;
		currentTrajectoryStep = 0;
		totalTrajectorySteps = fileTrajectorySteps;
		return;
	}
	/* Each proposal closes the one before it, so a rejected one is gone. */
	const double bound = trajectoryAcceptanceBound(currentTime);
	long proposals = 0;
	double originSize;
	do {
		if (++proposals > MAX_TRAJECTORY_PROPOSALS) {
			fprintf(stderr, "Error: %ld sweep trajectories were proposed and "
				"none was accepted. A trajectory is accepted with probability "
				"N(origin) / N(max): the size of population 0 at the allele's "
				"origin over its largest size (%g x N) at or before the start "
				"of the sweep. Make the size changes of population 0 during "
				"and before the sweep smaller.\n",
				MAX_TRAJECTORY_PROPOSALS, bound);
			exit(1);
		}
		originSize = proposeTrajectory(currentEventNumber, currentTrajectory,
			sizeRatio, sweepMode, initialFreq, finalFreq, alpha, f0,
			currentTime);
		assert(originSize <= bound);
	} while (ranf() > originSize / bound);
	mmapAcceptedTrajectory(proposedFd, totalTrajectorySteps);
	proposedFd = -1;
}

/* -S: read a sweep trajectory from a text file with one point per line:
 * the time before the present (4N0 generations) and the frequency of the
 * selected allele, separated by white space. Blank lines and lines that
 * start with '#' are skipped. Times must increase down the file; every
 * frequency but the last must be in (0, 1], and the last must be 0, the
 * allele's origin, which ends the sweep. Adds the sweep event at the first
 * time and sets the sweep up as a stochastic one (sweepMode 's').
 * getParameters() makes the sweep partial if the first frequency is below
 * 1, and builds the trajectory the sampler reads. Returns 0, or -1 after
 * printing an error. */
int setSweepTrajectoryFile(const char *path){
	FILE *in = fopen(path, "r");
	if (in == NULL) {
		fprintf(stderr, "Error: cannot open the trajectory file %s\n", path);
		return -1;
	}
	struct trajPoint *points = NULL;
	int rows = 0, capacity = 0, lineNumber = 0;
	const char *error = NULL;
	char line[1024], c;
	while (error == NULL && fgets(line, sizeof(line), in) != NULL) {
		lineNumber++;
		if (sscanf(line, " %c", &c) != 1 || c == '#')
			continue;
		double t, x;
		if (sscanf(line, "%lf %lf %c", &t, &x, &c) != 2)
			error = "expected a time and a frequency";
		else if (!isfinite(t) || t < 0.0 || (rows > 0 && t <= points[rows-1].time))
			error = "times must be >= 0 and increase down the file";
		else if (!(x >= 0.0 && x <= 1.0) || (rows > 0 && points[rows-1].freq == 0.0))
			error = "frequencies must be in [0, 1], and only the last can be 0";
		else {
			if (rows == capacity) {
				capacity = capacity ? 2 * capacity : 64;
				struct trajPoint *grown = realloc(points, sizeof(*points) * capacity);
				if (grown == NULL) {
					error = "out of memory";
					break;
				}
				points = grown;
			}
			points[rows].time = t;
			points[rows].freq = x;
			rows++;
		}
	}
	fclose(in);
	if (error != NULL)
		fprintf(stderr, "Error: %s, line %d: %s\n", path, lineNumber, error);
	else if (rows < 2 || points[rows-1].freq != 0.0)
		fprintf(stderr, "Error: %s: a trajectory needs at least two rows, and "
			"its last frequency must be 0 (the allele's origin)\n", path);
	else {
		for (int r = 0; r < rows; r++)
			points[r].time = userTimeToInternal(points[r].time);
		free(trajFilePoints);
		trajFilePoints = points;
		trajFileRows = rows;
		sweepMode = 's';
		ensureEventsCapacity();
		events[eventNumber].time = points[0].time;
		events[eventNumber].type = 's';
		eventNumber++;
		return 0;
	}
	free(points);
	return -1;
}

/* The -S trajectory on the sweep sampler's grid. The sampler reads value k
 * at sweep start + (k+1) dt, with dt = 1/(deltaTMod N); each value is the
 * file's trajectory interpolated linearly at that time, and the last is 0,
 * at or after the file's last time. Values are capped at 1 - 1/(2N), where
 * a simulated sweep starts: at frequency 1 the wild-type class has no
 * copies, and its coalescence rate is 0/0. The sampler ends the sweep at
 * the first value at or below 1/(2N), as for a simulated sweep. The
 * sweep's start, N and dt are the same in every replicate, so
 * getParameters() builds this once, and sampleSweepTrajectory() hands it
 * to the sampler. */
void buildSweepTrajectoryFromFile(void){
	const double dt = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	const double start = trajFilePoints[0].time;
	const double span = trajFilePoints[trajFileRows-1].time - start;
	if (span / dt >= 500000000.0) {
		fprintf(stderr, "Error: the trajectory file spans more than "
			"500000000 sweep time steps; use a smaller -N or -i\n");
		exit(1);
	}
	long int steps = (long int)ceil(span / dt);
	if (steps < 1) steps = 1;
	float *grid = malloc(sizeof(float) * steps);
	if (grid == NULL) {
		fprintf(stderr, "Error: out of memory for the trajectory file\n");
		exit(1);
	}
	/* nextafterf: with a large N, 1 - 1/(2N) rounds to 1 as a float. */
	const float cap = fminf((float)(1.0 - 1.0 / (2.0 * EFFECTIVE_POPN_SIZE)),
	                        nextafterf(1.0f, 0.0f));
	int r = 0;
	for (long int k = 0; k < steps; k++) {
		double t = start + (k + 1) * dt;
		while (r < trajFileRows - 1 && trajFilePoints[r+1].time <= t) r++;
		if (r == trajFileRows - 1) {
			grid[k] = 0.0f;
		} else {
			const struct trajPoint *p = &trajFilePoints[r];
			double w = (t - p[0].time) / (p[1].time - p[0].time);
			grid[k] = fminf((float)(p[0].freq + w * (p[1].freq - p[0].freq)), cap);
		}
	}
	grid[steps - 1] = 0.0f;
	if (grid[0] <= 1.0 / (2.0 * EFFECTIVE_POPN_SIZE)) {
		fprintf(stderr, "Error: the trajectory file falls to 1/(2N) or below "
			"within one sweep time step (1/(%g N) in 2N0 generations), so the "
			"sweep would have no effect\n", deltaTMod);
		exit(1);
	}
	free(fileTrajectory);
	fileTrajectory = grid;
	fileTrajectorySteps = steps;
}

/*proposeTrajectory-- this function creates a sweep trajectory and deals with
complications like changing population size, or soft sweeps, etc.
Returns the size of population 0 (relative to N) at the trajectory's end,
the allele's origin. sampleSweepTrajectory() accepts or rejects the
trajectory. */
double proposeTrajectory(int currentEventNumber, float *currentTrajectory, double *sizeRatio, char sweepMode, \
double initialFreq, double *finalFreq, double alpha, double f0, double currentTime)
{
	double tInc, tIncOrig, minF,ttau, N;
	double N_0 = (double) EFFECTIVE_POPN_SIZE;
	double localNextTime, currentSizeRatio;
	int i, insweepphase;
	long int j;
	double x;

	/* Integrated alpha_eff for the closed-form general detSweepFreq under
	 * non-constant shapes. Reset to 0 at function entry; updated in the
	 * inner loop's case 'd' when the dispatch chooses the general path. */
	double A_now = 0.0;
	double A_prev = 0.0;

	/* Save popShape state; we'll mutate it as we walk events forward to track
	 * 'n' and 'g' events, then restore at function exit so the caller's view
	 * is unchanged. */
	Shape saved_popShape[MAXPOPS];
	memcpy(saved_popShape, popShape, sizeof(saved_popShape));

	/* The trajectory goes to an anonymous temporary file (tmpfile): it has
	 * no name in the file system, so it disappears when it is closed or the
	 * process ends, however it ends. The previous proposal, not accepted,
	 * is closed here. */
	if (proposedFd != -1) {
		close(proposedFd);
		proposedFd = -1;
	}
	FILE *trajFile = tmpfile();
	if (!trajFile) {
		perror("Failed to create trajectory file");
		exit(1);
	}
	
	// Use buffered writes for efficiency
	float writeBuffer[1024];
	int bufferPos = 0;
	
	tIncOrig = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	j=0;
	/* Population 0's size at the sweep's start (sizeRatio[0] is not
	 * updated for growth). */
	currentSizeRatio = sizeAt(0, currentTime);
	N = (double) floor(N_0 * currentSizeRatio);
	x = initialFreq;
	minF = f0;
	const int fromStandingVariation = (f0 > 1.0/(2.*N));
	sweepFromStandingVariation = fromStandingVariation;
	insweepphase = 1;
	ttau=0.0;
	//go through each event until last one
	for(i=currentEventNumber;i<eventNumber;i++){
		if(i == eventNumber - 1){
			localNextTime = MAXTIME;
			}
		else{
			localNextTime = events[i+1].time;
		}
		if(events[i].type == 'n'){
			popShape[events[i].popID].type = SHAPE_CONSTANT;
			popShape[events[i].popID].anchor_value = events[i].popnSize;
			popShape[events[i].popID].rate_param = 0.0;
			popShape[events[i].popID].anchor_time = events[i].time;
			/* Only population 0 sweeps: the size of another population
			 * must not set the trajectory's N, its end or its acceptance. */
			if(events[i].popID == 0){
				currentSizeRatio = events[i].popnSize;
				N = floor(N_0 * events[i].popnSize);
			}
		}
		if(events[i].type == 'g'){
			popShape[events[i].popID].type = SHAPE_EXPONENTIAL;
			popShape[events[i].popID].anchor_value = sizeAt(events[i].popID, events[i].time);
			popShape[events[i].popID].rate_param = events[i].popnSize;
			popShape[events[i].popID].anchor_time = events[i].time;
		}
		if(events[i].type == 'l'){
			popShape[events[i].popID].type = SHAPE_LINEAR;
			popShape[events[i].popID].anchor_value = sizeAt(events[i].popID, events[i].time);
			popShape[events[i].popID].rate_param = events[i].popnSize;
			popShape[events[i].popID].anchor_time = events[i].time;
		}
		/* The frequency below which the trajectory is the allele's neutral
		 * history: f0 for a sweep from standing variation, and one copy,
		 * 1/(2N) with the N of this epoch, for a new mutation. During the
		 * sweep phase it follows N in both directions: if N is larger
		 * further back, the selected allele is still selected down to one
		 * copy of that larger population. In the neutral phase it is the
		 * ceiling of the neutral history, which must stay above x: it can
		 * only go up there. */
		if (insweepphase)
			minF = fromStandingVariation ? fmax(f0, 1.0/(2.*N)) : 1.0/(2.*N);
		else if (minF < 1.0/(2.*N))
			minF = 1.0/(2.*N);
		tInc = 1.0 / (deltaTMod * N);
		/* Shapes change only at the events applied above, not in this
		 * interval. */
		const int constantShapes = allShapesConstant();
		//iterate until epoch time or sweep freq
		while( x > trajectoryEnd(fromStandingVariation, N) && (currentTime+ttau) < localNextTime){
			ttau += tIncOrig;
			double sr_now = sizeAt(0, currentTime + ttau);
			N = floor(N_0 * sr_now);
			tInc = 1.0 / (deltaTMod * N);
			if(x > minF && insweepphase){
				//get next sweep allele freq
				switch(sweepMode){
					case 'd':
					if (constantShapes) {
						x = detSweepFreq(ttau, alpha * sr_now);
					} else {
						/* Time-varying alpha_eff: use the closed-form general formula with the
						 * incrementally-tracked integrated alpha. A_now = A_prev + alpha *
						 * integratedSizeRatio over the most recent dt step. */
						A_now = A_prev + alpha * integratedSizeRatio(0, currentTime + ttau - tIncOrig, tIncOrig);
						x = detSweepFreqGeneral(alpha, A_now);
						A_prev = A_now;
					}
					break;
					case 's':
					x = 1.0 - genicSelectionStochasticForwardsOptimized(tInc, (1.0 - x), alpha * sr_now);
					break;
					case 'N':
					x = neutralStochasticOptimized(tInc, x);
					break;
				}
			}
			else{
				/* Below f0: the neutral history of the allele before it
				 * became beneficial, back to its origin. */
				insweepphase = 0;
				x = neutralStochasticBelowCeiling(tInc, x, minF);
			}
			
			// Check trajectory size to prevent runaway
			if (j >= TRAJSTEPSTART) {
				fprintf(stderr, "Error: the sweep trajectory is too long (%ld steps); check the sweep parameters\n", j);
				fclose(trajFile);
				exit(1);
			}
			
			// Write to buffer
			writeBuffer[bufferPos++] = (float)x;
			if (bufferPos >= 1024) {
				// Flush buffer to file
				fwrite(writeBuffer, sizeof(float), bufferPos, trajFile);
				bufferPos = 0;
			}
			j++;
		}
		/* The trajectory has ended: later events are not part of it, and a
		 * later size change must not set its acceptance probability. */
		if (x <= trajectoryEnd(fromStandingVariation, N)) break;
	}
	
	// Flush any remaining data in buffer
	if (bufferPos > 0) {
		fwrite(writeBuffer, sizeof(float), bufferPos, trajFile);
	}
	/* Keep the file open for mapping (a second descriptor keeps it alive
	 * after the stream is closed). A write error (a full disk) stops the
	 * run here: a short file would fail later, when it is read through
	 * the mapping. */
	if (fflush(trajFile) != 0 || ferror(trajFile)
	    || (proposedFd = dup(fileno(trajFile))) == -1) {
		perror("Failed to write trajectory file");
		exit(1);
	}
	fclose(trajFile);
	currentTrajectoryStep = 0;
	totalTrajectorySteps = j;
	
	// Note: We don't mmap here because this function may be called multiple times
	// during rejection sampling. The accepted trajectory will be mmap'd later.

	/* Restore popShape state — caller expects it unchanged. */
	memcpy(popShape, saved_popShape, sizeof(saved_popShape));

	return currentSizeRatio;

}

/* The number of ancient sample lineages of population popn not yet added
 * (dormant: they join at their sampling time). */
static int dormantLineagesIn(int popn){
	compactNodes();
	int i, count = 0;
	for(i = 0; i < alleleNumber; i++)
		if(nodes[i]->population == (popn + 1) * -1)
			count++;
	return count;
}

/*sweepPhaseEventsConditionalTrajectory-- does sweep phase with trajectory created externally */
double sweepPhaseEventsConditionalTrajectory(double startTime, double endTime, double sweepSite,\
double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
double *sizeRatio, char sweepMode,double f0, double uA)
{

	double totRate;
	double  ttau, x, tIncOrig;
	double pCoalB, pCoalb, pRecB, pRecb, r, sum,eventRand,eventProb, pLeftRecB, pLeftRecb;
	double pRecurMut, pGCB, pGCb;
	double sweepPopTotRate,cRate[npops], rRate[npops], gcRate[npops];
	double totCRate, totRRate, totGCRate, eSum, r2;
	double N = (double) EFFECTIVE_POPN_SIZE;
	double minF;
	double cTime = startTime;
	int i;

	//initialize stuff
	pCoalB = pCoalb = pRecB = pRecb = totRate = pRecurMut = pLeftRecB = pLeftRecb = totGCRate = totCRate = totRRate = 0;
	sweepPopTotRate = pGCB = pGCb = 0;

	N = (double) floor(N * sizeAt(0, cTime));
	ttau = 0.0;
	x = initialFreq;
	minF = f0;
	if(minF < 1.0/(2.*N))
		minF = 1.0/(2.*N);

	if(!*stillSweeping){
		startSweep(sweepSite);
		*stillSweeping = 1;
	}
	
	//assume that sweep always happens in popn 0!!!
	//using popnSize global to manage bookkeeping
	sweepPopnSizes[1] = nodePopnSweepSize(0,1);
	sweepPopnSizes[0] = nodePopnSweepSize(0,0);

	//set time increment
	tIncOrig = 1.0 / (deltaTMod * EFFECTIVE_POPN_SIZE);
	
	//go for epoch time, sweep freq, or root
	/* The sweep ends at trajectoryEnd() for population 0's size at each
	 * step, as in proposeTrajectory(), so the two end on the same step. */
	const int fromStandingVariation = sweepFromStandingVariation;
	double xEnd = trajectoryEnd(fromStandingVariation, N);
	/* Ancient samples of population 0 that join later keep the sweep going
	 * even when fewer than two present-day lineages are left in it, or none
	 * in the beneficial class: they must see the sweep when they join. */
	const int ancientToCome = dormantLineagesIn(0);
	int trajectoryDone = 0;
	while( x > xEnd && (cTime+ttau) < endTime && (popnSizes[0] > 1 || ancientToCome > 0)){ 
		//rejection algorithm of Braverman et al. 1995
		eventRand = ranf();
		eventProb = 1.0;
		/* The recombination weights of the two classes change only at
		 * events, not while waiting for one. */
		const double weightB = (double) sweepClassWeight(0, 1);
		const double weightb = (double) sweepClassWeight(0, 0);
		//wait for something
		while(eventProb > eventRand && x > xEnd && (cTime+ttau) < endTime ){
			ttau += tIncOrig;

			/* The trajectory ended: the sweep ends at its last step. xEnd
			 * follows N as the trajectory's end test does, but the time here
			 * restarts after each event and can round differently, so this
			 * guards against reading past the trajectory. */
			if (currentTrajectoryStep >= totalTrajectorySteps) {
				ttau -= tIncOrig;
				trajectoryDone = 1;
				break;
			}
		x = currentTrajectory[currentTrajectoryStep++];
		const double size0 = sizeAt(0, cTime + ttau);
		xEnd = trajectoryEnd(fromStandingVariation, floor(EFFECTIVE_POPN_SIZE * size0));
		/* The origin of an allele from standing variation: no rates at
		 * frequency 0. (A sweep from a new mutation ends above 0.) */
		if (x <= 0.0) break;

			//calculate event probs
			//first 4 events are probs of events in population 0
			pCoalB = ((sweepPopnSizes[1] * (sweepPopnSizes[1] - 1) ) * 0.5)/x*tIncOrig / size0;
			pCoalb = ((sweepPopnSizes[0] * (sweepPopnSizes[0] - 1) ) * 0.5)/(1-x)*tIncOrig / size0;
			pRecB = rho * 0.5 * weightB / nSites * tIncOrig;
			pRecb = rho * 0.5 * weightb / nSites * tIncOrig;
			pGCB = my_gamma * 0.5 * gcWeight((int64_t) weightB, sweepPopnSizes[1]) / nSites * tIncOrig;
			pGCb = my_gamma * 0.5 * gcWeight((int64_t) weightb, sweepPopnSizes[0]) / nSites * tIncOrig;
			pRecurMut = (uA * sweepPopnSizes[1]*0.5 *tIncOrig)/x;///sizeRatio[0];
			if (sweepSite < 0.0){
				pLeftRecB = leftRho * sweepPopnSizes[1]*0.5 * tIncOrig * (1-x);
				pLeftRecb = leftRho * sweepPopnSizes[0]*0.5 * tIncOrig * x;
                        }
			sweepPopTotRate = pCoalB + pCoalb + pRecB + pRecb + pGCB + pGCb + pRecurMut + pLeftRecB + pLeftRecb;
			//now two events in population 1
			totRate = 0.0;
			totCRate = 0.0;
			totRRate = 0.0;
			totGCRate = 0.0;
			totRate = sweepPopTotRate;

			for(i=1;i<npops;i++){
				cRate[i] = pairCoalRate(popnSizes[i], sizeAt(i, cTime + ttau)) * tIncOrig;
				lineageEventRates(i, &rRate[i], &gcRate[i]);
				rRate[i] *= tIncOrig;
				gcRate[i] *= tIncOrig;
				totCRate += cRate[i];
				totRRate += rRate[i];
				totGCRate += gcRate[i];
				totRate += cRate[i] + rRate[i] + gcRate[i];
			}

			eventProb *= 1-totRate;

		}
		if(cTime+ttau >= endTime){
			/* The next event comes first (for example ancient samples, which
			 * need the allele's frequency now); the sweep goes on after it. */
			*finalFreq = x;
			return(cTime+ttau);
		}
		/* The trajectory reached the allele's origin, or its end, with no
		 * event. */
		if (x <= 0.0 || trajectoryDone) break;

		//Decide which event took place
		//first was it in population 0 (sweep) or not
		double eventTime = cTime + ttau;
		if(ranf() < sweepPopTotRate / totRate){
			/* An event in population 0, chosen in proportion to its rate among
			 * the sweep events (sweepPopTotRate). */
			r = ranf();
			if (r < (sum = pCoalB) / sweepPopTotRate)
				coalesceAtTimePopnSweep(eventTime, 0, 1);
			else if (r < (sum += pCoalb) / sweepPopTotRate)
				coalesceAtTimePopnSweep(eventTime, 0, 0);
			else if (r < (sum += pRecb) / sweepPopTotRate) {
				recombineAtTimePopnSweep(eventTime, 0, 0, sweepSite, (1.0-x));
			}
			else if (r < (sum += pRecB) / sweepPopTotRate) {
				recombineAtTimePopnSweep(eventTime, 0, 1, sweepSite, x);
			}
			else if (r < (sum += pGCB) / sweepPopTotRate)
				geneConversionAtTimePopnSweep(eventTime, 0, 1, sweepSite, x);
			else if (r < (sum += pGCb) / sweepPopTotRate)
				geneConversionAtTimePopnSweep(eventTime, 0, 0, sweepSite, 1.0 - x);
			else if (r < (sum += pLeftRecb) / sweepPopTotRate)
				recombineToLeftPopnSweep(0, 0, x);
			else if (r < (sum += pLeftRecB) / sweepPopTotRate)
				recombineToLeftPopnSweep(0, 1, x);
			else
				/* recurrent adaptive mutation: a lineage leaves the
				 * beneficial class */
				recurrentMutAtTime(eventTime, 0, 1);
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
						recombineAtTimePopn(eventTime,i);
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
	if(sweepPopnSizes[1]==0 && ancientToCome == 0) *stillSweeping = 0;

	*finalFreq = x;
	return(cTime+(ttau));
}

/*recurrentSweepPhaseGeneralPopNumber--coalescent, recombination, gc events, and sweeps! until
specified time. returns endTime. can handle multiple popns*/
double recurrentSweepPhaseGeneralPopNumber(double startTime, double endTime, double *finalFreq, double alpha, char sweepMode, double *sizeRatio){
	double cTime, cRate[npops], rRate[npops], gcRate[npops], mRate[npops],totRate, waitTime, r, r2;
	double totCRate, totRRate, totGCRate,totMRate, eSum,curSweepSite, initFreq;
	int  i,j;

	if(startTime == endTime){
		return(endTime);
	}
	/* The last interval has no end: stop with an error if the lineages can
	 * never meet, as the neutral sampler does, rather than run for ever. */
	int finalInterval = (isinf(endTime) && npops > 1 && activeSites > 0);
	if (finalInterval) {
		computeMigReach(startTime);
		requireLineagesCanMeet(startTime);
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
		for(i=0;i<npops;i++){
			cRate[i] = pairCoalRate(popnSizes[i], sizeAt(i, cTime));
			lineageEventRates(i, &rRate[i], &gcRate[i]);
			
			for(j=0;j<npops;j++) mRate[i]+=migAt(i, j, cTime);
			mRate[i] *= popnSizes[i] * MIG_LINEAGE_WEIGHT;
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
				recombineAtTimePopn(cTime,i);
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
						//pick dest popn
						eSum = migAt(i, 0, cTime) * popnSizes[i] * MIG_LINEAGE_WEIGHT;
						r2 = ranf();
						while(eSum/mRate[i] < r2){
							eSum += migAt(i, ++j, cTime) * popnSizes[i] * MIG_LINEAGE_WEIGHT;
							
						} 
						migrateAtTime(cTime,i,j);
						if (finalInterval) requireLineagesCanMeet(cTime);
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
							initFreq = sweepStartFreq(cTime);
							sampleSweepTrajectory(currentEventNumber, sizeRatio, sweepMode,
								initFreq, finalFreq, alpha, f0, cTime);
							cTime= sweepPhaseEventsConditionalTrajectory(cTime, endTime, curSweepSite,\
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
/* The total recombination weight of the lineages of population popn in
 * sweep class sp (lineageWeight()). */
static int64_t sweepClassWeight(int popn, int sp){
	int64_t w = 0;
	PopulationNodeList *list = &popLists[popn];
	for (int k = 0; k < list->count; k++)
		if (list->nodes[k]->sweepPopn == sp) w += fenwickGet(&list->weights, k);
	return w;
}

/* As drawLineageAndPosition(), among the lineages of sweep class sp: one
 * scan of the population's list (the sweep phase picks lineages by class
 * in O(n) per event). */
static rootedNode *drawSweepLineageAndPosition(int popn, int sp, int *position){
	PopulationNodeList *list = &popLists[popn];
	int64_t total = sweepClassWeight(popn, sp);
	if (total <= 0) return NULL;
	int64_t at = (int64_t) (ranf() * (double) total);
	if (at >= total) at = total - 1;
	for (int k = 0; k < list->count; k++) {
		if (list->nodes[k]->sweepPopn != sp) continue;
		int64_t w = fenwickGet(&list->weights, k);
		if (at < w) {
			*position = list->nodes[k]->lLim + 1 + (int) at;
			return list->nodes[k];
		}
		at -= w;
	}
	return NULL;
}

int recombineAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq){
	rootedNode *aNode, *lParent, *rParent;
	int xOver;
	double r;

		/* The breakpoint is in the lineage's span (drawn by weight). */
		aNode = drawSweepLineageAndPosition(popn, sp, &xOver);
		if (aNode != NULL){
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
			
			//determine sweep popn affinity
			//left side?
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
			
			splitLineage(aNode, lParent, rParent, xOver);
			
			//add in the nodes
			addNode(lParent);
			addNode(rParent);
			
			// Record edges in tskit for recombination during sweep (only in full mode)
			tsk_id_t child_tsk_id = get_tskit_node_id(aNode);
			tsk_id_t lparent_tsk_id = get_tskit_node_id(lParent);
			tsk_id_t rparent_tsk_id = get_tskit_node_id(rParent);
			
			recordRecombinationEdges(aNode, lParent, rParent, child_tsk_id, lparent_tsk_id, rparent_tsk_id);
			/* As in recombineAtTimePopn(): the child is not needed. */
			lParent->leftChild = rParent->leftChild = NULL;
			freeLineage(aNode);
			return xOver;
		}
		return -1;  /* no lineage of this class has a span */
}

int recombineToLeftPopnSweep(int popn, int sp, double popnFreq){
	rootedNode *aNode;

	aNode = pickNodePopnSweep(popn, sp);
	sweepPopnSizes[aNode->sweepPopn]--;

	aNode->sweepPopn = (sp == 0) ? 1:0;
	sweepPopnSizes[aNode->sweepPopn]++;
	return 0;
}

/* A gene conversion event in sweep class sp of population popn, as
 * geneConversionAtTimePopn() among the lineages of that class; the
 * parents' classes follow the sweep site (convertTract()). */
void geneConversionAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq){
	const double internal = (double) sweepClassWeight(popn, sp);
	const double total = gcWeight((int64_t) internal, sweepPopnSizes[sp]);
	rootedNode *aNode;
	int start, end;

	if (total <= 0) return;
	if (ranf() * total < internal) {
		aNode = drawSweepLineageAndPosition(popn, sp, &start);
		if (aNode == NULL) return;
		end = start + drawTractLength();
	} else {
		int64_t w;
		aNode = leftTractLineage((int64_t) internal, sweepPopnSizes[sp],
			pickNodePopnSweep, drawSweepLineageAndPosition, popn, sp, &w);
		if (aNode == NULL) return;
		drawLeftTract(aNode, w, &start, &end);
	}
	convertTract(cTime, popn, aNode, start, end, sp, sweepSite, popnFreq);
}

void coalesceAtTimePopnSweep(double cTime, int popn, int sp){
	rootedNode *temp, *lChild, *rChild;

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

	/* Minimal mode merges first, so that the parent enters its population
	 * list with its span and its weight is set once. The full ARG adds the
	 * parent first: the merge records it as the edges' parent. */
	if (minimalTreeSeq) {
		coalesceLineagesMinimal(temp, lChild, rChild);
		addNode(temp);
		return;
	}
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

	// Record edges using the new convenience function
	recordCoalescentEdges(temp, lChild, rChild, parent_tsk_id, lchild_tsk_id, rchild_tsk_id);
	
	/* The children's edges are recorded and the parent has its own copy of
	 * their segments, so nothing reads the children again. */
	temp->leftChild = temp->rightChild = NULL;
	freeLineage(lChild);
	freeLineage(rChild);
	
	//update active anc. material
	updateActiveMaterial(temp);
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
	migShape[popnSrc][popnDest].type = SHAPE_CONSTANT;
	migShape[popnSrc][popnDest].anchor_value = 0.0;
	migShape[popnSrc][popnDest].rate_param = 0.0;
	migShape[popnSrc][popnDest].anchor_time = currentTime;
	migShape[popnDest][popnSrc].type = SHAPE_CONSTANT;
	migShape[popnDest][popnSrc].anchor_value = 0.0;
	migShape[popnDest][popnSrc].rate_param = 0.0;
	migShape[popnDest][popnSrc].anchor_time = currentTime;
	
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
/* Add the dormant ancient samples of population popnDest at their
 * sampling time. If population 0 is sweeping then, each ancient lineage of
 * population 0 carries the beneficial allele (sweepPopn 1) with probability
 * currentFreq, the allele's frequency at that time, and is in the
 * wild-type class (sweepPopn 0) otherwise. */
void addAncientSample(int lineageNumber, int popnDest, double addTime, int stillSweeping, double currentFreq){
	compactNodes();
	int i;
	int count = 0;
	for(i=0; i < alleleNumber && count < lineageNumber; i++){
		if(nodes[i]->population == (popnDest+1) * -1){
			nodes[i]->population = popnDest;
			nodes[i]->time = addTime;
			// Add to destination population list
			addNodeToPopList(nodes[i], popnDest);
			if(stillSweeping == 1 && popnDest == 0){
				nodes[i]->sweepPopn = (ranf() < currentFreq) ? 1 : 0;
				/* An ancient carrier of a partial sweep's allele carries
				 * its selected SNP too. */
				if(partialSweepMode == 1 && nodes[i]->sweepPopn == 1 &&
				   partialSweepSite >= 0.0)
					addPartialSweepCarrier(nodes[i]);
			}
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
	ensureNodesCapacity(nodesUsed + 1);
	aNode->nodesIndex = nodesUsed;
	nodes[nodesUsed++] = aNode;
	
	/* The full ARG records a tskit node for every node. Minimal mode, as
	 * msprime, records none for recombination parents, and one for a
	 * coalescence parent only where its children's material overlaps: the
	 * merge (coalesceLineagesMinimal()) sets it before the parent is added,
	 * and a node from newRootedNode() has none. */
	if (!minimalTreeSeq) {
		tsk_id_t tsk_id = tskit_add_node(aNode->time, aNode->population, 0);  // is_sample=0
		set_tskit_node_id(aNode, tsk_id);
	}
	
	// Add to population list
	addNodeToPopList(aNode, aNode->population);
	
	alleleNumber += 1;
	popnSizes[aNode->population]+=1;
	if(aNode->population==0 && aNode->sweepPopn >= 0)
		sweepPopnSizes[aNode->sweepPopn]+=1;
}

/* Close the empty slots that removeNode() leaves in nodes[], keeping the
 * order of the lineages: after it, nodes[0 .. alleleNumber-1] are the
 * active lineages in the order in which they were added, as when each
 * removal shifted the later ones down. Code that walks nodes[] calls it
 * first; that order decides which lineage a random number picks, and
 * that dormant ancient samples of earlier events come first. */
static void compactNodes(void){
	if (nodesUsed == alleleNumber) return;
	int k = 0;
	for (int i = 0; i < nodesUsed; i++) {
		if (nodes[i] != NULL) {
			nodes[k] = nodes[i];
			nodes[k]->nodesIndex = k;
			k++;
		}
	}
	for (int i = k; i < nodesUsed; i++) nodes[i] = NULL;
	nodesUsed = k;
}

/* removeNode -- removes a lineage from nodes[] in O(1): its slot becomes
 * empty. The array is compacted when the empty slots are more than four
 * times the lineages plus 1024, so a compaction, O(slots), comes once in
 * many removals: the amortized cost stays O(1) and small. The slots hold
 * only pointers. */
void removeNode(rootedNode *aNode){
	popnSizes[aNode->population]-=1;
	if (aNode->population==0 && aNode->sweepPopn >= 0)
		sweepPopnSizes[aNode->sweepPopn]-=1;
	
	// Remove from population list
	removeNodeFromPopList(aNode);
	
	nodes[aNode->nodesIndex] = NULL;
	aNode->nodesIndex = -1;
	alleleNumber -= 1;
	if (nodesUsed - alleleNumber > 4 * alleleNumber + 1024) compactNodes();
	/* The event that removed it frees it (freeLineage()) once it has made
	 * the parents and recorded the edges. */
}

/*pickNodePopnSweep-- picks an allele at random from active nodes from specific popn and sweepPopn
*/
rootedNode *pickNodePopnSweep(int popn,int sp){
	compactNodes();
	int i, popnSize, indexArray[alleleNumber], index;

	//get popnSize -- may later decide to keep this is some sort of variable that is passed
	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
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

/*nodePopnSize-- returns popnSize of popn from
	nodes array */
int nodePopnSize(int popn){
	compactNodes();
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
	compactNodes();
	int i, popnSize;

	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popn && nodes[i]->sweepPopn == sp){
			popnSize++;
		}
	}
	return(popnSize);
}

/********** event stuff */
/* Stable insertion sort by time. Events at the same time keep the order
 * they were given in, as in ms: for example "-eg t ... -en t ..." ends with
 * no growth, because -en resets growth. qsort does not promise a stable
 * order (newer glibc and other C libraries reorder ties), and the demes
 * importer relies on its emission order for tied events. Event lists are
 * short enough that O(n^2) is cheap. */
void sortEventArray(struct event *eArray, int eNumber){
	for (int i = 1; i < eNumber; i++) {
		struct event key = eArray[i];
		int j = i - 1;
		while (j >= 0 && key.time < eArray[j].time) {
			eArray[j + 1] = eArray[j];
			j--;
		}
		eArray[j + 1] = key;
	}
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
/* Print the genotype matrix in ms format: one line per sample, one
 * character per site. tskit decodes the genotypes site by site, and ms
 * prints them sample by sample, so they are kept as one byte each
 * ('0' + genotype) in rows that end in a newline, and each row is written
 * with one fwrite(): one printf() per genotype took most of the time of
 * runs with many sites. The sites are decoded in blocks of SITE_BLOCK into
 * a small buffer, one line per site, and each sample's part of a block is
 * copied into its row in one piece: a write per genotype into rows far
 * apart missed the cache at every genotype.
 * A genotype outside 0..9 (more than one digit, or missing) prints as
 * printf("%d") does; the rows are then printed that way. Returns 0, or
 * -1 after an error message. */
#define SITE_BLOCK 64
static int printGenotypeRows(tsk_variant_t *variant, tsk_size_t num_samples,
                             tsk_size_t num_sites) {
	size_t rowLength = (size_t) num_sites + 1;
	char *rows = malloc((size_t) num_samples * rowLength);
	char *block = malloc((size_t) SITE_BLOCK * num_samples);
	int32_t *wide = NULL;  /* the genotypes as ints, if one is not a digit */
	int notDigit = 0;
	if (rows == NULL || block == NULL) {
		fprintf(stderr, "Error: Failed to allocate genotype matrix\n");
		free(rows);
		free(block);
		return -1;
	}
	for (tsk_size_t first = 0; first < num_sites && !notDigit; first += SITE_BLOCK) {
		tsk_size_t width = num_sites - first < SITE_BLOCK ? num_sites - first : SITE_BLOCK;
		for (tsk_size_t j = 0; j < width && !notDigit; j++) {
			int ret = tsk_variant_decode(variant, (tsk_id_t) (first + j), 0);
			if (ret != 0) {
				fprintf(stderr, "Error: Failed to decode variant at site %d: %s\n",
					(int) (first + j), tsk_strerror(ret));
				free(rows);
				free(block);
				return -1;
			}
			char *line = block + j * num_samples;
			for (tsk_size_t i = 0; i < num_samples; i++) {
				int32_t g = variant->genotypes[i];
				if (g < 0 || g > 9) {
					notDigit = 1;
					break;
				}
				line[i] = (char) ('0' + g);
			}
		}
		if (notDigit) break;
		for (tsk_size_t i = 0; i < num_samples; i++) {
			char *dest = rows + i * rowLength + first;
			for (tsk_size_t j = 0; j < width; j++)
				dest[j] = block[j * num_samples + i];
		}
	}
	free(block);
	if (notDigit) {
		/* Decode again, keeping every genotype as an int. */
		wide = malloc((size_t) num_samples * num_sites * sizeof(int32_t));
		if (wide == NULL) {
			fprintf(stderr, "Error: Failed to allocate genotype matrix\n");
			free(rows);
			return -1;
		}
		for (tsk_id_t site_id = 0; site_id < (tsk_id_t) num_sites; site_id++) {
			int ret = tsk_variant_decode(variant, site_id, 0);
			if (ret != 0) {
				fprintf(stderr, "Error: Failed to decode variant at site %d: %s\n",
					site_id, tsk_strerror(ret));
				free(rows);
				free(wide);
				return -1;
			}
			for (tsk_size_t i = 0; i < num_samples; i++)
				wide[i * num_sites + site_id] = variant->genotypes[i];
		}
	}
	for (tsk_size_t i = 0; i < num_samples; i++) {
		if (wide == NULL) {
			rows[i * rowLength + num_sites] = '\n';
			fwrite(rows + i * rowLength, 1, rowLength, stdout);
		} else {
			for (tsk_size_t j = 0; j < num_sites; j++)
				printf("%d", wide[i * num_sites + j]);
			printf("\n");
		}
	}
	free(rows);
	free(wide);
	return 0;
}

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
		if (printGenotypeRows(&variant, num_samples, num_sites) != 0) {
			tsk_variant_free(&variant);
			tsk_treeseq_free(&ts);
			return;
		}
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

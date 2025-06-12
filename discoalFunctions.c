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
#include "ancestryVerify.h"
#include "ancestryWrapper.h"
#include "activeSegment.h"
#include <time.h>
#include "discoal.h"
#include "discoalFunctions.h"
#include "ranlib.h"
#include "alleleTraj.h"


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
		// Optional: Print growth information for debugging
		// fprintf(stderr, "Debug: Grew breakPoints capacity to %d\n", breakPointsCapacity);
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

void initialize(){
	int i,j,p, count=0;
	int leafID=0;
	int tmpCount = 0;
	
	
	/* Initialize the arrays */
	totChunkNumber = 0;
	initializeBreakPoints();
	initializeNodeArrays();
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			// Initialize ancestry segment tree for leaf node
			nodes[count]->ancestryRoot = newSegment(0, nSites, NULL, NULL);

			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;
			allNodes[count] = nodes[count];
			count += 1;
		}
	}


	breakNumber = 0;
	// Initialize segment-based active material (all sites start as active)
	initializeActiveMaterial(&activeMaterialSegments, nSites);
	
	//set initial counts
	alleleNumber = sampleSize;
	totNodeNumber = sampleSize;
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
	//mask mode?
	if (mask){
//		getMasking(mFile);
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
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			// Initialize ancestry segment tree for leaf node
			nodes[count]->ancestryRoot = newSegment(0, nSites, NULL, NULL);
			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;
			//do stuff for leafs containers
			//nodes->leafs = calloc(sizeof(int) * sampleSize);
			//nodes->leafs[nodes[count]->id] = 1;
			allNodes[count] = nodes[count];
			count += 1;
		}
	}
	
	breakNumber = 0;
	// Initialize segment-based active material (all sites start as active)
	initializeActiveMaterial(&activeMaterialSegments, nSites);
	
	//set initial counts
	alleleNumber = sampleSize;
	totNodeNumber = sampleSize;
	activeSites = nSites;
	if (npops>1){
		if(tDiv==666 && mig[0] == 0 && mig[1] == 0){
			fprintf(stderr,"tDiv not set in population split model\n");
			exit(1);
		}
		
		eventFlag = 0;
	}
	//mask mode?
	if (mask){
//		getMasking(mFile);
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


rootedNode *newRootedNode(double cTime, int popn) {
	rootedNode *temp;

	temp =  malloc(sizeof(rootedNode));
	if(temp==NULL){
		printf("fucked: out o' memory from rootedNode alloc\n");
	}
	temp->rightParent = NULL;
	temp->leftParent = NULL;
	temp->rightChild = NULL;
	temp->leftChild = NULL; 
	temp->time = cTime;
	temp->branchLength=0.0;
	temp->mutationNumber = 0;
	temp->population = popn;
	temp->sweepPopn = -1;
	
	// Initialize muts array dynamically
	initializeMuts(temp, 10);  // Start small, grow as needed
//	temp->leafs = malloc(sizeof(int) * sampleSize);
//	for(i=0;i<nSites;i++)temp->ancSites[i]=1;

	// Initialize ancestry segment tree to NULL (will be set during initialization)
	temp->ancestryRoot = NULL;

	return temp;
}

////////  Migration Stuff ///////////

//migrateAtTime -- this only switches between two populations
void migrateAtTime(double cTime,int srcPopn, int destPopn){
	rootedNode *temp;
	
	temp = pickNodePopn(srcPopn);
	temp->population = destPopn;

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
	rootedNode *temp;
	int doMig=0;
	temp = pickNodePopn(srcPopn);
	
	if(isAncestralHere(temp,site) == 0 || ranf() < scalar) doMig = 1; //migration of site is 10% rate of other sites
	if(doMig == 1)
	{
		temp->population = destPopn;
		popnSizes[srcPopn]-=1;
		popnSizes[destPopn]+=1;
		if (srcPopn==0)
			sweepPopnSizes[temp->sweepPopn]-=1;
	}
	
	
}


////////

void coalesceAtTimePopn(double cTime, int popn){
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
	// Merge ancestry segment trees
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot);
	
	// Update stats from tree when in tree-only mode
	updateAncestryStatsFromTree(temp);
	
	// Verify consistency in debug mode
	#ifdef DEBUG_ANCESTRY
	if (!verifyAncestryConsistency(temp, nSites)) {
		fprintf(stderr, "Ancestry verification failed after coalescence\n");
		exit(1);
	}
	#endif
	
	//printNode(temp);
	//printf("coal-> lChild: %u rChild: %u parent: %u time: %f\n",lChild,rChild,temp,cTime);
	addNode(temp);
	//update active anc. material
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
	
	//printNode(temp);
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
		lParent->ancestryRoot = splitLeft(aNode->ancestryRoot, xOver);
		rParent->ancestryRoot = splitRight(aNode->ancestryRoot, xOver);
		
		// Update stats from tree when in tree-only mode
		updateAncestryStatsFromTree(lParent);
		updateAncestryStatsFromTree(rParent);
		
		addNode(lParent);
		addNode(rParent);
		
	//	printf("reco-> lParent: %u rParent:%u child: %u time:%f xover:%d\n",lParent,rParent,aNode,cTime,xOver);
		//printNode(aNode);
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
		} else {
			lParent->ancestryRoot = NULL;
			rParent->ancestryRoot = NULL;
		}
		
		// Update stats from the ancestry trees
		updateAncestryStatsFromTree(lParent);
		updateAncestryStatsFromTree(rParent);
		
		addNode(lParent);
		addNode(rParent);
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
					
					
					}
				}
			}
		}
	}
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
						if(isAncestralHere(nodes[i],sweepSite) && hidePartialSNP == 0)
							addMutation(nodes[i],sweepSite);
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
		if(ranf() < sweepPopTotRate / totRate){
			//event is in sweep pop; choose event
			r = ranf();
			sum = pCoalB;
			//coalescent in B?
			if(r < sum / sweepPopTotRate){
				coalesceAtTimePopnSweep(cTime+(ttau), 0,1);
			}
			else{
				sum += pCoalb;
				//coalescent in b?
				if(r < sum / sweepPopTotRate){
					coalesceAtTimePopnSweep(cTime+(ttau),0, 0);
				}
				else{
					sum += pRecb;
					//recombination in b, also need bookkeeping
					if(r < sum / sweepPopTotRate){
						bp = recombineAtTimePopnSweep(cTime + (ttau),0, 0, sweepSite, (1.0-x));
						if(bp != 666){
							addBreakPoint(bp);
							if(bp >= lSpot && bp < rSpot){
								condRecMet = 1;
							}
						}
					}
					else{
						sum += pRecB;
						if( r < sum / sweepPopTotRate){
							//recombination in B and bookkeeping
							bp = recombineAtTimePopnSweep(cTime + (ttau),0, 1, sweepSite, x);
							if(bp != 666){
								addBreakPoint(bp);
								if(bp >= lSpot && bp < rSpot){
									condRecMet = 1;
								}
							}
						}
						else{
							sum+= pGCB;
							if( r < sum / sweepPopTotRate){
								geneConversionAtTimePopnSweep(cTime + (ttau),0, 1, sweepSite, x);
							}
							else{
								sum+= pGCb;
								if( r < sum / sweepPopTotRate){
									geneConversionAtTimePopnSweep(cTime + (ttau),0, 0, sweepSite, x);
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
											//fprintf(stderr,"recurrent mutation at time %f; freq=%f\n", cTime+(ttau),x);
											//fprintf(stderr,"recurMut prob: %g; uA: %f, sweepPopnSizes[1]:%d; x:%f; tInc: %g; sizeRatio: %f\npopnSizes[0]: %d; popnSizes[1]: %d\n", pRecurMut,uA,sweepPopnSizes[1],x,tInc,sizeRatio,popnSizes[0],popnSizes[1]);
											recurrentMutAtTime(cTime+(ttau),0, 1);
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
				bp = recombineAtTimePopn(cTime,i);
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
					geneConversionAtTimePopn(cTime,i);
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
					coalesceAtTimePopn(cTime,i);
					
					
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
						if(isAncestralHere(nodes[i],sweepSite) && hidePartialSNP == 0)
							addMutation(nodes[i],sweepSite);
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
		if(ranf() < sweepPopTotRate / totRate){
			//event is in sweep pop; choose event
			r = ranf();
			sum = pCoalB;
			//coalescent in B?
			if(r < sum / sweepPopTotRate){
				coalesceAtTimePopnSweep(cTime+(ttau), 0,1);
			}
			else{
				sum += pCoalb;
				//coalescent in b?
				if(r < sum / sweepPopTotRate){
					coalesceAtTimePopnSweep(cTime+(ttau),0, 0);
				}
				else{
					sum += pRecb;
					//recombination in b, also need bookkeeping
					if(r < sum / sweepPopTotRate){
						bp = recombineAtTimePopnSweep(cTime + (ttau),0, 0, sweepSite, (1.0-x));
						if(bp != 666){
							addBreakPoint(bp);
							if(bp >= lSpot && bp < rSpot){
								condRecMet = 1;
							}
						}
					}
					else{
						sum += pRecB;
						if( r < sum / sweepPopTotRate){
							//recombination in B and bookkeeping
							bp = recombineAtTimePopnSweep(cTime + (ttau),0, 1, sweepSite, x);
							if(bp != 666){
								addBreakPoint(bp);
								if(bp >= lSpot && bp < rSpot){
									condRecMet = 1;
								}
							}
						}
						else{
							sum+= pGCB;
							if( r < sum / sweepPopTotRate){
								geneConversionAtTimePopnSweep(cTime + (ttau),0, 1, sweepSite, x);
							}
							else{
								sum+= pGCb;
								if( r < sum / sweepPopTotRate){
									geneConversionAtTimePopnSweep(cTime + (ttau),0, 0, sweepSite, x);
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
											//fprintf(stderr,"recurrent mutation at time %f; freq=%f\n", cTime+(ttau),x);
											//fprintf(stderr,"recurMut prob: %g; uA: %f, sweepPopnSizes[1]:%d; x:%f; tInc: %g; sizeRatio: %f\npopnSizes[0]: %d; popnSizes[1]: %d\n", pRecurMut,uA,sweepPopnSizes[1],x,tInc,sizeRatio,popnSizes[0],popnSizes[1]);
											recurrentMutAtTime(cTime+(ttau),0, 1);
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
						bp = recombineAtTimePopn(cTime,i);
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
							geneConversionAtTimePopn(cTime,i);
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
							coalesceAtTimePopn(cTime,i);
					
					
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
			
			// Split ancestry segment tree at crossover point
			lParent->ancestryRoot = splitLeft(aNode->ancestryRoot, xOver);
			rParent->ancestryRoot = splitRight(aNode->ancestryRoot, xOver);
			
			// Update stats from ancestry trees
			updateAncestryStatsFromTree(lParent);
			updateAncestryStatsFromTree(rParent);
			
			//add in the nodes
			addNode(lParent);
			addNode(rParent);
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
		if(sweepSite >= (float) xOver / nSites  || sweepSite < (float) (xOver+tractL) / nSites ){ //lParent has sweep site set to same as child node,rParent is then random
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
			//add in the nodes
		addNode(lParent);
		addNode(rParent);
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
	// Merge ancestry segment trees
	temp->ancestryRoot = mergeAncestryTrees(lChild->ancestryRoot, rChild->ancestryRoot);
	// Update stats from tree when in tree-only mode
	updateAncestryStatsFromTree(temp);
	//printNode(temp);
	addNode(temp);
	//update active anc. material
	updateActiveMaterial(temp);
	//popnSizes[popn]--; //decrese popnSize
	//sweepPopnSizes[sp]--;	
}



/*******************************************************/
void dropMutations(){
	int i, j, m;
	double p;
	double mutSite, error;
	
	//get time and set probs
	coaltime = totalTimeInTree();
	//printf("%f\n",coaltime);
	for(i=0;i<totNodeNumber;i++){
		//add Mutations
		p = allNodes[i]->blProb * theta * 0.5;
		if(p>0.0){
			m = ignpoi(p);
			while(m>0){
				mutSite = genunf((float)allNodes[i]->lLim / nSites, (float) allNodes[i]->rLim / nSites);
				while(isAncestralHere(allNodes[i],mutSite) != 1){
				//	printf("in here\n");
					if (allNodes[i]->lLim == allNodes[i]->rLim){
						if (mutSite*nSites < (float)allNodes[i]->lLim)
							error = (float)allNodes[i]->lLim-(mutSite*nSites);
						else
							error = 0.0;
						p = allNodes[i]->rLim + (1.0/nSites);
						mutSite = genunf(((double)allNodes[i]->lLim + error) / nSites, (double)(p + error) / nSites);
					}
					else
						mutSite = genunf((double)allNodes[i]->lLim / nSites, (double) allNodes[i]->rLim / nSites);
				}
				//printf("new mut allele: %d mut: %f\n",i,mutSite);
				addMutation(allNodes[i],mutSite);
				m--;
			}
		}
	}
//	printf("coaltime: %f totM: %d\n",coaltime,tm);
	//push mutations down tree
	for(i=totNodeNumber-1;i>=0;i--){
		for(j=0;j<allNodes[i]->mutationNumber;j++){
			mutSite = allNodes[i]->muts[j];
			if(allNodes[i]->leftChild != NULL && isAncestralHere(allNodes[i]->leftChild,mutSite))
				addMutation(allNodes[i]->leftChild,mutSite);
			if(allNodes[i]->rightChild != NULL && isAncestralHere(allNodes[i]->rightChild,mutSite))
				addMutation(allNodes[i]->rightChild,mutSite);
				
		}
	}
}

void dropMutationsRecurse(){
	int i, j, m;
	double p;
	float mutSite, error;
	
	//get time and set probs
	coaltime = totalTimeInTree();
	//printf("%f\n",coaltime);
	for(i=0;i<totNodeNumber;i++){
		//add Mutations
		p = allNodes[i]->blProb * theta * 0.5;
		if(p>0.0){
			m = ignpoi(p);
			while(m>0){
				mutSite = genunf((float)allNodes[i]->lLim / nSites, (float) allNodes[i]->rLim / nSites);
				while(isAncestralHere(allNodes[i],mutSite) != 1){
				//	printf("in here\n");
					if (allNodes[i]->lLim == allNodes[i]->rLim){
						if (mutSite*nSites < (float)allNodes[i]->lLim)
							error = (float)allNodes[i]->lLim-(mutSite*nSites);
						else
							error = 0.0;
						p = allNodes[i]->rLim + (1.0/nSites);
						mutSite = genunf(((double)allNodes[i]->lLim + error) / nSites, (double)(p + error) / nSites);
					}
					else
						mutSite = genunf((double)allNodes[i]->lLim / nSites, (double) allNodes[i]->rLim / nSites);
				}
		//	printf("mut: %f\n",mutSite);
				addMutation(allNodes[i],mutSite);
				m--;
			}
		}
	}
//	printf("coaltime: %f totM: %d\n",coaltime,tm);
	//push mutations down tree
	for(i=totNodeNumber-1;i>=0;i--){
		for(j=0;j<allNodes[i]->mutationNumber;j++){
			mutSite = allNodes[i]->muts[j];
			recurseTreePushMutation(allNodes[i],mutSite);
		}
	}
}
//calculates the total time in the tree, then set blProbs for each node
double totalTimeInTree(){
	double tTime, siteLength;
	int i;
	

	tTime=0.0;
	siteLength=1.0/nSites;
	for(i=0;i<totNodeNumber;i++){
		allNodes[i]->blProb = siteLength * allNodes[i]->nancSites * allNodes[i]->branchLength;
		tTime += allNodes[i]->blProb;
	}
//	for(i=0;i<totNodeNumber;i++)allNodes[i]->blProb /= tTime;
	return tTime;
}

void recurseTreePushMutation(rootedNode *aNode, float site){
	if(aNode->leftChild != NULL && isAncestralHere(aNode->leftChild,site))
		recurseTreePushMutation(aNode->leftChild,site);
	if(aNode->rightChild != NULL && isAncestralHere(aNode->rightChild,site))
		recurseTreePushMutation(aNode->rightChild,site);
	if( isLeaf(aNode) && isAncestralHere(aNode,site)){
		if(hasMutation(aNode,site)==0) addMutation(aNode, site);
	}
}
void addMutation(rootedNode *aNode, double site){
	ensureMutsCapacity(aNode, aNode->mutationNumber + 1);
	aNode->muts[aNode->mutationNumber] = site;
	aNode->mutationNumber += 1;

}

/* Sort mutations on a single node */
void sortNodeMutations(rootedNode *node) {
	if (node != NULL && node->mutationNumber > 1) {
		qsort(node->muts, node->mutationNumber, sizeof(double), compare_doubles);
	}
}

/* Sort mutations on all nodes after mutation placement */
void sortAllMutations() {
	int i;
	for (i = 0; i < totNodeNumber; i++) {
		if (allNodes[i] != NULL) {
			sortNodeMutations(allNodes[i]);
		}
	}
}

/* Binary search implementation for sorted mutations */
static int hasMutationBinary(rootedNode *aNode, double site) {
	int left = 0;
	int right = aNode->mutationNumber - 1;
	
	while (left <= right) {
		int mid = left + (right - left) / 2;
		if (aNode->muts[mid] == site) {
			return 1;
		} else if (aNode->muts[mid] < site) {
			left = mid + 1;
		} else {
			right = mid - 1;
		}
	}
	return 0;
}

int hasMutation(rootedNode *aNode, double site){
	if (aNode->mutationNumber == 0) return 0;
	
	/* Use linear search for small arrays (faster due to cache locality) */
	if (aNode->mutationNumber < 10) {
		int i;
		for (i = 0; i < aNode->mutationNumber; i++){
			if (aNode->muts[i] == site) return 1;
		}
		return 0;
	}
	
	/* Use binary search for larger arrays */
	return hasMutationBinary(aNode, site);
}

//findRootAtSite-- returns the index of the node that is the root at a given site
int findRootAtSite(float site){
	int j;
	j = 0;
	while(nAncestorsHere(allNodes[j], site) != sampleSize){
		j++;
	}
	return(j);
}

//printTreeAtSite-- prints a newick tree at a given site
void printTreeAtSite(float site){
	int rootIdx;
    float tPtr;
	
    tPtr = 0;
	rootIdx = findRootAtSite(site);
	newickRecurse(allNodes[rootIdx],site,tPtr);
	printf(";\n");

}

void newickRecurse(rootedNode *aNode, float site, float tempTime){
	//printf("site: %f tempTime: %f\n",site,tempTime);
	//printNode(aNode);
    if(isCoalNode(aNode)){
		
		if(hasMaterialHere(aNode->leftChild,site) && hasMaterialHere(aNode->rightChild,site)){	
			printf("(");
			newickRecurse(aNode->leftChild,site,0.0);
			printf(",");
			newickRecurse(aNode->rightChild,site,0.0);
			printf(")");
			if(nAncestorsHere(aNode, site) != sampleSize){
				printf(":%f",(aNode->branchLength + tempTime)*0.5);
			}
			
		}
		else{
            if(hasMaterialHere(aNode->leftChild,site)){
                tempTime += aNode->branchLength;
                newickRecurse(aNode->leftChild,site, \
                        tempTime);
            }
			else if(hasMaterialHere(aNode->rightChild,site)){
                tempTime += aNode->branchLength;
                newickRecurse(aNode->rightChild,site, \
                        tempTime);
            }
		}
	}
	else{
		if(isLeaf(aNode)){
			printf("%d:%f",aNode->id, \
                    (aNode->branchLength +tempTime)*0.5);
		}
		else{ //recombination node
			if(hasMaterialHere(aNode->leftChild,site) && \
                    hasMaterialHere(aNode,site)){
               tempTime+= aNode->branchLength; 
               newickRecurse(aNode->leftChild,site, \
                        tempTime);
            }
		}
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

/*makeGametesMS-- MS style sample output */
void makeGametesMS(int argc,const char *argv[]){
	int i,j, size, mutNumber;
	double allMuts[MAXMUTS];
	MutHashEntry *hashTable[MUTATION_HASH_SIZE];
	MutHashEntry *entry, *newEntry;

	/* Sort all mutations before output generation for binary search */
	sortAllMutations();

	/* Initialize hash table */
	for (i = 0; i < MUTATION_HASH_SIZE; i++) {
		hashTable[i] = NULL;
	}

	/* get unique list of muts using hash table for O(n) duplicate detection */
	size = 0;
	for (i = 0; i < sampleSize; i++){
		for (j = 0; j < allNodes[i]->mutationNumber; j++){
			double mut = allNodes[i]->muts[j];
			unsigned int hash = hashMutation(mut);
			
			/* Check if mutation already exists in hash table */
			int found = 0;
			for (entry = hashTable[hash]; entry != NULL; entry = entry->next) {
				if (entry->mutation == mut) {
					found = 1;
					break;
				}
			}
			
			if (!found) {
				assert(size < MAXMUTS);
				allMuts[size] = mut;
				size++;
				
				/* Add to hash table */
				newEntry = (MutHashEntry*)malloc(sizeof(MutHashEntry));
				newEntry->mutation = mut;
				newEntry->next = hashTable[hash];
				hashTable[hash] = newEntry;
			}
		}
	}
	
	/* Clean up hash table */
	for (i = 0; i < MUTATION_HASH_SIZE; i++) {
		entry = hashTable[i];
		while (entry != NULL) {
			MutHashEntry *temp = entry;
			entry = entry->next;
			free(temp);
		}
	}
	
	mutNumber = size;
	qsort(allMuts, size, sizeof(allMuts[0]), compare_doubles);
	printf("\n//\nsegsites: %d",mutNumber);
	if(mutNumber > 0) printf("\npositions: ");
	for(i = 0; i < mutNumber; i++)
		fprintf(stdout,"%6.6lf ",allMuts[i] );
	fprintf(stdout,"\n");

/* Phase 3 optimization: Pre-compute presence matrix */
	char *presenceMatrix = NULL;
	if (mutNumber > 0) {
		presenceMatrix = (char*)calloc(sampleSize * mutNumber, sizeof(char));
		if (presenceMatrix == NULL) {
			fprintf(stderr, "Error: Failed to allocate presence matrix\n");
			exit(1);
		}
		
		/* Pre-compute all ancestry and mutation information */
		for (i = 0; i < sampleSize; i++) {
			for (j = 0; j < mutNumber; j++) {
				int idx = i * mutNumber + j;
				if (isAncestralHere(allNodes[i], allMuts[j])) {
					if (hasMutation(allNodes[i], allMuts[j])) {
						presenceMatrix[idx] = '1';
					} else {
						presenceMatrix[idx] = '0';
					}
				} else {
					presenceMatrix[idx] = 'N';
				}
			}
		}
	}
	
	/* Output using pre-computed matrix */
	for (i = 0; i < sampleSize; i++) {
		for (j = 0; j < mutNumber; j++) {
			printf("%c", presenceMatrix[i * mutNumber + j]);
		}
		printf("\n");
	}
	
	/* Clean up */
	if (presenceMatrix) {
		free(presenceMatrix);
	}  
}

void errorCheckMutations(){
	int i,j;
	
	for (i = 0; i < sampleSize; i++){
		printf("allNodes[%d]:\n",i);
		for (j = 0; j < allNodes[i]->mutationNumber; j++){
			printf("muts[%d]=%lf\n",j,allNodes[i]->muts[j]);
		}
	}
}

//dropMutationsUntilTime-- places mutationson the tree iff they occur before time t
//this will not return the "correct" number of mutations conditional on theta
void dropMutationsUntilTime(double t){
	int i, j, m;
	double mutSite,p;
	//get time and set probs
	coaltime = totalTimeInTreeUntilTime(t);
	//printf("%f\n",coaltime);
	for(i=0;i<totNodeNumber;i++){
		//add Mutations
		p = allNodes[i]->blProb * theta * 0.5;
		if(p>0.0){
		  m = ignpoi(p);
		  while(m>0){
		  mutSite = genunf((float)allNodes[i]->lLim / nSites, (float) allNodes[i]->rLim / nSites);
		  while(isAncestralHere(allNodes[i],mutSite) != 1){
				  if (allNodes[i]->lLim == allNodes[i]->rLim){
				    p = allNodes[i]->rLim + (1.0/nSites);
				    mutSite = genunf((float)allNodes[i]->lLim / nSites, p / nSites);
				  }
				  else
				    mutSite = genunf((float)allNodes[i]->lLim / nSites, (float) allNodes[i]->rLim / nSites);
			  }
		  //	printf("mut: %f\n",mutSite);
			addMutation(allNodes[i],mutSite);
			 m--;
		  }
		}
	}
//	printf("coaltime: %f totM: %d\n",coaltime,tm);
	//push mutations down tree
	for(i=totNodeNumber-1;i>=0;i--){
		for(j=0;j<allNodes[i]->mutationNumber;j++){
			mutSite = allNodes[i]->muts[j];
			if(allNodes[i]->leftChild != NULL && isAncestralHere(allNodes[i]->leftChild,mutSite))
				addMutation(allNodes[i]->leftChild,mutSite);
			if(allNodes[i]->rightChild != NULL && isAncestralHere(allNodes[i]->rightChild,mutSite))
				addMutation(allNodes[i]->rightChild,mutSite);
				
		}
	}
}

//calculates the total time in the tree, then set blProbs for each node
double totalTimeInTreeUntilTime(double t){
	double tTime, siteLength;
	int i;
	

	tTime=0.0;
	siteLength=1.0/nSites;
	for(i=0;i<totNodeNumber;i++){
		//censor times
		if(allNodes[i]->leftParent == NULL || allNodes[i]->leftParent->time > t){
			allNodes[i]->branchLength = MAX(t - allNodes[i]->time,0);
		}
		allNodes[i]->blProb = siteLength * allNodes[i]->nancSites * allNodes[i]->branchLength;
		tTime += allNodes[i]->blProb;
	}
//	for(i=0;i<totNodeNumber;i++)allNodes[i]->blProb /= tTime;
	return tTime;
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
	int i, popnSize, indexArray[alleleNumber], index;

	//get popnSize -- may later decide to keep this is some sort of variable that is passed
	popnSize = 0;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popn){
			indexArray[popnSize] = i;  //store indexes as we go
			popnSize++;
		}
	}
	if(popnSize == 0){
		fprintf(stderr,"error encountered in pickNodePopn\n");
		fprintf(stderr,"tried to pick allele from popn %d, but popnSize is %d! Rho=%f\n",popn,popnSize,rho);
		exit(1);
	}
	//choose random index
	index = indexArray[ignuin(0, popnSize - 1)];

	return(nodes[index]);
}

void mergePopns(int popnSrc, int popnDest){
	int i;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popnSrc){
			nodes[i]->population = popnDest;
			popnSizes[popnDest]++;
			popnSizes[popnSrc]--;
		//	sweepPopnSizes[popnDest]++;
		//	sweepPopnSizes[popnSrc]--;
			
		}
	}
	//set migration rates to zero
	migMat[popnSrc][popnDest] = 0.0;
	migMat[popnDest][popnSrc] = 0.0;
	
}

void admixPopns(int popnSrc, int popnDest1, int popnDest2, double admixProp){
	int i;
	double rn;
	for(i = 0 ; i < alleleNumber; i++){
		if(nodes[i]->population == popnSrc){
			rn = ranf();
			if(rn < admixProp){
				nodes[i]->population = popnDest1;
				popnSizes[popnDest1]++;
				popnSizes[popnSrc]--;
			}
			else{
				nodes[i]->population = popnDest2;
				popnSizes[popnDest2]++;
				popnSizes[popnSrc]--;
			}
		}
	}
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
	int initialCapacity = 1000;  // Start small
	
	nodesCapacity = initialCapacity;
	allNodesCapacity = initialCapacity;
	
	nodes = malloc(sizeof(rootedNode*) * nodesCapacity);
	allNodes = malloc(sizeof(rootedNode*) * allNodesCapacity);
	
	if (nodes == NULL || allNodes == NULL) {
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
		
		nodes = newNodes;
		nodesCapacity = newCapacity;
	}
}

void ensureAllNodesCapacity(int requiredSize) {
	if (requiredSize >= allNodesCapacity) {
		int newCapacity = allNodesCapacity;
		while (newCapacity <= requiredSize) {
			newCapacity *= 2;
		}
		
		rootedNode **newAllNodes = realloc(allNodes, sizeof(rootedNode*) * newCapacity);
		if (newAllNodes == NULL) {
			fprintf(stderr, "Error: Failed to reallocate allNodes array (requested: %d pointers)\n", newCapacity);
			exit(1);
		}
		
		allNodes = newAllNodes;
		allNodesCapacity = newCapacity;
	}
}

void addNode(rootedNode *aNode){
	ensureNodesCapacity(alleleNumber + 1);
	ensureAllNodesCapacity(totNodeNumber + 1);
	
	nodes[alleleNumber] = aNode;
	allNodes[totNodeNumber] = aNode;
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

	for (i = index ; i < alleleNumber; i++){
		temp =  nodes[i + 1];
		nodes[i] = temp;
	}
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
	removeNodeAt(i);
}

/* addNodeAtIndex - variation on the theme here. adds
a rootedNode at a specific index. useful for adding pops
	or outgroups to current pops */
void addNodeAtIndex(rootedNode *aNode, int anIndex){
	nodes[anIndex] = aNode;
	allNodes[anIndex] = aNode;
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

	for(i = totNodeNumber - 1; i >= 0; i--){
		temp =  allNodes[i];
		allNodes[i+offset] = temp;
	}
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
		cleanupMuts(allNodes[i]);      // Free muts array
		// Free ancestry segment tree
		if (allNodes[i]->ancestryRoot) {
			freeSegmentTree(allNodes[i]->ancestryRoot);
			allNodes[i]->ancestryRoot = NULL;
		}
		free(allNodes[i]);
		allNodes[i] = NULL;
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

void printAllNodes(){
	int i;
	for(i = 0 ; i < totNodeNumber; i++){
		printNode(allNodes[i]);
	}
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
void initializeMuts(rootedNode *node, int capacity) {
	if (capacity <= 0) {
		capacity = 10;  // Start with small capacity, will grow as needed
	}
	
	node->mutsCapacity = capacity;
	node->muts = malloc(sizeof(double) * capacity);
	
	if (node->muts == NULL) {
		fprintf(stderr, "Error: Failed to allocate memory for muts array (capacity: %d)\n", capacity);
		exit(1);
	}
}

void ensureMutsCapacity(rootedNode *node, int requiredSize) {
	if (requiredSize > node->mutsCapacity) {
		int newCapacity = node->mutsCapacity;
		
		// Double capacity until we have enough
		while (newCapacity < requiredSize) {
			newCapacity *= 2;
		}
		
		double *newMuts = realloc(node->muts, sizeof(double) * newCapacity);
		if (newMuts == NULL) {
			fprintf(stderr, "Error: Failed to reallocate memory for muts array (capacity: %d -> %d)\n", 
					node->mutsCapacity, newCapacity);
			exit(1);
		}
		
		node->muts = newMuts;
		node->mutsCapacity = newCapacity;
	}
}

void cleanupMuts(rootedNode *node) {
	if (node->muts != NULL) {
		free(node->muts);
		node->muts = NULL;
		node->mutsCapacity = 0;
	}
}

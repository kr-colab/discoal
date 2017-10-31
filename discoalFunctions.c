
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <assert.h>
#include "discoal.h"
#include "discoalFunctions.h"
#include "ranlib.h"
#include "alleleTraj.h"



void initialize(){
	int i,j,p, count=0;
	int leafID=0;
	int tmpCount = 0;
	
	
	/* Initialize the arrays */
	totChunkNumber = 0;
	breakPoints[0] = 666;
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			for(j=0;j<nSites;j++)nodes[count]->ancSites[j]=1;

			if(p>0)nodes[count]->sweepPopn = 0;
			nodes[count]->id=leafID++;
			allNodes[count] = nodes[count];
			count += 1;
		}
	}


	breakNumber = 0;
	//setup active material array
	for(i=0;i<nSites;i++) activeMaterial[i] = 1;
	
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
	
	maxTrajSteps = TRAJSTEPSTART;
	currentTrajectory = malloc(sizeof(float) * TRAJSTEPSTART);
	
}

void initializeTwoSite(){
	int i,j,p, count=0;
	int leafID=0;
	/* initialize the arrays */
	totChunkNumber = 0;
	breakPoints[0] = 666;
	for(p=0;p<npops;p++){
		popnSizes[p]=sampleSizes[p];
		for( i = 0; i < sampleSizes[p]; i++){
			nodes[count] = newRootedNode(0,p);
			nodes[count]->nancSites = nSites;
			nodes[count]->lLim=0;
			nodes[count]->rLim=nSites-1;
			for(j=0;j<nSites;j++)nodes[count]->ancSites[j]=1;
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
	//setup active material array
	for(i=0;i<nSites;i++) activeMaterial[i] = 1;
	
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
//	temp->sweepPopn =1;
//	temp->leafs = malloc(sizeof(int) * sampleSize);
//	for(i=0;i<nSites;i++)temp->ancSites[i]=1;

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
	for(i=0;i<nSites;i++){
		temp->ancSites[i] = lChild->ancSites[i] + rChild->ancSites[i];
		if(temp->ancSites[i] > 0 && temp->ancSites[i] < sampleSize){
			temp->nancSites += 1;
			temp->rLim=i;
			temp->lLim=MIN(temp->lLim,i);
			
		}
	}
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
	for(i=0;i<nSites;i++){
		temp->ancSites[i] = lChild->ancSites[i] + rChild->ancSites[i];
		if(temp->ancSites[i] > 0 && temp->ancSites[i] < sampleSize){
			temp->nancSites += 1;
			temp->rLim=i;
			temp->lLim=MIN(temp->lLim,i);
			
		}
	}
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


/*updateActiveMaterial- does bookkeeping on the material that has found it's MRCA-- added for efficiency */
void updateActiveMaterial(rootedNode *aNode){
	int i;
	activeSites = 0;
	for(i=0;i<nSites;i++){
		if(aNode->ancSites[i]==sampleSize){
			activeMaterial[i]=0;
		}
		//set global activeSites
		activeSites+=activeMaterial[i];	
	}
}

/* asks whether a potential breakpoint is in material that hasn't found MRCA yet*/
int isActive(int site){
	return(activeMaterial[site]);
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
	if(aNode->ancSites[bp] > 0 && aNode->ancSites[bp] < sampleSize)
		return(1);
	else
		return(0);
}

int hasMaterialHere(rootedNode *aNode, float site){
	int bp;
	bp = floor(site * nSites);
	if(aNode->ancSites[bp] > 0)
		return(1);
	else
		return(0);
}

int nAncestorsHere(rootedNode *aNode, float site){
	int bp;
	bp = floor(site * nSites);
	return(aNode->ancSites[bp]);
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
//	printf("picked: %p\n",aNode);
	xOver = ignuin(0,nSites-1);
//	printf("xo: %d test1: %d test2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver) );
//	printNode(aNode);
//	for(i=0;i<nSites;i++)printf("%d",activeMaterial[i]);
//	printf("\n");
//	printf("there xover: %d t1: %d t2: %d\n",xOver,siteBetweenChunks(aNode, xOver),isActive(xOver));
		
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
		for(i=0;i<nSites;i++){
			if(i<xOver){
				lParent->ancSites[i] = aNode->ancSites[i];
				rParent->ancSites[i] = 0;
			}
			else{
				lParent->ancSites[i] = 0;
				rParent->ancSites[i] = aNode->ancSites[i];
			}

			if(lParent->ancSites[i] > 0 && lParent->ancSites[i] < sampleSize){
				lParent->nancSites += 1;
				lParent->rLim=i;
				lParent->lLim=MIN(lParent->lLim,i);
			}

			if(rParent->ancSites[i] > 0 && rParent->ancSites[i] < sampleSize){
				rParent->nancSites += 1;
				rParent->rLim=i;
				rParent->lLim=MIN(rParent->lLim,i);
			}
		}
	//	if(lParent->lLim < lParent->rLim)
			addNode(lParent);
	//	else
	//		free(lParent);
	//	if(rParent->lLim < rParent->rLim)
			addNode(rParent);
	//	else
	//		free(rParent)
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
//	printf("picked: %p\n",aNode);
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
		for(i=0;i<nSites;i++){
			if(i<xOver || i >= xOver+tractL){
				lParent->ancSites[i] = 0;
				rParent->ancSites[i] = aNode->ancSites[i];
			}
			else{
				lParent->ancSites[i] = aNode->ancSites[i];
				rParent->ancSites[i] = 0;
			}

			if(lParent->ancSites[i] > 0 && lParent->ancSites[i] < sampleSize){
				lParent->nancSites += 1;
				lParent->rLim=i;
				lParent->lLim=MIN(lParent->lLim,i);
			}

			if(rParent->ancSites[i] > 0 && rParent->ancSites[i] < sampleSize){
				rParent->nancSites += 1;
				rParent->rLim=i;
				rParent->lLim=MIN(rParent->lLim,i);
			}
		}
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
								bpArray[breakNumber] = bp;
								breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
							bpArray[breakNumber] = bp;
							breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
							bpArray[breakNumber] = bp;
							breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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

/*proposeTrajectory-- this function creates a sweep trajectory and deals with
complications like changing population size, or soft sweeps, etc 
returns the acceptance probability of the trajectory */
double proposeTrajectory(int currentEventNumber, float *currentTrajectory, double *sizeRatio, char sweepMode, \
double initialFreq, double *finalFreq, double alpha, double f0, double currentTime)
{	
	double tInc, tIncOrig, minF,ttau, N;
	double N_0 = (double) EFFECTIVE_POPN_SIZE;
	double Nmax, localNextTime,localCurrentTime, currentSizeRatio;
	int i,j, insweepphase;
	float x;
	
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
					x = 1.0 - genicSelectionStochasticForwards(tInc, (1.0 - x), alpha * currentSizeRatio);
					break;
					case 'N':
					x = neutralStochastic(tInc, x);
					break;
				}
			}
			else{
				insweepphase = 0;
				tInc = 1.0 / (deltaTMod * N );
				x = neutralStochastic(tInc, x);
			}
			if(j>=maxTrajSteps){
				maxTrajSteps += 1000000;
				currentTrajectory = (float *) realloc(currentTrajectory,sizeof(float)*maxTrajSteps);
			}
			currentTrajectory[j++]=x;
		}
	}
	currentTrajectoryStep=0;
	totalTrajectorySteps=j;
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
					x = 1.0 - genicSelectionStochasticForwards(tInc, (1.0 - x), alpha * sizeRatio[0]);
				//	printf("x here:%f ttau: %f alpha*sizeRatio: %f\n",x,ttau,alpha*sizeRatio[0]);
					break;
					case 'N':
				//	printf("here\n");
					x = neutralStochastic(tInc, x);
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
							bpArray[breakNumber] = bp;
							breakNumber += 1; 
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
								bpArray[breakNumber] = bp;
								breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
							bpArray[breakNumber] = bp;
							breakNumber += 1; 
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
								bpArray[breakNumber] = bp;
								breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
					bpArray[breakNumber] = bp;
					breakNumber += 1; 
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
								printf("probAccept: %lf\n",probAccept);
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
			for(i=0;i<nSites;i++){
				if(i<xOver){
					lParent->ancSites[i] = aNode->ancSites[i];
					rParent->ancSites[i] = 0;
				}
				else{
					lParent->ancSites[i] = 0;
					rParent->ancSites[i] = aNode->ancSites[i];
				}

				if(lParent->ancSites[i] > 0 && lParent->ancSites[i] < sampleSize){
					lParent->nancSites += 1;
					lParent->rLim=i;
					lParent->lLim=MIN(lParent->lLim,i);
				}

				if(rParent->ancSites[i] > 0 && rParent->ancSites[i] < sampleSize){
					rParent->nancSites += 1;
					rParent->rLim=i;
					rParent->lLim=MIN(rParent->lLim,i);
				}
			}
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
		for(i=0;i<nSites;i++){
			if(i<xOver || i >= xOver+tractL){
				lParent->ancSites[i] = 0;
				rParent->ancSites[i] = aNode->ancSites[i];
			}
			else{
				lParent->ancSites[i] = aNode->ancSites[i];
				rParent->ancSites[i] = 0;
			}

			if(lParent->ancSites[i] > 0 && lParent->ancSites[i] < sampleSize){
				lParent->nancSites += 1;
				lParent->rLim=i;
				lParent->lLim=MIN(lParent->lLim,i);
			}

			if(rParent->ancSites[i] > 0 && rParent->ancSites[i] < sampleSize){
				rParent->nancSites += 1;
				rParent->rLim=i;
				rParent->lLim=MIN(rParent->lLim,i);
			}
		}
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
	for(i=0;i<nSites;i++){
		temp->ancSites[i] = lChild->ancSites[i] + rChild->ancSites[i];
		if(temp->ancSites[i] > 0 && temp->ancSites[i] < sampleSize){
			temp->nancSites += 1;
			temp->rLim=i;
			temp->lLim=MIN(temp->lLim,i);
			
		}
	}
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
				  if (allNodes[i]->lLim == allNodes[i]->rLim){
		  		    if (mutSite*nSites < (float)allNodes[i]->lLim)
					  error = (float)allNodes[i]->lLim-(mutSite*nSites);
		  		    else
					  error = 0.0;
				    p = allNodes[i]->rLim + (1.0/nSites);
				    mutSite = genunf(((float)allNodes[i]->lLim + error) / nSites, (p + error) / nSites);
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


void addMutation(rootedNode *aNode, float site){
	aNode->muts[aNode->mutationNumber] = site;
	aNode->mutationNumber += 1;

}

int hasMutation(rootedNode *aNode, float site){
	int i;
	for (i = 0; i < aNode->mutationNumber; i++){
		if (aNode->muts[i] == site) return 1;
	}
	return 0;
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
	

	rootIdx = findRootAtSite(site);
	newickRecurse(allNodes[rootIdx],site,0);
	printf(";\n");

}

void newickRecurse(rootedNode *aNode, float site, float tempTime){
	
	if(isCoalNode(aNode)){
		
		if(hasMaterialHere(aNode->leftChild,site) && hasMaterialHere(aNode->rightChild,site)){	
			printf("(");
			newickRecurse(aNode->leftChild,site,tempTime);
			printf(",");
			newickRecurse(aNode->rightChild,site,tempTime);
			printf(")");
			if(nAncestorsHere(aNode, site) != sampleSize){
				printf(":%f",(aNode->leftParent->time - aNode->time)*0.5);
			}
			
		}
		else{
			if(hasMaterialHere(aNode->leftChild,site)) newickRecurse(aNode->leftChild,site,tempTime+(aNode->leftParent->time - aNode->time));
			else if(hasMaterialHere(aNode->rightChild,site)) newickRecurse(aNode->rightChild,site,tempTime+(aNode->leftParent->time - aNode->time));
		}
	}
	else{
		if(isLeaf(aNode)){
			printf("%d:%f",aNode->id,(aNode->leftParent->time + tempTime - aNode->time)*0.5);
		}
		else{ //recombination node
			if(hasMaterialHere(aNode->leftChild,site)) newickRecurse(aNode->leftChild,site, tempTime + (aNode->leftParent->time - aNode->time));
		}
	}
}

/*makeGametesMS-- MS style sample output */
void makeGametesMS(int argc,const char *argv[]){
	int i,j, size, count, k, mutNumber;
	float allMuts[MAXMUTS];

	/* get unique list of muts */
	size = 0;
	for (i = 0; i < sampleSize; i++){
		for (j = 0; j < allNodes[i]->mutationNumber; j++){
			count = 0;
			for( k = 0; k < size; k++){
				if (allMuts[k] == allNodes[i]->muts[j]){
					count += 1;
				}
			}
			if (count == 0){
                                assert(size < MAXMUTS);
				allMuts[size] = allNodes[i]->muts[j];
				size++;
			}
		}
	}
	
	mutNumber = size;
	qsort(allMuts, size, sizeof(allMuts[0]), compare_floats);
	printf("\n//\nsegsites: %d",mutNumber);
	if(mutNumber > 0) printf("\npositions: ");
	for(i = 0; i < mutNumber; i++)
		fprintf(stdout,"%7.5lf ",allMuts[i] );
	fprintf(stdout,"\n");

/* go through sites, if there is a mut is allMuts print out segSite */
	for (i = 0; i < sampleSize; i++){
		for (j = 0; j < mutNumber; j++){
			if(isAncestralHere(allNodes[i],allMuts[j] )){
				if (hasMutation(allNodes[i],allMuts[ j])){
					printf("1");
				}
				else{
					printf("0");
				}
			}
			else{
				printf("N");
			}
		}
		printf("\n");
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
void addAncientSample(int lineageNumber, int popnDest, double addTime){
	int i;
	int count = 0;
	for(i=0; i < alleleNumber && count < lineageNumber; i++){
		if(nodes[i]->population == (popnDest+1) * -1){
			nodes[i]->population = popnDest;
			nodes[i]->time = addTime;
			//printf("time for %d: %f\n", i, nodes[i]->time);
			popnSizes[popnDest]++;
			count++;
		}
	}
	
}


void addNode(rootedNode *aNode){
        assert(alleleNumber < MAXNODES);
        assert(totNodeNumber < MAXNODES);
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
	int i;
	printf("node: %p time: %f lLim: %d rLim: %d nancSites: %d popn: %d\n nsites:\n",aNode, aNode->time,aNode->lLim,\
		aNode->rLim, aNode->nancSites, aNode->population );
	for(i=0;i<nSites;i++)printf("%d",aNode->ancSites[i]);
	printf("\n");
}

void freeTree(rootedNode *aNode){
	int i;
	//printf("final nodeNumber = %d\n",totNodeNumber);
	//cleanup nodes
	for (i = 0; i < totNodeNumber; i++){
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

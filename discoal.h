// discoal.h
// defs, headers, etc for discoal.c

#ifndef __DISCOAL_GLOBALS__
#define __DISCOAL_GLOBALS__

/******************************************************************************/
/* here are just some defines to allocate initial global arrays and stuff     */
/* like that                                                                  */

#define MAXNODES 20000000
#define MAXSITES 220020
#define SMALLCHUNKS 100000
#define MAXBREAKS 1000000
#define MAXMUTS 40000
#define MAXTIME 100000.0
#define MAXPOPS 6

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) > (b)) ? (b) : (a))
#define MAXEVENTS 200
#define MAXLEAFS 200

/******************************************************************************/

/******************************************************************************/
/* A "rootedNode" is just a node in a coalescent graph. Since recombination   */
/* is possible, each node has a left and a right parent as well as 2 children */
/* Nodes also keep track of their time in the tree and their mutations (which */
/* they pass on to their children) as well as the number of descendents       */
/* at each site                                                               */
/* Nodes are used to build the coalescent tree and place mutations            */

typedef struct rootedNode
{
	struct rootedNode *leftParent, *rightParent, *leftChild, *rightChild;
	double time, branchLength, blProb;
	float muts[MAXMUTS];
	#ifdef BIG
	uint16_t ancSites[MAXSITES];
	#else
	uint8_t ancSites[MAXSITES];
	#endif
	int nancSites, lLim, rLim;
	int id, mutationNumber, population, sweepPopn;
	int ndes[2],*leafs;
	double times[2];

}
rootedNode;

/******************************************************************************/

/******************************************************************************/
/* Here is the event object used to keep track of demographic changes         */


typedef struct event
{
	double time, popnSize;
	char type;
	int popID, popID2, popID3; //some events like splits involve two populations, admixture is 3
	double admixProp; // admixture proportion from population popID2
}
event;




/******************************************************************************/

/******************************************************************************/
/* here are all the global variables that are needed in the program. Most of  */
/* these have to do with storage for the coalescent graph and the chunks of   */
/* ancestral dna. There are also parameters controlling the coalescent        */
/* process                                                                    */

rootedNode  *nodes[MAXNODES], *allNodes[MAXNODES];

int activeMaterial[MAXSITES];

int sampleSize, sampleNumber, breakNumber, segSites,alleleNumber, \
	totNodeNumber, totChunkNumber, npops, eventFlag, nSites, activeSites,\
	mask, finiteOutputFlag, outputStyle, effectiveSampleSize, runMode, gcMean,\
	breakPoints[MAXBREAKS],sampleSizes[MAXPOPS];

double  leftRho, rho, theta, tDiv, alpha, sweepSite, tau, my_gamma, \
 	lambda, timeRecovery,bottleNeckRatio, bottleNeckDuration, \
 	ancestralSizeRatio, sweepLeft, sweepRight, thetas[MAXPOPS];
double mig[MAXPOPS];

int sampleS, sampleFD, sampleHaps, rejectCount, sampleRMin, offset, winNumber, \
	popnSizes[MAXPOPS],sweepPopnSizes[MAXPOPS];


const char *mFile;

char sweepMode, windowMode;

double coaltime, currentTime, pAccept;
int mn, eventNumber, migFlag;


double SweepStartingFrequency, f0;

double uA;

double gammaCoRatio, gammaCoRatioMode;
int priorTheta, priorRho, priorAlpha, priorTau, priorX, priorF0, priorE1, priorE2, priorUA, priorC;
double gammaCoRatioMode, gammaCoRatio;
double pThetaUp, pThetaLow,pRhoMean,pRhoUp,pRhoLow,pAlphaUp,pAlphaLow,pTauUp,pTauLow,pXUp,pXLow,pF0Up,pF0Low,pUALow,pUAUp,pCUp,pCLow;
double pE2TLow,pE1TLow, pE2THigh, pE1THigh, pE1SLow, pE1SHigh, pE2SLow,pE2SHigh;
double migMat[MAXPOPS][MAXPOPS], migMatConst[MAXPOPS][MAXPOPS];
double recurSweepRate;

int EFFECTIVE_POPN_SIZE;


struct event events[MAXEVENTS];

int lSpot, rSpot, condRecMode;
int condRecMet;
int activeSweepFlag;
int recurSweepMode;
int partialSweepMode,softSweepMode;
double partialSweepFinalFreq;

double deltaTMod;
int treeOutputMode;


#endif

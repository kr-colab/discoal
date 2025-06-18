// discoal.h
// defs, headers, etc for discoal.c

#ifndef __DISCOAL_GLOBALS__
#define __DISCOAL_GLOBALS__

#include <stdint.h>
#include "ancestrySegment.h"
#include "activeSegment.h"

#include <tskit.h>

/******************************************************************************/
/* Global constants and limits                                                */
/*                                                                            */

/* Still needed for various static arrays and limits */
#define MAXSITES 100000000   /* Maximum number of sites - used for input validation */
#define MAXMUTS 400000        /* Maximum mutations for output formatting */
#define MAXTIME 100000.0     /* Sentinel value representing "infinite" time */
#define MAXPOPS 121          /* Maximum number of populations */

/* No longer needed after dynamic memory optimizations:
   - MAXNODES: nodes array is now dynamic
   - SMALLCHUNKS: was never used  
   - MAXBREAKS: breakPoints array is now dynamic
   - MAXLEAFS: was never used
   - MAXEVENTS: events array is now dynamic
*/

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) > (b)) ? (b) : (a))

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
	// Essential pointers for tree structure
	struct rootedNode *leftParent, *rightParent, *leftChild, *rightChild;
	
	// Essential time and branch information
	double time, branchLength;
	
	// Ancestry information
	int nancSites, lLim, rLim;  // Calculated from ancestry tree
	int id, population, sweepPopn;
	
	// Ancestry segment tree for tracking which sites this node is ancestral to
	AncestrySegment *ancestryRoot;
	
	// Tskit node ID for tree sequence recording
	tsk_id_t tskit_node_id;
	
	// Node state tracking for safe memory management
	unsigned char parentsRecorded;  // Bit flags: left parent (bit 0), right parent (bit 1)
	unsigned char isFullyRecorded;  // 1 when all genealogical info is in tskit
	unsigned char inActiveSet;      // 1 if still in nodes[0..alleleNumber-1]
	unsigned char carriesSweepMutation;  // 1 if this node carries the beneficial mutation
	
	// DEPRECATED FIELDS - temporarily kept for compilation, will be removed
	int ndes[2];             // DEPRECATED - never used
	int *leafs;              // DEPRECATED - never used
	double times[2];         // DEPRECATED - never used
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
	int lineageNumber; //some events involve a specified number of lineages, like ancient samples
	double admixProp; // admixture proportion from population popID2
}
event;




/******************************************************************************/

/******************************************************************************/
/* here are all the global variables that are needed in the program. Most of  */
/* these have to do with storage for the coalescent graph and the chunks of   */
/* ancestral dna. There are also parameters controlling the coalescent        */
/* process                                                                    */

// Track only active nodes and sample node IDs (tskit mode)
rootedNode  **nodes;
int nodesCapacity;
// Array to track tskit IDs for sample nodes
tsk_id_t *sample_node_ids;
int sample_node_count;
int sample_node_capacity;


// int activeMaterial[MAXSITES];  // DEPRECATED - replaced by segment structure
ActiveMaterial activeMaterialSegments;  // New segment-based structure

int sampleSize, sampleNumber, breakNumber, segSites,alleleNumber, \
	totNodeNumber, totChunkNumber, npops, eventFlag, nSites, activeSites,\
	mask, finiteOutputFlag, outputStyle, effectiveSampleSize, runMode, gcMean,\
	*breakPoints, breakPointsCapacity, sampleSizes[MAXPOPS];

double  leftRho, rho, theta, tDiv, alpha, sweepSite, tau, my_gamma, \
 	lambda, timeRecovery,bottleNeckRatio, bottleNeckDuration, \
 	ancestralSizeRatio, sweepLeft, sweepRight, thetas[MAXPOPS];
double mig[MAXPOPS];

int sampleS, sampleFD, sampleHaps, rejectCount, sampleRMin, offset, winNumber, \
	popnSizes[MAXPOPS],sweepPopnSizes[MAXPOPS];


const char *mFile;

char sweepMode, windowMode;

double coaltime, currentTime, pAccept;
int mn, eventNumber, migFlag, currentEventNumber;


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

// Trajectory support
#define TRAJSTEPSTART 500000000
#define TRAJ_GROWTH_FACTOR 2
long int  maxTrajSteps;
long int  trajectoryCapacity;
float *currentTrajectory;
long int currentTrajectoryStep, totalTrajectorySteps;

// Node memory management statistics
int freedNodeCount;

// Memory-mapped trajectory support
char trajectoryFilename[256];  // Current trajectory file
int trajectoryFd;              // File descriptor for mmap
size_t trajectoryFileSize;     // Size of mmap'd region

struct event *events;          /* Dynamic array of demographic events */
int eventsCapacity;            /* Allocated capacity for events array */

int lSpot, rSpot, condRecMode;
int condRecMet;
int activeSweepFlag;
int recurSweepMode;
int partialSweepMode,softSweepMode;
double partialSweepFinalFreq;

double deltaTMod;

int ancSampleSize, ancPopID, ancSampleFlag;
int ancSampleTime;

int hidePartialSNP;

// Tree sequence recording
int tskitOutputMode;
char tskitOutputFilename[1024];
int minimalTreeSeq;  // Flag for minimal tree sequence recording (coalescence only)


#endif

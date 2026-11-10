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
	int id;
	int16_t population;      // MAXPOPS is 121, so 16-bit is sufficient
	int16_t sweepPopn;       // MAXPOPS is 121, so 16-bit is sufficient
	
	// Ancestry segment tree for tracking which sites this node is ancestral to
	AncestrySegment *ancestryRoot;
	
	// Tskit node ID for tree sequence recording
	tsk_id_t tskit_node_id;
	
	// Node state tracking for safe memory management
	unsigned char parentsRecorded;  // Bit flags: left parent (bit 0), right parent (bit 1)
	unsigned char isFullyRecorded;  // 1 when all genealogical info is in tskit
	unsigned char inActiveSet;      // 1 if still in nodes[0..alleleNumber-1]
	
	// Population list tracking for O(1) node selection
	int popListIndex;        // Index in population's node array (-1 if not in list)
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
	/* For an 'n' event whose time and size are drawn for each replicate:
	 * 1 + its index in priorSizeChanges. 0 for every other event. */
	int priorSlot;
}
event;

/******************************************************************************/
/* Shape state for time-varying demographic parameters                         */
/******************************************************************************/

/* Command-line and YAML times are in units of 4N0 generations, as in ms.
 * Internal time is in units of 2N0 generations. A rate per unit of time
 * (-eg alpha, -el gamma) converts by the inverse factor, so the product
 * rate * time is the same in both scales. */
#define USER_TO_INTERNAL_TIME 2.0
static inline double userTimeToInternal(double t) { return t * USER_TO_INTERNAL_TIME; }
static inline double userRateToInternal(double r) { return r / USER_TO_INTERNAL_TIME; }
static inline double internalTimeToUser(double t) { return t / USER_TO_INTERNAL_TIME; }

typedef enum {
    SHAPE_CONSTANT = 0,
    SHAPE_EXPONENTIAL = 1,
    SHAPE_LINEAR = 2
} ShapeType;

typedef struct {
    ShapeType type;
    double anchor_value;   /* size or rate at anchor_time */
    double rate_param;     /* alpha for EXP, gamma for LIN, unused for CONST.
                            * Forward-time rate of change per unit of internal
                            * time (2N0 generations); LIN is in the value's own
                            * units (relative size, or 4Nm for migration). */
    double anchor_time;    /* time at which anchor_value applies */
} Shape;

/******************************************************************************/
/* Population node list structure for O(1) node selection                     */
/******************************************************************************/

typedef struct {
    rootedNode **nodes;      // Array of node pointers
    int count;              // Number of active nodes in this population
    int capacity;           // Current array capacity
} PopulationNodeList;




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

int sampleSize, sampleNumber, breakNumber, alleleNumber, \
	totNodeNumber, totChunkNumber, npops, eventFlag, nSites, activeSites,\
	mask, outputStyle, effectiveSampleSize, gcMean,\
	*breakPoints, breakPointsCapacity, sampleSizes[MAXPOPS];

double  leftRho, rho, theta, tDiv, alpha, sweepSite, tau, my_gamma, \
 	lambda, timeRecovery,bottleNeckRatio, bottleNeckDuration, \
 	ancestralSizeRatio, sweepLeft, sweepRight, thetas[MAXPOPS];
double mig[MAXPOPS];

int sampleS, sampleFD, sampleHaps, rejectCount, sampleRMin, offset, winNumber, \
	popnSizes[MAXPOPS],sweepPopnSizes[MAXPOPS];

// Population node lists for O(1) node selection
PopulationNodeList popLists[MAXPOPS];


const char *mFile;

char sweepMode, windowMode;

double coaltime, currentTime, pAccept;
int mn, eventNumber, migFlag, currentEventNumber;


double SweepStartingFrequency, f0;

double uA;

double gammaCoRatio, gammaCoRatioMode;
int priorTheta, priorRho, priorAlpha, priorTau, priorX, priorF0, priorUA, priorC;
/* -S file: the sweep's trajectory read from a file, as points of time
 * before the present (internal units, increasing) and frequency of the
 * selected allele; trajFileRows is 0 without -S. */
struct trajPoint { double time, freq; };
int trajFileRows;
struct trajPoint *trajFilePoints;
double gammaCoRatioMode, gammaCoRatio;
double pThetaUp, pThetaLow,pRhoMean,pRhoUp,pRhoLow,pAlphaUp,pAlphaLow,pTauUp,pTauLow,pXUp,pXLow,pF0Up,pF0Low,pUALow,pUAUp,pCUp,pCLow;
/* Population size changes with uniform priors on time (internal units)
 * and size (relative to N0): -Pe1, -Pe2 and YAML
 * priors.population_size_changes. Each has one 'n' event in the event
 * list, marked with priorSlot. */
#define MAX_PRIOR_SIZE_CHANGES 64
struct priorSizeChange { double timeLow, timeHigh, sizeLow, sizeHigh; };
struct priorSizeChange priorSizeChanges[MAX_PRIOR_SIZE_CHANGES];
int numPriorSizeChanges;
double migMat[MAXPOPS][MAXPOPS], migMatConst[MAXPOPS][MAXPOPS];
extern double *currentSize;
double currentSizeConst[MAXPOPS];  /* snapshot of initial currentSize; restored per replicate */
Shape popShape[MAXPOPS];
Shape migShape[MAXPOPS][MAXPOPS];
double recurSweepRate;

int EFFECTIVE_POPN_SIZE;

/* Set by -v. Progress and diagnostic messages go to stderr only when it is
 * set, so a normal run writes nothing to stderr unless there is a problem. */
int verboseFlag;
#define VERBOSE_LOG(...) do { if (verboseFlag) fprintf(stderr, __VA_ARGS__); } while (0)

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

int activeSweepFlag;
int recurSweepMode;
int partialSweepMode,softSweepMode;
double partialSweepFinalFreq;

double deltaTMod;
extern int leftRhoFlag;   /* -ls/-ld/-ln; defined in discoal_multipop.c */
/* -L (or YAML recurrent_sweeps_left_of_locus): recurrent sweeps outside the
 * locus; defined in discoal_multipop.c. */
extern int recurrentLeftOfLocus;
/* A sweep position was given (-x, or YAML sweep_position); defined in
 * discoal_multipop.c. */
extern int sweepPositionGiven;


int hidePartialSNP;

// Tree sequence recording
int tskitOutputMode;
char tskitOutputFilename[1024];
int minimalTreeSeq;  // Flag for minimal tree sequence recording (coalescence only)


#endif

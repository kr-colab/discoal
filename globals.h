#ifndef __COAL_GLOBALS__
#define __COAL_GLOBALS__

/******************************************************************************/
/* here are just some defines to allocate initial global arrays and stuff     */
/* like that                                                                  */

#define MAXNODES 1000000
#define MAXCHUNKS 1000000
#define SMALLCHUNKS 100000
#define MAXBREAKS 1000000
#define MAXMUTS 50000
#define MAXTIME 1000.0
#define MAXPOPS 2
#define EFFECTIVE_POPN_SIZE 1000000
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) > (b)) ? (b) : (a))

/******************************************************************************/

/******************************************************************************/
/* A "chunk" is a piece of contiguous ancestral DNA material, with a left and */
/* a right end, with 0 <= leftEnd < rightEnd <= 1                             */
/*                                                                            */
/* chunks are used to keep track of portions of recombining DNA               */

typedef struct chunk
{
	double leftEnd, rightEnd;
}
chunk; 

/******************************************************************************/

/******************************************************************************/
/* A "rootedNode" is just a node in a coalescent graph. Since recombination   */
/* is possible, each node has a left and a right parent as well as 2 children */
/* Nodes also keep track of their time in the tree and their mutations (which */
/* they pass on to their children in the infinite alleles model) as well as   */
/* which "chunks" of ancestral DNA they contain                               */
/*                                                                            */
/* Nodes are used to build the coalescent tree and place mutations            */

typedef struct rootedNode
{
	struct rootedNode *leftParent, *rightParent, *leftChild, *rightChild;
	double time, *muts;
	struct chunk **chunks;
	int id, chunkNumber, mutationNumber, population, descNumber;
	int ndes[2];
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

chunk *allChunks[MAXCHUNKS], *activeMaterial[MAXCHUNKS];

int sampleSize, sampleNumber, breakNumber, segSites,alleleNumber, \
	totNodeNumber, totChunkNumber, npops, eventFlag, nSites,activeChunks,\
	mask, finiteOutputFlag, outputStyle, effectiveSampleSize, runMode;

double breakPoints[MAXBREAKS], rho, theta, tDiv, alpha, sweepSite, tau, \
 	lambda, timeRecovery,bottleNeckRatio, bottleNeckDuration, \
 	ancestralSizeRatio, sweepLeft, sweepRight ;

int sampleS, sampleFD, sampleHaps, rejectCount, sampleRMin, offset, winNumber, \
	popnSizes[MAXPOPS];

double samplePi, tolerance, sampleTajD, winSize, *samplePiWins;

const char *mFile;

char sweepMode, windowMode;

double coaltime, currentTime, pAccept;
int mn;

double SweepStartingFrequency;

/******************************************************************************/

#endif

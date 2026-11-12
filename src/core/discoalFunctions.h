#ifndef __DISCOAL_H__
#define __DISCOAL_H__

void initialize();
void ensureEventsCapacity();

// Population list management functions for O(1) node selection
void initializePopLists();
void cleanupPopLists();
void addNodeToPopList(rootedNode *node, int popn);
void removeNodeFromPopList(rootedNode *node);
rootedNode *pickNodePopnFast(int popn);

// Sample node tracking functions for tskit mode
void initializeSampleNodeIds();
void ensureSampleNodeCapacity(int required_size);
void addSampleNodeId(tsk_id_t tskit_node_id);
void cleanupSampleNodeIds();
void initializeMuts(rootedNode *node, int capacity);
void ensureMutsCapacity(rootedNode *node, int requiredSize);
void cleanupMuts(rootedNode *node);
rootedNode *newRootedNode(double cTime, int popn);
void releaseFreeLineages(void);

// Node memory management
void markLeafNodeRecorded(rootedNode *node);
void cleanupRemainingNodes();

void coalesceAtTimePopn(double cTime, int popn);
void coalesceAtTimePopnSweep(double cTime, int popn, int sp);
void migrateAtTime(double cTime,int srcPopn, int destPopn);

int recombineAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);
int recombineToLeftPopnSweep(int popn, int sp, double popnFreq);
void geneConversionAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);


void updateActiveMaterial(rootedNode *aNode);
void updateAncestryStatsFromTree(rootedNode *node);
int isAncestralHere(rootedNode *aNode, float site);

void addMutation(rootedNode *aNode, double site);
int hasMutation(rootedNode *aNode, double site);
void sortNodeMutations(rootedNode *node);
void sortAllMutations();
void makeGametesMS(int argc,const char *argv[]);


void dropMutationsRecurse();
void recurseTreePushMutation(rootedNode *aNode, float site);
void errorCheckMutations();

void mergePopns(int popnSrc, int popnDest);
void admixPopns(int popnSrc, int popnDest1, int popnDest2, double admixProp);
void addAncientSample(int lineageNumber, int popnDest, double addTime, int stillSweeping, double currentFreq);


void recurrentMutAtTime(double cTime,int srcPopn, int sp);

void ensureTrajectoryCapacity(long int requiredSize);

// Memory-mapped trajectory functions
void releaseTrajectory(void);
void mmapAcceptedTrajectory(int fd, long int numSteps);
void initializeNodeArrays();
void ensureNodesCapacity(int requiredSize);
void ensureAllNodesCapacity(int requiredSize);

			
			
double recurrentSweepPhaseGeneralPopNumber(double startTime, double endTime, double *finalFreq, double alpha, char sweepMode, double *sizeRatio);
		
double proposeTrajectory(int currentEventNumber, float *currentTrajectory, double *sizeRatio, char sweepMode, \
	double initialFreq, double *finalFreq, double alpha, double f0, double currentTime);
double sweepPhaseEventsConditionalTrajectory(double startTime, double endTime, double sweepSite,\
	double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
	double *sizeRatio, char sweepMode,double f0, double uA);


		
double totalTimeInTree();

void lineageEventRates(int popn, double *recombination, double *geneConversion);
int recombineAtTimePopn(double cTime, int popn);
void geneConversionAtTimePopn(double cTime, int popn);

double neutralPhaseGeneralPopNumber(double startTime, double endTime, double *sizeRatio);

rootedNode *pickNodePopn(int popn);
void addNode(rootedNode *aNode);
void removeNode(rootedNode *aNode);
void printNode(rootedNode *aNode);
int nodePopnSize(int popn);
int nodePopnSweepSize(int popn, int sp);
rootedNode *pickNodePopnSweep(int popn,int sp);

void sortEventArray(struct event *eArray, int eNumber);
int addPriorSizeChange(const char *what, double timeLow, double timeHigh,
                       double sizeLow, double sizeHigh, int popID);
void drawPriors();
void sampleSweepTrajectory(int currentEventNumber, double *sizeRatio,
	char sweepMode, double initialFreq, double *finalFreq, double alpha,
	double f0, double currentTime);
int setSweepTrajectoryFile(const char *path);
void buildSweepTrajectoryFromFile(void);
double sweepStartFreq(double t);

int findRootAtSite(float site);
void newickRecurse(rootedNode *aNode, float site,float tempTime);
void printTreeAtSite(float site);
void printAllNodes();

unsigned int devrand(void);
int compare_doubles(const void *a,const void *b);

#endif

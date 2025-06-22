#ifndef __DISCOAL_H__
#define __DISCOAL_H__

void initialize();
void initializeBreakPoints();
void ensureBreakPointsCapacity();
void cleanupBreakPoints();

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
void addBreakPoint(int bp);
void initializeMuts(rootedNode *node, int capacity);
void ensureMutsCapacity(rootedNode *node, int requiredSize);
void cleanupMuts(rootedNode *node);
void cleanupNodeArrays();
void cleanupRemovedNodes();
rootedNode *newRootedNode(double cTime, int popn);

// Node memory management
void markParentRecorded(rootedNode *child, rootedNode *parent);
void markLeafNodeRecorded(rootedNode *node);
void tryFreeNode(rootedNode *node);
void cleanupRemainingNodes();

void coalesceAtTimePopn(double cTime, int popn);
void coalesceAtTimePopnSweep(double cTime, int popn, int sp);
void migrateAtTime(double cTime,int srcPopn, int destPopn);
void migrateExceptSite(double site, double scalar, int srcPopn, int destPopn);

int recombineAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);
int recombineToLeftPopnSweep(int popn, int sp, double popnFreq);
void geneConversionAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);


void updateActiveMaterial(rootedNode *aNode);
void updateAncestryStatsFromTree(rootedNode *node);
int isActive(int site);
int isAncestralHere(rootedNode *aNode, float site);
int nAncestorsHere(rootedNode *aNode, float site);

int siteBetweenChunks(rootedNode *aNode, int xOverSite);
void dropMutations();
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
void mmapAcceptedTrajectory(const char *filename, long int numSteps);
void cleanupRejectedTrajectory(const char *filename);
void initializeNodeArrays();
void ensureNodesCapacity(int requiredSize);
void ensureAllNodesCapacity(int requiredSize);

double sweepPhaseEventsGeneralPopNumber(int *bpArray, double startTime, double endTime, double sweepSite,\
			double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
			double *sizeRatio, char sweepMode,double f0, double uA);
			
			
double recurrentSweepPhaseGeneralPopNumber(int *bpArray,double startTime, double endTime, double *finalFreq, double alpha, char sweepMode, double *sizeRatio);
		
double proposeTrajectory(int currentEventNumber, float *currentTrajectory, double *sizeRatio, char sweepMode, \
	double initialFreq, double *finalFreq, double alpha, double f0, double currentTime);
double sweepPhaseEventsConditionalTrajectory(int *bpArray, double startTime, double endTime, double sweepSite,\
	double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
	double *sizeRatio, char sweepMode,double f0, double uA);


		
double totalTimeInTree();

int recombineAtTimePopn(double cTime, int popn);
void geneConversionAtTimePopn(double cTime, int popn);

double neutralPhase(int *bpArray,double startTime, double endTime, double sizeRatio);
double neutralPhaseMig(int *bpArray,double startTime, double endTime, double sizeRatio);
double neutralPhaseMigExclude(int *bpArray,double startTime, double endTime, double sizeRatio, double selSite, double migScale);
double neutralPhaseGeneralPopNumber(int *bpArray,double startTime, double endTime, double *sizeRatio);

rootedNode *pickNodePopn(int popn);
void addNode(rootedNode *aNode);
void removeNodeAt(int index);
void removeNode(rootedNode *aNode);
void addNodeAtIndex(rootedNode *aNode, int anIndex);
void shiftNodes(int offset);
void printNode(rootedNode *aNode);
void freeTree(rootedNode *aNode);
int nodePopnSize(int popn);
int nodePopnSweepSize(int popn, int sp);
rootedNode *pickNodePopnSweep(int popn,int sp);

int compare_events(const void *a,const void *b);
void sortEventArray(struct event *eArray, int eNumber);

int findRootAtSite(float site);
int hasMaterialHere(rootedNode *aNode, float site);
int isLeaf(rootedNode *aNode);
int isCoalNode(rootedNode *aNode);
void newickRecurse(rootedNode *aNode, float site,float tempTime);
void printTreeAtSite(float site);
void printAllNodes();
void printAllActiveNodes();

unsigned int devrand(void);
int compare_doubles(const void *a,const void *b);
int compare_floats(const void *a,const void *b);

#endif

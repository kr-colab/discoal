#ifndef __DISCOAL_H__
#define __DISCOAL_H__

void initialize();
rootedNode *newRootedNode(double cTime, int popn);

void coalesceAtTimePopn(double cTime, int popn);
void coalesceAtTimePopnSweep(double cTime, int popn, int sp);
void migrateAtTime(double cTime,int srcPopn, int destPopn);
void migrateExceptSite(double site, double scalar, int srcPopn, int destPopn);

int recombineAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);
void geneConversionAtTimePopnSweep(double cTime, int popn, int sp, double sweepSite, double popnFreq);


void updateActiveMaterial(rootedNode *aNode);
int isActive(int site);
int isAncestralHere(rootedNode *aNode, float site);
int siteBetweenChunks(rootedNode *aNode, int xOverSite);
void dropMutations();
void addMutation(rootedNode *aNode, float site);
int hasMutation(rootedNode *aNode, float site);
void makeGametesMS(int argc,const char *argv[]);

void mergePopns(int popnSrc, int popnDest);
void admixPopns(int popnSrc, int popnDest1, int popnDest2, double admixProp);


void recurrentMutAtTime(double cTime,int srcPopn, int sp);
double sweepPhaseEvents(int *bpArray, double startTime, double endTime, double sweepSite,\
				  double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
				double sizeRatio, char sweepMode,double f0);
double sweepPhaseEventsWithRecurrentMut(int *bpArray, double startTime, double endTime, double sweepSite,\
                                  double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
                                double sizeRatio, char sweepMode,double f0,double uA);
double sweepPhaseEventsGeneralPopNumber(int *bpArray, double startTime, double endTime, double sweepSite,\
			double initialFreq, double *finalFreq, int *stillSweeping, double alpha,\
			double *sizeRatio, char sweepMode,double f0, double uA);
			
double totalTimeInTree();
void dropMutationsUntilTime(double t);
double totalTimeInTreeUntilTime(double t);

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

unsigned int devrand(void);
int compare_doubles(const void *a,const void *b);
int compare_floats(const void *a,const void *b);

#endif

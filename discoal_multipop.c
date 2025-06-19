// discoal_multipop.c 
//
// a discrete sequence coalescent simulation
// A. Kern 10/3/14
//




#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <signal.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "ranlib.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "alleleTraj.h"
#include "tskitInterface.h"



int locusNumber; 
int leftRhoFlag=0;
int untilMode = 0;
const char *fileName;
double uTime;
double *currentSize;
long seed1, seed2;
double nextTime, currentFreq;
//float *currentTrajectory;

void getParameters(int argc,const char **argv);
void usage();
void cleanup_and_exit(int sig);
void initializeTrajectoryVariables(double N);

// Signal handler for cleanup
void cleanup_and_exit(int sig) {
	// Clean up trajectory storage if it exists
	if (trajectoryFd != -1) {
		if (currentTrajectory && currentTrajectory != MAP_FAILED) {
			munmap(currentTrajectory, trajectoryFileSize);
		}
		close(trajectoryFd);
		if (trajectoryFilename[0] != '\0') {
			unlink(trajectoryFilename);
		}
	}
	exit(sig);
}

// Initialize all trajectory-related variables
void initializeTrajectoryVariables(double N) {
	// Initialize time-related variables
	currentTime = 0;
	nextTime = 999;
	
	// Initialize population size
	currentSize[0] = 1.0;
	
	// Initialize frequency (just to initialize the value)
	currentFreq = 1.0 - (1.0 / (2.0 * N * currentSize[0]));
	
	// Initialize trajectory steps
	maxTrajSteps = trajectoryCapacity;
	
	// Initialize sweep flag
	activeSweepFlag = 0;
}

// Helper function to ensure events array has enough capacity
void ensureEventsCapacity() {
	if (eventNumber >= eventsCapacity) {
		int newCapacity = eventsCapacity * 2;
		struct event *newEvents = (struct event*) realloc(events, newCapacity * sizeof(struct event));
		if (newEvents == NULL) {
			fprintf(stderr, "Error: Failed to reallocate events array\n");
			exit(1);
		}
		events = newEvents;
		eventsCapacity = newCapacity;
		// Initialize new elements to zero
		memset(&events[eventNumber], 0, (newCapacity - eventNumber) * sizeof(struct event));
	}
}

int main(int argc, const char * argv[]){
	int i,j,k, totalSimCount;
	float tempSite;
	int lastBreak;
	double probAccept;
	
	
	
	getParameters(argc,argv);
	double N = EFFECTIVE_POPN_SIZE; // effective population size
	setall(seed1, seed2 );
	
	// Register signal handlers for cleanup
	signal(SIGINT, cleanup_and_exit);
	signal(SIGTERM, cleanup_and_exit);
	signal(SIGSEGV, cleanup_and_exit);

	//Hudson style header
	for(i=0;i<argc;i++)printf("%s ",argv[i]);
	printf("\n%ld %ld\n", seed1, seed2);
	
	i = 0;
        totalSimCount = 0;
	trajectoryCapacity = TRAJSTEPSTART;
	trajectoryFd = -1;  // Initialize to invalid
	trajectoryFilename[0] = '\0';  // Empty filename
	currentTrajectory = NULL;  // Will be mmap'd when needed
	
	// Initialize global trajectory generator (lazy approach)

	while(i < sampleNumber){
		// Initialize all trajectory-related variables
		initializeTrajectoryVariables(N);
		
		// Always initialize tskit for each replicate (must be before initialize())
		// Debug output commented out for production
		// if (tskitOutputMode && i == 0) {
		//	fprintf(stderr, "Initializing tskit with sequence length %d (replicate %d)\n", nSites, i+1);
		// }
		int ret = tskit_initialize((double)nSites);
		if (ret != 0) {
			fprintf(stderr, "Error initializing tskit: %s\n", tsk_strerror(ret));
			exit(1);
		}
		
		
		initialize();

		j=0;
		activeSweepFlag = 0;
		for(j=0;j<eventNumber && alleleNumber > 1;j++){
			currentEventNumber=j; //need this annoying global for trajectory generation
			if(j == eventNumber - 1){
				nextTime = MAXTIME;
				}
			else{
				nextTime = events[j+1].time;
			}
			// printf("type: %c popID: %d size: %f currentTime: %f nextTime: %f alleleNumber: %d\n",events[j].type,events[j].popID,events[j].popnSize,
			//	currentTime,nextTime,alleleNumber);	
			switch(events[j].type){
				case 'n':
				currentTime = events[j].time;
				currentSize[events[j].popID] = events[j].popnSize;
				//for(i=0;i<npops;i++)
				//	for(j=0;j<npops;j++) printf("%f\n",migMat[i][j]);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
			//	printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],
			//							sweepPopnSizes[0]);
				break;
				case 's':
				assert(activeSweepFlag == 0);
				currentTime = events[j].time;
				if (partialSweepMode == 1){
					currentFreq = MIN(partialSweepFinalFreq,1.0 - (1.0 / (2.0 * N * currentSize[0])));
				}
				else{
					currentFreq = 1.0 - (1.0 / (2.0 * N * currentSize[0]));
				}
			//	printf("event%d currentTime: %f nextTime: %f popnSize: %f\n",j,currentTime,nextTime,currentSize);

				//generate a proposed trajectory
				char previousTrajectoryFile[256] = "";
				probAccept = proposeTrajectory(currentEventNumber, currentTrajectory, currentSize, sweepMode, currentFreq, &currentFreq, alpha, f0, currentTime);
				while(ranf()>probAccept){
					// Clean up rejected trajectory
					if (previousTrajectoryFile[0] != '\0') {
						cleanupRejectedTrajectory(previousTrajectoryFile);
					}
					strcpy(previousTrajectoryFile, trajectoryFilename);
					
					probAccept = proposeTrajectory(currentEventNumber, currentTrajectory, currentSize, sweepMode, currentFreq, &currentFreq, alpha, f0, currentTime);
					//printf("probAccept: %lf\n",probAccept);
				}
				
				// Clean up any remaining rejected trajectory
				if (previousTrajectoryFile[0] != '\0' && strcmp(previousTrajectoryFile, trajectoryFilename) != 0) {
					cleanupRejectedTrajectory(previousTrajectoryFile);
				}
				
				// Now mmap the accepted trajectory
				mmapAcceptedTrajectory(trajectoryFilename, totalTrajectorySteps);
				
				currentTime = sweepPhaseEventsConditionalTrajectory(&breakPoints[0], currentTime, nextTime, sweepSite, \
					 currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
				//printf("currentFreqAfter: %f alleleNumber:%d currentTime:%f\n",currentFreq,alleleNumber,currentTime);
				//printf("pn0:%d pn1:%d alleleNumber: %d sp1: %d sp2: %d \n", popnSizes[0],popnSizes[1], alleleNumber,sweepPopnSizes[1],
				//			sweepPopnSizes[0]);
				if (currentTime < nextTime)
                                        currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
						
				break;
				case 'p': //merging populations
				currentTime = events[j].time;
				//printf("here at P flag time=%f\n",currentTime);
				mergePopns(events[j].popID, events[j].popID2);

				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                        		currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				break;
				case 'a':
				currentTime = events[j].time;
				admixPopns(events[j].popID, events[j].popID2, events[j].popID3, events[j].admixProp);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				break;
				case 'A':
				currentTime = events[j].time;
				//printAllActiveNodes();
				addAncientSample(events[j].lineageNumber, events[j].popID, events[j].time, activeSweepFlag, currentFreq);
				//printAllActiveNodes();
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode,currentSize);
					}
				}
				break;
			}
			
		}
		//finish up the coalescing action!
		if(alleleNumber > 1){
			currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, MAXTIME, currentSize);
		}
		
		// Sort and simplify the tree sequence after ancestry simulation is complete
		// This removes all non-ancestral nodes and edges before mutation placement
		if (tsk_tables != NULL) {
			// Sort tables (required before simplification)
			int ret = tsk_table_collection_sort(tsk_tables, NULL, 0);
			if (ret != 0) {
				fprintf(stderr, "Error sorting tables after ancestry simulation: %s\n", tsk_strerror(ret));
				exit(1);
			}
			
			// Get the sample node IDs before simplification
			tsk_id_t *samples = malloc(sample_node_count * sizeof(tsk_id_t));
			if (samples == NULL) {
				fprintf(stderr, "Error: Failed to allocate memory for samples\n");
				exit(1);
			}
			for (int i = 0; i < sample_node_count; i++) {
				samples[i] = sample_node_ids[i];
			}
			
			// Simplify to remove non-ancestral material
			tsk_id_t *node_map = malloc(tsk_tables->nodes.num_rows * sizeof(tsk_id_t));
			if (node_map == NULL) {
				fprintf(stderr, "Error: Failed to allocate memory for node map\n");
				free(samples);
				exit(1);
			}
			
			ret = tsk_table_collection_simplify(tsk_tables, samples, sample_node_count, 0, node_map);
			if (ret != 0) {
				fprintf(stderr, "Error simplifying tables after ancestry simulation: %s\n", tsk_strerror(ret));
				free(samples);
				free(node_map);
				exit(1);
			}
			
			// Update the sample node IDs to reflect the new node numbers after simplification
			for (int i = 0; i < sample_node_count; i++) {
				if (samples[i] != TSK_NULL && node_map[samples[i]] != TSK_NULL) {
					sample_node_ids[i] = node_map[samples[i]];
				}
			}
			
			free(samples);
			free(node_map);
		}
		
		//assign root
	//	root = nodes[0];
		//add Mutations
		// First record any sweep mutations
		if (sweepMode != 'n' && sweepSite >= 0.0) {
			if (tskit_record_sweep_mutations(sweepSite) < 0) {
				fprintf(stderr, "Error: Failed to record sweep mutations in tskit\n");
				exit(1);
			}
		}
		
		// Always use tskit for mutation placement
		extern double theta;
		if (untilMode==0) {
			// Use edge-based algorithm for tskit mutations
			if (tskit_place_mutations_edge_based(theta) < 0) {
				fprintf(stderr, "Error: Failed to place mutations in tskit\n");
				exit(1);
			}
			
			// Only populate discoal mutation arrays if we need ms output
			// (i.e., not in conditional recombination mode)
			if (condRecMode == 0) {
				if (tskit_populate_discoal_mutations() < 0) {
					fprintf(stderr, "Error: Failed to populate discoal mutations from tskit\n");
					exit(1);
				}
			}
		} else {
			// For untilMode, need to implement tskit-based solution
			// TODO: Implement tskit-based mutation placement until time
			fprintf(stderr, "Error: -u option not yet supported with tskit-only mode\n");
			exit(1);
		}

		if(condRecMode == 0){
			// Hudson style output
			//errorCheckMutations();
			// Only generate MS output if not using tskit-only output mode
			if (!tskitOutputMode) {
				makeGametesMS(argc,argv);
			}
			//printf("rep: %d\n",i);
			i++;
		}
		else{
			if(condRecMet == 1){
                                i++;
				// Only generate MS output if not using tskit-only output mode
				if (!tskitOutputMode) {
					makeGametesMS(argc,argv);
				}
				condRecMet = 0;
			}
	
		}
		
		// Debug: Print memory statistics (commented out for production)
		// if (tskitOutputMode || minimalTreeSeq) {
		// 	fprintf(stderr, "Replicate %d: Total nodes created: %d, Nodes freed during sim: %d, Net nodes: %d\n", 
		// 	        i, totNodeNumber, freedNodeCount, totNodeNumber - freedNodeCount);
		// }
		
		// Clean up nodes for next replicate (tskit-only mode)
		// In tskit-only mode, we need minimal cleanup since tskit handles most memory
		for (int cleanup_i = 0; cleanup_i < alleleNumber && cleanup_i < 1000; cleanup_i++) {
			if (nodes[cleanup_i] != NULL) {
				// Free ancestry segments if they exist
				if (nodes[cleanup_i]->ancestryRoot) {
					freeSegmentTree(nodes[cleanup_i]->ancestryRoot);
					nodes[cleanup_i]->ancestryRoot = NULL;
				}
				// Free the node itself
				free(nodes[cleanup_i]);
				nodes[cleanup_i] = NULL;
			}
		}
		cleanupBreakPoints();
		
		// Clean up trajectory after each simulation
		if (trajectoryFd != -1) {
			if (currentTrajectory && currentTrajectory != MAP_FAILED) {
				munmap(currentTrajectory, trajectoryFileSize);
			}
			close(trajectoryFd);
			if (trajectoryFilename[0] != '\0') {
				unlink(trajectoryFilename);
			}
			trajectoryFd = -1;
			trajectoryFilename[0] = '\0';
			currentTrajectory = NULL;
		}
		
		// Finalize and write tskit tree sequence for this replicate
		if (tskitOutputMode) {
			// Create replicate-specific filename with larger buffer to prevent overflow
			char replicate_filename[1200];  // Extra space for _repN suffix
			if (sampleNumber == 1) {
				// Single replicate - use original filename
				strncpy(replicate_filename, tskitOutputFilename, sizeof(replicate_filename) - 1);
				replicate_filename[sizeof(replicate_filename) - 1] = '\0';
			} else {
				// Multiple replicates - append replicate number
				// i was already incremented, so it gives us the correct 1-based replicate number
				int replicate_num = i;
				char *dot = strrchr(tskitOutputFilename, '.');
				if (dot) {
					// Insert replicate number before extension
					int prefix_len = dot - tskitOutputFilename;
					int ret = snprintf(replicate_filename, sizeof(replicate_filename), 
						"%.*s_rep%d%s", prefix_len, tskitOutputFilename, replicate_num, dot);
					if (ret >= sizeof(replicate_filename)) {
						fprintf(stderr, "Error: Output filename too long\\n");
						exit(1);
					}
				} else {
					// No extension - append replicate number
					int ret = snprintf(replicate_filename, sizeof(replicate_filename), 
						"%s_rep%d", tskitOutputFilename, replicate_num);
					if (ret >= sizeof(replicate_filename)) {
						fprintf(stderr, "Error: Output filename too long\\n");
						exit(1);
					}
				}
			}
			
			if (i == 0) {
				fprintf(stderr, "Writing tree sequence to %s\n", replicate_filename);
			}
			int ret = tskit_finalize(replicate_filename);
			if (ret != 0) {
				fprintf(stderr, "Error writing tree sequence: %s\n", tsk_strerror(ret));
			}
		} else {
			// No file output - just need to sort tables for internal consistency
			// (tskit_finalize with NULL filename would fail)
		}
		
		// Clean up tskit for this replicate
		tskit_cleanup();
		
                totalSimCount += 1;
	}
        if(condRecMode == 1)
        {
            fprintf(stderr, "Needed run %d simulations to get %d with a recombination event within the specified bounds.\n", totalSimCount, i);
        }
	// Clean up any remaining trajectory storage
	if (trajectoryFd != -1) {
		if (currentTrajectory && currentTrajectory != MAP_FAILED) {
			munmap(currentTrajectory, trajectoryFileSize);
		}
		close(trajectoryFd);
		if (trajectoryFilename[0] != '\0') {
			unlink(trajectoryFilename);
		}
	}
	
	// tskit finalization now handled per-replicate
	
	free(currentSize);
		free(events);
	
	// Clean up node arrays
	cleanupNodeArrays();
	
	// Debug: print freed node count
	if (freedNodeCount > 0) {
		fprintf(stderr, "Debug: Freed %d nodes during simulation\n", freedNodeCount);
	}
	
	return(0);
}



void getParameters(int argc,const char **argv){
	int args;
	int i,j;
	double migR;
	int selCheck,nChangeCheck;
	
	if( argc < 3){
		usage();
	}

	sampleSize = atoi(argv[1]);
	if(sampleSize > 65535){
		printf("Error: sampleSize > 65535. This exceeds the maximum supported by uint16_t ancestry counts.\n");
		exit(666);
	}
	sampleNumber = atoi(argv[2]);
	nSites = atoi(argv[3]);
	if(nSites>MAXSITES){
		printf("Error: number of sites set higher than current compilation limit. Please reduce the number of sites or change the MAXSITES define and recompile\n");
		exit(666);
	}
	args = 4;

	npops = 1;
	popnSizes[0]=sampleSize;
	popnSizes[1]=0;
	sampleSizes[0]=sampleSize;
	sampleSizes[1]=0;
	leftRho = 0.0;
	rho = 0.0;
	my_gamma = 0.0;
	gcMean = 0;
	theta = 0.0;
	alpha = 0.0;
	lambda = 0.0;
	tau = 0;
	ancestralSizeRatio = 1.0;
	f0=0.0;
        uA=0.0;

	seed1 = (long) (devrand() % 2147483399);
	seed2 = (long) (devrand() % 2147483399);
	

	EFFECTIVE_POPN_SIZE = 1000000;
	sweepSite = 0.5;
	tDiv=666;
        gammaCoRatioMode = 0;
	priorTheta=priorRho=priorAlpha=priorTau=priorX=priorF0=priorUA=priorC=0;

	eventFlag = 1;
	effectiveSampleSize = sampleSize;
	finiteOutputFlag = 0;
	outputStyle = 'h';
	mask = 0;
	migFlag = 0;
	deltaTMod = 40;
	recurSweepMode = 0;
	partialSweepMode = 0;
	softSweepMode = 0;
	ancSampleFlag = 0;
	ancSampleSize = 0;
	hidePartialSNP = 0;
	tskitOutputMode = 0;
	tskitOutputFilename[0] = '\0';
	minimalTreeSeq = 1;  // Default to minimal tree sequence mode
	
	// Initialize events array with initial capacity
	eventsCapacity = 50;  // Start with reasonable capacity
	events = (struct event*) calloc(eventsCapacity, sizeof(struct event));
	if (events == NULL) {
		fprintf(stderr, "Error: Failed to allocate events array\n");
		exit(1);
	}
	
	//set up first bogus event
	eventNumber = 0;
	events[eventNumber].time = 0.0;
	events[eventNumber].popID = 0;
	events[eventNumber].popnSize = 1.0;
	events[eventNumber].type = 'n';
	eventNumber++;
	currentSize = malloc(sizeof(double) * MAXPOPS);

	condRecMode= 0;
	while(args < argc){
		switch(argv[args][1]){
			case 'F' :
			if (!tskitOutputMode) {
				fprintf(stderr, "Error: -F flag requires -ts to be specified first\n");
				fprintf(stderr, "Usage: %s [options] -ts <filename.trees> -F\n", argv[0]);
				exit(1);
			}
			minimalTreeSeq = 0;  // Disable minimal mode, use full recording
			fprintf(stderr, "Note: Full tree sequence mode enabled (recording all edges including recombination)\n");
			break;
			case 'S' :
			runMode = 'S';
			fileName = argv[++args];
			break;
			case 's' :
			segSites =  atoi(argv[++args]);
			break;
			case 't' :
			if (argv[args][2] == 's') {  // -ts flag for tree sequence
				tskitOutputMode = 1;
				args++;
				if (args >= argc || argv[args] == NULL || argv[args][0] == '-' || strlen(argv[args]) == 0) {
					fprintf(stderr, "Error: -ts flag requires a filename argument\n");
					fprintf(stderr, "Usage: %s [options] -ts <filename.trees>\n", argv[0]);
					fprintf(stderr, "Example: %s 10 1 1000 -t 5 -ts output.trees\n", argv[0]);
					if (args < argc && argv[args] != NULL && argv[args][0] == '-') {
						fprintf(stderr, "Note: '%s' looks like another flag, not a filename\n", argv[args]);
					}
					exit(1);
				}
				strcpy(tskitOutputFilename, argv[args]);
			}
			else {  // -t flag for theta
				theta = atof(argv[++args]);
			}
			break;
			case 'i' :
			deltaTMod = atof(argv[++args]);
			break;
			case 'r' :
			rho = atof (argv[++args]);
			break;
			case 'g' :
                        if (argv[args][2] == 'r')
                        {
			  gammaCoRatio = atof (argv[++args]);
			  gcMean = atoi (argv[++args]);
                          gammaCoRatioMode = 1;
                        }
                        else
                        {
			  my_gamma = atof (argv[++args]);
			  gcMean = atoi (argv[++args]);
                        }
			break;
			case 'a' :
			alpha = atof(argv[++args]);
			break;
			case 'x' :
			sweepSite = atof(argv[++args]);
			break;
			case 'M' :
			if(npops==1){
				fprintf(stderr,"Error: attempting to set migration but only one population! Be sure that 'm' flags are specified after 'p' flag\n");
				exit(1);
			}
			migR = atof(argv[++args]);
			for(i=0;i<npops;i++){
				for(j=0;j<npops;j++){
					if(i!=j){
						migMatConst[i][j]=migR;
					}
					else{
						migMatConst[i][j]= 0.0;
					}
				}
			}
			migFlag = 1;
			break;
			case 'm' :
			// -m flag for migration
				if(npops==1){
					fprintf(stderr,"Error: attempting to set migration but only one population! Be sure that 'm' flags are specified after 'p' flag\n");
					exit(1);
				}
				i = atoi(argv[++args]);
				j = atoi(argv[++args]);
				migR = atof(argv[++args]);
				migMatConst[i][j]=migR;
				migFlag = 1;
			break;
			case 'p' :
			npops = atoi(argv[++args]);
			if(npops > MAXPOPS){
				fprintf(stderr,"Error: too many populations defined. Current maximum number = %d. Change MAXPOPS define in discoal.h and recompile... if you dare\n",MAXPOPS);
				exit(1);
			}
			for(i=0;i<npops;i++){
				sampleSizes[i]=atoi(argv[++args]);
				currentSize[i] = 1.0;
			}
			
			break;
			case 'e' :
				switch(argv[args][2]){
					case 'n':
					ensureEventsCapacity();
					events[eventNumber].time = atof(argv[++args]) * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popnSize = atof(argv[++args]);
					events[eventNumber].type = 'n'; //size change
					eventNumber++;
					break;
					case 'd' :
					tDiv =  atof(argv[++args]);
						ensureEventsCapacity();
					events[eventNumber].time = tDiv * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].type = 'p'; //pop split
					eventNumber++;
					break;
					case 'a' :
						ensureEventsCapacity();
					events[eventNumber].time = atof(argv[++args]) * 2.0;
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].popID3 = atoi(argv[++args]);
					events[eventNumber].admixProp = atof(argv[++args]);
					events[eventNumber].type = 'a'; //admix split
					eventNumber++;
					break;
				}
			break;
			case 'w':
				switch(argv[args][2]){
					case 'd':
					sweepMode = 'd';
					break;
					case 's':
					sweepMode = 's';
					break;
					case 'n':
					sweepMode = 'N';
					break;
					}
				tau = atof(argv[++args]) * 2.0;
					ensureEventsCapacity();
				events[eventNumber].time = tau;
				events[eventNumber].type = 's'; //sweep event
				eventNumber++;
				break;
			case 'l':
				switch(argv[args][2]){
					case 'd':
                	                sweepMode = 'd';
                        	        break;
                                	case 's':
	                                sweepMode = 's';
        	                        break;
                	                case 'n':
                        	        sweepMode = 'N';
                                	break;
				}
				sweepSite = -1.0;
				tau = atof(argv[++args]) * 2.0;
				leftRho = atof(argv[++args]) * 2.0;
				leftRhoFlag=1;
					ensureEventsCapacity();
				events[eventNumber].time = tau;
				events[eventNumber].type = 's'; //sweep event
				eventNumber++;
					break;
			case 'f':
			f0 = atof(argv[++args]);
			softSweepMode = 1;
			break;
                        case 'u':
                        uA = atof(argv[++args]);
                        break;
			case 'P' :
			  switch(argv[args][2]){
				case 't':
				  priorTheta = 1;
				  pThetaLow=atof(argv[++args]);
				  pThetaUp=atof(argv[++args]);
				break;
				case 'c':
				  priorC = 1;
			          partialSweepMode = 1;
				  pCLow=atof(argv[++args]);
				  pCUp=atof(argv[++args]);
				break;
                                case 'r':
                                  if (strlen(argv[args]) == 4 && argv[args][3] == 'e'){
                                        priorRho = 2;
                                        pRhoMean=atof(argv[++args]);
                                        pRhoUp=atof(argv[++args]);
                                  }
                                  else{
                                        priorRho = 1;
                                        pRhoLow=atof(argv[++args]);
                                        pRhoUp=atof(argv[++args]);
                                  }
                                break;
				case 'a':
				  priorAlpha = 1;
				  pAlphaLow=atof(argv[++args]);
				  pAlphaUp=atof(argv[++args]);
				break;
				case 'u':
                                	if (strlen(argv[args]) == 4 && argv[args][3] == 'A'){
					  priorUA = 1;
					  pUALow=atof(argv[++args]);
					  pUAUp=atof(argv[++args]);
					}
                                	else{
					  priorTau = 1;
					  pTauLow=atof(argv[++args]) * 2.0;
					  pTauUp=atof(argv[++args]) * 2.0;
					}
				break;
				case 'x':
				  priorX = 1;
				  pXLow=atof(argv[++args]);
				  pXUp=atof(argv[++args]);
				break;
				case 'f':
				  priorF0 = 1;
				  pF0Low=atof(argv[++args]);
				  pF0Up=atof(argv[++args]);
				break;
				case 'e':
				switch(argv[args][3]){
					case '1':
					priorE1 = 1;
					pE1TLow=atof(argv[++args])*2.0;
					pE1THigh=atof(argv[++args])*2.0;
					pE1SLow=atof(argv[++args]);
					pE1SHigh=atof(argv[++args]);
						ensureEventsCapacity();
					events[eventNumber].type = 'n';
					eventNumber++;
					break;
					case '2':
					priorE2 = 1;
					pE2TLow=atof(argv[++args])*2.0;
					pE2THigh=atof(argv[++args])*2.0;
					pE2SLow=atof(argv[++args]);
					pE2SHigh=atof(argv[++args]);
						ensureEventsCapacity();
					events[eventNumber].type = 'n';
					eventNumber++;
					break;
				}
				break;
			}
                        break;
			case 'U' :
			untilMode= 1;
			uTime=atof(argv[++args])*2.0;
			break; 
			case 'd' :
			seed1=atoi(argv[++args]);
			seed2=atoi(argv[++args]);
			break;
			case 'N' :
			EFFECTIVE_POPN_SIZE=atoi(argv[++args]);
			break;
			case 'C' :
			condRecMode= 1;
			condRecMet = 0;
			lSpot=atoi(argv[++args]);
			rSpot=atoi(argv[++args]);
			break;
			case 'R' :
			recurSweepMode = 1;
			sweepMode = 's';
			recurSweepRate = atof(argv[++args]);
			break;
			case 'L' :
			recurSweepMode = 1;
			sweepMode = 's';
			sweepSite = -1.0;
			recurSweepRate = atof(argv[++args]);
            if (recurSweepRate <= 0)
            {
                fprintf(stderr,"recurSweepRate must be > 0\n");
                exit(0);
            }
			break;
			case 'c' :
			partialSweepMode = 1;
			sweepMode = 's';
			partialSweepFinalFreq = atof(argv[++args]);
            if (partialSweepFinalFreq <= 0.0 || partialSweepFinalFreq >= 1.0)
            {
                fprintf(stderr,"partialSweepFinalFreq must be > 0 and < 1.0\n");
                exit(1);
            }
			break;
			case 'h' :
			hidePartialSNP = 1;
			break;
			case 'A' :
				ensureEventsCapacity();
			events[eventNumber].lineageNumber = atoi(argv[++args]);
			events[eventNumber].popID = atoi(argv[++args]);	 
			events[eventNumber].time = atof(argv[++args]) * 2.0; 
			ancSampleSize += events[eventNumber].lineageNumber;
			events[eventNumber].type = 'A'; //ancient sample
			ancSampleFlag = 1;
			eventNumber++;
			assert(events[eventNumber-1].lineageNumber < sampleSize);
			break;
			case 'T' :
			fprintf(stderr, "Error: -T tree output mode has been removed.\\n");
			fprintf(stderr, "Use -ts <filename.trees> for tree sequence output instead.\\n");
			exit(1);
			break;
			 
		}
		args++;
	}
	sortEventArray(events,eventNumber);

	//make sure events are kosher
	selCheck = 0;
	nChangeCheck=0;
	for(i=0;i<eventNumber;i++){
		//printf("event %d: type is %c\n", i, events[i].type);
	 		if(events[i].type == 's'){
				selCheck = 1;
	 		}
			if(events[i].type == 'n'){
				nChangeCheck += 1;
	 		}
	}
	if(selCheck == 1){
		if(recurSweepMode == 1){
			printf("Error with event specification: a single sweep event has been found but recurrentSweep mode has been specified\n");
			exit(666);
		}
		if(nChangeCheck > 1 && sweepMode=='d'){
			printf("Error with event specification: you chose 1 or more population size changes with a deterministic sweep. Please us -ws flag instead\n");
			exit(666);
		}
		if(softSweepMode == 1 && partialSweepMode == 1){
			if(f0 >= partialSweepFinalFreq){
				printf("Error with event specification: you specified a partial soft sweep but final frequency of partial sweep <= f_0\n");
				exit(666);
			}
		}
	}
	if(leftRhoFlag && sweepSite >= 0.0){
		printf("Error with event specification: you chose leftRho mode but the sweep site is within the locus\n");
		exit(666);
	}
	if(softSweepMode == 1 && recurSweepMode == 1){
		printf("Error with event specification: currently recurrent soft sweeps are not implemented. this will be a future addition\n");
		exit(666);
	}
	
}
		


void usage(){
	fprintf(stderr,"usage: discoal sampleSize numReplicates nSites -ws tau\n");
	fprintf(stderr,"parameters: \n");

	fprintf(stderr,"\t -t theta \n");
	fprintf(stderr,"\t -r rho (=zero if not specified)\n");
	fprintf(stderr,"\t -g conversionRate tractLengthMean (gene conversion)\n");
	fprintf(stderr,"\t -gr conversionToCrossoverRatio tractLengthMean (gene conversion where initiation rate = rho*conversionToCrossoverRatio)\n");
	fprintf(stderr,"\t -p npops sampleSize1 sampleSize2 etc.\n");
	fprintf(stderr,"\t -en time popnID size (changes size of popID)\n");	
	fprintf(stderr,"\t -ed time popnID1 popnID2 (joins popnID1 into popnID2)\n");
	fprintf(stderr,"\t -ea time daughterPopnID founderPopnID1 founderPopnID2 admixProp (admixture-- back in time daughterPopnID into two founders)\n");
	
	fprintf(stderr,"\t -ws tau (sweep happend tau generations ago- stochastic sweep)\n");  
	fprintf(stderr,"\t -wd tau (sweep happend tau generations ago- deterministic sweep)\n"); 
	fprintf(stderr,"\t -wn tau (sweep happend tau generations ago- neutral sweep)\n");
	fprintf(stderr,"\t -ls tau leftRho (stochastic sweep some genetic distance to the left of the simulated window--specified by leftRho=4Nr)\n");
	fprintf(stderr,"\t\t similarly, ld and ln simulate deterministic and neutral sweeps to the left of the window, respectively\n");
	fprintf(stderr,"\t -f first frequency at which selection acts on allele (F0; sweep models only)\n");
	fprintf(stderr,"\t -uA rate at which adaptive mutation recurs during the sweep phase (sweep models only)\n");
	fprintf(stderr,"\t -N sweepEffectivePopnSize (sweep models only)\n");	
	fprintf(stderr,"\t -a alpha (=2Ns)\n");
	fprintf(stderr,"\t -x sweepSite (0-1)\n");
	fprintf(stderr,"\t -c partialSweepFinalFrequency (partial sweeps)\n");	
	fprintf(stderr,"\t -i dt (sweep time increment scalar; default 400 -> 1/400N)\n");
	
	fprintf(stderr,"\t -M migRate (sets all rates to migRate)\n");
	fprintf(stderr,"\t -m popnID1 popnID2 migRate (sets migRate from popnID1 to popnID2)\n");
	fprintf(stderr,"\t -A sampleSize popnID time (ancient sample from popnID at specified time)\n");
	
	fprintf(stderr,"\t -Pt low high (prior on theta)\n");
	fprintf(stderr,"\t -Pr low high (prior on rho)\n");
        fprintf(stderr,"\t -Pre mean upperBound (prior on rho -- exponentially distributed but truncated at an upper bound)\n");
	fprintf(stderr,"\t -Pa low high (prior on alpha)\n");
	fprintf(stderr,"\t -Pu low high (prior on tau; sweep models only; still must use \"-ws tau\" and \"tau\" will be ignored)\n");
	fprintf(stderr,"\t -PuA low high (prior on uA; sweep models only)\n");
	fprintf(stderr,"\t -Px low high (prior on sweepSite; sweep models only)\n");
	fprintf(stderr,"\t -Pf low high (prior on F0; sweep models only)\n");
	fprintf(stderr,"\t -Pc low high (prior on partialSweepFinalFreq; sweep models only)\n");
	fprintf(stderr,"\t -Pe1 lowTime highTime lowSize highSize (priors on first demographic move time and size)\n");
	fprintf(stderr,"\t -Pe2 lowTime highTime lowSize highSize (priors on second demographic move time and size)\n");
	//fprintf(stderr,"\t -U time (only record mutations back to specified time)\n");
	fprintf(stderr,"\t -R rhhRate (recurrent hitch hiking mode at the locus; rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -L rhhRate (recurrent hitch hiking mode to the side of locus; leftRho is ~Unif(0,4Ns); rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -h (hide selected SNP in partial sweep mode)\n");
	fprintf(stderr,"\t -ts filename (tree sequence output mode - saves to tskit file)\n");
	fprintf(stderr,"\t -F (full ARG mode - record all edges including recombination; use with -ts)\n");
	fprintf(stderr,"\t -d seed1 seed2 (set random number generator seeds)\n");
	
	exit(1);
}


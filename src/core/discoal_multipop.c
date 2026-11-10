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
#include <limits.h>
#include <errno.h>
#include "ranlib.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "shapes.h"
#include "alleleTraj.h"
#include "tskitInterface.h"
#include <tskit/tables.h>
#include "version.h"
#include "segmentPool.h"
#include "demesInterface.h"
#include "configInterface.h"



int locusNumber;
int leftRhoFlag=0;
int recurrentLeftOfLocus=0;
double *currentSize;
long seed1, seed2;
double nextTime, currentFreq;
//float *currentTrajectory;

void getParameters(int argc,const char **argv);
void usage();
void cleanup_and_exit(int sig);
void initializeTrajectoryVariables(double N);

// Remove the sweep trajectory file, if any. Registered with atexit, so
// every exit path (error exits included) removes it.
static void removeTrajectoryFile(void) {
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
	}
}

// Signal handler for cleanup
void cleanup_and_exit(int sig) {
	removeTrajectoryFile();
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
	
	getParameters(argc,argv);
	atexit(removeTrajectoryFile);
	/* Snapshot for per-replicate restore in initialize(); see currentSizeConst. */
	for (i = 0; i < MAXPOPS; i++) currentSizeConst[i] = currentSize[i];
	double N = EFFECTIVE_POPN_SIZE; // effective population size
	// fprintf(stderr, "DEBUG: About to call setall() with seeds: %ld, %ld\n", seed1, seed2);
	setall(seed1, seed2 );
	
	// Store command line for tskit provenance after parsing (now we know if -ts was used)
	if (tskitOutputMode) {
		tskit_store_command_line(argc, argv);
	}
	
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

	// Initialize segment pool with size based on parameters
	size_t initial_pool_size = calculateInitialPoolSize(sampleSize, rho, nSites);
	initSegmentPool(initial_pool_size);

	// fprintf(stderr, "Starting simulation with %d replicates\n", sampleNumber);
	while(i < sampleNumber){
		// fprintf(stderr, "Replicate %d/%d\n", i+1, sampleNumber);
		// fprintf(stderr, "Loop start: i=%d, sampleNumber=%d, continuing=%d\n", i, sampleNumber, i < sampleNumber);
		// fprintf(stderr, "REP %d\n", i+1);
		// fprintf(stderr, "DEBUG: Starting replicate %d of %d\n", i+1, sampleNumber);
		// Clean up from previous replicate if not the first one
		if (i > 0) {
			// fprintf(stderr, "DEBUG: Cleaning up after replicate %d, nodesCapacity=%d\n", i, nodesCapacity);
			// fprintf(stderr, "Starting cleanup for rep %d\n", i+1);
			// fprintf(stderr, "Starting cleanup for rep %d, totNodeNumber=%d\n", i+1, totNodeNumber);
			// fprintf(stderr, "Debug: Cleaning up %d total nodes from previous replicate\n", totNodeNumber);
			// fprintf(stderr, "Starting node cleanup, totNodeNumber=%d\n", totNodeNumber);
			// Clean all nodes in the array, not just up to totNodeNumber
			// This ensures we don't leave any nodes from previous replicates
			int non_null_count = 0;
			for (int cleanup_i = 0; cleanup_i < nodesCapacity; cleanup_i++) {
				if (nodes[cleanup_i] != NULL) non_null_count++;
				// if (cleanup_i % 500 == 0) {
				// 	fprintf(stderr, "Cleaning node %d/%d\n", cleanup_i, cleanup_limit);
				// }
				if (nodes[cleanup_i] != NULL) {
					// Free ancestry segments if they exist
					if (nodes[cleanup_i]->ancestryRoot) {
						// if (cleanup_i >= 1000 && cleanup_i < 1010) {
						// 	fprintf(stderr, "  Freeing ancestry tree for node %d\n", cleanup_i);
						// }
						freeSegmentTree(nodes[cleanup_i]->ancestryRoot);
						nodes[cleanup_i]->ancestryRoot = NULL;
					}
					// Free the node itself
					free(nodes[cleanup_i]);
					nodes[cleanup_i] = NULL;
				}
			}
			// fprintf(stderr, "DEBUG: Found %d non-null nodes to clean\n", non_null_count);
			// fprintf(stderr, "Node cleanup completed\n");
			// fprintf(stderr, "Node cleanup done\n");
			
			// Clean up removed nodes list
			cleanupRemovedNodes();
			
			// Fast O(1) reset of segment pool between replicates
			resetSegmentPool();
			
			cleanupBreakPoints();
			// fprintf(stderr, "BreakPoints cleanup done\n");
			cleanupPopLists();
			// fprintf(stderr, "PopLists cleanup done\n");
			cleanupSampleNodeIds();  // Clean up sample node IDs between replicates
			// fprintf(stderr, "SampleNodeIds cleanup done\n");
			// fprintf(stderr, "Cleanup done for rep %d\n", i+1);
		}
		
		// fprintf(stderr, "About to initialize trajectory variables\n");
		// Initialize all trajectory-related variables
		initializeTrajectoryVariables(N);
		
		// Always initialize tskit for each replicate (must be before initialize())
		// Debug output commented out for production
		// if (tskitOutputMode && i == 0) {
		//	fprintf(stderr, "Initializing tskit with sequence length %d (replicate %d)\n", nSites, i+1);
		// }
		// fprintf(stderr, "DEBUG: About to call tskit_initialize\n");
		int ret = tskit_initialize((double)nSites);
		if (ret != 0) {
			fprintf(stderr, "Error initializing tskit: %s\n", tsk_strerror(ret));
			exit(1);
		}
		// fprintf(stderr, "DEBUG: tskit_initialize returned %d\n", ret);
		
		// fprintf(stderr, "Calling initialize() for rep %d\n", i+1);
		initialize();

		j=0;
		activeSweepFlag = 0;
		currentTime = 0;  // Reset current time for new replicate
		// fprintf(stderr, "Before event loop: i=%d, eventNumber=%d, alleleNumber=%d\n", i, eventNumber, alleleNumber);
		// fprintf(stderr, "DEBUG: Before event loop, eventNumber=%d, alleleNumber=%d\n", eventNumber, alleleNumber);
		for(j=0;j<eventNumber && alleleNumber > 1;j++){
			currentEventNumber=j; //need this annoying global for trajectory generation
			if(j == eventNumber - 1){
				/* No event follows, so the last interval has no end. */
				nextTime = INFINITY;
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
				popShape[events[j].popID].type = SHAPE_CONSTANT;
				popShape[events[j].popID].anchor_value = events[j].popnSize;
				popShape[events[j].popID].rate_param = 0.0;
				popShape[events[j].popID].anchor_time = events[j].time;
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
				case 'g':
				currentTime = events[j].time;
				popShape[events[j].popID].type = SHAPE_EXPONENTIAL;
				popShape[events[j].popID].anchor_value = sizeAt(events[j].popID, currentTime);
				popShape[events[j].popID].rate_param = events[j].popnSize;  /* alpha is stored in popnSize field */
				popShape[events[j].popID].anchor_time = currentTime;
				/* Now run the inter-event interval the same way 'n' does. */
				if(activeSweepFlag == 0){
					if(recurSweepMode == 0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode == 0){
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
				case 'l':
				currentTime = events[j].time;
				popShape[events[j].popID].type = SHAPE_LINEAR;
				popShape[events[j].popID].anchor_value = sizeAt(events[j].popID, currentTime);
				popShape[events[j].popID].rate_param = events[j].popnSize;
				popShape[events[j].popID].anchor_time = currentTime;
				/* Run the inter-event interval as 'n' / 'g' does. */
				if(activeSweepFlag == 0){
					if(recurSweepMode == 0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode == 0){
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
				case 'm':
				currentTime = events[j].time;
				migShape[events[j].popID2][events[j].popID].type = SHAPE_CONSTANT;
				migShape[events[j].popID2][events[j].popID].anchor_value = events[j].popnSize;
				migShape[events[j].popID2][events[j].popID].rate_param = 0.0;
				migShape[events[j].popID2][events[j].popID].anchor_time = currentTime;
				/* Run the inter-event interval as 'n' / 'g' / 'l' do. */
				if(activeSweepFlag == 0){
					if(recurSweepMode == 0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode == 0){
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
				case 's':
				assert(activeSweepFlag == 0);
				currentTime = events[j].time;
				currentFreq = sweepStartFreq(currentTime);
			//	printf("event%d currentTime: %f nextTime: %f popnSize: %f\n",j,currentTime,nextTime,currentSize);

				sampleSweepTrajectory(currentEventNumber, currentSize, sweepMode,
					currentFreq, &currentFreq, alpha, f0, currentTime);
				
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
			currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, INFINITY, currentSize);
		}
		/* Never write a genealogy without a common ancestor at every site.
		 * The last interval has no end, so this happens only if no event can
		 * occur: for example a population that grows without limit going
		 * back in time, so that its lineages may never coalesce. */
		if (activeSites > 0) {
			fprintf(stderr, "Error: in replicate %d, %d sites have no common "
			        "ancestor and no further event can occur; the lineages "
			        "cannot coalesce under this model.\n", i + 1, activeSites);
			exit(1);
		}
		
		// Sort and simplify the tree sequence after ancestry simulation is complete
		// This removes all non-ancestral nodes and edges before mutation placement
		if (tsk_tables != NULL) {
			/* Edges are buffered and written at the next node or flush.
			 * The last events of a replicate can come from any sampler
			 * (neutral, sweep or recurrent sweep); write their edges now,
			 * or simplify drops the nodes they need. */
			if (tskit_flush_edges() != 0) {
				fprintf(stderr, "Error: failed to write buffered tree sequence edges\n");
				exit(1);
			}
			// fprintf(stderr, "Debug: Sorting tskit tables\n");
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
			// fprintf(stderr, "Debug: Starting simplification with %d samples\n", sample_node_count);
			tsk_id_t *node_map = malloc(tsk_tables->nodes.num_rows * sizeof(tsk_id_t));
			if (node_map == NULL) {
				fprintf(stderr, "Error: Failed to allocate memory for node map\n");
				free(samples);
				exit(1);
			}
			
			// In full ARG mode, keep unary nodes (recombination nodes)
			tsk_flags_t simplify_flags = minimalTreeSeq ? 0 : TSK_SIMPLIFY_KEEP_UNARY;
			ret = tsk_table_collection_simplify(tsk_tables, samples, sample_node_count, simplify_flags, node_map);
			if (ret != 0) {
				fprintf(stderr, "Error simplifying tables after ancestry simulation: %s\n", tsk_strerror(ret));
				free(samples);
				free(node_map);
				exit(1);
			}
			// fprintf(stderr, "Debug: Simplification completed successfully\n");
			
			// Update the sample node IDs to reflect the new node numbers after simplification
			for (int i = 0; i < sample_node_count; i++) {
				if (samples[i] != TSK_NULL && node_map[samples[i]] != TSK_NULL) {
					sample_node_ids[i] = node_map[samples[i]];
				}
			}

			free(samples);
			free(node_map);
		}
		
		// fprintf(stderr, "After event loop: i=%d, alleleNumber=%d\n", i, alleleNumber);
		// fprintf(stderr, "DEBUG: After event loop, alleleNumber=%d\n", alleleNumber);
		
		//assign root
	//	root = nodes[0];
		// Place mutations using tskit edge-based algorithm
		extern double theta;
		if (tskit_place_mutations_edge_based(theta) < 0) {
			fprintf(stderr, "Error: Failed to place mutations in tskit\n");
			exit(1);
		}

		// Hudson style output
		//errorCheckMutations();
		// Only generate MS output if not using tskit-only output mode
		if (!tskitOutputMode) {
			// fprintf(stderr, "About to call makeGametesMS for rep %d\n", i+1);
			makeGametesMS(argc,argv);
		}
		//printf("rep: %d\n",i);
		// fprintf(stderr, "DEBUG: Incrementing i from %d to %d\n", i, i+1);
		i++;
		// fprintf(stderr, "End of rep %d, i=%d, sampleNumber=%d\n", i, i, sampleNumber);
		
		// Debug: Print memory statistics (commented out for production)
		// if (tskitOutputMode || minimalTreeSeq) {
		// 	fprintf(stderr, "Replicate %d: Total nodes created: %d, Nodes freed during sim: %d, Net nodes: %d\n", 
		// 	        i, totNodeNumber, freedNodeCount, totNodeNumber - freedNodeCount);
		// }
		
		// Cleanup is now done at the beginning of each replicate to prevent overwrites
		// No cleanup needed here between replicates
		
		// Clean up trajectory after each simulation
		if (trajectoryFd != -1) {
			// fprintf(stderr, "Cleaning up trajectory\n");
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
		
		// fprintf(stderr, "After tskit_cleanup: i=%d\n", i);
		// fprintf(stderr, "Debug: Completed replicate %d (i=%d, totalSimCount=%d, sampleNumber=%d)\n", 
		//         totalSimCount + 1, i, totalSimCount, sampleNumber);
                totalSimCount += 1;
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
	
	// Final cleanup after all replicates
	// fprintf(stderr, "Debug: Final cleanup of nodes\n");
	// int nodes_freed = 0;
	// Clean all nodes in the array
	for (int cleanup_i = 0; cleanup_i < nodesCapacity; cleanup_i++) {
		if (nodes[cleanup_i] != NULL) {
			// nodes_freed++;
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
	// fprintf(stderr, "Debug: Freed %d nodes in final cleanup\n", nodes_freed);
	cleanupBreakPoints();
	cleanupPopLists();
	
	// Free the node arrays themselves (but not the nodes, we already did that)
	if (nodes != NULL) {
		free(nodes);
		nodes = NULL;
	}
	cleanupSampleNodeIds();
	
	// Debug: print freed node count
	// if (freedNodeCount > 0) {
	//	fprintf(stderr, "Debug: Freed %d nodes during simulation\n", freedNodeCount);
	// }
	
	// Clean up segment pool
	destroySegmentPool();
	
	return(0);
}

// Helper function to safely get next argument with bounds checking
static const char* getNextArg(int argc, const char **argv, int *args, const char *option) {
	(*args)++;
	if (*args >= argc) {
		fprintf(stderr, "Error: Option %s requires an argument\n", option);
		usage();
		exit(1);
	}
	return argv[*args];
}

// Helper function to safely parse integer argument
static int parseIntArg(int argc, const char **argv, int *args, const char *option) {
	const char *arg = getNextArg(argc, argv, args, option);
	char *endptr;
	long val = strtol(arg, &endptr, 10);
	
	if (*endptr != '\0' || val > INT_MAX || val < INT_MIN) {
		fprintf(stderr, "Error: Invalid integer argument for %s: '%s'\n", option, arg);
		exit(1);
	}
	
	return (int)val;
}

// Helper function to safely parse double argument
static double parseDoubleArg(int argc, const char **argv, int *args, const char *option) {
	const char *arg = getNextArg(argc, argv, args, option);
	char *endptr;
	double val = strtod(arg, &endptr);
	
	if (*endptr != '\0') {
		fprintf(stderr, "Error: Invalid numeric argument for %s: '%s'\n", option, arg);
		exit(1);
	}
	
	return val;
}


void getParameters(int argc,const char **argv){
	int args;
	int i,j;
	double migR;
	int selCheck,nChangeCheck;
	int yaml_config_loaded = 0;
	
	// Initialize all variables first (before any YAML or command line processing)
	npops = 1;
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
	// fprintf(stderr, "DEBUG: Initial random seeds set: %ld, %ld\n", seed1, seed2);

	EFFECTIVE_POPN_SIZE = 1000000;
	sweepSite = 0.5;
	tDiv=666;
        gammaCoRatioMode = 0;
	priorTheta=priorRho=priorAlpha=priorTau=priorX=priorF0=priorUA=priorC=0;
	numPriorSizeChanges = 0;

	eventFlag = 1;
	outputStyle = 'h';
	mask = 0;
	migFlag = 0;
	deltaTMod = 40;
	recurSweepMode = 0;
	recurrentLeftOfLocus = 0;
	trajFileRows = 0;
	partialSweepMode = 0;
	softSweepMode = 0;
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
	/* calloc, not malloc: the -p handler treats `currentSize[i] == 0.0` as
	 * "not yet set by demes". glibc malloc on a fresh process happens to
	 * return zeroed pages so the sentinel works; ASan's allocator returns
	 * a non-zero poison pattern, which made -p skip the initialization and
	 * left currentSize at garbage (#46). */
	currentSize = calloc(MAXPOPS, sizeof(double));

	// Locate -Y first so we can validate "-Y xor positional args" before
	// doing any work. The positional form "sampleSize numReplicates nSites"
	// and -Y are mutually exclusive: loading YAML and then running the
	// positional parser would feed "-Y" or another flag into atoi() (when -Y
	// comes before the positional slots) or overwrite the per-pop state YAML
	// just installed (when positional slots come first).
	// -v must be known before the YAML file (and any demes file it names)
	// is loaded, because loading prints its diagnostics then.
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) verboseFlag = 1;
	}

	int yaml_arg_idx = -1;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-' && argv[i][1] == 'Y') {
			yaml_arg_idx = i;
			break;
		}
	}
	if (yaml_arg_idx != -1 && argc >= 2 && argv[1][0] != '-') {
		fprintf(stderr,
			"Error: positional arguments (sampleSize numReplicates nSites) "
			"cannot be combined with -Y. Use either '-Y <config.yaml>' "
			"(optionally with flag overrides) or the traditional positional "
			"form without -Y.\n");
		exit(1);
	}

	if (yaml_arg_idx != -1) {
		i = yaml_arg_idx;
		if (i + 1 >= argc || argv[i + 1] == NULL || argv[i + 1][0] == '-' || strlen(argv[i + 1]) == 0) {
			fprintf(stderr, "Error: -Y flag requires a YAML configuration filename argument\n");
			fprintf(stderr, "Usage: %s [options] -Y <config_file.yaml>\n", argv[0]);
			exit(1);
		}
		const char *configFile = argv[i + 1];

		struct discoal_config *config = NULL;
		int ret = load_yaml_config(configFile, &config);
		if (ret != EXIT_SUCCESS) {
			fprintf(stderr, "Error: Failed to load YAML configuration file '%s'\n", configFile);
			exit(1);
		}

		ret = apply_yaml_config(config);
		if (ret != EXIT_SUCCESS) {
			fprintf(stderr, "Error: Failed to apply configuration from '%s'\n", configFile);
			exit(1);
		}

		VERBOSE_LOG("Loaded configuration from YAML file '%s'\n", configFile);
		yaml_config_loaded = 1;
	}

	if (!yaml_config_loaded && argc < 4) {
		usage();
	}

	if (!yaml_config_loaded) {
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

		popnSizes[0]=sampleSize;
		popnSizes[1]=0;
		sampleSizes[0]=sampleSize;
		sampleSizes[1]=0;
		effectiveSampleSize = sampleSize;
	} else {
		// sampleSize/sampleSizes[]/popnSizes[] already set via apply_yaml_config;
		// skip positional parsing and start the flag loop at argv[1].
		args = 1;
		if (sampleSize > 0) {
			effectiveSampleSize = sampleSize;
		}
	}

	while(args < argc){
		// Skip YAML config arguments if already processed
		if (argv[args][0] == '-' && argv[args][1] == 'Y') {
			args += 2;  // Skip -Y and its filename argument
			continue;
		}
		
		// Check if argument starts with '-'
		if (argv[args][0] != '-') {
			fprintf(stderr, "Error: Unexpected argument '%s'\n", argv[args]);
			fprintf(stderr, "All options must start with '-'. Try '%s' without arguments for usage.\n", argv[0]);
			exit(1);
		}
		
		// Check for empty option after '-'
		if (argv[args][1] == '\0') {
			fprintf(stderr, "Error: Empty option '-'\n");
			exit(1);
		}

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
			if (setSweepTrajectoryFile(getNextArg(argc, argv, &args, "-S")) != 0)
				exit(1);
			break;
			case 's' :
			/* ms's -s is not supported; reject it, not ignore it. */
			fprintf(stderr, "Error: -s (a fixed number of segregating sites) "
				"is not supported; place mutations with tskit on tree "
				"sequence output (-ts) instead\n");
			exit(1);
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
				theta = parseDoubleArg(argc, argv, &args, "-t");
			}
			break;
			case 'i' :
			deltaTMod = parseDoubleArg(argc, argv, &args, "-i");
			break;
			case 'r' :
			rho = parseDoubleArg(argc, argv, &args, "-r");
			break;
			case 'g' :
                        if (argv[args][2] == 'r')
                        {
			  gammaCoRatio = parseDoubleArg(argc, argv, &args, "-gr");
			  gcMean = parseIntArg(argc, argv, &args, "-gr");
                          gammaCoRatioMode = 1;
                        }
                        else
                        {
			  my_gamma = parseDoubleArg(argc, argv, &args, "-g");
			  gcMean = parseIntArg(argc, argv, &args, "-g");
                        }
			break;
			case 'a' :
			alpha = parseDoubleArg(argc, argv, &args, "-a");
			break;
			case 'D' :
			fprintf(stderr,
				"Error: the -D command-line flag for loading demes files is no "
				"longer supported. Use a YAML config (-Y) with "
				"demography.demes_filename instead. See "
				"config_examples/demes_example.yaml.\n");
			exit(1);
			break;
			case 'x' :
			sweepSite = parseDoubleArg(argc, argv, &args, "-x");
			break;
			case 'M' :
			if(npops==1){
				fprintf(stderr,"Error: attempting to set migration but only one population! Be sure that 'm' flags are specified after 'p' flag\n");
				exit(1);
			}
			migR = parseDoubleArg(argc, argv, &args, "-M");
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
				i = parseIntArg(argc, argv, &args, "-m");
				j = parseIntArg(argc, argv, &args, "-m");
				migR = parseDoubleArg(argc, argv, &args, "-m");
				migMatConst[i][j]=migR;
				migFlag = 1;
			break;
			case 'p' :
			npops = parseIntArg(argc, argv, &args, "-p");
			if(npops > MAXPOPS){
				fprintf(stderr,"Error: too many populations defined. Current maximum number = %d. Change MAXPOPS define in discoal.h and recompile... if you dare\n",MAXPOPS);
				exit(1);
			}
			for(i=0;i<npops;i++){
				sampleSizes[i] = parseIntArg(argc, argv, &args, "-p");
				// Don't overwrite currentSize if it was already set (e.g., by demes)
				// Only set to 1.0 if it's still 0 (uninitialized)
				if (currentSize[i] == 0.0) {
					currentSize[i] = 1.0;
				}
			}
			
			break;
			case 'e' :
				switch(argv[args][2]){
					case 'n':
					ensureEventsCapacity();
					events[eventNumber].time = userTimeToInternal(atof(argv[++args]));
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popnSize = atof(argv[++args]);
					events[eventNumber].type = 'n'; //size change
					eventNumber++;
					break;
					case 'd' :
					case 'j' :  // -ej is equivalent to -ed (for ms compatibility)
					tDiv =  atof(argv[++args]);
						ensureEventsCapacity();
					events[eventNumber].time = userTimeToInternal(tDiv);
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].type = 'p'; //pop split
					eventNumber++;
					break;
					case 'a' :
						ensureEventsCapacity();
					events[eventNumber].time = userTimeToInternal(atof(argv[++args]));
					events[eventNumber].popID = atoi(argv[++args]);
					events[eventNumber].popID2 = atoi(argv[++args]);
					events[eventNumber].popID3 = atoi(argv[++args]);
					events[eventNumber].admixProp = atof(argv[++args]);
					events[eventNumber].type = 'a'; //admix split
					eventNumber++;
					break;
					/* -eg/-eG alpha and -el/-eL gamma are rates per 4N0
					 * generations, as in ms; userRateToInternal converts them
					 * like the times. */
					case 'g':
						/* -eg time popID alpha */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(atof(argv[++args]));
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = userRateToInternal(atof(argv[++args]));  /* alpha stored in popnSize field */
						events[eventNumber].type = 'g';
						eventNumber++;
						break;
					case 'G':
						/* -eG time alpha -- applies to all populations */
						{
							double t = userTimeToInternal(atof(argv[++args]));
							double alpha_val = userRateToInternal(atof(argv[++args]));
							for (int p = 0; p < npops; p++) {
								ensureEventsCapacity();
								events[eventNumber].time = t;
								events[eventNumber].popID = p;
								events[eventNumber].popnSize = alpha_val;
								events[eventNumber].type = 'g';
								eventNumber++;
							}
						}
						break;
					case 'l':
						/* -el time popID gamma: gamma is in relative size (size / N0) */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(atof(argv[++args]));
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = userRateToInternal(atof(argv[++args]));  /* gamma stored in popnSize field */
						events[eventNumber].type = 'l';
						eventNumber++;
						break;
					case 'L':
						/* -eL time gamma — applies to all populations */
						{
							double t = userTimeToInternal(atof(argv[++args]));
							double gamma_val = userRateToInternal(atof(argv[++args]));
							for (int p = 0; p < npops; p++) {
								ensureEventsCapacity();
								events[eventNumber].time = t;
								events[eventNumber].popID = p;
								events[eventNumber].popnSize = gamma_val;
								events[eventNumber].type = 'l';
								eventNumber++;
							}
						}
						break;
					case 'm':
						/* -em time srcPop dstPop rate — change migration rate at this time */
						{
							ensureEventsCapacity();
							events[eventNumber].time = userTimeToInternal(atof(argv[++args]));
							events[eventNumber].popID2 = atoi(argv[++args]);  /* source pop */
							events[eventNumber].popID = atoi(argv[++args]);   /* destination pop */
							events[eventNumber].popnSize = atof(argv[++args]); /* new rate stored in popnSize field */
							events[eventNumber].type = 'm';  /* lowercase 'm' for time-varying migration event */
							/* Migration that starts later still counts as
							 * migration for initialize()'s check. */
							if (events[eventNumber].popnSize > 0) migFlag = 1;
							eventNumber++;
						}
						break;
					case 'M':
						/* -eM time rate — set all off-diagonal pairs to the given rate at this time */
						{
							double t = userTimeToInternal(atof(argv[++args]));
							double rate_val = atof(argv[++args]);
							for (int src = 0; src < npops; src++) {
								for (int dst = 0; dst < npops; dst++) {
									if (src == dst) continue;
									ensureEventsCapacity();
									events[eventNumber].time = t;
									events[eventNumber].popID2 = src;
									events[eventNumber].popID = dst;
									events[eventNumber].popnSize = rate_val;
									events[eventNumber].type = 'm';
									eventNumber++;
								}
							}
							/* As for -em: migration that starts later still
							 * counts for initialize()'s check. */
							if (rate_val > 0) migFlag = 1;
						}
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
				tau = userTimeToInternal(atof(argv[++args]));
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
				tau = userTimeToInternal(atof(argv[++args]));
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
					  pTauLow=userTimeToInternal(atof(argv[++args]));
					  pTauUp=userTimeToInternal(atof(argv[++args]));
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
				/* -Pe1 / -Pe2 lowTime highTime lowSize highSize: a size
				 * change of population 0 with uniform priors on time and
				 * size. Both work the same way; the digit only names them. */
				if (argv[args][3] == '1' || argv[args][3] == '2') {
					const char *opt = argv[args];
					double tLow = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					double tHigh = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					double sLow = parseDoubleArg(argc, argv, &args, opt);
					double sHigh = parseDoubleArg(argc, argv, &args, opt);
					if (addPriorSizeChange(opt, tLow, tHigh, sLow, sHigh, 0) != 0)
						exit(1);
				}
				break;
			}
                        break;
			case 'd' :
			// fprintf(stderr, "DEBUG: Processing -d flag, overriding seeds\n");
			seed1 = parseIntArg(argc, argv, &args, "-d");
			seed2 = parseIntArg(argc, argv, &args, "-d");
			// fprintf(stderr, "DEBUG: Seeds set by -d flag: %ld, %ld\n", seed1, seed2);
			break;
			case 'N' :
			EFFECTIVE_POPN_SIZE = parseIntArg(argc, argv, &args, "-N");
			break;
			case 'R' :
			recurSweepMode = 1;
			sweepMode = 's';
			recurSweepRate = atof(argv[++args]);
			break;
			case 'L' :
			recurSweepMode = 1;
			recurrentLeftOfLocus = 1;
			sweepMode = 's';
			sweepSite = -1.0;
			recurSweepRate = atof(argv[++args]);
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
			case 'v' :
			// verboseFlag is set before argument parsing; see above.
			break;
			case 'h' :
			hidePartialSNP = 1;
			break;
			case 'A' :
				ensureEventsCapacity();
			events[eventNumber].lineageNumber = atoi(argv[++args]);
			events[eventNumber].popID = atoi(argv[++args]);	 
			events[eventNumber].time = userTimeToInternal(atof(argv[++args])); 
			events[eventNumber].type = 'A'; //ancient sample
			eventNumber++;
			break;
			case 'T' :
			fprintf(stderr, "Error: -T tree output mode has been removed.\n");
			fprintf(stderr, "Use -ts <filename.trees> for tree sequence output instead.\n");
			exit(1);
			break;
			default:
			fprintf(stderr, "Error: Unknown option '-%c'\n", argv[args][1]);
			fprintf(stderr, "Try '%s' without arguments for usage information.\n", argv[0]);
			exit(1);
			break;
		}
		args++;
	}
	// fprintf(stderr, "DEBUG: Final seeds at end of getParameters: %ld, %ld\n", seed1, seed2);
	// fprintf(stderr, "DEBUG: Final theta at end of getParameters: %f\n", theta);
	sortEventArray(events,eventNumber);

	if (validateShapeTrajectories(events, eventNumber) != 0) {
		exit(1);
	}

	/* Priors (the -P flags, or YAML priors) follow the ranges of their YAML
	 * fields. A soft sweep that is partial needs every drawn f0 below every
	 * drawn final frequency; -Pf counts as a soft sweep. */
	{
		if(checkUniformPriorBounds() != 0)
			exit(1);
		/* -Pre mean upper: an exponential prior on rho, capped at upper. */
		if(priorRho == 2 && !(pRhoMean > 0.0 && pRhoUp > 0.0)){
			fprintf(stderr, "Error: -Pre %g %g: the mean and the upper bound must "
				"be > 0\n", pRhoMean, pRhoUp);
			exit(1);
		}
	}
	/* Fixed values, with the rules of their YAML fields. A negative rate
	 * could make a run hang (-r, -g, -R). A sweep outside the locus (-ls,
	 * -ld, -ln, -L) has no position in it (sweepSite -1). A value that a prior
	 * replaces (-t with -Pt, and so on) is not used, so it is not checked.
	 * Each test is written so that NaN fails it. */
	{
		/* A sweep time must be >= 0; with -Pu the time is replaced by a draw,
		 * but a NaN time still breaks the sorting of the events, so NaN is
		 * rejected either way. */
		const int sweepOutsideLocus = leftRhoFlag || recurrentLeftOfLocus;
		double tauBadValue = 0.0, admixBad = 0.0;
		int tauBad = 0;
		int admixOut = 0, sampleSum = 0, sampleNegative = 0, sampleBad = 0;
		for(i = 0; i < eventNumber; i++){
			if(events[i].type == 's' && (isnan(events[i].time) ||
			   (!priorTau && events[i].time < 0.0))){
				tauBad = 1;
				tauBadValue = events[i].time;
			}
			if(events[i].type == 'a' && !(events[i].admixProp >= 0.0 && events[i].admixProp <= 1.0)){
				admixOut = 1;
				admixBad = events[i].admixProp;
			}
		}
		for(i = 0; i < npops; i++){
			sampleSum += sampleSizes[i];
			if(sampleSizes[i] < 0){
				sampleNegative = 1;
				sampleBad = sampleSizes[i];
			}
		}
		const struct { int bad; const char *what; double value; const char *rule; } fixed[] = {
			{!priorTheta && !(theta >= 0.0), "-t", theta, ">= 0"},
			{!priorRho && !(rho >= 0.0), "-r", rho, ">= 0"},
			{!(my_gamma >= 0.0), "-g (rate)", my_gamma, ">= 0"},
			{!(gammaCoRatio >= 0.0), "-gr (ratio)", gammaCoRatio, ">= 0"},
			{gcMean < 0, "-g or -gr (tract length)", (double)gcMean, ">= 0"},
			{!priorAlpha && !(alpha >= 0.0), "-a", alpha, ">= 0"},
			{!priorUA && !(uA >= 0.0), "-uA", uA, ">= 0"},
			{!(leftRho >= 0.0), "-ls, -ld or -ln (leftRho)", leftRho / 2.0, ">= 0"},
			{!priorX && !sweepOutsideLocus && !(sweepSite >= 0.0 && sweepSite <= 1.0),
				"-x", sweepSite, "in [0, 1]"},
			{!priorF0 && !(f0 >= 0.0 && f0 < 1.0), "-f", f0, "in [0, 1)"},
			{tauBad, "a sweep time (-ws, -wd, -wn, -ls, -ld, -ln)",
				tauBadValue / USER_TO_INTERNAL_TIME, ">= 0"},
			{!(deltaTMod > 0.0), "-i", deltaTMod, "> 0"},
			{EFFECTIVE_POPN_SIZE < 1, "-N", (double)EFFECTIVE_POPN_SIZE, ">= 1"},
			{recurSweepMode && !(recurSweepRate > 0.0), "-R or -L (rate)",
				recurSweepRate, "> 0"},
			{admixOut, "-ea (proportion)", admixBad, "in [0, 1]"},
			{sampleNegative, "a -p sample size", (double)sampleBad, ">= 0"},
			{sampleSum != sampleSize, "the sum of the -p sample sizes",
				(double)sampleSum, "equal to the sample size"},
		};
		for(size_t k = 0; k < sizeof(fixed) / sizeof(fixed[0]); k++){
			if(fixed[k].bad){
				fprintf(stderr, "Error: %s (%g) must be %s\n", fixed[k].what,
					fixed[k].value, fixed[k].rule);
				exit(1);
			}
		}
	}
	/* Ancient samples are taken from a population's sample count, so the
	 * ancient samples of one population must not outnumber its samples.
	 * initialize() relies on this when it makes sample nodes ancient. A
	 * negative time would put a sample after the present. */
	{
		int ancientPerPop[MAXPOPS] = {0};
		for (i = 0; i < eventNumber; i++) {
			if (events[i].type != 'A') continue;
			int p = events[i].popID;
			if (p < 0 || p >= npops || events[i].lineageNumber < 1) {
				fprintf(stderr, "Error: ancient sample of %d lineages from "
					"population %d; the count must be >= 1 and the population "
					"between 0 and %d\n", events[i].lineageNumber, p, npops - 1);
				exit(1);
			}
			if (events[i].time < 0) {
				fprintf(stderr, "Error: ancient sample of %d lineages from "
					"population %d at time %g; the time must be >= 0\n",
					events[i].lineageNumber, p,
					internalTimeToUser(events[i].time));
				exit(1);
			}
			ancientPerPop[p] += events[i].lineageNumber;
			if (ancientPerPop[p] > sampleSizes[p]) {
				fprintf(stderr, "Error: population %d has %d ancient samples "
					"but only %d samples in all; ancient samples are part of "
					"the population's sample count\n",
					p, ancientPerPop[p], sampleSizes[p]);
				exit(1);
			}
		}
	}

	//make sure events are kosher
	selCheck = 0;
	nChangeCheck=0;
	int sweepEvents = 0;
	for(i=0;i<eventNumber;i++){
		//printf("event %d: type is %c\n", i, events[i].type);
	 		if(events[i].type == 's'){
				selCheck = 1;
				sweepEvents++;
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
	}
	const int softSweep = softSweepMode || priorF0;
	if(softSweep && partialSweepMode && trajFileRows == 0){
		const double maxF0 = priorF0 ? pF0Up : f0;
		const double minC = priorC ? pCLow : partialSweepFinalFreq;
		if(maxF0 >= minC){
			fprintf(stderr, "Error: the initial frequency of a soft sweep (-f, "
				"-Pf, or YAML initial_frequency; %g) must be below the final "
				"frequency of a partial sweep (-c, -Pc, or YAML final_frequency; "
				"%g), for every value of their priors\n", maxF0, minC);
			exit(1);
		}
	}
	/* Priors that contradict the sweep model, from the CLI or YAML: a
	 * drawn position puts the sweep inside the locus, and a drawn tau
	 * needs a single sweep. */
	/* A sweep outside the locus is set by its options, not by the value
	 * of sweepSite: with -Px, -x is a placeholder that may be anything. */
	if(priorX && (leftRhoFlag || recurrentLeftOfLocus)){
		fprintf(stderr, "Error: a prior on the sweep position (-Px, or YAML "
			"priors.sweep_position) cannot be combined with a sweep outside "
			"the locus (-ls, -ld, -ln, -L)\n");
		exit(1);
	}
	if(priorTau && recurSweepMode){
		fprintf(stderr, "Error: a prior on the sweep time (-Pu, or YAML "
			"priors.fixation_time_ago) needs a single sweep, not recurrent "
			"sweeps (-R, -L)\n");
		exit(1);
	}
	if((leftRhoFlag || recurrentLeftOfLocus) && !priorX && sweepSite >= 0.0){
		fprintf(stderr, "Error: -x (%g) cannot be combined with a sweep outside "
			"the locus (-ls, -ld, -ln, -L)\n", sweepSite);
		exit(1);
	}
	if(softSweep && recurSweepMode == 1){
		fprintf(stderr, "Error: recurrent soft sweeps (-R or -L with -f or "
			"-Pf) are not implemented\n");
		exit(1);
	}
	/* Sweep trajectories are not modeled under growth: proposeTrajectory()
	 * accepts a trajectory by a ratio of sizes that follows only size
	 * changes (-en). So population 0 may not grow or shrink continuously
	 * (-eg, -el, -eG, -eL, or an exponential epoch of a demes model) at or
	 * after the earliest start of a sweep, going back in time. Growth in
	 * other populations, and growth of population 0 that surely ends before
	 * the sweep, are allowed. Recurrent sweeps draw coalescence in every
	 * population at constant rates between sweeps, so with them no
	 * population may grow at any time. A trajectory read from a file (-S)
	 * does not use that ratio, and its sampler follows the size at each
	 * step, so -S is not restricted. The events are sorted here, with ties
	 * in command-line order. */
	if((selCheck && trajFileRows == 0) || recurSweepMode){
		int growthAtSweep = 0;
		int sweepIndex = -1;
		double sweepStart = 0.0;
		for(i = 0; i < eventNumber; i++){
			if(events[i].type == 's' && sweepIndex < 0){
				sweepIndex = i;  /* the earliest sweep */
				sweepStart = events[i].time;
			}
			if(recurSweepMode && (events[i].type == 'g' || events[i].type == 'l') &&
			   events[i].popnSize != 0.0)
				growthAtSweep = 1;
		}
		if(priorTau) sweepStart = pTauLow;
		for(i = 0; !recurSweepMode && i < eventNumber; i++){
			if(events[i].popID != 0 || (events[i].type != 'g' && events[i].type != 'l') ||
			   events[i].popnSize == 0.0)
				continue;
			/* Growth of population 0 from events[i]: it must start before
			 * the sweep and surely end before it, with a later size change
			 * of population 0 (an -en, or growth at rate 0). */
			const double tg = events[i].time;
			if(tg > sweepStart || (tg == sweepStart && !priorTau && i > sweepIndex)){
				growthAtSweep = 1;
				break;
			}
			int ended = 0;
			for(j = i + 1; j < eventNumber && !ended; j++){
				if(events[j].popID != 0) continue;
				const char type = events[j].type;
				if(!(type == 'n' || ((type == 'g' || type == 'l') && events[j].popnSize == 0.0)))
					continue;
				/* A size change with a prior on its time can fall anywhere
				 * in [timeLow, timeHigh]. */
				double lo = events[j].time, hi = events[j].time;
				if(type == 'n' && events[j].priorSlot > 0){
					lo = priorSizeChanges[events[j].priorSlot - 1].timeLow;
					hi = priorSizeChanges[events[j].priorSlot - 1].timeHigh;
				}
				const int afterGrowth = (lo > tg) || (lo == tg && lo == hi);
				const int beforeSweep = priorTau ? (hi <= sweepStart)
					: (hi < sweepStart || (hi == sweepStart && lo == hi && j < sweepIndex));
				ended = afterGrowth && beforeSweep;
			}
			if(!ended){
				growthAtSweep = 1;
				break;
			}
		}
		if(growthAtSweep){
			fprintf(stderr, "Error: population 0 grows or shrinks (-eg, -el, "
				"-eG, -eL, or an exponential demes epoch) at or after the start "
				"of a sweep, or a population grows with recurrent sweeps (-R, -L); "
				"sweeps are not modeled under growth. Use size changes (-en) "
				"instead\n");
			exit(1);
		}
	}
	/* A sweep needs at least one individual in population 0 at its start.
	 * Without growth at the sweep, the size there is one of population 0's
	 * sizes (its present size, an -en, or the low bound of a size prior), so
	 * check them all now rather than during the run. A size of 0 is not a
	 * population (the present size before it is set, or a deme that does
	 * not exist then). */
	if(selCheck || recurSweepMode){
		double smallest = INFINITY;
		if(currentSize[0] > 0.0) smallest = currentSize[0];
		for(i = 0; i < eventNumber; i++){
			if(events[i].popID != 0 || events[i].type != 'n') continue;
			double size = events[i].popnSize;
			if(events[i].priorSlot > 0)
				size = priorSizeChanges[events[i].priorSlot - 1].sizeLow;
			if(size > 0.0 && size < smallest) smallest = size;
		}
		if(isfinite(smallest) && floor(EFFECTIVE_POPN_SIZE * smallest) < 1.0){
			fprintf(stderr, "Error: with a sweep, population 0 must have at "
				"least one individual at every size (N x size = %g; N from -N)\n",
				EFFECTIVE_POPN_SIZE * smallest);
			exit(1);
		}
	}
	/* A partial sweep that ends below one copy of the allele has no
	 * trajectory: the sweep ends when the frequency is 1/(2N). A prior
	 * needs only its high bound above that; a draw below one copy gives a
	 * sweep that ends at once. */
	const double oneCopy = 1.0 / (2.0 * EFFECTIVE_POPN_SIZE);
	if(partialSweepMode && trajFileRows == 0 &&
	   (priorC ? pCUp : partialSweepFinalFreq) <= oneCopy){
		fprintf(stderr, "Error: the final frequency of a partial sweep (-c, YAML "
			"final_frequency, or the high bound of a prior on it) must be above "
			"1/(2N) = %g, one copy of the allele (N from -N)\n", oneCopy);
		exit(1);
	}
	/* -S: the file gives the sweep's time, its frequencies, and so the
	 * strength of selection, so options that set these are rejected. A
	 * first frequency below 1 makes the sweep partial, as -c does: each
	 * lineage carries the allele with that probability, and its carriers
	 * carry the selected SNP. */
	if(trajFileRows > 0){
		if(sweepEvents > 1 || recurSweepMode || leftRhoFlag || softSweepMode ||
		   partialSweepMode || alpha != 0.0 || priorAlpha || priorTau ||
		   priorF0 || priorC){
			fprintf(stderr, "Error: a trajectory file (-S, or YAML sweep_mode "
				"'trajectory_file') gives the sweep's time and frequencies; it "
				"cannot be combined with -ws, -wd, -wn, -ls, -ld, -ln, -R, -L, "
				"-a, -f, -c, -Pa, -Pu, -Pf or -Pc\n");
			exit(1);
		}
		if(trajFilePoints[0].freq < 1.0){
			partialSweepMode = 1;
			partialSweepFinalFreq = trajFilePoints[0].freq;
		}
		buildSweepTrajectoryFromFile();
	}
	
}
		


void usage(){
	fprintf(stderr, "▗▄▄▄ ▗▄▄▄▖ ▗▄▄▖ ▗▄▄▖ ▗▄▖  ▗▄▖ ▗▖   \n");
	fprintf(stderr, "▐▌  █  █  ▐▌   ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌   \n");
	fprintf(stderr, "▐▌  █  █   ▝▀▚▖▐▌   ▐▌ ▐▌▐▛▀▜▌▐▌   \n");
	fprintf(stderr, "▐▙▄▄▀▗▄█▄▖▗▄▄▞▘▝▚▄▄▖▝▚▄▞▘▐▌ ▐▌▐▙▄▄\n\n\n");
	
	fprintf(stderr,"discoal version %s\n", DISCOAL_VERSION);
	fprintf(stderr,"usage: discoal sampleSize numReplicates nSites -ws tau\n");
	fprintf(stderr,"parameters: \n");

	fprintf(stderr,"\t -t theta \n");
	fprintf(stderr,"\t -r rho (=zero if not specified)\n");
	fprintf(stderr,"\t -g conversionRate tractLengthMean (gene conversion)\n");
	fprintf(stderr,"\t -gr conversionToCrossoverRatio tractLengthMean (gene conversion where initiation rate = rho*conversionToCrossoverRatio)\n");
	fprintf(stderr,"\t -p npops sampleSize1 sampleSize2 etc.\n");
	fprintf(stderr,"\t -D demesFile.yaml (load demographic model from demes format file)\n");
	fprintf(stderr,"\t -en time popnID size (changes size of popnID)\n");
	fprintf(stderr,"\t -eg time popnID alpha (sets popnID to exponential growth at rate alpha per 4N0 generations, as in ms; positive alpha = past was smaller)\n");
	fprintf(stderr,"\t -eG time alpha (-eg applied to all populations)\n");
	fprintf(stderr,"\t -el time popnID gamma (sets popnID to linear growth; gamma is the change in relative size per 4N0 generations; positive gamma = past was smaller)\n");
	fprintf(stderr,"\t -eL time gamma (-el applied to all populations)\n");
	fprintf(stderr,"\t -ed time popnID1 popnID2 (joins popnID1 into popnID2; -ej is an alias)\n");
	fprintf(stderr,"\t -ea time daughterPopnID founderPopnID1 founderPopnID2 admixProp (admixture-- back in time daughterPopnID into two founders)\n");
	fprintf(stderr,"\t -em time popnID1 popnID2 migRate (changes migration rate from popnID1 to popnID2 starting at time)\n");
	fprintf(stderr,"\t -eM time migRate (sets all off-diagonal migration rates to migRate starting at time)\n");
	
	fprintf(stderr,"\t -ws tau (sweep happend tau generations ago- stochastic sweep)\n");  
	fprintf(stderr,"\t -wd tau (sweep happend tau generations ago- deterministic sweep)\n"); 
	fprintf(stderr,"\t -wn tau (sweep happend tau generations ago- neutral sweep)\n");
	fprintf(stderr,"\t -ls tau leftRho (stochastic sweep some genetic distance to the left of the simulated window--at 4Nr = 2 x leftRho)\n");
	fprintf(stderr,"\t\t similarly, ld and ln simulate deterministic and neutral sweeps to the left of the window, respectively\n");
	fprintf(stderr,"\t -f first frequency at which selection acts on allele (F0; sweep models only)\n");
	fprintf(stderr,"\t -uA rate at which adaptive mutation recurs during the sweep phase (sweep models only)\n");
	fprintf(stderr,"\t -N sweepEffectivePopnSize (sweep models only)\n");	
	fprintf(stderr,"\t -a alpha (=2Ns)\n");
	fprintf(stderr,"\t -x sweepSite (0-1)\n");
	fprintf(stderr,"\t -c partialSweepFinalFrequency (partial sweeps)\n");	
	fprintf(stderr,"\t -i dt (sweep time increment scalar; default 40 -> 1/40N)\n");
	fprintf(stderr,"\t -S file (sweep trajectory from a file: lines of time (4N0 generations, increasing) and allele frequency, ending at frequency 0)\n");
	
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
	fprintf(stderr,"\t -R rhhRate (recurrent hitch hiking mode at the locus; rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -L rhhRate (recurrent hitch hiking mode to the side of locus; leftRho is ~Unif(0,4Ns); rhh is rate per 2N individuals / generation)\n");
	fprintf(stderr,"\t -h (hide selected SNP in partial sweep mode)\n");
	fprintf(stderr,"\t -v (verbose: print progress and demes conversion details to stderr)\n");
	fprintf(stderr,"\t -ts filename (tree sequence output mode - saves to tskit file)\n");
	fprintf(stderr,"\t -F (full ARG mode - record all edges including recombination; use with -ts)\n");
	fprintf(stderr,"\t -d seed1 seed2 (set random number generator seeds)\n");
	
	exit(1);
}


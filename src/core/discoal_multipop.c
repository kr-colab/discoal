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

/* The value of a macro as a string literal, e.g. for an error message. */
#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
#include "segmentPool.h"
#include "demesInterface.h"
#include "configInterface.h"

int locusNumber;
int leftRhoFlag=0;
int recurrentLeftOfLocus=0;
int sweepPositionGiven=0;
int demesModelLoaded=0;
double *currentSize;
long seed1, seed2;
double nextTime, currentFreq;

void getParameters(int argc,const char **argv);
static void printUsage(FILE *out);
void usage();
void initializeTrajectoryVariables(double N);

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

/* At exit, flush stdout: a write error (a full disk, for example) must
 * not end as success with cut-off output. */
static void checkStdout(void){
	if (fflush(stdout) != 0 || ferror(stdout)) {
		fprintf(stderr, "Error: cannot write the output: %s\n", strerror(errno));
		_exit(1);
	}
}

int main(int argc, const char * argv[]){
	int i,j;
	
	atexit(checkStdout);
	getParameters(argc,argv);
	/* Snapshot for per-replicate restore in initialize(); see currentSizeConst. */
	for (i = 0; i < MAXPOPS; i++) currentSizeConst[i] = currentSize[i];
	double N = EFFECTIVE_POPN_SIZE; // effective population size
	setall(seed1, seed2 );
	
	// Store command line for tskit provenance after parsing (now we know if -ts was used)
	if (tskitOutputMode) {
		tskit_store_command_line(argc, argv);
	}

	//Hudson style header
	for(i=0;i<argc;i++)printf("%s ",argv[i]);
	printf("\n%ld %ld\n", seed1, seed2);
	
	i = 0;
	trajectoryCapacity = TRAJSTEPSTART;
	trajectoryFd = -1;  // Initialize to invalid
	currentTrajectory = NULL;  // Will be mmap'd when needed
	
	// Initialize global trajectory generator (lazy approach)

	// Initialize segment pool with size based on parameters
	size_t initial_pool_size = calculateInitialPoolSize(sampleSize, rho);
	initSegmentPool(initial_pool_size);

	while(i < sampleNumber){
		// Clean up from previous replicate if not the first one
		if (i > 0) {
			// Free every node left in the array by the previous replicate
			for (int cleanup_i = 0; cleanup_i < nodesCapacity; cleanup_i++) {
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
			
			// Fast O(1) reset of segment pool between replicates
			resetSegmentPool();
			
			cleanupPopLists();
			cleanupSampleNodeIds();  // Clean up sample node IDs between replicates
		}
		
		// Initialize all trajectory-related variables
		initializeTrajectoryVariables(N);
		
		// Always initialize tskit for each replicate (must be before initialize())
		int ret = tskit_initialize((double)nSites);
		if (ret != 0) {
			fprintf(stderr, "Error initializing tskit: %s\n", tsk_strerror(ret));
			exit(1);
		}
		
		initialize();

		j=0;
		activeSweepFlag = 0;
		currentTime = 0;  // Reset current time for new replicate
		for(j=0;j<eventNumber && alleleNumber > 1;j++){
			currentEventNumber=j; //need this annoying global for trajectory generation
			if(j == eventNumber - 1){
				/* No event follows, so the last interval has no end. */
				nextTime = INFINITY;
				}
			else{
				nextTime = events[j+1].time;
			}
			switch(events[j].type){
				case 'n':
				currentTime = events[j].time;
				currentSize[events[j].popID] = events[j].popnSize;
				popShape[events[j].popID].type = SHAPE_CONSTANT;
				popShape[events[j].popID].anchor_value = events[j].popnSize;
				popShape[events[j].popID].rate_param = 0.0;
				popShape[events[j].popID].anchor_time = events[j].time;
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
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
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode == 0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
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
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode == 0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
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
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode == 0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
							currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				break;
				case 's':
				assert(activeSweepFlag == 0);
				currentTime = events[j].time;
				currentFreq = sweepStartFreq(currentTime);

				sampleSweepTrajectory(currentEventNumber, sweepMode,
					currentFreq, alpha, f0, currentTime);
				
				currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
				if (currentTime < nextTime)
                                        currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
						
				break;
				case 'p': //merging populations
				currentTime = events[j].time;
				mergePopns(events[j].popID, events[j].popID2);

				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                        		currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				break;
				case 'a':
				currentTime = events[j].time;
				admixPopns(events[j].popID, events[j].popID2, events[j].popID3, events[j].admixProp);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				break;
				case 'A':
				currentTime = events[j].time;
				addAncientSample(events[j].lineageNumber, events[j].popID, events[j].time, activeSweepFlag, currentFreq);
				if(activeSweepFlag == 0){
					if(recurSweepMode ==0){
						currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				else{
					if(recurSweepMode ==0){
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
						 	currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
	                                        	currentTime = neutralPhaseGeneralPopNumber(currentTime, nextTime);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(currentTime, nextTime, sweepSite, \
					 		currentFreq, &currentFreq, &activeSweepFlag, f0, uA);
						if (currentTime < nextTime)
                                               		currentTime = recurrentSweepPhaseGeneralPopNumber(currentTime, nextTime, &currentFreq, alpha, sweepMode);
					}
				}
				break;
			}
			
		}
		//finish up the coalescing action!
		if(alleleNumber > 1){
			currentTime = neutralPhaseGeneralPopNumber(currentTime, INFINITY);
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
			
			// In full ARG mode, keep unary nodes (recombination nodes)
			tsk_flags_t simplify_flags = minimalTreeSeq ? 0 : TSK_SIMPLIFY_KEEP_UNARY;
			ret = tsk_table_collection_simplify(tsk_tables, samples, sample_node_count, simplify_flags, node_map);
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
		// Place mutations using tskit edge-based algorithm
		extern double theta;
		if (tskit_place_mutations_edge_based(theta) < 0) {
			fprintf(stderr, "Error: Failed to place mutations in tskit\n");
			exit(1);
		}

		// Hudson style output
		// Only generate MS output if not using tskit-only output mode
		if (!tskitOutputMode) {
			makeGametesMS();
		}
		i++;

		// Cleanup is now done at the beginning of each replicate to prevent overwrites
		// No cleanup needed here between replicates
		
		releaseTrajectory();
		
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
					if (ret < 0 || (size_t) ret >= sizeof(replicate_filename)) {
						fprintf(stderr, "Error: Output filename too long\n");
						exit(1);
					}
				} else {
					// No extension - append replicate number
					int ret = snprintf(replicate_filename, sizeof(replicate_filename), 
						"%s_rep%d", tskitOutputFilename, replicate_num);
					if (ret < 0 || (size_t) ret >= sizeof(replicate_filename)) {
						fprintf(stderr, "Error: Output filename too long\n");
						exit(1);
					}
				}
			}
			
			VERBOSE_LOG("Writing tree sequence to %s\n", replicate_filename);
			if (tskit_finalize(replicate_filename) != 0) exit(1);
		} else {
			// No file output - just need to sort tables for internal consistency
			// (tskit_finalize with NULL filename would fail)
		}
		
		// Clean up tskit for this replicate
		tskit_cleanup();
		
	}
	
	// tskit finalization now handled per-replicate
	
	free(currentSize);
		free(events);
	
	// Final cleanup after all replicates
	// Clean all nodes in the array
	for (int cleanup_i = 0; cleanup_i < nodesCapacity; cleanup_i++) {
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
	releaseFreeLineages();
	cleanupPopLists();
	
	// Free the node arrays themselves (but not the nodes, we already did that)
	if (nodes != NULL) {
		free(nodes);
		nodes = NULL;
	}
	cleanupSampleNodeIds();
	
	// Clean up segment pool
	destroySegmentPool();

	return(0);
}

/* -eG, -eL and -eM apply to every population, but -p, which sets the
 * number of populations, can come after them. They are parsed into one
 * event with allPopulations set, which is replaced here, after parsing,
 * by one event for each population (-eG, -eL) or for each ordered pair
 * (-eM). The new events take the place of the old one, so events at the
 * same time keep the order of the options. */
static void expandAllPopulationEvents(void) {
	int numAll = 0;
	for (int i = 0; i < eventNumber; i++)
		if (events[i].allPopulations) numAll++;
	if (numAll == 0) return;
	const int oldNumber = eventNumber;
	struct event *old = malloc(sizeof(struct event) * oldNumber);
	if (old == NULL) {
		fprintf(stderr, "Error: out of memory\n");
		exit(1);
	}
	memcpy(old, events, sizeof(struct event) * oldNumber);
	eventNumber = 0;
	for (int i = 0; i < oldNumber; i++) {
		if (!old[i].allPopulations) {
			ensureEventsCapacity();
			events[eventNumber++] = old[i];
			continue;
		}
		for (int a = 0; a < npops; a++) {
			for (int b = 0; b < npops; b++) {
				/* Growth: one event per population (b only counts once).
				 * Migration: one per ordered pair, source a, destination b. */
				if (old[i].type == 'm' ? a == b : b > 0) continue;
				ensureEventsCapacity();
				events[eventNumber] = old[i];
				events[eventNumber].allPopulations = 0;
				if (old[i].type == 'm') {
					events[eventNumber].popID2 = a;
					events[eventNumber].popID = b;
				} else {
					events[eventNumber].popID = a;
				}
				eventNumber++;
			}
		}
	}
	free(old);
}

/* Check the demographic events of the command line (and of YAML, which
 * checks them already, and of a demes model) before they are sorted: a NaN
 * time breaks the sort, a negative time puts an event after the present, a
 * population outside 0 .. npops-1 was ignored, and a NaN rate made the
 * simulation hang. Ancient samples and sweeps have their own checks. */
static void checkDemographicEvents(void) {
	for (int i = 0; i < eventNumber; i++) {
		const struct event *e = &events[i];
		const char *what;
		switch (e->type) {
			case 'n': what = "-en (YAML population_size_changes)"; break;
			case 'g': what = "-eg or -eG (YAML exponential_growth_changes)"; break;
			case 'l': what = "-el or -eL (YAML linear_growth_changes)"; break;
			case 'm': what = "-em or -eM (YAML migration_rate_changes)"; break;
			case 'p': what = "-ed or -ej (YAML population_splits)"; break;
			case 'a': what = "-ea (YAML admixtures)"; break;
			default: continue;
		}
		const double t = internalTimeToUser(e->time);
		if (!isfinite(t) || t < 0) {
			fprintf(stderr, "Error: %s: the time (%g) must be a finite "
				"number >= 0\n", what, t);
			exit(1);
		}
		/* Population fields in use: popID for all; popID2 for -em (the
		 * source), -ed (the ancestral population) and -ea (source 1); popID3
		 * for -ea (source 2). */
		int pops[3] = {e->popID, e->popID2, e->popID3};
		int numPops = (e->type == 'a') ? 3 : (e->type == 'm' || e->type == 'p') ? 2 : 1;
		for (int k = 0; k < numPops; k++) {
			if (pops[k] < 0 || pops[k] >= npops) {
				fprintf(stderr, "Error: %s at time %g: population %d must be "
					"between 0 and %d (the number of populations, -p, minus 1)\n",
					what, t, pops[k], npops - 1);
				exit(1);
			}
		}
		/* A merge of a population into itself never ends, and migration
		 * into the same population has no meaning. An -ea source can be the
		 * admixed population itself: the lineages that do not move stay,
		 * which makes a pulse (the demes importer uses this). */
		if ((e->type == 'm' || e->type == 'p') && e->popID == e->popID2) {
			fprintf(stderr, "Error: %s at time %g: the two populations must "
				"be different (population %d is given twice)\n",
				what, t, e->popID);
			exit(1);
		}
		/* The value: a size (-en), a growth rate (-eg, -el; stored per 2N0
		 * generations) or a migration rate (-em). Its range is checked in
		 * validateShapeTrajectories(); here only that it is a number. */
		double value = e->popnSize;
		if (e->type == 'g' || e->type == 'l') value *= USER_TO_INTERNAL_TIME;
		if ((e->type == 'n' || e->type == 'g' || e->type == 'l' || e->type == 'm')
		    && !isfinite(value)) {
			fprintf(stderr, "Error: %s at time %g: the value (%g) must be a "
				"finite number\n", what, t, value);
			exit(1);
		}
	}
}

/* The value of an option: the next argument. Stops with a message if
 * there is none. */
static const char* getNextArg(int argc, const char **argv, int *args, const char *option) {
	(*args)++;
	if (*args >= argc) {
		fprintf(stderr, "Error: option %s needs more values; try 'discoal --help'\n", option);
		exit(1);
	}
	return argv[*args];
}

/* An integer value of an option. The whole argument must be an integer in
 * the range of int: atoi() would read "abc" as 0 and "12x" as 12. */
static int parseIntArg(int argc, const char **argv, int *args, const char *option) {
	const char *arg = getNextArg(argc, argv, args, option);
	char *endptr;
	errno = 0;
	long val = strtol(arg, &endptr, 10);
	if (endptr == arg || *endptr != '\0' || errno == ERANGE || val > INT_MAX || val < INT_MIN) {
		fprintf(stderr, "Error: Invalid integer argument for %s: '%s'\n", option, arg);
		exit(1);
	}
	return (int)val;
}

/* A numeric value of an option; the whole argument must be a finite
 * number (strtod() also reads "nan", "inf" and, as infinity, numbers too
 * large for a double). */
static double parseDoubleArg(int argc, const char **argv, int *args, const char *option) {
	const char *arg = getNextArg(argc, argv, args, option);
	char *endptr;
	double val = strtod(arg, &endptr);
	if (endptr == arg || *endptr != '\0' || !isfinite(val)) {
		fprintf(stderr, "Error: Invalid numeric argument for %s: '%s'\n", option, arg);
		exit(1);
	}
	return val;
}

/* Every option name. An argument that is not one of them is an error:
 * the parser below looks only at the first letters, so without this check
 * -e, -P or -uX would be taken for another option or ignored. -s, -D and
 * -T are here so that the parser gives their own messages. A new option
 * goes here, in longOptions[], in the parser of getParameters(), and in
 * printUsage(). */
static int isOption(const char *arg) {
	static const char *const names[] = {
		"-t", "-r", "-g", "-gr", "-p", "-Y", "-en", "-eg", "-eG", "-el",
		"-eL", "-ed", "-ej", "-ea", "-em", "-eM", "-ws", "-wd", "-wn", "-ls",
		"-ld", "-ln", "-f", "-uA", "-N", "-a", "-x", "-c", "-i", "-S", "-M",
		"-m", "-A", "-Pt", "-Pr", "-Pre", "-Pa", "-Pu", "-PuA", "-Px", "-Pf",
		"-Pc", "-Pe1", "-Pe2", "-R", "-L", "-h", "-v", "-ts", "-F", "-d",
		"-s", "-D", "-T",
	};
	for (size_t k = 0; k < sizeof(names) / sizeof(names[0]); k++)
		if (strcmp(arg, names[k]) == 0) return 1;
	return 0;
}

/* The long name of each option and the option it stands for. A long name
 * is an exact synonym: the same values in the same order. The names
 * follow the YAML keys. printUsage() lists them too. */
static const struct { const char *longName, *option; } longOptions[] = {
		{"--theta", "-t"},
		{"--rho", "-r"},
		{"--gene-conversion", "-g"},
		{"--gene-conversion-ratio", "-gr"},
		{"--populations", "-p"},
		{"--migration", "-M"},
		{"--migration-rate", "-m"},
		{"--ancient-sample", "-A"},
		{"--size-change", "-en"},
		{"--exponential-growth", "-eg"},
		{"--exponential-growth-all", "-eG"},
		{"--linear-growth", "-el"},
		{"--linear-growth-all", "-eL"},
		{"--join", "-ed"},
		{"--admixture", "-ea"},
		{"--migration-change", "-em"},
		{"--migration-change-all", "-eM"},
		{"--sweep", "-ws"},
		{"--sweep-deterministic", "-wd"},
		{"--sweep-neutral", "-wn"},
		{"--linked-sweep", "-ls"},
		{"--linked-sweep-deterministic", "-ld"},
		{"--linked-sweep-neutral", "-ln"},
		{"--trajectory-file", "-S"},
		{"--recurrent-sweeps", "-R"},
		{"--recurrent-linked-sweeps", "-L"},
		{"--alpha", "-a"},
		{"--sweep-position", "-x"},
		{"--initial-frequency", "-f"},
		{"--beneficial-mutation-rate", "-uA"},
		{"--final-frequency", "-c"},
		{"--hide-selected-snp", "-h"},
		{"--sweep-population-size", "-N"},
		{"--sweep-time-step", "-i"},
		{"--prior-theta", "-Pt"},
		{"--prior-rho", "-Pr"},
		{"--prior-rho-exponential", "-Pre"},
		{"--prior-alpha", "-Pa"},
		{"--prior-sweep-time", "-Pu"},
		{"--prior-beneficial-mutation-rate", "-PuA"},
		{"--prior-sweep-position", "-Px"},
		{"--prior-initial-frequency", "-Pf"},
		{"--prior-final-frequency", "-Pc"},
		{"--prior-size-change", "-Pe1"},
		{"--config", "-Y"},
		{"--tree-sequence", "-ts"},
		{"--full-arg", "-F"},
		{"--seed", "-d"},
		{"--verbose", "-v"},
};

/* The arguments with each long option name replaced by its option, so
 * that the parser sees one name for each option. An argument that starts
 * with "--" and is not a long name is left as it is; the parser reports
 * it as an unknown option. */
static const char **withShortNames(int argc, const char **argv) {
	static const char **copy = NULL;
	free(copy);
	copy = malloc(sizeof(const char *) * (argc + 1));
	if (copy == NULL) {
		fprintf(stderr, "Error: out of memory\n");
		exit(1);
	}
	for (int i = 0; i < argc; i++) {
		copy[i] = argv[i];
		if (argv[i][0] != '-' || argv[i][1] != '-') continue;
		for (size_t k = 0; k < sizeof(longOptions) / sizeof(longOptions[0]); k++)
			if (strcmp(argv[i], longOptions[k].longName) == 0)
				copy[i] = longOptions[k].option;
	}
	copy[argc] = NULL;
	return copy;
}

void getParameters(int argc,const char **argv){
	int args;
	int i,j;
	double migR;
	int nChangeCheck;
	int yaml_config_loaded = 0;
	
	// Initialize all variables first (before any YAML or command line processing)
	npops = 1;
	leftRho = 0.0;
	rho = 0.0;
	my_gamma = 0.0;
	gcMean = 0;
	theta = 0.0;
	alpha = 0.0;
	tau = 0;
	f0=0.0;
        uA=0.0;

	/* setall() takes seeds in [1, 2147483562] and [1, 2147483398]. */
	seed1 = (long) (devrand() % 2147483398) + 1;
	seed2 = (long) (devrand() % 2147483398) + 1;

	EFFECTIVE_POPN_SIZE = 1000000;
	sweepSite = 0.5;
	tDiv=666;
        gammaCoRatioMode = 0;
	priorTheta=priorRho=priorAlpha=priorTau=priorX=priorF0=priorUA=priorC=0;
	numPriorSizeChanges = 0;

	migFlag = 0;
	deltaTMod = 40;
	recurSweepMode = 0;
	recurrentLeftOfLocus = 0;
	sweepPositionGiven = 0;
	demesModelLoaded = 0;
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

	/* --help (also -help, and -h as the only argument) and --version (-V)
	 * write to stdout and end with status 0, before anything else is read.
	 * In a run, -h hides the selected SNP. */
	const int onlyH = (argc == 2 && strcmp(argv[1], "-h") == 0);
	argv = withShortNames(argc, argv);
	for (i = 1; i < argc; i++) {
		if (onlyH || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0) {
			printUsage(stdout);
			exit(0);
		}
		if (strcmp(argv[i], "--version") == 0 || strcmp(argv[i], "-V") == 0) {
			printf("discoal %s\n", DISCOAL_VERSION);
			exit(0);
		}
	}

	// Locate -Y first so we can validate "-Y xor positional args" before
	// doing any work. The positional form "sampleSize numReplicates nSites"
	// and -Y are mutually exclusive: loading YAML and then running the
	// positional parser would read "-Y" or another flag as a number (when -Y
	// comes before the positional slots) or overwrite the per-pop state YAML
	// just installed (when positional slots come first).
	// -v must be known before the YAML file (and any demes file it names)
	// is loaded, because loading prints its diagnostics then.
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) verboseFlag = 1;
	}

	int yaml_arg_idx = -1;
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-Y") != 0) continue;
		if (yaml_arg_idx != -1) {
			fprintf(stderr, "Error: -Y can be given only once\n");
			exit(1);
		}
		yaml_arg_idx = i;
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
		/* The three positional arguments are argv[1 .. 3]; the options
		 * start at argv[4]. */
		args = 0;
		sampleSize = parseIntArg(argc, argv, &args, "sampleSize");
		sampleNumber = parseIntArg(argc, argv, &args, "numReplicates");
		nSites = parseIntArg(argc, argv, &args, "nSites");
		args++;

		popnSizes[0]=sampleSize;
		popnSizes[1]=0;
		sampleSizes[0]=sampleSize;
		sampleSizes[1]=0;
	} else {
		// sampleSize/sampleSizes[]/popnSizes[] already set via apply_yaml_config;
		// skip positional parsing and start the flag loop at argv[1].
		args = 1;
	}

	while(args < argc){
		// Skip YAML config arguments if already processed
		if (strcmp(argv[args], "-Y") == 0) {
			args += 2;  // Skip -Y and its filename argument
			continue;
		}
		
		// Check if argument starts with '-'
		if (argv[args][0] != '-') {
			fprintf(stderr, "Error: Unexpected argument '%s' (an option must start with '-'); try 'discoal --help'\n", argv[args]);
			exit(1);
		}
		
		if (!isOption(argv[args])) {
			fprintf(stderr, "Error: Unknown option '%s'; try 'discoal --help'\n", argv[args]);
			exit(1);
		}
		/* The option's name, for messages about its values. */
		const char *opt = argv[args];

		/* A demes model gives the demography, and the YAML its samples,
		 * checked against the demes' lifetimes (as the YAML allows only
		 * ancient_samples with demes_filename). Options that change the
		 * demography or the samples would bypass these checks. */
		if (demesModelLoaded &&
		    (strcmp(argv[args], "-p") == 0 || strcmp(argv[args], "-A") == 0 ||
		     strcmp(argv[args], "-m") == 0 || strcmp(argv[args], "-M") == 0 ||
		     strcmp(argv[args], "-Pe1") == 0 || strcmp(argv[args], "-Pe2") == 0 ||
		     argv[args][1] == 'e')) {
			fprintf(stderr, "Error: %s cannot be used with a demes model "
				"(demography.demes_filename): the demes file gives the "
				"demography, and demography.deme_sample_size and "
				"demographic_events.ancient_samples the samples\n", argv[args]);
			exit(1);
		}

		switch(argv[args][1]){
			case 'F' :
			if (!tskitOutputMode) {
				fprintf(stderr, "Error: -F needs -ts before it, as in '-ts file.trees -F'\n");
				exit(1);
			}
			minimalTreeSeq = 0;  // Disable minimal mode, use full recording
			break;
			case 'S' :
			if (setSweepTrajectoryFile(getNextArg(argc, argv, &args, opt)) != 0)
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
				const char *filename = getNextArg(argc, argv, &args, opt);
				if (filename[0] == '\0' || filename[0] == '-') {
					fprintf(stderr, "Error: -ts needs a filename, not '%s'; try 'discoal --help'\n", filename);
					exit(1);
				}
				if (strlen(filename) >= sizeof(tskitOutputFilename)) {
					fprintf(stderr, "Error: the -ts filename is too long (at most %zu characters)\n",
						sizeof(tskitOutputFilename) - 1);
					exit(1);
				}
				strcpy(tskitOutputFilename, filename);
			}
			else {  // -t flag for theta
				theta = parseDoubleArg(argc, argv, &args, opt);
			}
			break;
			case 'i' :
			deltaTMod = parseDoubleArg(argc, argv, &args, opt);
			break;
			case 'r' :
			rho = parseDoubleArg(argc, argv, &args, opt);
			break;
			case 'g' :
                        if (argv[args][2] == 'r')
                        {
			  gammaCoRatio = parseDoubleArg(argc, argv, &args, opt);
			  gcMean = parseIntArg(argc, argv, &args, opt);
                          gammaCoRatioMode = 1;
                        }
                        else
                        {
			  my_gamma = parseDoubleArg(argc, argv, &args, opt);
			  gcMean = parseIntArg(argc, argv, &args, opt);
                        }
			break;
			case 'a' :
			alpha = parseDoubleArg(argc, argv, &args, opt);
			break;
			case 'D' :
			fprintf(stderr,
				"Error: the -D command-line flag for loading demes files is no "
				"longer supported. Use a YAML config (-Y) with "
				"demography.demes_filename instead. See "
				"config_examples/demes_example.yaml.\n");
			exit(1);
			case 'x' :
			sweepSite = parseDoubleArg(argc, argv, &args, opt);
			sweepPositionGiven = 1;
			break;
			case 'M' :
			if(npops==1){
				fprintf(stderr,"Error: %s needs more than one population; give -p before %s\n", opt, opt);
				exit(1);
			}
			migR = parseDoubleArg(argc, argv, &args, opt);
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
					fprintf(stderr,"Error: %s needs more than one population; give -p before %s\n", opt, opt);
					exit(1);
				}
				i = parseIntArg(argc, argv, &args, opt);
				j = parseIntArg(argc, argv, &args, opt);
				migR = parseDoubleArg(argc, argv, &args, opt);
				if (i < 0 || i >= npops || j < 0 || j >= npops || i == j) {
					fprintf(stderr, "Error: -m %d %d: the two populations must be different and between 0 and %d (give -p before -m)\n",
						i, j, npops - 1);
					exit(1);
				}
				migMatConst[i][j]=migR;
				migFlag = 1;
			break;
			case 'p' :
			npops = parseIntArg(argc, argv, &args, opt);
			if(npops < 1){
				fprintf(stderr,"Error: -p needs at least 1 population (got %d)\n", npops);
				exit(1);
			}
			if(npops > MAXPOPS){
				fprintf(stderr,"Error: too many populations defined. Current maximum number = %d. Change MAXPOPS in discoal.h and recompile\n",MAXPOPS);
				exit(1);
			}
			for(i=0;i<npops;i++){
				sampleSizes[i] = parseIntArg(argc, argv, &args, opt);
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
					events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].popnSize = parseDoubleArg(argc, argv, &args, opt);
					events[eventNumber].type = 'n'; //size change
					eventNumber++;
					break;
					case 'd' :
					case 'j' :  // -ej is equivalent to -ed (for ms compatibility)
					tDiv =  parseDoubleArg(argc, argv, &args, opt);
						ensureEventsCapacity();
					events[eventNumber].time = userTimeToInternal(tDiv);
					events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].popID2 = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].type = 'p'; //pop split
					eventNumber++;
					break;
					case 'a' :
						ensureEventsCapacity();
					events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].popID2 = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].popID3 = parseIntArg(argc, argv, &args, opt);
					events[eventNumber].admixProp = parseDoubleArg(argc, argv, &args, opt);
					events[eventNumber].type = 'a'; //admix split
					eventNumber++;
					break;
					/* -eg/-eG alpha and -el/-eL gamma are rates per 4N0
					 * generations, as in ms; userRateToInternal converts them
					 * like the times. */
					case 'g':
						/* -eg time popID alpha */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);
						events[eventNumber].popnSize = userRateToInternal(parseDoubleArg(argc, argv, &args, opt));  /* alpha stored in popnSize field */
						events[eventNumber].type = 'g';
						eventNumber++;
						break;
					case 'G':
						/* -eG time alpha: -eg for every population. -p can
						 * come later, so expandAllPopulationEvents() makes
						 * the per-population events after parsing. */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].allPopulations = 1;
						events[eventNumber].popnSize = userRateToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].type = 'g';
						eventNumber++;
						break;
					case 'l':
						/* -el time popID gamma: gamma is in relative size (size / N0) */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);
						events[eventNumber].popnSize = userRateToInternal(parseDoubleArg(argc, argv, &args, opt));  /* gamma stored in popnSize field */
						events[eventNumber].type = 'l';
						eventNumber++;
						break;
					case 'L':
						/* -eL time gamma: -el for every population, expanded
						 * after parsing as for -eG. */
						ensureEventsCapacity();
						events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].allPopulations = 1;
						events[eventNumber].popnSize = userRateToInternal(parseDoubleArg(argc, argv, &args, opt));
						events[eventNumber].type = 'l';
						eventNumber++;
						break;
					case 'm':
						/* -em time srcPop dstPop rate — change migration rate at this time */
						{
							ensureEventsCapacity();
							events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
							events[eventNumber].popID2 = parseIntArg(argc, argv, &args, opt);  /* source pop */
							events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);   /* destination pop */
							events[eventNumber].popnSize = parseDoubleArg(argc, argv, &args, opt); /* new rate stored in popnSize field */
							events[eventNumber].type = 'm';  /* lowercase 'm' for time-varying migration event */
							/* Migration that starts later still counts as
							 * migration for initialize()'s check. */
							if (events[eventNumber].popnSize > 0) migFlag = 1;
							eventNumber++;
						}
						break;
					case 'M':
						/* -eM time rate: -em for every ordered pair of
						 * populations, expanded after parsing as for -eG. */
						{
							ensureEventsCapacity();
							events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
							events[eventNumber].allPopulations = 1;
							events[eventNumber].popnSize = parseDoubleArg(argc, argv, &args, opt);
							events[eventNumber].type = 'm';
							/* As for -em: migration that starts later still
							 * counts for initialize()'s check. */
							if (events[eventNumber].popnSize > 0) migFlag = 1;
							eventNumber++;
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
				tau = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
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
				tau = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
				leftRho = parseDoubleArg(argc, argv, &args, opt) * 2.0;
				leftRhoFlag=1;
					ensureEventsCapacity();
				events[eventNumber].time = tau;
				events[eventNumber].type = 's'; //sweep event
				eventNumber++;
					break;
			case 'f':
			f0 = parseDoubleArg(argc, argv, &args, opt);
			softSweepMode = 1;
			break;
                        case 'u':
                        uA = parseDoubleArg(argc, argv, &args, opt);
                        break;
			case 'P' :
			  switch(argv[args][2]){
				case 't':
				  priorTheta = 1;
				  pThetaLow=parseDoubleArg(argc, argv, &args, opt);
				  pThetaUp=parseDoubleArg(argc, argv, &args, opt);
				break;
				case 'c':
				  priorC = 1;
			          partialSweepMode = 1;
				  pCLow=parseDoubleArg(argc, argv, &args, opt);
				  pCUp=parseDoubleArg(argc, argv, &args, opt);
				break;
                                case 'r':
                                  if (strcmp(opt, "-Pre") == 0){
                                        priorRho = 2;
                                        pRhoMean=parseDoubleArg(argc, argv, &args, opt);
                                        pRhoUp=parseDoubleArg(argc, argv, &args, opt);
                                  }
                                  else{
                                        priorRho = 1;
                                        pRhoLow=parseDoubleArg(argc, argv, &args, opt);
                                        pRhoUp=parseDoubleArg(argc, argv, &args, opt);
                                  }
                                break;
				case 'a':
				  priorAlpha = 1;
				  pAlphaLow=parseDoubleArg(argc, argv, &args, opt);
				  pAlphaUp=parseDoubleArg(argc, argv, &args, opt);
				break;
				case 'u':
                                	if (strcmp(opt, "-PuA") == 0){
					  priorUA = 1;
					  pUALow=parseDoubleArg(argc, argv, &args, opt);
					  pUAUp=parseDoubleArg(argc, argv, &args, opt);
					}
                                	else{
					  priorTau = 1;
					  pTauLow=userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					  pTauUp=userTimeToInternal(parseDoubleArg(argc, argv, &args, opt));
					}
				break;
				case 'x':
				  priorX = 1;
				  pXLow=parseDoubleArg(argc, argv, &args, opt);
				  pXUp=parseDoubleArg(argc, argv, &args, opt);
				break;
				case 'f':
				  priorF0 = 1;
				  pF0Low=parseDoubleArg(argc, argv, &args, opt);
				  pF0Up=parseDoubleArg(argc, argv, &args, opt);
				break;
				case 'e':
				/* -Pe1 / -Pe2 lowTime highTime lowSize highSize: a size
				 * change of population 0 with uniform priors on time and
				 * size. Both work the same way; the digit only names them. */
				{
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
			seed1 = parseIntArg(argc, argv, &args, opt);
			seed2 = parseIntArg(argc, argv, &args, opt);
			break;
			case 'N' :
			EFFECTIVE_POPN_SIZE = parseIntArg(argc, argv, &args, opt);
			break;
			case 'R' :
			recurSweepMode = 1;
			sweepMode = 's';
			recurSweepRate = parseDoubleArg(argc, argv, &args, opt);
			break;
			case 'L' :
			recurSweepMode = 1;
			recurrentLeftOfLocus = 1;
			sweepMode = 's';
			sweepSite = -1.0;
			recurSweepRate = parseDoubleArg(argc, argv, &args, opt);
			break;
			case 'c' :
			partialSweepMode = 1;
			sweepMode = 's';
			partialSweepFinalFreq = parseDoubleArg(argc, argv, &args, opt);
            if (partialSweepFinalFreq <= 0.0 || partialSweepFinalFreq >= 1.0)
            {
                fprintf(stderr,"Error: -c (%g) must be above 0 and below 1\n", partialSweepFinalFreq);
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
			events[eventNumber].lineageNumber = parseIntArg(argc, argv, &args, opt);
			events[eventNumber].popID = parseIntArg(argc, argv, &args, opt);	 
			events[eventNumber].time = userTimeToInternal(parseDoubleArg(argc, argv, &args, opt)); 
			events[eventNumber].type = 'A'; //ancient sample
			eventNumber++;
			break;
			case 'T' :
			fprintf(stderr, "Error: -T tree output mode has been removed.\n");
			fprintf(stderr, "Use -ts <filename.trees> for tree sequence output instead.\n");
			exit(1);
			default:
			/* isOption() passed a name that no case handles. */
			fprintf(stderr, "Error: internal error: option %s is not handled\n", opt);
			exit(1);
		}
		args++;
	}
	expandAllPopulationEvents();
	checkDemographicEvents();
	sortEventArray(events,eventNumber);

	if (validateShapeTrajectories(events, eventNumber) != 0) {
		exit(1);
	}

	/* A sweep outside the locus is set by its options (-ls, -ld, -ln, -L, or
	 * their YAML fields), not by the value of sweepSite, which the code
	 * sets to -1 for it: a user -x value could be -1 too. */
	const int sweepOutsideLocus = leftRhoFlag || recurrentLeftOfLocus;
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
			/* The tract length is geometric on 1, 2, ... with mean gcMean,
			 * so a mean below 1 has no such distribution. */
			{gcMean < 0 || ((my_gamma > 0.0 || gammaCoRatio > 0.0) && gcMean < 1),
				"the -g or -gr tract length (YAML gene_conversion_tract_length)",
				(double)gcMean, ">= 1 when the gene conversion rate is above 0, and >= 0"},
			/* Ancestry counts are uint16_t. */
			{sampleSize > 65535, "the sample size", (double)sampleSize, "<= 65535"},
			{sampleSize < 2, "the sample size", (double)sampleSize, ">= 2"},
			{sampleNumber < 1, "the number of replicates", (double)sampleNumber, ">= 1"},
			{nSites < 1, "the number of sites", (double)nSites, ">= 1"},
			{nSites > MAXSITES, "the number of sites", (double)nSites,
				"<= " STRINGIFY(MAXSITES) " (MAXSITES, set at compile time)"},
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
			/* The bounds of the generator's setall(). */
			{seed1 < 1 || seed1 > 2147483562, "the first seed (-d or YAML seed)",
				(double)seed1, "in [1, 2147483562]"},
			{seed2 < 1 || seed2 > 2147483398, "the second seed (-d or YAML seed)",
				(double)seed2, "in [1, 2147483398]"},
		};
		for(size_t k = 0; k < sizeof(fixed) / sizeof(fixed[0]); k++){
			if(fixed[k].bad){
				fprintf(stderr, "Error: %s (%.10g) must be %s\n", fixed[k].what,
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
	nChangeCheck=0;
	/* The single sweep, if any: its event index and its start time. */
	int sweepEvents = 0, sweepIndex = -1;
	double sweepStart = 0.0;
	for(i=0;i<eventNumber;i++){
	 		if(events[i].type == 's'){
				sweepEvents++;
				sweepIndex = i;
				sweepStart = events[i].time;
	 		}
			if(events[i].type == 'n'){
				nChangeCheck += 1;
	 		}
	}
	/* The sweep sampler follows one sweep at a time; recurrent sweeps (-R,
	 * -L) give many. */
	if(sweepEvents > 1){
		fprintf(stderr, "Error: only one single sweep can be given (-ws, -wd, "
			"-wn, -ls, -ld, -ln, -S, or a single sweep in a YAML selection "
			"block); for many sweeps use recurrent sweeps (-R, -L)\n");
		exit(1);
	}
	if(sweepEvents > 0){
		if(recurSweepMode == 1){
			fprintf(stderr, "Error: a single sweep (-ws, -wd, -wn, -ls, -ld, -ln, -S) cannot be combined with recurrent sweeps (-R, -L)\n");
			exit(1);
		}
		if(nChangeCheck > 1 && sweepMode=='d'){
			fprintf(stderr, "Error: a deterministic sweep (-wd) cannot be combined with more than one population size change (-en); use a stochastic sweep (-ws)\n");
			exit(1);
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
	/* Options that contradict the sweep model, from the CLI or YAML: a
	 * drawn position puts the sweep inside the locus, and a drawn tau needs
	 * a single sweep. A sweep position (-x or YAML sweep_position) is
	 * rejected with a sweep outside the locus whatever its value, also when
	 * one comes from the command line and the other from YAML. */
	if(priorX && sweepOutsideLocus){
		fprintf(stderr, "Error: a prior on the sweep position (-Px, or YAML "
			"priors.sweep_position) cannot be combined with a sweep outside "
			"the locus (-ls, -ld, -ln, -L, or YAML left_rho or "
			"recurrent_sweeps_left_of_locus)\n");
		exit(1);
	}
	if(priorTau && recurSweepMode){
		fprintf(stderr, "Error: a prior on the sweep time (-Pu, or YAML "
			"priors.fixation_time_ago) needs a single sweep, not recurrent "
			"sweeps (-R, -L)\n");
		exit(1);
	}
	if(sweepOutsideLocus && sweepPositionGiven){
		fprintf(stderr, "Error: -x (or YAML sweep_position) cannot be combined "
			"with a sweep outside the locus (-ls, -ld, -ln, -L, or YAML left_rho "
			"or recurrent_sweeps_left_of_locus)\n");
		exit(1);
	}
	if(softSweep && recurSweepMode == 1){
		fprintf(stderr, "Error: recurrent soft sweeps (-R or -L with -f or "
			"-Pf) are not implemented\n");
		exit(1);
	}
	/* Sweep trajectories are not modeled under growth: a trajectory is
	 * drawn, and accepted (sampleSweepTrajectory()), with sizes that follow
	 * only size changes (-en). So population 0 may not grow or shrink continuously
	 * (-eg, -el, -eG, -eL, or an exponential epoch of a demes model) at or
	 * after the earliest start of the sweep (with -Pu, its low bound), going
	 * back in time. Growth in other populations, and growth of population 0
	 * that surely ends before the sweep, are allowed. Recurrent sweeps draw coalescence in every
	 * population at constant rates between sweeps, so with them no
	 * population may grow at any time. A trajectory read from a file (-S)
	 * does not use that ratio, and its sampler follows the size at each
	 * step, so -S is not restricted. The events are sorted here, with ties
	 * in command-line order. */
	if((sweepEvents > 0 && trajFileRows == 0) || recurSweepMode){
		int growthAtSweep = 0;
		for(i = 0; recurSweepMode && i < eventNumber; i++)
			if((events[i].type == 'g' || events[i].type == 'l') &&
			   events[i].popnSize != 0.0)
				growthAtSweep = 1;
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
	if(sweepEvents > 0 || recurSweepMode){
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
		if(recurSweepMode || softSweepMode || partialSweepMode || alpha != 0.0 ||
		   priorAlpha || priorTau || priorF0 || priorC){
			fprintf(stderr, "Error: a trajectory file (-S, or YAML sweep_mode "
				"'trajectory_file') gives the sweep's time and frequencies; it "
				"cannot be combined with -R, -L, -a, -f, -c, -Pa, -Pu, -Pf or "
				"-Pc\n");
			exit(1);
		}
		if(trajFilePoints[0].freq < 1.0){
			partialSweepMode = 1;
			partialSweepFinalFreq = trajFilePoints[0].freq;
		}
		buildSweepTrajectoryFromFile();
	}
	
}

/* The usage text. --help writes it to stdout; a call with too few
 * arguments writes it to stderr (usage()). */
static void printUsage(FILE *out){
	fprintf(out, "▗▄▄▄ ▗▄▄▄▖ ▗▄▄▖ ▗▄▄▖ ▗▄▖  ▗▄▖ ▗▖   \n");
	fprintf(out, "▐▌  █  █  ▐▌   ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌   \n");
	fprintf(out, "▐▌  █  █   ▝▀▚▖▐▌   ▐▌ ▐▌▐▛▀▜▌▐▌   \n");
	fprintf(out, "▐▙▄▄▀▗▄█▄▖▗▄▄▞▘▝▚▄▄▖▝▚▄▞▘▐▌ ▐▌▐▙▄▄\n\n\n");
	
	fprintf(out,
		"discoal version %s\n"
		"A coalescent simulator with recombination, population structure and\n"
		"selective sweeps.\n"
		"\n"
		"usage: discoal sampleSize numReplicates nSites [options]\n"
		"       discoal -Y config.yaml [options]\n"
		"\n"
		"examples:\n"
		"  discoal 10 5 10000 -t 10 -r 10\n"
		"      neutral, with recombination\n"
		"  discoal 10 5 10000 -t 10 -r 10 -ws 0.01 -a 1000 -x 0.5\n"
		"      a hard sweep at the middle of the locus\n"
		"  discoal 20 1 10000 -t 10 -p 2 10 10 -ed 1.0 1 0\n"
		"      two populations that join\n"
		"  discoal 10 3 10000 -t 10 -r 10 -ts out.trees\n"
		"      tree sequences (tskit)\n"
		"\n"
		"Each option has a short name, as in ms, and a long name; they do the same.\n"
		"\n"
		"Times are in units of 4N0 generations, sizes in units of N0, and theta,\n"
		"rho and migration rates are 4N0 times the rate (as in ms). Populations are\n"
		"numbered from 0.\n"
		"\n"
		"mutation and recombination:\n"
		"  -t, --theta theta       4N0 u for the whole locus\n"
		"  -r, --rho rho           4N0 r for the whole locus (default 0)\n"
		"  -g, --gene-conversion gamma tractLength\n"
		"                          gene conversion: 4N0 g for the whole locus, and the\n"
		"                          mean tract length in sites\n"
		"  -gr, --gene-conversion-ratio ratio tractLength\n"
		"                          gene conversion at ratio x rho\n"
		"\n"
		"populations and samples:\n"
		"  -p, --populations npops n1 n2 ...\n"
		"                          number of populations and the sample size of each\n"
		"  -M, --migration migRate\n"
		"                          migration between all populations at migRate\n"
		"                          (give -p first)\n"
		"  -m, --migration-rate pop1 pop2 migRate\n"
		"                          migration from pop1 to pop2 (give -p first)\n"
		"  -A, --ancient-sample n pop time\n"
		"                          n of the samples of pop are ancient, taken at time\n"
		"\n"
		"demographic events:\n"
		"  -en, --size-change time pop size\n"
		"                          size of pop changes to size\n"
		"  -eg, --exponential-growth time pop alpha\n"
		"                          exponential growth of pop at rate alpha per 4N0\n"
		"                          generations (positive: smaller in the past)\n"
		"  -eG, --exponential-growth-all time alpha\n"
		"                          -eg for all populations\n"
		"  -el, --linear-growth time pop gamma\n"
		"                          linear growth of pop; gamma is the change in size\n"
		"                          per 4N0 generations (positive: smaller in the past)\n"
		"  -eL, --linear-growth-all time gamma\n"
		"                          -el for all populations\n"
		"  -ed, --join time pop1 pop2\n"
		"                          pop1 joins pop2, backward in time (-ej is the same)\n"
		"  -ea, --admixture time pop src1 src2 prop\n"
		"                          admixture: backward in time a lineage of pop goes\n"
		"                          to src1 with probability prop, else to src2\n"
		"  -em, --migration-change time pop1 pop2 migRate\n"
		"                          migration from pop1 to pop2 changes to migRate\n"
		"  -eM, --migration-change-all time migRate\n"
		"                          all migration rates change to migRate\n"
		"\n"
		"sweeps:\n"
		"  -ws, --sweep tau        stochastic sweep that ended tau ago\n"
		"  -wd, --sweep-deterministic tau\n"
		"                          deterministic sweep that ended tau ago\n"
		"  -wn, --sweep-neutral tau\n"
		"                          neutral fixation that ended tau ago\n"
		"  -ls, --linked-sweep tau leftRho\n"
		"                          stochastic sweep to the left of the locus, at a\n"
		"                          distance of 4N0 r = 2 x leftRho\n"
		"  -ld, --linked-sweep-deterministic tau leftRho\n"
		"                          as -ls, deterministic\n"
		"  -ln, --linked-sweep-neutral tau leftRho\n"
		"                          as -ls, neutral fixation\n"
		"  -S, --trajectory-file file\n"
		"                          sweep whose trajectory is in file: lines of time\n"
		"                          (increasing) and allele frequency, ending at 0\n"
		"  -R, --recurrent-sweeps rate\n"
		"                          recurrent sweeps in the locus; rate per 2N\n"
		"                          individuals per generation\n"
		"  -L, --recurrent-linked-sweeps rate\n"
		"                          recurrent sweeps to the left of the locus, at a\n"
		"                          distance drawn uniformly from (0, 4Ns)\n"
		"  -a, --alpha alpha       strength of selection, 2Ns\n"
		"  -x, --sweep-position position\n"
		"                          position of the selected site in the locus, 0 to 1\n"
		"  -f, --initial-frequency f0\n"
		"                          soft sweep: selection starts at frequency f0\n"
		"  -uA, --beneficial-mutation-rate rate\n"
		"                          soft sweep: rate of recurrent beneficial mutation\n"
		"  -c, --final-frequency frequency\n"
		"                          partial sweep: final frequency of the allele\n"
		"  -h, --hide-selected-snp\n"
		"                          hide the selected SNP of a partial sweep\n"
		"  -N, --sweep-population-size size\n"
		"                          effective population size during a sweep\n"
		"                          (default 1000000)\n"
		"  -i, --sweep-time-step dt\n"
		"                          time step of sweep trajectories, 1/(dt N)\n"
		"                          (default 40)\n"
		"\n"
		"priors (a value is drawn for each replicate):\n"
		"  -Pt, --prior-theta low high\n"
		"                          theta, uniform\n"
		"  -Pr, --prior-rho low high\n"
		"                          rho, uniform\n"
		"  -Pre, --prior-rho-exponential mean upperBound\n"
		"                          rho, exponential, cut at upperBound\n"
		"  -Pa, --prior-alpha low high\n"
		"                          alpha, uniform\n"
		"  -Pu, --prior-sweep-time low high\n"
		"                          tau, uniform (give a sweep option too; its tau is\n"
		"                          not used)\n"
		"  -PuA, --prior-beneficial-mutation-rate low high\n"
		"                          uA, uniform\n"
		"  -Px, --prior-sweep-position low high\n"
		"                          position of the selected site, uniform\n"
		"  -Pf, --prior-initial-frequency low high\n"
		"                          f0, uniform\n"
		"  -Pc, --prior-final-frequency low high\n"
		"                          final frequency of a partial sweep, uniform\n"
		"  -Pe1, --prior-size-change lowTime highTime lowSize highSize\n"
		"                          a size change of population 0: time and size,\n"
		"                          uniform\n"
		"  -Pe2, --prior-size-change lowTime highTime lowSize highSize\n"
		"                          a second such size change\n"
		"\n"
		"input and output:\n"
		"  -Y, --config file       read the model from a YAML file; a demes model\n"
		"                          goes in demography.demes_filename\n"
		"  -ts, --tree-sequence file\n"
		"                          write tree sequences (tskit) to file; with several\n"
		"                          replicates, to file_rep1, file_rep2, ...\n"
		"  -F, --full-arg          with -ts (and after it): the full ARG, with the\n"
		"                          nodes and edges of recombination\n"
		"  -d, --seed seed1 seed2  seeds of the random number generator\n"
		"  -v, --verbose           print progress and demes details to stderr\n"
		"  --help                  print this text\n"
		"  --version               print the version\n"
		"\n"
		"Documentation: https://discoal.readthedocs.io\n",
		DISCOAL_VERSION);
}

void usage(){
	printUsage(stderr);
	exit(1);
}

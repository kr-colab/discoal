#include "unity.h"
#include <stdio.h>
#include "../../discoal.h"
#include "../../discoalFunctions.h"
#include "../../ranlib.h"
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <time.h>

// External declarations for test fixtures
extern rootedNode* testNode;
extern event testEvent;
extern AncestrySegment* testSegment;
extern ActiveMaterial testActiveMaterial;

// Forward declarations for all test functions
// From test_node.c
void test_node_initialization(void);
void test_node_set_properties(void);
void test_newRootedNode_creation(void);

// From test_event.c
void test_event_initialization(void);
void test_event_set_properties(void);

// From test_node_operations.c
void test_newRootedNode(void);
void test_addRemoveNode(void);
void test_pickNodePopn(void);
void test_nodePopnSize(void);

// From test_mutations.c
void test_basicNodeCreation(void);
void test_mutationArrayAccess(void);
void test_simpleManualMutation(void);

// From test_ancestry_segment.c
void test_newSegment_creates_valid_segment(void);
void test_newSegment_handles_invalid_range(void);
void test_reference_counting_retain_release(void);
void test_shallowCopySegment_shares_references(void);
void test_getAncestryCount_single_segment(void);
void test_hasAncestry_queries(void);
void test_getAncestryCount_multiple_segments(void);
void test_copySegmentTree_creates_independent_copy(void);
void test_null_safety(void);
void test_mergeAncestryTrees_non_overlapping(void);
void test_splitLeft_basic(void);
void test_splitRight_basic(void);
void test_split_edge_cases(void);

// From test_active_segment.c
void test_initializeActiveMaterial_all_sites_active(void);
void test_isActiveSite_queries(void);
void test_getActiveSiteCount(void);
void test_removeFixedRegion_single_region(void);
void test_removeFixedRegion_multiple_regions(void);
void test_removeFixedRegion_edge_cases(void);
void test_newActiveSegment(void);
void test_segment_coalescing(void);
void test_segment_no_coalescing(void);
void test_updateActiveMaterialFromAncestry_basic(void);
void test_active_segment_null_safety(void);
void test_verifyActiveMaterial(void);

// From test_trajectory.c
void test_ensureTrajectoryCapacity_within_limit(void);
void test_ensureTrajectoryCapacity_exceeds_limit(void);
void test_cleanupRejectedTrajectory_existing_file(void);
void test_cleanupRejectedTrajectory_nonexistent_file(void);
void test_cleanupRejectedTrajectory_null_filename(void);
void test_mmapAcceptedTrajectory_basic(void);
void test_mmapAcceptedTrajectory_replaces_previous(void);
void test_mmapAcceptedTrajectory_large_file(void);
void test_mmapAcceptedTrajectory_read_only(void);
void test_trajectory_filename_generation(void);
void test_multiple_trajectory_cleanup(void);
void test_trajectory_file_persistence(void);

// From test_coalescence_recombination.c
void test_coalesceAtTimePopn_basic(void);
void test_coalesceAtTimePopn_ancestry_merge(void);
void test_coalesceAtTimePopn_overlapping_ancestry(void);
void test_recombineAtTimePopn_basic(void);
void test_recombineAtTimePopn_ancestry_split(void);
void test_geneConversionAtTimePopn_basic(void);
void test_makeGametesMS_mutation_collection(void);
void test_updateActiveMaterial_basic(void);
void test_siteBetweenChunks_basic(void);
void test_pickNodePopn_basic(void);
void test_nodePopnSize_count(void);

// From test_memory_management.c
void test_initializeBreakPoints_basic(void);
void test_cleanupBreakPoints_basic(void);
void test_ensureBreakPointsCapacity_growth(void);
void test_addBreakPoint_multiple(void);
void test_initializeNodeArrays_basic(void);
void test_cleanupNodeArrays_basic(void);
void test_ensureNodesCapacity_growth(void);
void test_ensureAllNodesCapacity_growth(void);
void test_addNode_with_growth(void);
void test_initializeMuts_basic(void);
void test_initializeMuts_zero_capacity(void);
void test_cleanupMuts_basic(void);
void test_ensureMutsCapacity_growth(void);
void test_mutations_stress_test(void);
void test_reinitialize_breakPoints(void);
void test_cleanup_null_safety(void);
void test_integrated_memory_usage(void);

// Per-suite setup/teardown functions
void setUp_node(void) {
    testNode = (rootedNode*)malloc(sizeof(rootedNode));
    testNode->leftParent = NULL;
    testNode->rightParent = NULL;
    testNode->leftChild = NULL;
    testNode->rightChild = NULL;
    testNode->time = 0.0;
    testNode->branchLength = 0.0;
    testNode->blProb = 0.0;
    testNode->nancSites = 0;
    testNode->lLim = 0;
    testNode->rLim = 0;
    testNode->id = 0;
    testNode->mutationNumber = 0;
    testNode->population = 0;
    testNode->sweepPopn = 0;
    testNode->ndes[0] = 0;
    testNode->ndes[1] = 0;
    testNode->leafs = NULL;
    testNode->times[0] = 0.0;
    testNode->times[1] = 0.0;
}

void tearDown_node(void) {
    free(testNode);
}

void setUp_event(void) {
    testEvent.time = 0.0;
    testEvent.popnSize = 0.0;
    testEvent.type = 'N';
    testEvent.popID = 0;
    testEvent.popID2 = 0;
    testEvent.popID3 = 0;
    testEvent.lineageNumber = 0;
    testEvent.admixProp = 0.0;
}

void tearDown_event(void) {
    // Nothing to tear down for event
}

void setUp_node_operations(void) {
    // Initialize any test setup
    initialize();
}

void tearDown_node_operations(void) {
    // Nothing specific to tear down
}

void setUp_mutations(void) {
    // Minimal setup
    sampleSize = 10;
    nSites = 100;
}

void tearDown_mutations(void) {
    // Nothing specific to tear down
}

void setUp_ancestry_segment(void) {
    // testSegment is initialized to NULL
}

void tearDown_ancestry_segment(void) {
    // Cleanup handled in individual tests
}

void setUp_active_segment(void) {
    nSites = 100;  // Set global for tests
    testActiveMaterial.segments = NULL;
    testActiveMaterial.avlTree = NULL;
    testActiveMaterial.totalActive = 0;
}

void tearDown_active_segment(void) {
    freeActiveMaterial(&testActiveMaterial);
}

// External test data for trajectory tests
extern char testTrajFilename[256];
extern float *testData;
extern int testDataSize;

void setUp_trajectory(void) {
    // Initialize globals
    trajectoryFd = -1;
    currentTrajectory = NULL;
    trajectoryFileSize = 0;
    
    // Initialize test-specific data
    testData = NULL;
    testDataSize = 0;
    
    // Create a unique test filename
    snprintf(testTrajFilename, 256, 
             "/tmp/test_traj_%d_%ld.tmp", getpid(), time(NULL));
}

void tearDown_trajectory(void) {
    // Clean up any open trajectory
    if (trajectoryFd != -1) {
        if (currentTrajectory && currentTrajectory != MAP_FAILED) {
            munmap(currentTrajectory, trajectoryFileSize);
        }
        close(trajectoryFd);
        trajectoryFd = -1;
        currentTrajectory = NULL;
    }
    
    // Remove test file if it exists
    unlink(testTrajFilename);
    
    // Free test data
    if (testData) {
        free(testData);
        testData = NULL;
    }
}

// External for coalescence/recombination tests
extern rootedNode *testNode1, *testNode2, *testNode3;
extern int originalSampleSize;
extern int originalNSites;
extern int originalActiveSites;
extern int originalBreakNumber;

void setUp_coalescence_recombination(void) {
    // Save original globals
    originalSampleSize = sampleSize;
    originalNSites = nSites;
    originalActiveSites = activeSites;
    originalBreakNumber = breakNumber;
    
    // Initialize globals for testing
    sampleSize = 10;
    nSites = 100;
    activeSites = 100;
    breakNumber = 0;
    npops = 2;  // Set number of populations
    alleleNumber = 0;  // Reset active node count
    totNodeNumber = 0;  // Reset total node count
    
    // Initialize population sizes
    for (int i = 0; i < MAXPOPS; i++) {
        popnSizes[i] = 0;
        sweepPopnSizes[i] = 0;
    }
    
    // Initialize required arrays
    initializeNodeArrays();
    initializeBreakPoints();
    initializeActiveMaterial(&activeMaterialSegments, nSites);
    
    // Initialize random number generator
    setall(12345, 67890);
    
    // Initialize test nodes
    testNode1 = NULL;
    testNode2 = NULL;
    testNode3 = NULL;
}

void tearDown_coalescence_recombination(void) {
    // Clean up test nodes
    if (testNode1) {
        cleanupMuts(testNode1);
        free(testNode1);
    }
    if (testNode2) {
        cleanupMuts(testNode2);
        free(testNode2);
    }
    if (testNode3) {
        cleanupMuts(testNode3);
        free(testNode3);
    }
    
    // Clean up node arrays
    cleanupNodeArrays();
    cleanupBreakPoints();
    freeActiveMaterial(&activeMaterialSegments);
    
    // Restore original globals
    sampleSize = originalSampleSize;
    nSites = originalNSites;
    activeSites = originalActiveSites;
    breakNumber = originalBreakNumber;
}

// External for memory management tests
extern int originalBreakPointsCapacity;
extern int originalNodesCapacity;
extern int originalAllNodesCapacity;
rootedNode *testNodeMem = NULL;  // Local for memory management

void setUp_memory_management(void) {
    // Save original values
    originalBreakPointsCapacity = breakPointsCapacity;
    originalNodesCapacity = nodesCapacity;
    originalAllNodesCapacity = allNodesCapacity;
    
    // Reset globals
    breakPoints = NULL;
    breakPointsCapacity = 0;
    breakNumber = 0;
    nodes = NULL;
    nodesCapacity = 0;
    allNodes = NULL;
    allNodesCapacity = 0;
    alleleNumber = 0;
    totNodeNumber = 0;
    
    testNodeMem = NULL;
}

void tearDown_memory_management(void) {
    // Clean up any allocations
    cleanupBreakPoints();
    cleanupNodeArrays();
    
    if (testNodeMem) {
        cleanupMuts(testNodeMem);
        free(testNodeMem);
        testNodeMem = NULL;
    }
    
    // Restore original values
    breakPointsCapacity = originalBreakPointsCapacity;
    nodesCapacity = originalNodesCapacity;
    allNodesCapacity = originalAllNodesCapacity;
}

// Global setUp and tearDown that dispatch to appropriate suite functions
void (*current_setUp)(void) = NULL;
void (*current_tearDown)(void) = NULL;

void setUp(void) {
    if (current_setUp) {
        current_setUp();
    }
}

void tearDown(void) {
    if (current_tearDown) {
        current_tearDown();
    }
}

int main(void) {
    UNITY_BEGIN();
    
    printf("\n========== Running Node Tests ==========\n");
    current_setUp = setUp_node;
    current_tearDown = tearDown_node;
    RUN_TEST(test_node_initialization);
    RUN_TEST(test_node_set_properties);
    RUN_TEST(test_newRootedNode_creation);
    
    printf("\n========== Running Event Tests ==========\n");
    current_setUp = setUp_event;
    current_tearDown = tearDown_event;
    RUN_TEST(test_event_initialization);
    RUN_TEST(test_event_set_properties);
    
    printf("\n========== Running Node Operations Tests ==========\n");
    current_setUp = setUp_node_operations;
    current_tearDown = tearDown_node_operations;
    RUN_TEST(test_newRootedNode);
    RUN_TEST(test_addRemoveNode);
    RUN_TEST(test_pickNodePopn);
    RUN_TEST(test_nodePopnSize);
    
    printf("\n========== Running Mutation Tests ==========\n");
    current_setUp = setUp_mutations;
    current_tearDown = tearDown_mutations;
    RUN_TEST(test_basicNodeCreation);
    RUN_TEST(test_mutationArrayAccess);
    RUN_TEST(test_simpleManualMutation);
    
    printf("\n========== Running Ancestry Segment Tests ==========\n");
    current_setUp = setUp_ancestry_segment;
    current_tearDown = tearDown_ancestry_segment;
    RUN_TEST(test_newSegment_creates_valid_segment);
    RUN_TEST(test_newSegment_handles_invalid_range);
    RUN_TEST(test_reference_counting_retain_release);
    RUN_TEST(test_shallowCopySegment_shares_references);
    RUN_TEST(test_getAncestryCount_single_segment);
    RUN_TEST(test_hasAncestry_queries);
    RUN_TEST(test_getAncestryCount_multiple_segments);
    RUN_TEST(test_copySegmentTree_creates_independent_copy);
    RUN_TEST(test_null_safety);
    RUN_TEST(test_mergeAncestryTrees_non_overlapping);
    RUN_TEST(test_splitLeft_basic);
    RUN_TEST(test_splitRight_basic);
    RUN_TEST(test_split_edge_cases);
    
    printf("\n========== Running Active Segment Tests ==========\n");
    current_setUp = setUp_active_segment;
    current_tearDown = tearDown_active_segment;
    RUN_TEST(test_initializeActiveMaterial_all_sites_active);
    RUN_TEST(test_isActiveSite_queries);
    RUN_TEST(test_getActiveSiteCount);
    RUN_TEST(test_removeFixedRegion_single_region);
    RUN_TEST(test_removeFixedRegion_multiple_regions);
    RUN_TEST(test_removeFixedRegion_edge_cases);
    RUN_TEST(test_newActiveSegment);
    RUN_TEST(test_segment_coalescing);
    RUN_TEST(test_segment_no_coalescing);
    RUN_TEST(test_updateActiveMaterialFromAncestry_basic);
    RUN_TEST(test_active_segment_null_safety);
    RUN_TEST(test_verifyActiveMaterial);
    
    printf("\n========== Running Trajectory Tests ==========\n");
    current_setUp = setUp_trajectory;
    current_tearDown = tearDown_trajectory;
    RUN_TEST(test_ensureTrajectoryCapacity_within_limit);
    RUN_TEST(test_ensureTrajectoryCapacity_exceeds_limit);
    RUN_TEST(test_cleanupRejectedTrajectory_existing_file);
    RUN_TEST(test_cleanupRejectedTrajectory_nonexistent_file);
    RUN_TEST(test_cleanupRejectedTrajectory_null_filename);
    RUN_TEST(test_mmapAcceptedTrajectory_basic);
    RUN_TEST(test_mmapAcceptedTrajectory_replaces_previous);
    RUN_TEST(test_mmapAcceptedTrajectory_large_file);
    RUN_TEST(test_mmapAcceptedTrajectory_read_only);
    RUN_TEST(test_trajectory_filename_generation);
    RUN_TEST(test_multiple_trajectory_cleanup);
    RUN_TEST(test_trajectory_file_persistence);
    
    printf("\n========== Running Coalescence/Recombination Tests ==========\n");
    current_setUp = setUp_coalescence_recombination;
    current_tearDown = tearDown_coalescence_recombination;
    RUN_TEST(test_coalesceAtTimePopn_basic);
    RUN_TEST(test_coalesceAtTimePopn_ancestry_merge);
    RUN_TEST(test_coalesceAtTimePopn_overlapping_ancestry);
    RUN_TEST(test_recombineAtTimePopn_basic);
    RUN_TEST(test_recombineAtTimePopn_ancestry_split);
    RUN_TEST(test_geneConversionAtTimePopn_basic);
    RUN_TEST(test_makeGametesMS_mutation_collection);
    RUN_TEST(test_updateActiveMaterial_basic);
    RUN_TEST(test_siteBetweenChunks_basic);
    RUN_TEST(test_pickNodePopn_basic);
    RUN_TEST(test_nodePopnSize_count);
    
    printf("\n========== Running Memory Management Tests ==========\n");
    current_setUp = setUp_memory_management;
    current_tearDown = tearDown_memory_management;
    RUN_TEST(test_initializeBreakPoints_basic);
    RUN_TEST(test_cleanupBreakPoints_basic);
    RUN_TEST(test_ensureBreakPointsCapacity_growth);
    RUN_TEST(test_addBreakPoint_multiple);
    RUN_TEST(test_initializeNodeArrays_basic);
    RUN_TEST(test_cleanupNodeArrays_basic);
    RUN_TEST(test_ensureNodesCapacity_growth);
    RUN_TEST(test_ensureAllNodesCapacity_growth);
    RUN_TEST(test_addNode_with_growth);
    RUN_TEST(test_initializeMuts_basic);
    RUN_TEST(test_initializeMuts_zero_capacity);
    RUN_TEST(test_cleanupMuts_basic);
    RUN_TEST(test_ensureMutsCapacity_growth);
    RUN_TEST(test_mutations_stress_test);
    RUN_TEST(test_reinitialize_breakPoints);
    RUN_TEST(test_cleanup_null_safety);
    RUN_TEST(test_integrated_memory_usage);
    
    return UNITY_END();
}
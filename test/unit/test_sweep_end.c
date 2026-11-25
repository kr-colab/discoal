/*
 * test_sweep_end.c
 *
 * The end of the sweep phase of a sweep from a new mutation.
 *
 * The sweep phase waits for the next event one trajectory step at a time.
 * It stops waiting when an event is due, or when the trajectory has ended
 * (the selected allele is down to one copy). When it stops because the
 * trajectory ended, with no event due, it must not go on to carry out an
 * event: that would be a coalescence or recombination that did not happen.
 *
 * The allele began as a single copy, so when the sweep is over, the
 * lineages still carrying it all descend from that copy, and they coalesce
 * into one.
 *
 * These tests build three lineages by hand, give the sweep phase a short
 * trajectory by hand, and look at the lineages afterwards. They need no
 * replication: with no recombination, the only event that can be drawn in
 * the last step is a coalescence, and the chance of a real event in so
 * short a step (deltaTMod is huge) is negligible whatever the seed.
 *
 * discoalFunctions.c is included here, not linked, so that the tests can
 * set its static state (whether the sweep is from standing variation).
 */

#include "unity.h"
#include "discoal.h"
#include "../../src/core/discoalFunctions.c"
#include "tskitInterface.h"
#include "segmentPool.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

extern void setall(long iseed1, long iseed2);

static float trajectory[8];

/* Three lineages in population 0, a locus of 100 sites, a population of
 * N = 100, and no recombination or gene conversion. */
void setUp(void) {
    setall(11, 12);

    sampleSize = 3;
    sampleNumber = 1;
    nSites = 100;
    npops = 1;
    memset(sampleSizes, 0, sizeof(sampleSizes));
    sampleSizes[0] = 3;
    EFFECTIVE_POPN_SIZE = 100;
    deltaTMod = 1.0e9;      /* a step is 1e-11 of a time unit */
    rho = 0.0;
    leftRho = 0.0;
    my_gamma = 0.0;
    gcMean = 0;
    sweepFromStandingVariation = 0;

    free(currentSize);
    currentSize = (double *) calloc(MAXPOPS, sizeof(double));
    currentSizeConst[0] = 1.0;
    memset(popnSizes, 0, sizeof(popnSizes));
    eventNumber = 0;

    TEST_ASSERT_EQUAL(0, tskit_initialize((double) nSites));
    initialize();
}

void tearDown(void) {
    for (int i = 0; i < nodesCapacity; i++) {
        if (nodes[i] != NULL) {
            if (nodes[i]->ancestryRoot) freeSegmentTree(nodes[i]->ancestryRoot);
            free(nodes[i]);
            nodes[i] = NULL;
        }
    }
    resetSegmentPool();
    cleanupPopLists();
    cleanupSampleNodeIds();
    tskit_cleanup();
    currentTrajectory = NULL;
}

/* The first nBeneficial of the three lineages are in the beneficial class
 * and the rest in the ancestral class. The given trajectory is run through
 * the sweep phase from frequency 0.99. The sweep is already under way
 * (stillSweeping is 1), so the classes set here stand. */
static void runSweep(int nBeneficial, const float *steps, int nSteps) {
    double finalFreq;
    int stillSweeping = 1;

    for (int k = 0; k < 3; k++) nodes[k]->sweepPopn = (k < nBeneficial) ? 1 : 0;
    for (int i = 0; i < nSteps; i++) trajectory[i] = steps[i];
    currentTrajectory = trajectory;
    currentTrajectoryStep = 0;
    totalTrajectorySteps = nSteps;

    sweepPhaseEventsConditionalTrajectory(0.0, 10.0, 0.5, 0.99, &finalFreq,
                                          &stillSweeping, 0.0, 0.0);
}

/* One lineage is left, and the coalescences were recorded in an order that
 * a tree sequence accepts: a parent strictly older than its children.
 * tskit_finalize() sorts the tables and reports it if not. */
static void expectOneLineageAndValidTree(void) {
    char path[] = "/tmp/test_sweep_end_XXXXXX";
    int fd = mkstemp(path);
    TEST_ASSERT_TRUE(fd >= 0);
    close(fd);

    TEST_ASSERT_EQUAL(1, nodePopnSize(0));
    TEST_ASSERT_EQUAL(0, tskit_finalize(path));
    unlink(path);
}

/* The trajectory ends in its last step, at a frequency below one copy
 * (1 / 2N = 0.005), with no event due. The two ancestral lineages could
 * coalesce, and nothing else can happen, so an event carried out at the end
 * would join them: all three lineages must still be there. */
void test_no_event_when_the_trajectory_ends(void) {
    const float steps[] = {0.5f, 0.2f, 0.004f};

    runSweep(1, steps, 3);

    TEST_ASSERT_EQUAL(3, nodePopnSize(0));
    TEST_ASSERT_EQUAL(1, nodePopnSweepSize(0, 1));
    TEST_ASSERT_EQUAL(2, nodePopnSweepSize(0, 0));
}

/* Three beneficial lineages at the end of the sweep coalesce, all three,
 * into one. */
void test_beneficial_lineages_coalesce_at_the_origin(void) {
    const float steps[] = {0.5f, 0.2f, 0.004f};

    runSweep(3, steps, 3);

    expectOneLineageAndValidTree();
}

/* The same when the trajectory ends at frequency 0, as a trajectory read
 * from a file does. */
void test_beneficial_lineages_coalesce_when_the_trajectory_ends_at_zero(void) {
    const float steps[] = {0.5f, 0.2f, 0.0f};

    runSweep(3, steps, 3);

    expectOneLineageAndValidTree();
}

/* The same when an event is due in the last step. A step is long
 * (deltaTMod small), so that the rate of coalescence in the beneficial
 * class is more than 1 in the one step of the trajectory, at frequency
 * 0.004: the event is certain, and it joins one pair. The third lineage
 * must still join them at the origin. */
void test_origin_coalescence_after_an_event_in_the_last_step(void) {
    const float steps[] = {0.004f};
    deltaTMod = 1.0 / (1.5 * EFFECTIVE_POPN_SIZE);   /* a step is 1.5 time units */

    runSweep(3, steps, 1);

    expectOneLineageAndValidTree();
}

/* A sweep from standing variation ends at frequency 0, and nothing is
 * carried out there: it does not coalesce the beneficial class either. */
void test_standing_variation_sweep_ends_without_an_event(void) {
    const float steps[] = {0.5f, 0.2f, 0.004f, 0.0f};
    sweepFromStandingVariation = 1;

    runSweep(1, steps, 4);

    TEST_ASSERT_EQUAL(3, nodePopnSize(0));
}

/* An event that is due in the last step still happens: the rate of
 * coalescence in the ancestral class is more than 1 in the one step, and
 * the two ancestral lineages coalesce. */
void test_event_due_in_the_last_step_still_happens(void) {
    const float steps[] = {0.004f};
    deltaTMod = 1.0 / (1.5 * EFFECTIVE_POPN_SIZE);

    runSweep(1, steps, 1);

    TEST_ASSERT_EQUAL(2, nodePopnSize(0));
    TEST_ASSERT_EQUAL(1, nodePopnSweepSize(0, 1));
    TEST_ASSERT_EQUAL(1, nodePopnSweepSize(0, 0));
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_no_event_when_the_trajectory_ends);
    RUN_TEST(test_beneficial_lineages_coalesce_at_the_origin);
    RUN_TEST(test_beneficial_lineages_coalesce_when_the_trajectory_ends_at_zero);
    RUN_TEST(test_origin_coalescence_after_an_event_in_the_last_step);
    RUN_TEST(test_standing_variation_sweep_ends_without_an_event);
    RUN_TEST(test_event_due_in_the_last_step_still_happens);
    return UNITY_END();
}
#endif

/*
 * test_sweep_switch.c
 *
 * The redraw of a lineage's background (beneficial or ancestral) by a
 * crossover that falls between the selected site and the lineage's
 * material.
 *
 * During a sweep, a lineage that has material on only one side of the
 * selected site is still moved to a new background by a crossover between
 * the site and its material: all of the material goes to the parent that
 * did not carry the site, and that parent's background is drawn from the
 * population. The ordinary recombination weight counts only breakpoints
 * inside the lineage's own span, so these crossovers need their own weight,
 * nodeSwitchMass().
 *
 * Each test below runs a table of cases. The first looks at the weight as a
 * function, and uses no random numbers. The second runs the sweep phase on a
 * hand-made trajectory (see test_sweep_end.c), with a step so long that an
 * event is certain in it, and with lineages chosen so that the only event
 * with a rate is the one under test.
 *
 * discoalFunctions.c is included here, not linked, to reach its static
 * functions.
 */

#include "unity.h"
#include "discoal.h"
#include "../../src/core/discoalFunctions.c"
#include "tskitInterface.h"
#include "segmentPool.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void setall(long iseed1, long iseed2);

static float trajectory[8];

/* A failure message that names the row of a table. */
static char rowName[32];
static const char *row(int i) {
    snprintf(rowName, sizeof(rowName), "row %d", i);
    return rowName;
}

/* Two lineages in population 0, a locus of 100 sites, a population of
 * N = 100, no recombination or gene conversion unless a test turns it on. */
void setUp(void) {
    setall(11, 12);

    sampleSize = 2;
    sampleNumber = 1;
    nSites = 100;
    npops = 1;
    memset(sampleSizes, 0, sizeof(sampleSizes));
    sampleSizes[0] = 2;
    EFFECTIVE_POPN_SIZE = 100;
    deltaTMod = 1.0 / (1.5 * EFFECTIVE_POPN_SIZE);   /* a step is 1.5 time units */
    rho = 0.0;
    leftRho = 0.0;
    my_gamma = 0.0;
    gcMean = 0;

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

/* A lineage record with only the fields the weights read. */
static rootedNode spanOnly(int first, int last, int nanc) {
    rootedNode n;
    memset(&n, 0, sizeof(n));
    n.lLim = first;
    n.rLim = last;
    n.nancSites = nanc;
    return n;
}

/* A lineage with material only at sites first .. last. */
static void giveMaterial(rootedNode *node, int first, int last) {
    freeSegmentTree(node->ancestryRoot);
    node->ancestryRoot = newSegment(first, last + 1, get_tskit_node_id(node), NULL, NULL);
    SegmentSpan span = {last - first + 1, first, last};
    setLineageSpan(node, &span);
}

/* ---- The weights, as functions ---- */

/* The number of breakpoints from the selected site to the nearer end of
 * the span; none if the span contains the site. A span of a single site has
 * no breakpoints of its own, but crossovers between it and the site count.
 * The site is a fraction of the locus, and which breakpoints are to its left
 * follows the test the recombination code uses, sweepSite < (float) x / nSites:
 * with the site at 0.3 and 1000 sites, 0.3 is less than the float 300 / 1000,
 * so the site is left of breakpoint 300 and a span starting at 300 has one
 * breakpoint between it and the site. */
void test_switch_mass_of_a_span(void) {
    static const struct { int sites, first, last, nanc; double site; int mass; } rows[] = {
        /* sites first last nanc site  mass */
        {   100,   60,   70,  11, 0.5,   10 },   /* right of the site */
        {   100,   10,   20,  11, 0.5,   30 },   /* left of it */
        {   100,   40,   60,  21, 0.5,    0 },   /* containing it */
        {   100,   60,   70,   0, 0.5,    0 },   /* no material */
        {   100,   60,   60,   1, 0.5,   10 },   /* a single site, right */
        {   100,   20,   20,   1, 0.5,   30 },   /* a single site, left */
        {   100,    5,    9,   5, 0.0,    5 },   /* the site at the left end */
        {   100,    0,    9,  10, 0.0,    0 },   /* a span that holds site 0 */
        {   100,   10,   20,  11, 0.99,  78 },   /* near the right end */
        {  1000,  300,  310,  11, 0.3,    1 },   /* the site test, above */
    };
    for (int i = 0; i < (int) (sizeof(rows) / sizeof(rows[0])); i++) {
        nSites = rows[i].sites;
        rootedNode n = spanOnly(rows[i].first, rows[i].last, rows[i].nanc);
        TEST_ASSERT_EQUAL_INT_MESSAGE(rows[i].mass, (int) nodeSwitchMass(&n, rows[i].site), row(i));
    }
    nSites = 100;
}

/* ---- The events, through the sweep phase ---- */

/* Lineage 0 has a single site of material at the given position, in class
 * cls; lineage 1 has a single site at the selected site (site 50) in the
 * other class. Neither has breakpoints of its own, and no two lineages share
 * a class, so nothing else can happen. The sweep is from standing
 * variation, so that its end does not coalesce the beneficial lineages (see
 * test_sweep_end.c), which would hide the move. The trajectory has one step,
 * and the step is so long that the rate of the event under test is above 1
 * in it: the event is certain. The redraw always moves the lineage. With no
 * recombination, nothing happens. */
void test_a_redraw_moves_a_lineage_through_the_sweep_phase(void) {
    static const struct { double rho; int position, cls, moves; } rows[] = {
        /* rho position cls moves */
        { 100,   70,     1,  1 },                /* a crossover, lineage right of the site */
        { 100,   30,     0,  1 },                /* left of it */
        {   0,   70,     1,  0 },                /* no recombination: nothing moves */
    };
    for (int i = 0; i < (int) (sizeof(rows) / sizeof(rows[0])); i++) {
        double finalFreq;
        int stillSweeping = 1;
        const float x = 0.5f;

        rho = rows[i].rho;
        sweepFromStandingVariation = 1;
        giveMaterial(nodes[0], rows[i].position, rows[i].position);
        giveMaterial(nodes[1], 50, 50);
        nodes[0]->sweepPopn = rows[i].cls;
        nodes[1]->sweepPopn = 1 - rows[i].cls;
        trajectory[0] = x;
        currentTrajectory = trajectory;
        currentTrajectoryStep = 0;
        totalTrajectorySteps = 1;

        sweepPhaseEventsConditionalTrajectory(0.0, 10.0, 0.5, 0.99, &finalFreq,
                                              &stillSweeping, 0.0, 0.0);

        TEST_ASSERT_EQUAL_INT_MESSAGE(2, nodePopnSize(0), row(i));
        TEST_ASSERT_EQUAL_INT_MESSAGE(rows[i].moves ? 1 - rows[i].cls : rows[i].cls,
                                      nodes[0]->sweepPopn, row(i));

        /* The next row starts from a fresh pair of lineages. */
        tearDown();
        setUp();
    }
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_switch_mass_of_a_span);
    RUN_TEST(test_a_redraw_moves_a_lineage_through_the_sweep_phase);
    return UNITY_END();
}
#endif

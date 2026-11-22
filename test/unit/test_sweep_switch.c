/*
 * test_sweep_switch.c
 *
 * The redraw of a lineage's background (beneficial or ancestral) by a
 * crossover or a gene conversion tract that falls between the selected site
 * and the lineage's material.
 *
 * During a sweep, a lineage that has material on only one side of the
 * selected site is still moved to a new background by a crossover between
 * the site and its material: all of the material goes to the parent that
 * did not carry the site, and that parent's background is drawn from the
 * population. The ordinary recombination weight counts only breakpoints
 * inside the lineage's own span, so these crossovers need their own weight
 * (nodeSwitchMass()), and gene conversion tracts that separate the site from
 * the material need theirs (nodeGcSwitchWeight()).
 *
 * Each test below runs a table of cases. The first two look at the weights
 * as functions, and the third at convertTract(); none uses random numbers
 * that matter. The fourth runs the sweep phase on a hand-made trajectory
 * (see test_sweep_end.c), with a step so long that an event is certain in
 * it, and with lineages chosen so that the only event with a rate is the one
 * under test.
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

/* A lineage with material at sites a1 .. a2 and b1 .. b2, and a gap
 * between them. */
static void giveMaterialWithGap(rootedNode *node, int a1, int a2, int b1, int b2) {
    AncestrySegment *left = newSegment(a1, a2 + 1, get_tskit_node_id(node), NULL, NULL);
    AncestrySegment *right = newSegment(b1, b2 + 1, get_tskit_node_id(node), NULL, NULL);
    freeSegmentTree(node->ancestryRoot);
    node->ancestryRoot = mergeAncestryTrees(left, right, get_tskit_node_id(node));
    SegmentSpan span = {(a2 - a1 + 1) + (b2 - b1 + 1), a1, b2};
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

/* The weight of a lineage for gene conversion tracts, which are geometric
 * with mean gcMean: gcMean (1 - pc^d) (1 + pc^w), pc = 1 - 1 / gcMean, d the
 * crossover weight above and w the number of breakpoints in the span. A
 * mean of 10, d = 2, w = 3 gives 10 x (1 - 0.81) x (1 + 0.729) = 3.2851. With
 * tracts of one site it is 2 for a single site of material and 1 for a
 * longer span. For tracts much longer than the locus it tends to twice the
 * crossover weight: such tracts act as crossovers at twice the conversion
 * rate, which testing/gc_crossover_check.py relies on. */
void test_gc_switch_weight(void) {   /* gene conversion */
    static const struct { int first, last, nanc, mean; double weight, tol; } rows[] = {
        /* first last nanc mean        weight  tolerance */
        {    60,   60,   1,  1,          2.0,  1e-12 },    /* tracts of one site */
        {    60,   65,   6,  1,          1.0,  1e-12 },
        {    40,   60,  21,  1,          0.0,  1e-12 },    /* containing the site */
        {    52,   55,   4, 10,       3.2851,  1e-9  },
        {    60,   65,   6, 1000000000, 20.0,  1e-3  },    /* 2 x 10 */
    };
    for (int i = 0; i < (int) (sizeof(rows) / sizeof(rows[0])); i++) {
        gcMean = rows[i].mean;
        rootedNode n = spanOnly(rows[i].first, rows[i].last, rows[i].nanc);
        TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(rows[i].tol, rows[i].weight,
                                          nodeGcSwitchWeight(&n, 0.5), row(i));
    }
}

/* ---- A gene conversion tract that does not cut the material ---- */

/* A lineage in the beneficial class, with material at sites 0 .. 49 or, for
 * the last two rows, at 0 .. 10 and 80 .. 90 with a gap between. A tract
 * [start, end) and the selected site are given. A tract that covers the site
 * and none of the material, or all of the material and not the site,
 * redraws the background; the class frequency is passed as 0, so the
 * lineage then always moves to the other class. A tract that covers both, or
 * neither, changes nothing. */
void test_tract_that_does_not_cut_the_material(void) {   /* gene conversion */
    static const struct { int gap, start, end; double site; int classAfter; } rows[] = {
        /* gap start end  site  class after */
        {   0,   60,  70, 0.65,  0 },   /* the site, no material */
        {   0,    0,  60, 0.8,   0 },   /* all the material, not the site */
        {   0,    0,  90, 0.8,   1 },   /* both */
        {   0,   60,  70, 0.8,   1 },   /* neither */
        {   1,   40,  60, 0.5,   0 },   /* the site, in the gap */
        {   1,   40,  60, 0.95,  1 },   /* a tract in the gap, not over the site */
    };
    for (int i = 0; i < (int) (sizeof(rows) / sizeof(rows[0])); i++) {
        if (rows[i].gap) giveMaterialWithGap(nodes[0], 0, 10, 80, 90);
        else giveMaterial(nodes[0], 0, 49);
        nodes[0]->sweepPopn = 1;
        nodes[1]->sweepPopn = 0;
        sweepPopnSizes[1] = nodePopnSweepSize(0, 1);
        sweepPopnSizes[0] = nodePopnSweepSize(0, 0);

        convertTract(1.0, 0, nodes[0], rows[i].start, rows[i].end, 1, rows[i].site, 0.0);

        TEST_ASSERT_EQUAL_INT_MESSAGE(rows[i].classAfter, nodes[0]->sweepPopn, row(i));
        TEST_ASSERT_EQUAL_INT_MESSAGE(nodePopnSweepSize(0, 1), sweepPopnSizes[1], row(i));
        TEST_ASSERT_EQUAL_INT_MESSAGE(nodePopnSweepSize(0, 0), sweepPopnSizes[0], row(i));
    }
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
 * recombination or conversion, nothing happens. */
void test_a_redraw_moves_a_lineage_through_the_sweep_phase(void) {
    static const struct { double rho, gamma; int mean, position, cls, moves; } rows[] = {
        /* rho gamma mean position cls moves */
        { 100,   0,   0,   70,     1,  1 },                /* a crossover, lineage right of the site */
        { 100,   0,   0,   30,     0,  1 },                /* left of it */
        {   0, 100,   5,   70,     0,  1 },                /* a conversion tract */    /* gene conversion */
        {   0, 100,   5,   30,     1,  1 },                                           /* gene conversion */
        {   0,   0,   0,   70,     1,  0 },                /* neither: nothing moves */
    };
    for (int i = 0; i < (int) (sizeof(rows) / sizeof(rows[0])); i++) {
        double finalFreq;
        int stillSweeping = 1;
        const float x = 0.5f;

        rho = rows[i].rho;
        my_gamma = rows[i].gamma;
        gcMean = rows[i].mean;
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
    RUN_TEST(test_gc_switch_weight);
    RUN_TEST(test_tract_that_does_not_cut_the_material);
    RUN_TEST(test_a_redraw_moves_a_lineage_through_the_sweep_phase);
    return UNITY_END();
}
#endif

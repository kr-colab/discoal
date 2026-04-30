/*
 * test_globals.c
 *
 * Provides minimal global variables required by discoal functions
 * when they are tested in isolation. These globals are normally
 * defined in discoal_multipop.c but need to be provided for unit tests.
 */

#include <stdlib.h>
#include <string.h>
#include "discoal.h"

// Random number generator seeds
long seed1 = 12345;
long seed2 = 67890;

// Current population size (used by some functions)
double *currentSize = NULL;

// Event management globals
event *events = NULL;
int eventNumber = 0;
int eventsCapacity = 0;

// Stub implementation of ensureEventsCapacity for testing
void ensureEventsCapacity() {
    if (events == NULL) {
        eventsCapacity = 100;
        events = (event *)malloc(eventsCapacity * sizeof(event));
    } else if (eventNumber >= eventsCapacity - 1) {
        eventsCapacity *= 2;
        events = (event *)realloc(events, eventsCapacity * sizeof(event));
    }
}

/*
 * reset_config_globals
 *
 * Zero every discoal global that the parse_* functions in
 * src/core/configInterface.c write into. Tests that call
 * apply_yaml_config invoke this from setUp() so each test starts
 * from a known state.
 *
 * The list below corresponds to the `extern` declarations inside
 * the parse_* functions in configInterface.c. If a parse_*
 * function gains a new global write, add the corresponding reset
 * here too — otherwise the next test's setUp will see the previous
 * test's value.
 *
 * Notes:
 *   - currentSize and events are heap-allocated and managed by the
 *     setUp / tearDown code in test_config_interface.c, not here.
 *   - sampleSizes, popnSizes, migMatConst, and tskitOutputFilename
 *     are fixed-size arrays declared in discoal.h, so sizeof(...)
 *     gives the array length in bytes (the -fcommon build setting
 *     keeps the declarations and definitions consistent).
 */
void reset_config_globals(void) {
    extern int sampleSize, sampleNumber, nSites;
    extern int sampleSizes[MAXPOPS];
    extern int popnSizes[MAXPOPS];
    extern int npops;
    extern int migFlag;
    extern int EFFECTIVE_POPN_SIZE;
    extern int gcMean, finiteOutputFlag;
    extern int recurSweepMode, partialSweepMode, softSweepMode;
    extern int tskitOutputMode, minimalTreeSeq, hidePartialSNP;
    extern double theta, rho;
    extern double gammaCoRatio, my_gamma, gammaCoRatioMode;
    extern double alpha, sweepSite, tau, f0, uA;
    extern double partialSweepFinalFreq, recurSweepRate;
    extern double tDiv;
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern char sweepMode;
    extern char tskitOutputFilename[1024];

    sampleSize = 0;
    sampleNumber = 0;
    nSites = 0;
    npops = 0;
    migFlag = 0;
    EFFECTIVE_POPN_SIZE = 0;
    gcMean = 0;
    finiteOutputFlag = 0;
    recurSweepMode = 0;
    partialSweepMode = 0;
    softSweepMode = 0;
    tskitOutputMode = 0;
    minimalTreeSeq = 0;
    hidePartialSNP = 0;

    theta = 0.0;
    rho = 0.0;
    gammaCoRatio = 0.0;
    my_gamma = 0.0;
    gammaCoRatioMode = 0.0;
    alpha = 0.0;
    sweepSite = 0.0;
    tau = 0.0;
    f0 = 0.0;
    uA = 0.0;
    partialSweepFinalFreq = 0.0;
    recurSweepRate = 0.0;
    tDiv = 0.0;

    sweepMode = '\0';

    seed1 = 0;
    seed2 = 0;

    memset(sampleSizes, 0, sizeof(sampleSizes));
    memset(popnSizes, 0, sizeof(popnSizes));
    memset(migMatConst, 0, sizeof(migMatConst));
    memset(tskitOutputFilename, 0, sizeof(tskitOutputFilename));
}

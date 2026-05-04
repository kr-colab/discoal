#ifndef SHAPES_H
#define SHAPES_H

#include "discoal.h"

/* Evaluate the size/rate of a shape at absolute time t. */
double sizeAt(int popID, double t);
double migAt(int srcPopID, int dstPopID, double t);

/* Evaluate the integrated hazard for k lineages from t0 over duration T.
 * Used in tests to verify drawWaitingTime against numerical quadrature. */
double integratedHazardSize(int popID, double t0, double T, int k);
double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k);

/* Draw a waiting time T such that integratedHazard(t0, T, k) == xi.
 * Returns -1.0 if xi exceeds the integrated hazard over [t0, infinity)
 * or if the linear shape would drive N(t) to zero before the draw resolves. */
double drawWaitingTimeSize(int popID, double t0, double xi, int k);
double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k);

/* Initialize popShape[] and migShape[][] from the current values in
 * currentSize[] and migMatConst[][]. Sets all shapes to SHAPE_CONSTANT
 * anchored at t=0. Called once per replicate from initialize(). */
void initializeShapesFromGlobals(void);

/* Returns 1 if every popShape[i] for i in [0, npops) and every migShape[i][j]
 * for i,j in [0, npops) has type SHAPE_CONSTANT. Returns 0 otherwise.
 * Used by the inner-loop sampler to dispatch between the constant-rate
 * Exp(total) path and the non-homogeneous Poisson process per-component
 * path. */
int allShapesConstant(void);

/* Integrated size ratio: returns int_{t0}^{t0+T} sizeAt(popID, s) ds.
 * Sister to integratedHazardSize but integrates sizeAt directly rather than
 * 1/sizeAt. Used by detSweepFreqGeneral to compute the integrated selection
 * coefficient A(tau) = alpha * integratedSizeRatio(0, sweep_start, tau). */
double integratedSizeRatio(int popID, double t0, double T);

/* Preflight check on the parsed events array. Walks events in time order,
 * tracks each population's size shape state and each (src,dst) pair's
 * migration shape state, and returns -1 if either trajectory would reach
 * an invalid value before the next shape-changing event for that
 * pop/pair:
 *   - size:      'n' with popnSize <= 0; or linear-growth shape whose
 *                zero crossing precedes the next event for that pop.
 *                (Sizes must be strictly positive — coalescent rate would
 *                blow up.)
 *   - migration: 'm' with rate < 0; initial migMatConst[src][dst] < 0;
 *                or linear migration shape whose zero crossing precedes
 *                the next migration event for that pair. (Rates must be
 *                non-negative — zero is fine, negative would corrupt the
 *                aggregated mRate / dest-pop sampling in the inner loop.)
 * Prints an explanatory error to stderr on rejection. Returns 0 on
 * success. */
struct event;  /* forward decl; full type in discoal.h */
int validateShapeTrajectories(struct event *events, int eventNumber);

#endif /* SHAPES_H */

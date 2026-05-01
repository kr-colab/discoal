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

#endif /* SHAPES_H */

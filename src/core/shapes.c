#include <math.h>
#include <stdlib.h>
#include "shapes.h"

/* sizeAt -- evaluate the size shape for population popID at absolute time t.
 *
 * Sign convention (msprime-aligned, applies uniformly across shape types):
 * rate_param is the per-generation FORWARD-time rate of change. The
 * coalescent simulator runs backward in time, so positive rate_param
 * means the past held a smaller value.
 *   SHAPE_CONSTANT:    N(t) = anchor_value (rate_param unused)
 *   SHAPE_EXPONENTIAL: N(t) = anchor_value * exp(-rate_param * (t - anchor_time))
 *   SHAPE_LINEAR:      N(t) = anchor_value - rate_param * (t - anchor_time)
 */
double sizeAt(int popID, double t) {
    Shape *s = &popShape[popID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        case SHAPE_EXPONENTIAL:
            return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:
            return s->anchor_value - s->rate_param * (t - s->anchor_time);
        default:
            return 0.0;  /* other shapes implemented in subsequent tasks */
    }
}

double migAt(int srcPopID, int dstPopID, double t) {
    (void)srcPopID; (void)dstPopID; (void)t;
    return 0.0;
}

double integratedHazardSize(int popID, double t0, double T, int k) {
    (void)popID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double drawWaitingTimeSize(int popID, double t0, double xi, int k) {
    (void)popID; (void)t0; (void)xi; (void)k;
    return -1.0;
}

double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)xi; (void)k;
    return -1.0;
}

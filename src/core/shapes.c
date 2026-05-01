#include <math.h>
#include <stdlib.h>
#include "shapes.h"

/* sizeAt -- evaluate the size shape for population popID at absolute time t.
 *
 * Sign convention: rate_param holds the FORWARD-time growth rate.
 * Coalescent simulation runs backward in time, so:
 *   SHAPE_EXPONENTIAL: N(t) = anchor_value * exp(-rate_param * (t - anchor_time))
 *     (forward growth alpha > 0 means backward decline)
 *   SHAPE_LINEAR:      N(t) = anchor_value + rate_param * (t - anchor_time)
 *     (gamma > 0 means N grows as t increases backward in time)
 *   SHAPE_CONSTANT:    rate_param is unused.
 */
double sizeAt(int popID, double t) {
    Shape *s = &popShape[popID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        case SHAPE_EXPONENTIAL:
            return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:
            return s->anchor_value + s->rate_param * (t - s->anchor_time);
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

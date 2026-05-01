#include <math.h>
#include <stdlib.h>
#include "shapes.h"

/* sizeAt / migAt -- evaluate a shape (size or migration rate) at absolute
 * time t. sizeAt indexes popShape[popID]; migAt indexes migShape[src][dst].
 *
 * Sign convention (msprime-aligned, applies uniformly across shape types
 * and to both sizes and migration rates): rate_param is the per-generation
 * FORWARD-time rate of change. The coalescent simulator runs backward in
 * time, so positive rate_param means the past held a smaller value.
 *   SHAPE_CONSTANT:    value(t) = anchor_value (rate_param unused)
 *   SHAPE_EXPONENTIAL: value(t) = anchor_value * exp(-rate_param * (t - anchor_time))
 *   SHAPE_LINEAR:      value(t) = anchor_value - rate_param * (t - anchor_time)
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
    Shape *s = &migShape[srcPopID][dstPopID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        case SHAPE_EXPONENTIAL:
            return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:
            return s->anchor_value - s->rate_param * (t - s->anchor_time);
        default:
            return 0.0;
    }
}

double integratedHazardSize(int popID, double t0, double T, int k) {
    Shape *s = &popShape[popID];
    double pairs = (double)k * (k - 1) / 2.0;
    if (pairs == 0.0) return 0.0;
    double N0 = sizeAt(popID, t0);
    switch (s->type) {
        case SHAPE_CONSTANT:
            return pairs * T / N0;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return pairs * T / N0;
            return pairs * expm1(a * T) / (N0 * a);
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            if (g == 0.0) return pairs * T / N0;
            double frac = g * T / N0;
            if (frac >= 1.0) return INFINITY;  /* end <= 0 case */
            return -pairs * log1p(-frac) / g;
        }
        default:
            return 0.0;
    }
}

double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k) {
    Shape *s = &migShape[srcPopID][dstPopID];
    double m0 = migAt(srcPopID, dstPopID, t0);
    if (k <= 0 || m0 == 0.0) return 0.0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return k * m0 * T;
        case SHAPE_EXPONENTIAL: {
            double b = s->rate_param;
            if (b == 0.0) return k * m0 * T;
            return k * m0 * -expm1(-b * T) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            return k * (m0 * T - 0.5 * d * T * T);
        }
        default:
            return 0.0;
    }
}

double drawWaitingTimeSize(int popID, double t0, double xi, int k) {
    if (k < 2) return -1.0;
    Shape *s = &popShape[popID];
    double pairs = (double)k * (k - 1) / 2.0;
    double N0 = sizeAt(popID, t0);
    if (N0 <= 0.0) return -1.0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return xi * N0 / pairs;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return xi * N0 / pairs;
            return log1p(N0 * a * xi / pairs) / a;
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            if (g == 0.0) return xi * N0 / pairs;
            double T = (N0 / g) * (1.0 - exp(-g * xi / pairs));
            if (g > 0.0) {
                double T_cross = N0 / g;
                if (T >= T_cross) T = nextafter(T_cross, 0.0);
            }
            return T;
        }
        default:
            return -1.0;
    }
}

double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)xi; (void)k;
    return -1.0;
}

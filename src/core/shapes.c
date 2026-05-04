#include <math.h>
#include <stdlib.h>
#include "shapes.h"

/* sizeAt / migAt -- evaluate a shape (size or migration rate) at absolute
 * time t. sizeAt indexes popShape[popID]; migAt indexes migShape[src][dst].
 *
 * Sign convention: rate_param is the per-generation
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
            return 0.0;
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
    if (k <= 0) return -1.0;
    Shape *s = &migShape[srcPopID][dstPopID];
    double m0 = migAt(srcPopID, dstPopID, t0);
    if (m0 <= 0.0) return -1.0;
    double km0 = k * m0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return xi / km0;
        case SHAPE_EXPONENTIAL: {
            double b = s->rate_param;
            if (b == 0.0) return xi / km0;
            double arg = b * xi / km0;
            if (arg >= 1.0) return -1.0;  /* unreachable */
            return -log1p(-arg) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            if (d == 0.0) return xi / km0;
            double disc = m0 * m0 - 2.0 * d * xi / k;
            if (disc < 0.0) return -1.0;
            double T = (m0 - sqrt(disc)) / d;
            if (T < 0.0) return -1.0;
            return T;
        }
        default:
            return -1.0;
    }
}

void initializeShapesFromGlobals(void) {
    extern double *currentSize;
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern int npops;

    for (int i = 0; i < npops; i++) {
        popShape[i].type = SHAPE_CONSTANT;
        popShape[i].anchor_value = currentSize[i];
        popShape[i].rate_param = 0.0;
        popShape[i].anchor_time = 0.0;
        for (int j = 0; j < npops; j++) {
            migShape[i][j].type = SHAPE_CONSTANT;
            migShape[i][j].anchor_value = migMatConst[i][j];
            migShape[i][j].rate_param = 0.0;
            migShape[i][j].anchor_time = 0.0;
        }
    }
}

int allShapesConstant(void) {
    extern int npops;
    for (int i = 0; i < npops; i++) {
        if (popShape[i].type != SHAPE_CONSTANT) return 0;
        for (int j = 0; j < npops; j++) {
            if (migShape[i][j].type != SHAPE_CONSTANT) return 0;
        }
    }
    return 1;
}

double integratedSizeRatio(int popID, double t0, double T) {
    Shape *s = &popShape[popID];
    double N0 = sizeAt(popID, t0);
    switch (s->type) {
        case SHAPE_CONSTANT:
            return N0 * T;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return N0 * T;
            return N0 * (1.0 - exp(-a * T)) / a;
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            return N0 * T - 0.5 * g * T * T;
        }
        default:
            return 0.0;
    }
}

static double sizeFromShape(Shape *s, double t) {
    switch (s->type) {
        case SHAPE_CONSTANT:    return s->anchor_value;
        case SHAPE_EXPONENTIAL: return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:      return s->anchor_value - s->rate_param * (t - s->anchor_time);
        default:                return 0.0;
    }
}

int validateShapeTrajectories(struct event *events, int eventNumber) {
    extern double *currentSize;
    extern double currentSizeConst[MAXPOPS];

    Shape state[MAXPOPS];
    int merged[MAXPOPS];
    for (int p = 0; p < MAXPOPS; p++) {
        double init = (currentSize != NULL && currentSize[p] > 0.0)
                          ? currentSize[p]
                          : (currentSizeConst[p] > 0.0 ? currentSizeConst[p] : 1.0);
        state[p].type = SHAPE_CONSTANT;
        state[p].anchor_value = init;
        state[p].rate_param = 0.0;
        state[p].anchor_time = 0.0;
        merged[p] = 0;
    }

    for (int i = 0; i < eventNumber; i++) {
        struct event *e = &events[i];
        char t = e->type;
        int p = e->popID;
        if (p < 0 || p >= MAXPOPS) continue;
        if (merged[p]) continue;

        /* Pop-p's prior linear shape must not have crossed zero by e->time. */
        if (state[p].type == SHAPE_LINEAR && state[p].rate_param > 0.0) {
            double t_cross = state[p].anchor_time + state[p].anchor_value / state[p].rate_param;
            if (t_cross < e->time) {
                fprintf(stderr,
                    "Error: linear-growth trajectory drives population %d size to zero "
                    "at internal time %g (before the next event for that population at "
                    "time %g). Adjust the linear rate or shorten the interval.\n",
                    p, t_cross, e->time);
                return -1;
            }
        }

        if (t == 'p') { merged[p] = 1; continue; }
        if (t != 'n' && t != 'g' && t != 'l') continue;

        double size_at_event = sizeFromShape(&state[p], e->time);
        if (size_at_event <= 0.0) {
            fprintf(stderr,
                "Error: population %d size becomes non-positive (%g) at internal time %g.\n",
                p, size_at_event, e->time);
            return -1;
        }

        if (t == 'n') {
            if (e->popnSize <= 0.0) {
                fprintf(stderr,
                    "Error: -en sets population %d size to %g at time %g; size must be "
                    "strictly positive.\n", p, e->popnSize, e->time);
                return -1;
            }
            state[p].type = SHAPE_CONSTANT;
            state[p].anchor_value = e->popnSize;
            state[p].rate_param = 0.0;
            state[p].anchor_time = e->time;
        } else if (t == 'g') {
            state[p].type = SHAPE_EXPONENTIAL;
            state[p].anchor_value = size_at_event;
            state[p].rate_param = e->popnSize;  /* alpha */
            state[p].anchor_time = e->time;
        } else if (t == 'l') {
            state[p].type = SHAPE_LINEAR;
            state[p].anchor_value = size_at_event;
            state[p].rate_param = e->popnSize;  /* gamma */
            state[p].anchor_time = e->time;
        }
    }
    return 0;
}

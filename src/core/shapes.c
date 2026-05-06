#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "shapes.h"

/* All current shapes (CONSTANT/EXPONENTIAL/LINEAR) are handled explicitly in
 * each switch below. The default branches are unreachable; if a new shape
 * type is added without updating one of these switches, abort loudly rather
 * than silently returning a sentinel that would propagate as wrong rates. */
static void unknownShape(const char *fn, int type) {
    fprintf(stderr, "discoal: internal error: unhandled shape type %d in %s\n",
            type, fn);
    exit(1);
}

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
    }
    unknownShape("sizeAt", s->type);
    return 0.0;
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
    }
    unknownShape("migAt", s->type);
    return 0.0;
}

/* Integrated coalescent hazard for k lineages over [t0, t0+T].
 *
 *   H(T) = integral_{0}^{T} (k choose 2) / N(t0 + s)  ds
 *
 * Closed forms below.  Let N0 = N(t0); recall the sign convention that
 * rate_param is the FORWARD-time per-generation rate of change, so going
 * backward in time N(t0+s) = N0 * exp(-a*s) (exponential) or
 * N(t0+s) = N0 - g*s (linear).
 *
 *   CONSTANT:    H(T) = pairs * T / N0
 *   EXPONENTIAL: integrand = pairs * exp(a*s) / N0
 *                H(T) = pairs * (exp(a*T) - 1) / (N0 * a)
 *   LINEAR:      integrand = pairs / (N0 - g*s)
 *                H(T) = -(pairs/g) * log(1 - g*T/N0)
 *                Diverges when g*T/N0 -> 1 (the linear shape would drive
 *                N(t) to zero); we report INFINITY so the caller treats it
 *                as "no event before the trajectory crashes".
 */
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
            /* expm1(a*T) = exp(a*T) - 1, accurate near a*T = 0. */
            return pairs * expm1(a * T) / (N0 * a);
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            if (g == 0.0) return pairs * T / N0;
            double frac = g * T / N0;
            if (frac >= 1.0) return INFINITY;  /* N(t0+T) <= 0 */
            /* log1p(-frac) = log(1 - frac), accurate near frac = 0. */
            return -pairs * log1p(-frac) / g;
        }
    }
    unknownShape("integratedHazardSize", s->type);
    return 0.0;
}

/* Integrated migration hazard for k lineages over [t0, t0+T].
 *
 *   H(T) = integral_{0}^{T} k * m(t0 + s)  ds
 *
 * Let m0 = m(t0).  Closed forms:
 *
 *   CONSTANT:    H(T) = k * m0 * T
 *   EXPONENTIAL: integrand = k * m0 * exp(-b*s)
 *                H(T) = k * m0 * (1 - exp(-b*T)) / b
 *   LINEAR:      integrand = k * (m0 - d*s) -- the antiderivative is
 *                quadratic in T.  H_raw(T) = k*(m0*T - 0.5*d*T^2).
 *                For d > 0 the migration rate hits zero at T_zero = m0/d
 *                and is physically clamped at zero for s > T_zero, so the
 *                integral plateaus at H(T_zero) = k*m0^2/(2*d).  Without
 *                the clamp H_raw would peak at T_zero and then *decrease*
 *                (eventually going negative at T = 2*m0/d), which is
 *                meaningless for an integrated hazard.  We clamp T to
 *                T_zero before evaluating.
 */
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
            /* -expm1(-b*T) = 1 - exp(-b*T), accurate near b*T = 0. */
            return k * m0 * -expm1(-b * T) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            double T_eff = T;
            if (d > 0.0) {
                double T_zero = m0 / d;
                if (T_eff > T_zero) T_eff = T_zero;
            }
            return k * (m0 * T_eff - 0.5 * d * T_eff * T_eff);
        }
    }
    unknownShape("integratedHazardMig", s->type);
    return 0.0;
}

/* Inversion of integratedHazardSize: given xi ~ Exp(1), solve H(T) = xi
 * for T.  Each branch is the algebraic inverse of the corresponding
 * integratedHazardSize case.
 *
 *   CONSTANT:    T = xi * N0 / pairs
 *   EXPONENTIAL: xi = pairs * (exp(a*T) - 1) / (N0 * a)
 *                =>  T = log1p(N0 * a * xi / pairs) / a
 *   LINEAR:      xi = -(pairs/g) * log(1 - g*T/N0)
 *                =>  T = (N0/g) * (1 - exp(-g*xi/pairs))
 *                For g > 0 the solution is bounded above by T_cross = N0/g
 *                (where N -> 0); cap T just shy of T_cross to avoid
 *                returning a value that the caller would interpret as a
 *                physical event past the trajectory's crash time.
 */
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
    }
    unknownShape("drawWaitingTimeSize", s->type);
    return -1.0;
}

/* Inversion of integratedHazardMig: solve H(T) = xi for T.  Returns -1.0
 * when xi exceeds the integrated hazard over [0, infinity) for the given
 * shape (i.e. no event ever, given the rate trajectory).
 *
 *   CONSTANT:    T = xi / (k*m0)
 *   EXPONENTIAL: xi = k*m0*(1 - exp(-b*T))/b
 *                =>  T = -log(1 - b*xi/(k*m0)) / b
 *                Unreachable when b*xi/(k*m0) >= 1 (the migration rate
 *                decays to 0 in less hazard than xi requires).
 *   LINEAR:      xi = k*(m0*T - 0.5*d*T^2)  =>  quadratic in T, take the
 *                smaller positive root: T = (m0 - sqrt(m0^2 - 2*d*xi/k))/d.
 *                Discriminant negative <=> xi > k*m0^2/(2*d), i.e. xi
 *                exceeds the plateau hazard accumulated by the time the
 *                rate hits zero.
 */
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
            if (arg >= 1.0) return -1.0;  /* xi exceeds k*m0/b, the asymptote */
            return -log1p(-arg) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            if (d == 0.0) return xi / km0;
            double disc = m0 * m0 - 2.0 * d * xi / k;
            if (disc < 0.0) return -1.0;  /* xi exceeds the plateau */
            double T = (m0 - sqrt(disc)) / d;
            if (T < 0.0) return -1.0;
            return T;
        }
    }
    unknownShape("drawWaitingTimeMig", s->type);
    return -1.0;
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

/* Integral of N(t) (not 1/N(t)) from t0 to t0+T.  Used by detSweepFreqGeneral
 * to track A(tau) = alpha * integratedSizeRatio over the sweep window.
 *
 *   CONSTANT:    N0 * T
 *   EXPONENTIAL: N0 * (1 - exp(-a*T)) / a
 *   LINEAR:      N0*T - 0.5*g*T^2  (antiderivative of N0 - g*s)
 *
 * The linear branch is left unguarded for now: callers track the sweep
 * trajectory in small dt steps where the polynomial is well-behaved, and
 * validateShapeTrajectories rejects any model where N(t) crosses zero
 * before the next event.
 */
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
    }
    unknownShape("integratedSizeRatio", s->type);
    return 0.0;
}

static double sizeFromShape(Shape *s, double t) {
    switch (s->type) {
        case SHAPE_CONSTANT:    return s->anchor_value;
        case SHAPE_EXPONENTIAL: return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:      return s->anchor_value - s->rate_param * (t - s->anchor_time);
    }
    unknownShape("sizeFromShape", s->type);
    return 0.0;
}

int validateShapeTrajectories(struct event *events, int eventNumber) {
    extern double *currentSize;
    extern double currentSizeConst[MAXPOPS];
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern int npops;

    Shape state[MAXPOPS];
    int merged[MAXPOPS];

    /* Per-pair migration shape state. MAXPOPS^2 Shapes is too large for
     * the stack (~580KB at MAXPOPS=121), so heap-allocate and index by
     * src * npops + dst. */
    int npops_eff = (npops > 0 && npops <= MAXPOPS) ? npops : MAXPOPS;
    Shape *migState = (Shape *)calloc((size_t)npops_eff * npops_eff, sizeof(Shape));
    if (migState == NULL) {
        fprintf(stderr, "discoal: out of memory in validateShapeTrajectories\n");
        return -1;
    }

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
    for (int src = 0; src < npops_eff; src++) {
        for (int dst = 0; dst < npops_eff; dst++) {
            double m_init = migMatConst[src][dst];
            if (m_init < 0.0) {
                fprintf(stderr,
                    "Error: initial migration rate from population %d to %d is %g; "
                    "rate must be non-negative.\n", src, dst, m_init);
                free(migState);
                return -1;
            }
            Shape *ms = &migState[src * npops_eff + dst];
            ms->type = SHAPE_CONSTANT;
            ms->anchor_value = m_init;
            ms->rate_param = 0.0;
            ms->anchor_time = 0.0;
        }
    }

    for (int i = 0; i < eventNumber; i++) {
        struct event *e = &events[i];
        char t = e->type;

        /* Migration events: validate trajectory and the new rate. */
        if (t == 'm') {
            int src = e->popID2;
            int dst = e->popID;
            if (src < 0 || src >= npops_eff || dst < 0 || dst >= npops_eff) continue;
            if (merged[src] || merged[dst]) continue;

            Shape *ms = &migState[src * npops_eff + dst];

            /* Prior linear shape must not have crossed zero by e->time. */
            if (ms->type == SHAPE_LINEAR && ms->rate_param > 0.0) {
                double t_cross = ms->anchor_time + ms->anchor_value / ms->rate_param;
                if (t_cross < e->time) {
                    fprintf(stderr,
                        "Error: linear migration trajectory from population %d to %d "
                        "drives the rate to zero at internal time %g (before the "
                        "next event for that pair at time %g).\n",
                        src, dst, t_cross, e->time);
                    free(migState);
                    return -1;
                }
            }

            double mig_at_event = sizeFromShape(ms, e->time);
            if (mig_at_event < 0.0) {
                fprintf(stderr,
                    "Error: migration rate from population %d to %d becomes negative "
                    "(%g) at internal time %g.\n",
                    src, dst, mig_at_event, e->time);
                free(migState);
                return -1;
            }

            if (e->popnSize < 0.0) {
                fprintf(stderr,
                    "Error: -em sets migration from population %d to %d to %g at time "
                    "%g; rate must be non-negative.\n",
                    src, dst, e->popnSize, e->time);
                free(migState);
                return -1;
            }

            ms->type = SHAPE_CONSTANT;
            ms->anchor_value = e->popnSize;
            ms->rate_param = 0.0;
            ms->anchor_time = e->time;
            continue;
        }

        /* Size events. */
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
                free(migState);
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
            free(migState);
            return -1;
        }

        if (t == 'n') {
            if (e->popnSize <= 0.0) {
                fprintf(stderr,
                    "Error: -en sets population %d size to %g at time %g; size must be "
                    "strictly positive.\n", p, e->popnSize, e->time);
                free(migState);
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
    free(migState);
    return 0;
}

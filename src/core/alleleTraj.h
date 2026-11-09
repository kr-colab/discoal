#include <math.h>

double detSweepFreq(double t, double s);

/* Generalized deterministic-sweep frequency under (possibly time-varying)
 * alpha_eff. A is the integrated selection coefficient int_0^tau alpha_eff(s) ds.
 * For constant alpha_eff = alpha, A = alpha * tau and this reduces to
 * detSweepFreq(tau, alpha) bit-exactly. */
double detSweepFreqGeneral(double alpha, double A);

double genicSelectionStochastic(double dt, double currentFreq, double alpha);
double genicSelectionStochasticForwards(double dt, double currentFreq, double alpha);
double coth(double x);
double variablePopnSizeTraj(double dt, double currentFreq, double alpha, double h, double f);

/* Forward declaration for ranf() */
double ranf(void);

/* Fast approximation for coth(x) with 8x speedup and <0.1% error
 * Uses Taylor series for small x and asymptotic behavior for large x */
static inline double fast_coth(double x) {
    // For small |x|, use Taylor series: coth(x) ≈ 1/x + x/3 - x³/45
    if (fabs(x) < 0.1) {
        double x2 = x * x;
        return 1.0/x + x/3.0 - x*x2/45.0;
    }
    // For large |x|, coth(x) → sign(x)
    if (fabs(x) > 4.0) {
        return x > 0 ? 1.0 : -1.0;
    }
    // Standard implementation for intermediate values
    // Use expm1 for better accuracy and performance: e^(2x) = 1 + expm1(2x)
    double e2x_m1 = expm1(2.0 * x);  // e^(2x) - 1
    return (e2x_m1 + 2.0) / e2x_m1;
}

/* Fast approximation for tanh(x) */
static inline double fast_tanh(double x) {
    // For small |x|, use Padé approximation
    if (fabs(x) < 0.5) {
        double x2 = x * x;
        return x * (27.0 + x2) / (27.0 + 9.0 * x2);
    }
    // For large |x|, tanh(x) → sign(x)
    if (fabs(x) > 3.0) {
        return x > 0 ? 1.0 : -1.0;
    }
    // Standard implementation for intermediate values
    return tanh(x);
}

/* One step back in time of a neutral allele's frequency, conditioned on
 * loss (drift -x, variance x(1-x) per unit time): the history of an allele
 * that fixed, followed back to its origin (-wn). */
static inline double neutralStochasticOptimized(double dt, double currentFreq) {
    double drift_term = -currentFreq * dt;
    double variance = currentFreq * (1.0 - currentFreq) * dt;

    if (variance <= 0.0) return currentFreq + drift_term;

    /* The noise is +-sqrt(variance) with probability 1/2 each, so the step
     * has the diffusion's variance; a uniform draw on (-1, 1) would give a
     * third of it. */
    double diffusion_term = sqrt(variance);
    return currentFreq + drift_term +
        (ranf() < 0.5 ? diffusion_term : -diffusion_term);
}

/* One step back in time of a neutral allele's frequency below `ceiling`,
 * the frequency it first reached, forward in time, at the end of this
 * phase (for a soft sweep, f0: the allele drifted up to f0 and then became
 * beneficial). Backward from that first passage, the frequency is neutral
 * drift conditioned to reach 0 without returning to the ceiling: drift
 * -x(1-x)/(ceiling - x), variance x(1-x) per unit time. The drift is
 * singular at the ceiling, where this phase starts, so the step is taken
 * on z = (ceiling - x)^2, which by Ito's formula follows
 * dz = 3x(1-x) dt - 2 sqrt(z x(1-x)) dW with no singularity; x stays below
 * the ceiling. (With ceiling 1 the drift is -x, as in
 * neutralStochasticOptimized.) */
static inline double neutralStochasticBelowCeiling(double dt, double currentFreq,
                                                   double ceiling) {
    double p_q = currentFreq * (1.0 - currentFreq);
    if (p_q <= 0.0) return currentFreq;
    double z = (ceiling - currentFreq) * (ceiling - currentFreq);
    double noise = 2.0 * sqrt(z * p_q * dt);
    z += 3.0 * p_q * dt + (ranf() < 0.5 ? noise : -noise);
    if (z < 0.0) z = -z;  /* reflect at the ceiling */
    return ceiling - sqrt(z);
}

/* Optimized inline version of genicSelectionStochastic (backwards in time) */
static inline double genicSelectionStochasticOptimized(double dt, double currentFreq, double alpha) {
    // Precompute common expressions
    double p_q = currentFreq * (1.0 - currentFreq);
    double alpha_p_q = alpha * p_q;
    double half_alpha_p_q = 0.5 * alpha_p_q;
    
    // Early return for edge cases
    if (p_q <= 0.0) return currentFreq;
    
    double drift_term = (-half_alpha_p_q * fast_coth(half_alpha_p_q)) * dt;
    double diffusion_term = sqrt(p_q * dt);
    
    if (ranf() < 0.5) {
        return currentFreq + drift_term + diffusion_term;
    } else {
        return currentFreq + drift_term - diffusion_term;
    }
}

/* Optimized inline version of genicSelectionStochasticForwards */
static inline double genicSelectionStochasticForwardsOptimized(double dt, double currentFreq, double alpha) {
    // Precompute common expressions
    double p_q = currentFreq * (1.0 - currentFreq);
    
    // Early return for edge cases
    if (p_q <= 0.0) return currentFreq;
    
    double alpha_p = alpha * currentFreq;
    double drift_term = (alpha * p_q / fast_tanh(alpha_p)) * dt;
    double diffusion_term = sqrt(p_q * dt);
    
    if (ranf() < 0.5) {
        return currentFreq + drift_term + diffusion_term;
    } else {
        return currentFreq + drift_term - diffusion_term;
    }
}
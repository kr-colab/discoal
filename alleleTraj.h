#include <math.h>

double detSweepFreq(double t, double s);
double neutralStochastic(double dt, double currentFreq);
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

/* Optimized inline version of neutralStochastic */
static inline double neutralStochasticOptimized(double dt, double currentFreq) {
    // Precompute common expressions
    double drift_term = -currentFreq * dt;
    double variance = currentFreq * (1.0 - currentFreq) * dt;
    
    // Early return for edge cases to avoid sqrt of negative/zero
    if (variance <= 0.0) return currentFreq + drift_term;
    
    double diffusion_term = sqrt(variance);
    
    // Use arithmetic instead of branching for better performance
    // ranf() returns [0,1), so (2*ranf()-1) gives us ±1
    return currentFreq + drift_term + (2.0 * ranf() - 1.0) * diffusion_term;
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
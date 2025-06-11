#include <math.h>

double detSweepFreq(double t, double s);
double neutralStochastic(double dt, double currentFreq);
double genicSelectionStochastic(double dt, double currentFreq, double alpha);
double genicSelectionStochasticForwards(double dt, double currentFreq, double alpha);
double coth(double x);
double variablePopnSizeTraj(double dt, double currentFreq, double alpha, double h, double f);

/* Forward declaration for ranf() */
double ranf(void);

/* Optimized inline version of neutralStochastic */
static inline double neutralStochasticOptimized(double dt, double currentFreq) {
    // Precompute common expressions
    double drift_term = -currentFreq * dt;
    double variance = currentFreq * (1.0 - currentFreq) * dt;
    
    // Early return for edge cases to avoid sqrt of negative/zero
    if (variance <= 0.0) return currentFreq + drift_term;
    
    double diffusion_term = sqrt(variance);
    
    // Keep original branching structure but optimize the calculation
    if (ranf() < 0.5) {
        return currentFreq + drift_term + diffusion_term;
    } else {
        return currentFreq + drift_term - diffusion_term;
    }
}

/* Optimized inline version of genicSelectionStochastic (backwards in time) */
static inline double genicSelectionStochasticOptimized(double dt, double currentFreq, double alpha) {
    // Precompute common expressions
    double p_q = currentFreq * (1.0 - currentFreq);
    double alpha_p_q = alpha * p_q;
    double half_alpha_p_q = 0.5 * alpha_p_q;
    
    // Early return for edge cases
    if (p_q <= 0.0) return currentFreq;
    
    double drift_term = (-half_alpha_p_q * coth(half_alpha_p_q)) * dt;
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
    double drift_term = (alpha * p_q / tanh(alpha_p)) * dt;
    double diffusion_term = sqrt(p_q * dt);
    
    if (ranf() < 0.5) {
        return currentFreq + drift_term + diffusion_term;
    } else {
        return currentFreq + drift_term - diffusion_term;
    }
}
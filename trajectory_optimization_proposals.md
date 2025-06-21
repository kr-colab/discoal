# Trajectory Simulation Optimization Proposals

## Current Bottlenecks

From profiling and analysis:
1. **Hyperbolic functions** (coth, tanh) are expensive (~10-20 cycles each)
2. **Repeated calculations** of p*(1-p) and sqrt terms
3. **Function call overhead** in tight loops
4. **Branch misprediction** from random sign selection

## Proposed Optimizations

### 1. Fast Approximations for Hyperbolic Functions

#### For `coth(x) = (e^(2x) + 1)/(e^(2x) - 1)`

```c
static inline double fast_coth(double x) {
    // For |x| < 0.01, use Taylor series: coth(x) ≈ 1/x + x/3
    if (fabs(x) < 0.01) {
        return 1.0/x + x/3.0;
    }
    // For |x| > 5, coth(x) ≈ sign(x)
    if (fabs(x) > 5.0) {
        return x > 0 ? 1.0 : -1.0;
    }
    // Otherwise use optimized version with single exp
    double e2x = exp(2.0 * x);
    return (e2x + 1.0) / (e2x - 1.0);
}
```

#### For `tanh(x)`

```c
static inline double fast_tanh(double x) {
    // For |x| < 0.01, tanh(x) ≈ x
    if (fabs(x) < 0.01) {
        return x;
    }
    // For |x| > 5, tanh(x) ≈ sign(x)
    if (fabs(x) > 5.0) {
        return x > 0 ? 1.0 : -1.0;
    }
    // Otherwise use standard implementation
    return tanh(x);
}
```

### 2. Lookup Table for Common Values

Since alpha is often fixed during a simulation run:

```c
typedef struct {
    double alpha;
    double coth_values[101];  // For p = 0.00, 0.01, ..., 1.00
    double tanh_values[101];
    int initialized;
} TrajectoryCache;

static TrajectoryCache traj_cache = {0};

void init_trajectory_cache(double alpha) {
    if (traj_cache.initialized && traj_cache.alpha == alpha) return;
    
    traj_cache.alpha = alpha;
    for (int i = 0; i <= 100; i++) {
        double p = i / 100.0;
        double arg = 0.5 * alpha * p * (1.0 - p);
        traj_cache.coth_values[i] = (arg > 0.001) ? coth(arg) : 1.0/arg;
        
        arg = alpha * p;
        traj_cache.tanh_values[i] = (arg > 0.001) ? tanh(arg) : arg;
    }
    traj_cache.initialized = 1;
}

// Use linear interpolation for lookup
static inline double lookup_coth(double p) {
    int idx = (int)(p * 100);
    double frac = p * 100 - idx;
    if (idx >= 100) return traj_cache.coth_values[100];
    return traj_cache.coth_values[idx] * (1-frac) + 
           traj_cache.coth_values[idx+1] * frac;
}
```

### 3. Vectorized Trajectory Generation

Generate multiple trajectories simultaneously using SIMD:

```c
#include <immintrin.h>  // For AVX2

void generate_trajectories_vectorized(double *frequencies, int n_traj, 
                                     double alpha, double dt, int steps) {
    // Process 4 trajectories at a time with AVX2
    for (int i = 0; i < n_traj; i += 4) {
        __m256d freq = _mm256_load_pd(&frequencies[i]);
        
        for (int step = 0; step < steps; step++) {
            // Calculate p*(1-p) for all 4 trajectories
            __m256d one = _mm256_set1_pd(1.0);
            __m256d p_q = _mm256_mul_pd(freq, _mm256_sub_pd(one, freq));
            
            // Generate 4 random numbers at once
            __m256d rand_vals = _mm256_set_pd(ranf(), ranf(), ranf(), ranf());
            
            // ... vectorized drift and diffusion calculations ...
        }
        
        _mm256_store_pd(&frequencies[i], freq);
    }
}
```

### 4. Adaptive Time Stepping

Reduce function calls by using larger time steps when possible:

```c
double adaptive_dt(double currentFreq, double alpha, double base_dt) {
    double p_q = currentFreq * (1.0 - currentFreq);
    
    // Near boundaries, use smaller steps
    if (p_q < 0.01) return base_dt;
    
    // In the middle, can use larger steps
    if (p_q > 0.20) return base_dt * 4.0;
    
    return base_dt * 2.0;
}
```

### 5. Optimized Random Sign Selection

Replace branching with arithmetic:

```c
// Instead of:
if (ranf() < 0.5) {
    return currentFreq + drift + diffusion;
} else {
    return currentFreq + drift - diffusion;
}

// Use:
double sign = 2.0 * ranf() - 1.0;  // Generates -1 or +1
return currentFreq + drift + sign * diffusion;
```

### 6. Expression Caching and Reordering

```c
static inline double genicSelectionStochasticOptimized2(
    double dt, double currentFreq, double alpha) {
    
    // Cache all repeated calculations
    double p = currentFreq;
    double q = 1.0 - p;
    double p_q = p * q;
    
    // Early exit for boundary cases
    if (p_q < 1e-10) return p;
    
    // Precompute expensive terms
    double alpha_p_q = alpha * p_q;
    double half_alpha_p_q = 0.5 * alpha_p_q;
    
    // Use fast approximation
    double coth_term = (half_alpha_p_q < 0.01) ? 
                       1.0/half_alpha_p_q + half_alpha_p_q/3.0 :
                       fast_coth(half_alpha_p_q);
    
    double drift = -half_alpha_p_q * coth_term * dt;
    double diffusion = sqrt(p_q * dt);
    
    // Optimized random sign
    return p + drift + (2.0 * ranf() - 1.0) * diffusion;
}
```

### 7. Batch Processing

Process multiple time steps with fewer function calls:

```c
void trajectory_batch(double *freq_array, int batch_size, 
                     double alpha, double dt) {
    // Generate all random numbers at once
    double randoms[batch_size];
    for (int i = 0; i < batch_size; i++) {
        randoms[i] = 2.0 * ranf() - 1.0;
    }
    
    // Process trajectory
    for (int i = 0; i < batch_size; i++) {
        double p = *freq_array;
        double p_q = p * (1.0 - p);
        
        // ... calculate drift and diffusion ...
        
        *freq_array = p + drift + randoms[i] * diffusion;
        freq_array++;
    }
}
```

## Expected Performance Gains

1. **Fast approximations**: 2-5x speedup for hyperbolic functions
2. **Lookup tables**: Near-constant time for fixed alpha values
3. **Vectorization**: 3-4x speedup for multiple trajectories
4. **Adaptive stepping**: 20-50% fewer function calls
5. **Expression optimization**: 10-20% improvement
6. **Batch processing**: Better cache utilization

## Implementation Priority

1. **High priority** (easy wins):
   - Fast approximations for coth/tanh
   - Expression caching in existing functions
   - Optimized random sign

2. **Medium priority** (moderate effort):
   - Lookup tables for common values
   - Adaptive time stepping
   - Batch processing

3. **Low priority** (requires architecture changes):
   - Full vectorization
   - GPU acceleration
   - Parallel trajectory generation

## Testing Requirements

- Verify approximation accuracy (< 0.1% error)
- Statistical validation of trajectory distributions
- Performance benchmarks across parameter ranges
- Edge case testing (p near 0 or 1)
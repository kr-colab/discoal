# DSC Parameters Refactoring Design

## Overview

This document outlines the proposed refactoring of discoal's parameter system from scattered global variables to a structured, extensible parameter system that supports both command-line and YAML configuration.

## Design Goals

1. **Centralized Configuration**: All parameters in a single hierarchical structure
2. **Type Safety**: Clear parameter types and validation
3. **Extensibility**: Easy to add new parameters
4. **Multiple Input Formats**: Support CLI, YAML, and potentially JSON
5. **Backward Compatibility**: Maintain existing CLI interface
6. **Validation**: Centralized parameter validation and dependency checking

## Proposed Structure

```c
// Main simulation parameters structure
typedef struct {
    SimulationCore      core;
    EvolutionaryForces  forces;
    Demographics        demographics;
    SelectionParams     selection;
    OutputConfig        output;
    PriorConfig         priors;
    DebugConfig         debug;
    RandomConfig        random;
} SimulationParams;

// Core simulation parameters
typedef struct {
    int     total_samples;      // Total number of samples across all populations
    int     num_replicates;     // Number of simulation replicates
    int     num_sites;          // Number of sites in sequence
    int     seg_sites;          // Number of segregating sites (if fixed)
    int     num_populations;    // Number of populations
    int     sample_sizes[MAXPOPS]; // Sample size per population
    double  Ne;                 // Effective population size
} SimulationCore;

// Evolutionary forces
typedef struct {
    // Mutation
    double  theta;              // 4Nu
    double  theta_site;         // Per-site theta
    
    // Recombination
    double  rho;                // 4Nr
    double  left_rho;           // Recombination to the left
    double  gene_conversion_ratio; // Fraction of recombination that is gene conversion
    int     gc_tract_mean;      // Mean gene conversion tract length
    
    // Migration
    double  migration_matrix[MAXPOPS][MAXPOPS];
    bool    symmetric_migration;
    double  migration_rate;     // For simple island model
} EvolutionaryForces;

// Demographic events and structure
typedef struct {
    PopulationSizes     pop_sizes;
    DemographicEvent    *events;
    int                 num_events;
    int                 events_capacity;
    AncientSamples      ancient_samples;
} Demographics;

// Population sizes
typedef struct {
    double  sizes[MAXPOPS];     // Current sizes
    double  initial_sizes[MAXPOPS]; // Starting sizes
    double  growth_rates[MAXPOPS];  // Exponential growth rates
} PopulationSizes;

// Improved demographic event structure
typedef struct {
    double  time;
    enum {
        EVENT_SIZE_CHANGE,
        EVENT_SPLIT,
        EVENT_JOIN,
        EVENT_ADMIX,
        EVENT_MIGRATION_CHANGE,
        EVENT_GROWTH_RATE
    } type;
    
    union {
        struct { int pop; double size; } size_change;
        struct { int from; int to1; int to2; double prop; } split;
        struct { int pop1; int pop2; int dest; } join;
        struct { int from; int to; double prop; } admix;
        struct { int from; int to; double rate; } migration;
        struct { int pop; double rate; } growth;
    } params;
} DemographicEvent;

// Selection parameters
typedef struct {
    double  alpha;              // 2Ns
    double  sweep_position;     // Location on chromosome (0-1)
    double  tau;                // Time of sweep
    double  f0;                 // Starting frequency
    double  adaptive_mutation_rate; // Rate of adaptive mutations (uA)
    
    enum {
        SWEEP_DETERMINISTIC = 'd',
        SWEEP_STOCHASTIC = 's',
        SWEEP_NEUTRAL = 'N',
        SWEEP_NEUT_STOCHASTIC = 'n'
    } sweep_mode;
    
    bool    partial_sweep;
    double  partial_sweep_final_freq;
    
    bool    soft_sweep;
    int     soft_sweep_instances;
    
    bool    recurrent_sweep;
} SelectionParams;

// Prior distributions
typedef struct {
    bool            enabled;
    PriorDist       *distributions;
    int             num_priors;
} PriorConfig;

typedef struct {
    char            *parameter_name;  // e.g., "theta", "rho"
    enum { UNIFORM, LOG_UNIFORM, NORMAL, EXPONENTIAL } type;
    double          lower_bound;
    double          upper_bound;
    double          mean;            // For normal/exponential
    double          stddev;          // For normal
} PriorDist;

// Output configuration
typedef struct {
    enum { OUTPUT_MS, OUTPUT_FASTA, OUTPUT_NEXUS } format;
    bool            tree_sequence_output;
    char            output_filename[PATH_MAX];
    bool            minimal_tree_seq;
    bool            hide_singletons;
    bool            hide_partial_snp;
} OutputConfig;

// Debug and diagnostic options
typedef struct {
    bool    verbose;
    bool    print_trajectory;
    bool    print_trees;
    int     debug_level;
} DebugConfig;

// Random number generation
typedef struct {
    unsigned long   seed1;
    unsigned long   seed2;
    bool            use_xoshiro;  // Use modern RNG
} RandomConfig;
```

## Parameter Loading Interface

```c
// Main parameter loading functions
SimulationParams* params_create(void);
void params_destroy(SimulationParams *params);

// Loading from different sources
int params_load_from_args(SimulationParams *params, int argc, char *argv[]);
int params_load_from_yaml(SimulationParams *params, const char *filename);
int params_load_from_json(SimulationParams *params, const char *filename);

// Validation
int params_validate(const SimulationParams *params, char *error_msg, size_t msg_size);

// Utility functions
void params_print_summary(const SimulationParams *params, FILE *output);
int params_save_to_yaml(const SimulationParams *params, const char *filename);
SimulationParams* params_copy(const SimulationParams *params);
```

## YAML Configuration Example

```yaml
# discoal simulation parameters
simulation:
  samples: 100
  replicates: 1000
  sites: 10000
  populations: 2
  sample_distribution: [50, 50]
  
evolution:
  mutation:
    theta: 100.0
  recombination:
    rho: 100.0
    gene_conversion:
      enabled: true
      ratio: 0.1
      tract_length: 500
  migration:
    symmetric: true
    rate: 1.0

selection:
  enabled: true
  coefficient: 100.0  # 2Ns
  position: 0.5
  mode: deterministic
  partial_sweep:
    enabled: true
    final_frequency: 0.8

demographics:
  events:
    - time: 0.01
      type: size_change
      population: 0
      size: 10000
    - time: 0.05
      type: split
      from: 0
      to: [1, 2]
      proportions: [0.3, 0.7]

output:
  format: ms
  tree_sequence:
    enabled: true
    filename: "simulation.trees"
    minimal: false

random:
  seed1: 12345
  seed2: 67890
  engine: xoshiro256++
```

## Implementation Plan

### Phase 1: Core Structure (Week 1)
- [ ] Define parameter structures in new header `params.h`
- [ ] Implement basic creation/destruction functions
- [ ] Create parameter validation framework

### Phase 2: Command Line Parser (Week 2)
- [ ] Refactor `getParameters()` to populate new structures
- [ ] Maintain backward compatibility with existing CLI
- [ ] Add parameter validation after parsing

### Phase 3: YAML Support (Week 3)
- [ ] Integrate libyaml or similar
- [ ] Implement YAML parser for parameters
- [ ] Add schema validation

### Phase 4: Migration (Week 4)
- [ ] Replace global variables with parameter structure
- [ ] Update all functions to take `SimulationParams*`
- [ ] Extensive testing of all parameter combinations

### Phase 5: Documentation and Testing (Week 5)
- [ ] Update documentation
- [ ] Add comprehensive parameter tests
- [ ] Create migration guide for users

## Benefits

1. **Cleaner Code**: All parameters in one place, easier to understand
2. **Better Testing**: Can easily create parameter sets for testing
3. **Configuration Files**: Users can save and share complex parameter sets
4. **Validation**: Centralized validation catches errors early
5. **Extensibility**: Adding new parameters is straightforward
6. **Reduced Globals**: Move away from global state

## Backward Compatibility

The existing command-line interface will be preserved. Users can continue using:
```bash
./discoal 100 1000 10000 -t 100 -r 100
```

But can also use:
```bash
./discoal --config simulation.yaml
./discoal --config simulation.yaml -t 200  # Override theta from config
```

## Open Questions

1. Should we support JSON in addition to YAML?
2. How to handle parameter precedence (CLI overrides config file)?
3. Should we version the parameter schema for future compatibility?
4. Memory management strategy for dynamic arrays in parameters?
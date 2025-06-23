/**
 * params.h - Centralized parameter structure for discoal simulations
 * 
 * This header defines the parameter structures used to configure
 * discoal simulations, replacing the previous global variable approach.
 */

#ifndef PARAMS_H
#define PARAMS_H

#include <stdbool.h>
#include <limits.h>

#define MAXPOPS 2000
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/* Forward declarations */
typedef struct SimulationParams SimulationParams;

/* Core simulation parameters */
typedef struct {
    int     total_samples;          /* Total number of samples across all populations */
    int     num_replicates;         /* Number of simulation replicates */
    int     num_sites;              /* Number of sites in sequence */
    int     seg_sites;              /* Number of segregating sites (if fixed) */
    int     num_populations;        /* Number of populations */
    int     sample_sizes[MAXPOPS];  /* Sample size per population */
    double  Ne;                     /* Effective population size */
} SimulationCore;

/* Evolutionary forces */
typedef struct {
    /* Mutation */
    double  theta;                  /* 4Nu */
    double  theta_site;             /* Per-site theta */
    
    /* Recombination */
    double  rho;                    /* 4Nr */
    double  left_rho;               /* Recombination to the left */
    double  gene_conversion_ratio;  /* Fraction of recombination that is gene conversion */
    int     gc_tract_mean;          /* Mean gene conversion tract length */
    
    /* Migration */
    double  migration_matrix[MAXPOPS][MAXPOPS];
    bool    symmetric_migration;
    double  migration_rate;         /* For simple island model */
} EvolutionaryForces;

/* Demographic event types */
typedef enum {
    EVENT_SIZE_CHANGE,
    EVENT_SPLIT,
    EVENT_JOIN,
    EVENT_ADMIX,
    EVENT_MIGRATION_CHANGE,
    EVENT_GROWTH_RATE
} EventType;

/* Demographic event structure */
typedef struct {
    double      time;
    EventType   type;
    
    union {
        struct { int pop; double size; } size_change;
        struct { int from; int to1; int to2; double prop; } split;
        struct { int pop1; int pop2; int dest; } join;
        struct { int from; int to; double prop; } admix;
        struct { int from; int to; double rate; } migration;
        struct { int pop; double rate; } growth;
    } params;
} DemographicEvent;

/* Population sizes */
typedef struct {
    double  sizes[MAXPOPS];         /* Current sizes */
    double  initial_sizes[MAXPOPS]; /* Starting sizes */
    /* Note: discoal doesn't support exponential growth, only discrete size changes */
} PopulationSizes;

/* Sample specification */
typedef struct {
    char    *population;    /* Population name (for Demes) or index */
    int     size;          /* Number of samples */
    double  time;          /* Sampling time (default 0 = present) */
} SampleSpec;

/* Ancient samples configuration */
typedef struct {
    bool    enabled;
    int     *sample_times;          /* Time for each ancient sample */
    int     *sample_pops;           /* Population for each ancient sample */
    int     num_ancient_samples;
} AncientSamples;

/* Demographics configuration */
typedef struct {
    PopulationSizes     pop_sizes;
    DemographicEvent    *events;
    int                 num_events;
    int                 events_capacity;
    AncientSamples      ancient_samples;
    SampleSpec          *sample_specs;
    int                 num_sample_specs;
} Demographics;

/* Sweep modes */
typedef enum {
    SWEEP_DETERMINISTIC = 'd',
    SWEEP_STOCHASTIC = 's',
    SWEEP_NEUTRAL = 'N',
    SWEEP_NEUT_STOCHASTIC = 'n'
} SweepMode;

/* Selection parameters */
typedef struct {
    double      alpha;                  /* 2Ns */
    double      sweep_position;         /* Location on chromosome (0-1) */
    double      tau;                    /* Time of sweep */
    double      f0;                     /* Starting frequency */
    double      adaptive_mutation_rate; /* Rate of adaptive mutations (uA) */
    SweepMode   sweep_mode;
    
    /* Partial sweep */
    bool        partial_sweep;
    double      partial_sweep_final_freq;
    
    /* Soft sweep */
    bool        soft_sweep;
    int         soft_sweep_instances;
    
    /* Recurrent sweep */
    bool        recurrent_sweep;
    double      recurrent_sweep_rate;   /* Rate of recurrent sweeps (for -R option) */
} SelectionParams;

/* Prior distribution types */
typedef enum {
    PRIOR_UNIFORM,
    PRIOR_LOG_UNIFORM,
    PRIOR_NORMAL,
    PRIOR_EXPONENTIAL
} PriorType;

/* Prior distribution specification */
typedef struct {
    char        *parameter_name;    /* e.g., "theta", "rho" */
    PriorType   type;
    double      lower_bound;
    double      upper_bound;
    double      mean;               /* For normal/exponential */
    double      stddev;             /* For normal */
} PriorDist;

/* Prior configuration */
typedef struct {
    bool        enabled;
    PriorDist   *distributions;
    int         num_priors;
} PriorConfig;

/* Output formats */
typedef enum {
    OUTPUT_MS,
    OUTPUT_FASTA,
    OUTPUT_NEXUS
} OutputFormat;

/* Output configuration */
typedef struct {
    OutputFormat    format;
    bool            tree_sequence_output;
    char            output_filename[PATH_MAX];
    bool            minimal_tree_seq;
    bool            hide_singletons;
    bool            hide_partial_snp;
    int             mask;               /* Masking parameter (reserved for future use) */
} OutputConfig;

/* Debug configuration */
typedef struct {
    bool    verbose;
    bool    print_trajectory;
    bool    print_trees;
    int     debug_level;
} DebugConfig;

/* Random number configuration */
typedef struct {
    unsigned long   seed1;
    unsigned long   seed2;
    bool            use_xoshiro;    /* Use modern RNG */
} RandomConfig;

/* Main simulation parameters structure */
struct SimulationParams {
    SimulationCore      core;
    EvolutionaryForces  forces;
    Demographics        demographics;
    SelectionParams     selection;
    OutputConfig        output;
    PriorConfig         priors;
    DebugConfig         debug;
    RandomConfig        random;
};

/* Function prototypes */

/* Creation and destruction */
SimulationParams* params_create(void);
void params_destroy(SimulationParams *params);

/* Loading from different sources */
int params_load_from_args(SimulationParams *params, int argc, char *argv[]);

/* YAML support using libyaml - implemented in yaml_loader.c */
int yaml_load_params(SimulationParams *params, const char *filename);
int yaml_save_params(const SimulationParams *params, const char *filename);

/* Validation */
int params_validate(const SimulationParams *params, char *error_msg, size_t msg_size);

/* Utility functions */
void params_print_summary(const SimulationParams *params, FILE *output);
SimulationParams* params_copy(const SimulationParams *params);

/* Helper functions for common parameter patterns */
void params_set_symmetric_migration(SimulationParams *params, double rate);
void params_add_demographic_event(SimulationParams *params, DemographicEvent *event);
void params_set_simple_sweep(SimulationParams *params, double s, double position);

#endif /* PARAMS_H */
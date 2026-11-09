#ifndef CONFIGINTERFACE_H
#define CONFIGINTERFACE_H

#include <stdlib.h>
#include <stdio.h>
#include <cyaml/cyaml.h>
#include "discoal.h"
#include "demesInterface.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/* To add new options to the YAML interface:
 *
 *    - If the option type is a string, create a enum for different
 *      possibilities and a cyaml_strval_t array mapping enum to string.
 *
 *    - If the option is a sequence, define the datatype of the sequence
 *      with a cyaml_schema_value_t.
 *
 *    - Add the option as a member to the appropriate subconfig struct
 *      (e.g. "genetics_config"). If it is a sequence, it should
 *      be a pointer.
 *
 *    - Add the field to the cyaml schema for the struct (e.g. 
 *      "genetics_fields_schema"), using CYAML_FLAG_OPTIONAL if it
 *      is optional (otherwise the YAML parser will expect the option and exit
 *      with an informative error).
 *
 *    - Handle initialization and parsing in the appropriate parse_* function
 *
 *    - Add the new option to `config_examples/all_options.yaml`, even if there
 *      are conflicts (i.e. if the option is incompatible with another option and
 *      will error out discoal)
 */


/* enumerate string-valued options */
enum output_types {
    OUTPUT_DEFAULT = 0,      /* output_type not given: haplotype */
    OUTPUT_HAPLOTYPE = 1,
    OUTPUT_SNP_ARRAY = 2,
    OUTPUT_TREE_SEQN = 3,
};

enum sweep_modes {
    SWEEP_STOCHASTIC = 1,
    SWEEP_DETERMINISTIC = 2,
    SWEEP_NEUTRAL = 3,
};


/* for YAML parsing, use a distinct struct for each type of demographic event
 * so as to enforce unique key-to-member mapping */
struct population_size_change {
    double time, size;
    int population;
};
struct migration_rate_change {
    double time, rate;
    int source, destination;
};
struct population_split {
    double time;
    int derived, ancestral;
};
/* -ea: backward in time, the lineages of `population` move to source_1
 * with probability `proportion` and to source_2 otherwise. */
struct admixture {
    double time;
    int population, source_1, source_2;
    double proportion;
};
struct ancient_sample {
    double time;
    int sample_size, population;
};
struct demographic_events {  
    struct population_size_change *population_size_changes;
    unsigned num_population_size_changes;
    struct migration_rate_change *migration_rate_changes;
    unsigned num_migration_rate_changes;
    struct population_split *population_splits;
    unsigned num_population_splits; 
    struct ancient_sample *ancient_samples;
    unsigned num_ancient_samples;
    struct admixture *admixtures;
    unsigned num_admixtures;
};


/* map migration matrix to a sequence of rows */
struct migration_matrix_row {
    double *rates;
    unsigned num_cols;
};


/* blocks of config */
struct simulation_config {   
    int sample_size;
    int num_replicates; 
    int num_sites;
    int *seed;          /* optional, two seeds expected */
};
struct genetics_config {
    double mutation_rate;             
    double recombination_rate;        
    double gene_conversion_rate;      /* optional */
    int gene_conversion_tract_length; /* optional */
    double crossover_ratio;           /* optional */
};
struct demography_config {
    int *deme_sample_size;
    unsigned num_demes;
    double effective_population_size;              /* optional, applied to all populations */
    struct demographic_events *demographic_events; /* optional, contains arrays of events  */
    struct migration_matrix_row *migration_matrix; /* optional, array of arrays of floats  */
    unsigned num_migration_matrix_rows;
    const char *demes_filename;                    /* optional, path to demes YAML */
};
struct output_config {
    enum output_types output_type;       /* optional, default haplotype */
    bool hide_partial_snp;               /* optional, as -h */
    bool unsimplified_tree_sequence;     /* optional, don't simplify tree sequence */
    const char *tree_sequence_filename;  /* required if output type is "tree_sequence" */
};
struct selection_config {
    enum sweep_modes sweep_mode;
    double selection_coefficient;    /* -a */
    double *sweep_position;          /* optional, -x; NULL if absent (0.5) */
    double *fixation_time_ago;       /* optional, -ws/-wd/-wn tau; NULL if absent */
    double initial_frequency;        /* optional, -f; triggers soft sweep */
    double final_frequency;          /* optional, -c; triggers partial sweep */
    double beneficial_mutation_rate; /* optional, -uA */
    double recurrent_sweep_rate;     /* optional, -R (or -L, below) */
    double *left_rho;                /* optional, -ls/-ld/-ln leftRho; NULL if absent */
    bool recurrent_sweeps_left_of_locus; /* optional, with recurrent_sweep_rate: -L */
    double *sweep_time_step;         /* optional, -i; NULL if absent */
};


/* top level struct for config */
/* Priors: each replaces a fixed value with a new draw for every replicate,
 * as the -P flags do. Uniform on [low, high]. */
struct uniform_prior {
    double low, high;
};
/* -Pre: exponential with this mean, capped (not truncated) at upper. */
struct exponential_prior {
    double mean, upper;
};
/* -Pe1 / -Pe2: a size change with uniform priors on its time and size. */
struct prior_size_change {
    double time_low, time_high, size_low, size_high;
    int population;                       /* optional, default 0 */
};
struct priors_config {
    struct uniform_prior *mutation_rate;              /* -Pt */
    struct uniform_prior *recombination_rate;         /* -Pr */
    struct exponential_prior *recombination_rate_exponential; /* -Pre */
    struct uniform_prior *selection_coefficient;      /* -Pa */
    struct uniform_prior *fixation_time_ago;          /* -Pu */
    struct uniform_prior *beneficial_mutation_rate;   /* -PuA */
    struct uniform_prior *sweep_position;             /* -Px */
    struct uniform_prior *initial_frequency;          /* -Pf */
    struct uniform_prior *final_frequency;            /* -Pc */
    struct prior_size_change *population_size_changes; /* -Pe1, -Pe2 */
    unsigned num_population_size_changes;
};
struct discoal_config {
    struct simulation_config *simulation;
    struct genetics_config *genetics;
    struct demography_config *demography; /* optional */
    struct output_config *output;         /* optional */
    struct selection_config *selection;   /* optional */
    struct priors_config *priors;         /* optional */
};


int parse_simulation_block(struct simulation_config*);
int parse_genetics_block(struct genetics_config*);
int parse_demography_block(struct demography_config*);
int parse_selection_block(struct selection_config*);
int parse_output_block(struct output_config*);
int parse_priors_block(struct priors_config*, int has_selection,
                       unsigned num_demes, const struct demographic_events*);
int apply_yaml_config(struct discoal_config*);
int load_yaml_config(const char*, struct discoal_config**);


#endif // CONFIGINTERFACE_H


/**
 * params.c - Implementation of centralized parameter management for discoal
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "params.h"
#include "yaml_loader.h"
#include "demes_loader.h"

/* Forward declaration for internal functions */
static int compare_events(const void *a, const void *b);

/**
 * Create a new SimulationParams structure with default values
 */
SimulationParams* params_create(void) {
    SimulationParams *params = calloc(1, sizeof(SimulationParams));
    if (!params) {
        fprintf(stderr, "Error: Failed to allocate memory for simulation parameters\n");
        return NULL;
    }
    
    /* Set default values */
    
    /* Core defaults */
    params->core.total_samples = 0;
    params->core.num_replicates = 1;
    params->core.num_sites = 0;
    params->core.seg_sites = -1;  /* -1 means use infinite sites model */
    params->core.num_populations = 1;
    params->core.Ne = 1.0;
    
    /* Initialize all population sizes to 0 */
    memset(params->core.sample_sizes, 0, sizeof(params->core.sample_sizes));
    
    /* Evolutionary forces defaults */
    params->forces.theta = 0.0;
    params->forces.theta_site = 0.0;
    params->forces.rho = 0.0;
    params->forces.left_rho = 0.0;
    params->forces.gene_conversion_ratio = 0.0;
    params->forces.gc_tract_mean = 500;  /* Default tract length */
    params->forces.symmetric_migration = false;
    params->forces.migration_rate = 0.0;
    
    /* Initialize migration matrix to 0 */
    memset(params->forces.migration_matrix, 0, sizeof(params->forces.migration_matrix));
    
    /* Demographics defaults */
    params->demographics.events = NULL;
    params->demographics.num_events = 0;
    params->demographics.events_capacity = 0;
    params->demographics.ancient_samples.enabled = false;
    params->demographics.ancient_samples.sample_times = NULL;
    params->demographics.ancient_samples.sample_pops = NULL;
    params->demographics.ancient_samples.num_ancient_samples = 0;
    params->demographics.sample_specs = NULL;
    params->demographics.num_sample_specs = 0;
    
    /* Initialize population sizes to 1.0 */
    for (int i = 0; i < MAXPOPS; i++) {
        params->demographics.pop_sizes.sizes[i] = 1.0;
        params->demographics.pop_sizes.initial_sizes[i] = 1.0;
    }
    
    /* Selection defaults */
    params->selection.alpha = 0.0;
    params->selection.sweep_position = 0.5;
    params->selection.tau = 0.0;
    params->selection.f0 = 1.0 / (2.0 * params->core.Ne);  /* Default starting frequency */
    params->selection.adaptive_mutation_rate = 0.0;
    params->selection.sweep_mode = SWEEP_NEUTRAL;
    params->selection.partial_sweep = false;
    params->selection.partial_sweep_final_freq = 1.0;
    params->selection.soft_sweep = false;
    params->selection.soft_sweep_instances = 1;
    params->selection.recurrent_sweep = false;
    params->selection.recurrent_sweep_rate = 0.0;
    
    /* Output defaults */
    params->output.format = OUTPUT_MS;
    params->output.tree_sequence_output = false;
    strcpy(params->output.output_filename, "discoal_output.trees");
    params->output.minimal_tree_seq = false;
    params->output.hide_singletons = false;
    params->output.hide_partial_snp = false;
    params->output.mask = 0;
    
    /* Prior defaults */
    params->priors.enabled = false;
    params->priors.distributions = NULL;
    params->priors.num_priors = 0;
    
    /* Debug defaults */
    params->debug.verbose = false;
    params->debug.print_trajectory = false;
    params->debug.print_trees = false;
    params->debug.debug_level = 0;
    
    /* Random defaults */
    params->random.seed1 = 0;  /* 0 means use time-based seed */
    params->random.seed2 = 0;
    params->random.use_xoshiro = false;
    
    return params;
}

/**
 * Destroy a SimulationParams structure and free all associated memory
 */
void params_destroy(SimulationParams *params) {
    if (!params) return;
    
    /* Free demographic events */
    if (params->demographics.events) {
        free(params->demographics.events);
    }
    
    /* Free ancient samples data */
    if (params->demographics.ancient_samples.sample_times) {
        free(params->demographics.ancient_samples.sample_times);
    }
    if (params->demographics.ancient_samples.sample_pops) {
        free(params->demographics.ancient_samples.sample_pops);
    }
    
    /* Free sample specifications */
    if (params->demographics.sample_specs) {
        for (int i = 0; i < params->demographics.num_sample_specs; i++) {
            if (params->demographics.sample_specs[i].population) {
                free(params->demographics.sample_specs[i].population);
            }
        }
        free(params->demographics.sample_specs);
    }
    
    /* Free prior distributions */
    if (params->priors.distributions) {
        for (int i = 0; i < params->priors.num_priors; i++) {
            if (params->priors.distributions[i].parameter_name) {
                free(params->priors.distributions[i].parameter_name);
            }
        }
        free(params->priors.distributions);
    }
    
    /* Free the main structure */
    free(params);
}

/**
 * Validate simulation parameters for consistency and correctness
 */
int params_validate(const SimulationParams *params, char *error_msg, size_t msg_size) {
    assert(params != NULL);
    assert(error_msg != NULL);
    
    /* Check core parameters */
    if (params->core.total_samples <= 0) {
        snprintf(error_msg, msg_size, "Total number of samples must be positive");
        return -1;
    }
    
    if (params->core.num_replicates <= 0) {
        snprintf(error_msg, msg_size, "Number of replicates must be positive");
        return -1;
    }
    
    if (params->core.num_sites <= 0) {
        snprintf(error_msg, msg_size, "Number of sites must be positive");
        return -1;
    }
    
    if (params->core.num_populations <= 0 || params->core.num_populations > MAXPOPS) {
        snprintf(error_msg, msg_size, "Number of populations must be between 1 and %d", MAXPOPS);
        return -1;
    }
    
    /* Verify sample sizes sum to total */
    int sample_sum = 0;
    for (int i = 0; i < params->core.num_populations; i++) {
        if (params->core.sample_sizes[i] < 0) {
            snprintf(error_msg, msg_size, "Sample size for population %d cannot be negative", i);
            return -1;
        }
        sample_sum += params->core.sample_sizes[i];
    }
    
    if (sample_sum != params->core.total_samples) {
        snprintf(error_msg, msg_size, "Sum of population sample sizes (%d) does not equal total samples (%d)", 
                 sample_sum, params->core.total_samples);
        return -1;
    }
    
    /* Check evolutionary parameters */
    if (params->forces.theta < 0) {
        snprintf(error_msg, msg_size, "Theta (4Nu) cannot be negative");
        return -1;
    }
    
    if (params->forces.rho < 0) {
        snprintf(error_msg, msg_size, "Rho (4Nr) cannot be negative");
        return -1;
    }
    
    if (params->forces.gene_conversion_ratio < 0 || params->forces.gene_conversion_ratio > 1) {
        snprintf(error_msg, msg_size, "Gene conversion ratio must be between 0 and 1");
        return -1;
    }
    
    /* Check selection parameters */
    if (params->selection.sweep_position < 0 || params->selection.sweep_position > 1) {
        snprintf(error_msg, msg_size, "Sweep position must be between 0 and 1");
        return -1;
    }
    
    if (params->selection.partial_sweep && 
        (params->selection.partial_sweep_final_freq <= 0 || 
         params->selection.partial_sweep_final_freq > 1)) {
        snprintf(error_msg, msg_size, "Partial sweep final frequency must be between 0 and 1");
        return -1;
    }
    
    /* Check demographic events are properly ordered by time */
    for (int i = 1; i < params->demographics.num_events; i++) {
        if (params->demographics.events[i].time < params->demographics.events[i-1].time) {
            snprintf(error_msg, msg_size, "Demographic events must be ordered by increasing time");
            return -1;
        }
    }
    
    /* Migration matrix validation */
    for (int i = 0; i < params->core.num_populations; i++) {
        for (int j = 0; j < params->core.num_populations; j++) {
            if (i != j && params->forces.migration_matrix[i][j] < 0) {
                snprintf(error_msg, msg_size, "Migration rate from pop %d to %d cannot be negative", i, j);
                return -1;
            }
            if (i == j && params->forces.migration_matrix[i][j] != 0) {
                snprintf(error_msg, msg_size, "Migration rate from pop %d to itself must be 0", i);
                return -1;
            }
        }
    }
    
    return 0;  /* Success */
}

/**
 * Print a summary of simulation parameters
 */
void params_print_summary(const SimulationParams *params, FILE *output) {
    assert(params != NULL);
    assert(output != NULL);
    
    fprintf(output, "=== Simulation Parameters Summary ===\n");
    
    /* Core parameters */
    fprintf(output, "Core:\n");
    fprintf(output, "  Samples: %d total across %d population(s)\n", 
            params->core.total_samples, params->core.num_populations);
    fprintf(output, "  Sample distribution: ");
    for (int i = 0; i < params->core.num_populations; i++) {
        fprintf(output, "%d ", params->core.sample_sizes[i]);
    }
    fprintf(output, "\n");
    fprintf(output, "  Replicates: %d\n", params->core.num_replicates);
    fprintf(output, "  Sites: %d\n", params->core.num_sites);
    if (params->core.seg_sites > 0) {
        fprintf(output, "  Segregating sites: %d (fixed)\n", params->core.seg_sites);
    }
    
    /* Evolutionary forces */
    fprintf(output, "\nEvolutionary Forces:\n");
    if (params->forces.theta > 0) {
        fprintf(output, "  Theta (4Nu): %.4f\n", params->forces.theta);
    }
    if (params->forces.rho > 0) {
        fprintf(output, "  Rho (4Nr): %.4f\n", params->forces.rho);
        if (params->forces.gene_conversion_ratio > 0) {
            fprintf(output, "  Gene conversion: %.2f%% of recombination events\n", 
                    params->forces.gene_conversion_ratio * 100);
            fprintf(output, "  Mean tract length: %d bp\n", params->forces.gc_tract_mean);
        }
    }
    
    /* Selection */
    if (params->selection.alpha != 0) {
        fprintf(output, "\nSelection:\n");
        fprintf(output, "  Selection coefficient (2Ns): %.4f\n", params->selection.alpha);
        fprintf(output, "  Sweep position: %.4f\n", params->selection.sweep_position);
        fprintf(output, "  Sweep mode: %c\n", params->selection.sweep_mode);
        if (params->selection.partial_sweep) {
            fprintf(output, "  Partial sweep to frequency: %.4f\n", 
                    params->selection.partial_sweep_final_freq);
        }
    }
    
    /* Demographics */
    if (params->demographics.num_events > 0) {
        fprintf(output, "\nDemographic Events: %d\n", params->demographics.num_events);
    }
    
    /* Output */
    fprintf(output, "\nOutput:\n");
    const char *format_names[] = {"MS", "FASTA", "NEXUS"};
    fprintf(output, "  Format: %s\n", format_names[params->output.format]);
    if (params->output.tree_sequence_output) {
        fprintf(output, "  Tree sequence output: %s\n", params->output.output_filename);
    }
    
    fprintf(output, "=====================================\n");
}

/**
 * Helper function to set symmetric migration
 */
void params_set_symmetric_migration(SimulationParams *params, double rate) {
    assert(params != NULL);
    
    params->forces.symmetric_migration = true;
    params->forces.migration_rate = rate;
    
    /* Set migration matrix */
    for (int i = 0; i < params->core.num_populations; i++) {
        for (int j = 0; j < params->core.num_populations; j++) {
            if (i != j) {
                params->forces.migration_matrix[i][j] = rate;
            } else {
                params->forces.migration_matrix[i][j] = 0.0;
            }
        }
    }
}

/**
 * Helper function to add a demographic event
 */
void params_add_demographic_event(SimulationParams *params, DemographicEvent *event) {
    assert(params != NULL);
    assert(event != NULL);
    
    /* Grow events array if needed */
    if (params->demographics.num_events >= params->demographics.events_capacity) {
        int new_capacity = params->demographics.events_capacity == 0 ? 10 : 
                          params->demographics.events_capacity * 2;
        DemographicEvent *new_events = realloc(params->demographics.events, 
                                              new_capacity * sizeof(DemographicEvent));
        if (!new_events) {
            fprintf(stderr, "Error: Failed to allocate memory for demographic events\n");
            return;
        }
        params->demographics.events = new_events;
        params->demographics.events_capacity = new_capacity;
    }
    
    /* Add the event */
    params->demographics.events[params->demographics.num_events] = *event;
    params->demographics.num_events++;
    
    /* Sort events by time */
    qsort(params->demographics.events, params->demographics.num_events, 
          sizeof(DemographicEvent), compare_events);
}

/**
 * Helper function to set up a simple selective sweep
 */
void params_set_simple_sweep(SimulationParams *params, double s, double position) {
    assert(params != NULL);
    
    params->selection.alpha = 2.0 * params->core.Ne * s;
    params->selection.sweep_position = position;
    params->selection.sweep_mode = SWEEP_DETERMINISTIC;
    params->selection.tau = 0.0;  /* Recent sweep */
    params->selection.f0 = 1.0 / (2.0 * params->core.Ne);
}

/**
 * Helper function to copy a SimulationParams structure
 */
SimulationParams* params_copy(const SimulationParams *params) {
    if (!params) return NULL;
    
    SimulationParams *copy = params_create();
    if (!copy) return NULL;
    
    /* Copy all fields */
    *copy = *params;
    
    /* Deep copy dynamic arrays */
    if (params->demographics.events) {
        copy->demographics.events = malloc(params->demographics.events_capacity * sizeof(DemographicEvent));
        if (!copy->demographics.events) {
            params_destroy(copy);
            return NULL;
        }
        memcpy(copy->demographics.events, params->demographics.events, 
               params->demographics.num_events * sizeof(DemographicEvent));
    }
    
    if (params->demographics.ancient_samples.sample_times) {
        int n = params->demographics.ancient_samples.num_ancient_samples;
        copy->demographics.ancient_samples.sample_times = malloc(n * sizeof(int));
        copy->demographics.ancient_samples.sample_pops = malloc(n * sizeof(int));
        if (!copy->demographics.ancient_samples.sample_times || 
            !copy->demographics.ancient_samples.sample_pops) {
            params_destroy(copy);
            return NULL;
        }
        memcpy(copy->demographics.ancient_samples.sample_times, 
               params->demographics.ancient_samples.sample_times, n * sizeof(int));
        memcpy(copy->demographics.ancient_samples.sample_pops, 
               params->demographics.ancient_samples.sample_pops, n * sizeof(int));
    }
    
    if (params->priors.distributions) {
        copy->priors.distributions = malloc(params->priors.num_priors * sizeof(PriorDist));
        if (!copy->priors.distributions) {
            params_destroy(copy);
            return NULL;
        }
        for (int i = 0; i < params->priors.num_priors; i++) {
            copy->priors.distributions[i] = params->priors.distributions[i];
            if (params->priors.distributions[i].parameter_name) {
                copy->priors.distributions[i].parameter_name = 
                    strdup(params->priors.distributions[i].parameter_name);
            }
        }
    }
    
    return copy;
}

/* Helper functions for command line parsing */
static const char* getNextArg(int argc, char **argv, int *args, const char *option) {
    (*args)++;
    if (*args >= argc) {
        fprintf(stderr, "Error: Option %s requires an argument\n", option);
        exit(1);
    }
    return argv[*args];
}

static int parseIntArg(int argc, char **argv, int *args, const char *option) {
    const char *arg = getNextArg(argc, argv, args, option);
    char *endptr;
    long val = strtol(arg, &endptr, 10);
    
    if (*endptr != '\0' || val > INT_MAX || val < INT_MIN) {
        fprintf(stderr, "Error: Invalid integer argument for %s: '%s'\n", option, arg);
        exit(1);
    }
    
    return (int)val;
}

static double parseDoubleArg(int argc, char **argv, int *args, const char *option) {
    const char *arg = getNextArg(argc, argv, args, option);
    char *endptr;
    double val = strtod(arg, &endptr);
    
    if (*endptr != '\0') {
        fprintf(stderr, "Error: Invalid numeric argument for %s: '%s'\n", option, arg);
        exit(1);
    }
    
    return val;
}

/* Sort demographic events by time */
static int compare_events(const void *a, const void *b) {
    const DemographicEvent *ea = (const DemographicEvent *)a;
    const DemographicEvent *eb = (const DemographicEvent *)b;
    
    if (ea->time < eb->time) return -1;
    if (ea->time > eb->time) return 1;
    return 0;
}

static void sort_events(SimulationParams *params) {
    if (params->demographics.num_events > 0) {
        qsort(params->demographics.events, params->demographics.num_events, 
              sizeof(DemographicEvent), compare_events);
    }
}

/**
 * Load parameters from command line arguments
 * Maintains backward compatibility with original discoal command line format
 */
int params_load_from_args(SimulationParams *params, int argc, char **argv) {
    assert(params != NULL);
    int args;
    int i, j;
    
    /* Check for configuration file options */
    if (argc >= 3) {
        /* Check for YAML config file */
        if (strcmp(argv[1], "-yaml") == 0 || strcmp(argv[1], "--yaml") == 0) {
            if (argc < 3) {
                fprintf(stderr, "Error: -yaml requires a filename\n");
                return -1;
            }
            /* Check if it's a hybrid config with Demes */
            int ret = yaml_load_params_with_demes(params, argv[2]);
            if (ret != 0) {
                fprintf(stderr, "Error loading YAML configuration from '%s'\n", argv[2]);
                return -1;
            }
            /* Process remaining command-line arguments starting from argv[3] */
            args = 3;
            goto process_options;
        }
        /* Check for Demes file */
        else if (strcmp(argv[1], "-demes") == 0 || strcmp(argv[1], "--demes") == 0) {
            if (argc < 6) {
                fprintf(stderr, "Error: -demes requires: filename samples replicates sites\n");
                return -1;
            }
            /* Load demographics from Demes file */
            int ret = demes_load_demographics(params, argv[2]);
            if (ret != 0) {
                fprintf(stderr, "Error loading Demes file '%s'\n", argv[2]);
                return -1;
            }
            /* Parse remaining positional arguments */
            params->core.total_samples = atoi(argv[3]);
            params->core.num_replicates = atoi(argv[4]);
            params->core.num_sites = atoi(argv[5]);
            
            if (params->core.total_samples > 65535) {
                fprintf(stderr, "Error: sampleSize > 65535. This exceeds the maximum supported.\n");
                return -1;
            }
            
            /* Distribute samples evenly across populations if not specified */
            if (params->core.num_populations > 1) {
                int samples_per_pop = params->core.total_samples / params->core.num_populations;
                int remainder = params->core.total_samples % params->core.num_populations;
                for (int i = 0; i < params->core.num_populations; i++) {
                    params->core.sample_sizes[i] = samples_per_pop;
                    if (i < remainder) {
                        params->core.sample_sizes[i]++;
                    }
                }
            } else {
                params->core.sample_sizes[0] = params->core.total_samples;
            }
            
            args = 6;
            goto process_options;
        }
    }
    
    if (argc < 4) {
        return -1;  /* Not enough arguments */
    }
    
    /* Parse positional arguments */
    params->core.total_samples = atoi(argv[1]);
    if (params->core.total_samples > 65535) {
        fprintf(stderr, "Error: sampleSize > 65535. This exceeds the maximum supported.\n");
        return -1;
    }
    
    params->core.num_replicates = atoi(argv[2]);
    params->core.num_sites = atoi(argv[3]);
    
    /* Set default single population */
    params->core.num_populations = 1;
    params->core.sample_sizes[0] = params->core.total_samples;
    
    args = 4;
    
process_options:
    /* Process optional arguments */
    while (args < argc) {
        if (argv[args][0] != '-') {
            fprintf(stderr, "Error: Unexpected argument '%s'\n", argv[args]);
            return -1;
        }
        
        switch(argv[args][1]) {
            case 't':
                if (argv[args][2] == 's') {  /* -ts for tree sequence */
                    params->output.tree_sequence_output = true;
                    args++;
                    if (args >= argc) {
                        fprintf(stderr, "Error: -ts requires filename\n");
                        return -1;
                    }
                    strncpy(params->output.output_filename, argv[args], PATH_MAX - 1);
                } else {  /* -t for theta */
                    params->forces.theta = parseDoubleArg(argc, argv, &args, "-t");
                }
                break;
                
            case 'r':
                params->forces.rho = parseDoubleArg(argc, argv, &args, "-r");
                break;
                
            case 'g':
                if (argv[args][2] == 'r') {  /* -gr for gene conversion ratio */
                    params->forces.gene_conversion_ratio = parseDoubleArg(argc, argv, &args, "-gr");
                    params->forces.gc_tract_mean = parseIntArg(argc, argv, &args, "-gr");
                } else {  /* -g for gene conversion */
                    double gc_rate = parseDoubleArg(argc, argv, &args, "-g");
                    params->forces.gc_tract_mean = parseIntArg(argc, argv, &args, "-g");
                    /* Convert absolute rate to ratio */
                    if (params->forces.rho > 0) {
                        params->forces.gene_conversion_ratio = gc_rate / params->forces.rho;
                    }
                }
                break;
                
            case 'a':
                params->selection.alpha = parseDoubleArg(argc, argv, &args, "-a");
                break;
                
            case 'x':
                params->selection.sweep_position = parseDoubleArg(argc, argv, &args, "-x");
                break;
                
            case 'N':
                params->core.Ne = (double)parseIntArg(argc, argv, &args, "-N");
                break;
                
            case 'p':  /* Multiple populations */
                params->core.num_populations = parseIntArg(argc, argv, &args, "-p");
                if (params->core.num_populations > MAXPOPS) {
                    fprintf(stderr, "Error: Too many populations (max %d)\n", MAXPOPS);
                    return -1;
                }
                
                params->core.total_samples = 0;
                for (i = 0; i < params->core.num_populations; i++) {
                    params->core.sample_sizes[i] = parseIntArg(argc, argv, &args, "-p");
                    params->core.total_samples += params->core.sample_sizes[i];
                }
                break;
                
            case 'M':  /* Symmetric migration */
                if (params->core.num_populations == 1) {
                    fprintf(stderr, "Error: -M requires multiple populations (use -p first)\n");
                    return -1;
                }
                params_set_symmetric_migration(params, parseDoubleArg(argc, argv, &args, "-M"));
                break;
                
            case 'm':  /* Pairwise migration */
                if (params->core.num_populations == 1) {
                    fprintf(stderr, "Error: -m requires multiple populations (use -p first)\n");
                    return -1;
                }
                i = parseIntArg(argc, argv, &args, "-m");
                j = parseIntArg(argc, argv, &args, "-m");
                params->forces.migration_matrix[i][j] = parseDoubleArg(argc, argv, &args, "-m");
                break;
                
            case 'e':  /* Demographic events */
                {
                    DemographicEvent event;
                    event.time = 0;
                    
                    switch(argv[args][2]) {
                        case 'n':  /* Population size change */
                            event.time = parseDoubleArg(argc, argv, &args, "-en") * 2.0;
                            event.type = EVENT_SIZE_CHANGE;
                            event.params.size_change.pop = parseIntArg(argc, argv, &args, "-en");
                            event.params.size_change.size = parseDoubleArg(argc, argv, &args, "-en");
                            break;
                            
                        case 'd':  /* Population split/join */
                        case 'j':  /* ms compatibility */
                            event.time = parseDoubleArg(argc, argv, &args, "-ed") * 2.0;
                            event.type = EVENT_JOIN;
                            event.params.join.pop1 = parseIntArg(argc, argv, &args, "-ed");
                            event.params.join.pop2 = parseIntArg(argc, argv, &args, "-ed");
                            event.params.join.dest = event.params.join.pop2;
                            break;
                            
                        case 'a':  /* Admixture */
                            event.time = parseDoubleArg(argc, argv, &args, "-ea") * 2.0;
                            event.type = EVENT_ADMIX;
                            event.params.admix.to = parseIntArg(argc, argv, &args, "-ea");
                            int pop2 = parseIntArg(argc, argv, &args, "-ea");
                            int pop3 = parseIntArg(argc, argv, &args, "-ea");
                            event.params.admix.from = pop2;  /* Main source */
                            event.params.admix.prop = parseDoubleArg(argc, argv, &args, "-ea");
                            /* Note: pop3 is the secondary source, prop goes to pop2 */
                            break;
                    }
                    
                    if (event.time >= 0) {
                        params_add_demographic_event(params, &event);
                    }
                }
                break;
                
            case 'w':  /* Sweep modes */
                switch(argv[args][2]) {
                    case 'd':
                        params->selection.sweep_mode = SWEEP_DETERMINISTIC;
                        break;
                    case 's':
                        params->selection.sweep_mode = SWEEP_STOCHASTIC;
                        break;
                    case 'n':
                        params->selection.sweep_mode = SWEEP_NEUTRAL;
                        break;
                }
                params->selection.tau = parseDoubleArg(argc, argv, &args, "-w") * 2.0;
                break;
                
            case 'f':  /* Starting frequency for soft sweep */
                params->selection.f0 = parseDoubleArg(argc, argv, &args, "-f");
                params->selection.soft_sweep = true;
                break;
                
            case 'u':  /* Adaptive mutation rate */
                params->selection.adaptive_mutation_rate = parseDoubleArg(argc, argv, &args, "-u");
                break;
                
            case 'c':  /* Partial sweep */
                params->selection.partial_sweep = true;
                params->selection.partial_sweep_final_freq = parseDoubleArg(argc, argv, &args, "-c");
                if (params->selection.sweep_mode != SWEEP_STOCHASTIC) {
                    params->selection.sweep_mode = SWEEP_STOCHASTIC;
                }
                break;
                
            case 'R':  /* Recurrent sweep */
                params->selection.recurrent_sweep = true;
                params->selection.sweep_mode = SWEEP_STOCHASTIC;
                params->selection.recurrent_sweep_rate = parseDoubleArg(argc, argv, &args, "-R");
                params->selection.adaptive_mutation_rate = params->selection.recurrent_sweep_rate;
                break;
                
            case 'F':  /* Full tree sequence mode */
                if (!params->output.tree_sequence_output) {
                    fprintf(stderr, "Error: -F requires -ts to be specified first\n");
                    return -1;
                }
                params->output.minimal_tree_seq = false;
                break;
                
            case 'h':  /* Hide partial sweep SNP */
                params->output.hide_partial_snp = true;
                break;
                
            case 'd':  /* Random seeds */
                params->random.seed1 = (unsigned long)parseIntArg(argc, argv, &args, "-d");
                params->random.seed2 = (unsigned long)parseIntArg(argc, argv, &args, "-d");
                break;
                
            /* TODO: Add prior parsing (-P flags) */
            
            default:
                fprintf(stderr, "Error: Unknown option '-%c'\n", argv[args][1]);
                return -1;
        }
        
        args++;
    }
    
    /* Sort events by time */
    sort_events(params);
    
    /* Validate parameters */
    char error_msg[256];
    if (params_validate(params, error_msg, sizeof(error_msg)) != 0) {
        fprintf(stderr, "Parameter validation error: %s\n", error_msg);
        return -1;
    }
    
    return 0;
}
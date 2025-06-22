/**
 * params.c - Implementation of centralized parameter management for discoal
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "params.h"

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
    
    /* Initialize population sizes to 1.0 */
    for (int i = 0; i < MAXPOPS; i++) {
        params->demographics.pop_sizes.sizes[i] = 1.0;
        params->demographics.pop_sizes.initial_sizes[i] = 1.0;
        params->demographics.pop_sizes.growth_rates[i] = 0.0;
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
    
    /* Output defaults */
    params->output.format = OUTPUT_MS;
    params->output.tree_sequence_output = false;
    strcpy(params->output.output_filename, "discoal_output.trees");
    params->output.minimal_tree_seq = false;
    params->output.hide_singletons = false;
    params->output.hide_partial_snp = false;
    
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
    
    /* TODO: Sort events by time */
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
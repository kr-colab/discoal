/**
 * demes_loader.c - Implementation of Demes integration for discoal
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "demes_loader.h"
#include "yaml_loader.h"

/* Helper to convert time from Demes units to coalescent units */
double demes_time_to_coalescent(double time, struct demes_graph *graph) {
    /* For discoal, we require time_units = "4N" */
    if (strcmp((char*)graph->time_units, "4N") == 0) {
        /* Already in coalescent units */
        return time;
    } else {
        fprintf(stderr, "Error: discoal requires time_units='4N' in Demes files\n");
        fprintf(stderr, "Found time_units='%s' instead\n", graph->time_units);
        return -1;
    }
}

/* Helper to find or create population index for a deme name */
static int get_or_create_population_index(SimulationParams *params, 
                                         const char *deme_name,
                                         char **deme_name_map,
                                         int *next_pop_idx) {
    /* Check if we already have this deme */
    for (int i = 0; i < *next_pop_idx; i++) {
        if (strcmp(deme_name_map[i], deme_name) == 0) {
            return i;
        }
    }
    
    /* New deme - add it */
    if (*next_pop_idx >= MAXPOPS) {
        fprintf(stderr, "Error: Too many populations (max %d)\n", MAXPOPS);
        return -1;
    }
    
    deme_name_map[*next_pop_idx] = strdup(deme_name);
    return (*next_pop_idx)++;
}

/* Convert Demes graph to discoal demographics */
int demes_convert_graph(SimulationParams *params, struct demes_graph *graph) {
    /* Map deme names to population indices */
    char *deme_name_map[MAXPOPS] = {0};
    int next_pop_idx = 0;
    
    /* First pass: count populations and create mapping */
    params->core.num_populations = graph->n_demes;
    if (params->core.num_populations > MAXPOPS) {
        fprintf(stderr, "Error: Too many demes (%zu, max %d)\n", 
                graph->n_demes, MAXPOPS);
        return -1;
    }
    
    /* Check time units */
    if (strcmp((char*)graph->time_units, "4N") != 0) {
        fprintf(stderr, "Error: discoal requires time_units='4N' in Demes files\n");
        fprintf(stderr, "Found time_units='%s' instead\n", graph->time_units);
        return -1;
    }
    
    /* Initialize population sizes to 0 */
    for (int i = 0; i < MAXPOPS; i++) {
        params->demographics.pop_sizes.sizes[i] = 0;
        params->demographics.pop_sizes.initial_sizes[i] = 0;
    }
    
    double reference_size = 0;  /* Will be set to first population's present-day size */
    
    /* Process each deme */
    for (size_t i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        int pop_idx = get_or_create_population_index(params, (char*)deme->name,
                                                     deme_name_map, &next_pop_idx);
        if (pop_idx < 0) {
            goto cleanup;
        }
        
        /* Use the most recent epoch (last in array) for current size */
        if (deme->n_epochs > 0) {
            /* In Demes, epochs are ordered from past to present */
            struct demes_epoch *recent_epoch = &deme->epochs[deme->n_epochs - 1];
            double present_day_size;
            
            /* For the current population size, use the end_size if exponential, else start_size */
            if (recent_epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL && 
                recent_epoch->end_time == 0.0) {
                /* Present-day size is end_size for exponential growth to present */
                present_day_size = recent_epoch->end_size;
            } else {
                present_day_size = recent_epoch->start_size;
            }
            
            /* Store absolute sizes for now - will convert to relative after */
            params->demographics.pop_sizes.sizes[pop_idx] = present_day_size;
            params->demographics.pop_sizes.initial_sizes[pop_idx] = present_day_size;
            
            /* Set reference size to first population's size (pop index 0) */
            if (pop_idx == 0) {
                reference_size = present_day_size;
            }
        }
    }
    
    /* Validate reference size */
    if (reference_size <= 0) {
        fprintf(stderr, "Error: No valid reference population size found\n");
        goto cleanup;
    }
    
    /* Convert all population sizes to be relative to reference population */
    for (int i = 0; i < params->core.num_populations; i++) {
        params->demographics.pop_sizes.sizes[i] /= reference_size;
        params->demographics.pop_sizes.initial_sizes[i] /= reference_size;
    }
    
    /* Convert demographic events */
    /* Count events first */
    int num_events = 0;
    
    /* Size changes from epochs */
    for (size_t i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        /* Each epoch transition is an event */
        if (deme->n_epochs > 1) {
            num_events += deme->n_epochs - 1;
        }
    }
    
    /* Population splits/joins from ancestry */
    for (size_t i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (deme->n_ancestors > 0) {
            num_events++; /* Split or admixture event */
        }
    }
    
    /* Migration events */
    num_events += graph->n_migrations;
    
    /* Pulse events */
    num_events += graph->n_pulses;
    
    /* Allocate events */
    if (num_events > 0) {
        params->demographics.events = calloc(num_events, sizeof(DemographicEvent));
        if (!params->demographics.events) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            goto cleanup;
        }
        params->demographics.events_capacity = num_events;
    }
    
    int event_idx = 0;
    
    /* Convert epoch transitions to size change events */
    for (size_t i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        int pop_idx = get_or_create_population_index(params, (char*)deme->name,
                                                     deme_name_map, &next_pop_idx);
        if (pop_idx < 0) continue;
        
        /* Process epoch boundaries where size changes */
        for (size_t e = 1; e < deme->n_epochs; e++) {
            struct demes_epoch *prev_epoch = &deme->epochs[e-1];
            struct demes_epoch *curr_epoch = &deme->epochs[e];
            
            /* Check for instantaneous size change at epoch boundary */
            /* In Demes, the transition happens at prev_epoch.end_time */
            if (curr_epoch->start_size != prev_epoch->end_size || 
                (prev_epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL &&
                 prev_epoch->end_size != curr_epoch->start_size)) {
                
                DemographicEvent *event = &params->demographics.events[event_idx++];
                event->time = demes_time_to_coalescent(prev_epoch->end_time, graph);
                event->type = EVENT_SIZE_CHANGE;
                event->params.size_change.pop = pop_idx;
                
                /* Size in discoal is relative to reference population */
                event->params.size_change.size = curr_epoch->start_size / reference_size;
            }
        }
    }
    
    /* TODO: Implement remaining conversions:
     * - Exponential growth to discrete size changes (discoal has no -eg)
     * - Ancestry relationships to split/join events
     * - Migrations to migration rate change events
     * - Pulses to admixture events
     * 
     * CRITICAL LIMITATION: discoal does not support exponential growth!
     * Only discrete population size changes via -en events are supported.
     * Exponential growth in Demes models must be approximated.
     */
    
    params->demographics.num_events = event_idx;
    
    /* Set up simple symmetric migration if specified */
    if (graph->n_migrations > 0 && params->core.num_populations > 1) {
        /* For simplicity, use the first migration rate for all pairs */
        double mig_rate = graph->migrations[0].rate;
        params->forces.symmetric_migration = true;
        params->forces.migration_rate = mig_rate;
        params_set_symmetric_migration(params, mig_rate);
    }
    
    /* Clean up */
    for (int i = 0; i < next_pop_idx; i++) {
        free(deme_name_map[i]);
    }
    
    return 0;

cleanup:
    for (int i = 0; i < next_pop_idx; i++) {
        free(deme_name_map[i]);
    }
    return -1;
}

/* Load demographics from a Demes file */
int demes_load_demographics(SimulationParams *params, const char *filename) {
    struct demes_graph *graph = NULL;
    int ret;
    
    /* Load the Demes file */
    ret = demes_graph_load(filename, &graph);
    if (ret != 0) {
        fprintf(stderr, "Error loading Demes file '%s': error code %d\n", 
                filename, ret);
        return -1;
    }
    
    /* Convert to discoal format */
    ret = demes_convert_graph(params, graph);
    
    /* Free the Demes graph */
    demes_graph_free(graph);
    
    return ret;
}

/* Extended YAML loading with Demes support */
int yaml_load_params_with_demes(SimulationParams *params, const char *config_file) {
    FILE *file = fopen(config_file, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open YAML file '%s'\n", config_file);
        return -1;
    }
    
    yaml_parser_t parser;
    yaml_document_t document;
    
    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Error: Failed to initialize YAML parser\n");
        fclose(file);
        return -1;
    }
    
    yaml_parser_set_input_file(&parser, file);
    
    /* First load standard parameters */
    int ret = yaml_parse_document(&parser, params);
    if (ret != 0) {
        yaml_parser_delete(&parser);
        fclose(file);
        return ret;
    }
    
    /* Now check for demographics section with demes_file */
    /* Reset parser to beginning of file */
    fseek(file, 0, SEEK_SET);
    yaml_parser_delete(&parser);
    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, file);
    
    if (!yaml_parser_load(&parser, &document)) {
        yaml_parser_delete(&parser);
        fclose(file);
        return -1;
    }
    
    /* Look for demographics->demes_file */
    yaml_node_t *root = yaml_document_get_root_node(&document);
    if (root && root->type == YAML_MAPPING_NODE) {
        yaml_node_pair_t *pair;
        for (pair = root->data.mapping.pairs.start; 
             pair < root->data.mapping.pairs.top; pair++) {
            yaml_node_t *key = yaml_document_get_node(&document, pair->key);
            yaml_node_t *value = yaml_document_get_node(&document, pair->value);
            
            if (key->type == YAML_SCALAR_NODE && 
                strcmp((char*)key->data.scalar.value, "demographics") == 0 &&
                value->type == YAML_MAPPING_NODE) {
                
                /* Found demographics section */
                yaml_node_pair_t *demo_pair;
                for (demo_pair = value->data.mapping.pairs.start;
                     demo_pair < value->data.mapping.pairs.top; demo_pair++) {
                    yaml_node_t *demo_key = yaml_document_get_node(&document, demo_pair->key);
                    yaml_node_t *demo_value = yaml_document_get_node(&document, demo_pair->value);
                    
                    if (demo_key->type == YAML_SCALAR_NODE &&
                        strcmp((char*)demo_key->data.scalar.value, "demes_file") == 0 &&
                        demo_value->type == YAML_SCALAR_NODE) {
                        
                        /* Load demographics from Demes file */
                        const char *demes_file = (char*)demo_value->data.scalar.value;
                        ret = demes_load_demographics(params, demes_file);
                        if (ret != 0) {
                            fprintf(stderr, "Error loading demographics from '%s'\n", 
                                    demes_file);
                        }
                        break;
                    }
                }
                break;
            }
        }
    }
    
    yaml_document_delete(&document);
    yaml_parser_delete(&parser);
    fclose(file);
    
    return ret;
}

/* Get population index for a deme name */
int demes_get_population_index(SimulationParams *params, const char *deme_name) {
    /* This would need to maintain a mapping between deme names and indices */
    /* For now, return -1 (not implemented) */
    return -1;
}
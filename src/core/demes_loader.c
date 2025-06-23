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

/* Helper to find population index for a deme name (no creation) */
static int find_population_index(const char *deme_name,
                                char **deme_name_map,
                                int num_pops) {
    for (int i = 0; i < num_pops; i++) {
        if (deme_name_map[i] && strcmp(deme_name_map[i], deme_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* Comparison function for sorting demographic events by time */
static int compare_demographic_events(const void *a, const void *b) {
    const DemographicEvent *event_a = (const DemographicEvent *)a;
    const DemographicEvent *event_b = (const DemographicEvent *)b;
    
    /* Sort by time (ascending order - earliest events first) */
    if (event_a->time < event_b->time) return -1;
    if (event_a->time > event_b->time) return 1;
    
    /* If times are equal, maintain a consistent order by event type */
    /* Priority: size changes < splits/joins < admixture < migration */
    const int type_priority[] = {
        ['n'] = 0,  /* SIZE_CHANGE - do these first */
        ['p'] = 1,  /* JOIN/SPLIT */
        ['a'] = 2,  /* ADMIX */
        ['m'] = 3   /* MIGRATION_CHANGE */
    };
    
    /* Map EVENT_* enums to characters for legacy compatibility */
    char type_a = 'n', type_b = 'n';
    switch (event_a->type) {
        case EVENT_SIZE_CHANGE: type_a = 'n'; break;
        case EVENT_JOIN: type_a = 'p'; break;
        case EVENT_ADMIX: type_a = 'a'; break;
        case EVENT_MIGRATION_CHANGE: type_a = 'm'; break;
        default: type_a = 'n';
    }
    switch (event_b->type) {
        case EVENT_SIZE_CHANGE: type_b = 'n'; break;
        case EVENT_JOIN: type_b = 'p'; break;
        case EVENT_ADMIX: type_b = 'a'; break;
        case EVENT_MIGRATION_CHANGE: type_b = 'm'; break;
        default: type_b = 'n';
    }
    
    int priority_a = (type_a >= 0 && type_a < 128) ? type_priority[(int)type_a] : 99;
    int priority_b = (type_b >= 0 && type_b < 128) ? type_priority[(int)type_b] : 99;
    
    return priority_a - priority_b;
}

/* Convert Demes graph to discoal demographics */
int demes_convert_graph(SimulationParams *params, struct demes_graph *graph) {
    return demes_convert_graph_with_names(params, graph, NULL, NULL);
}

/* Convert Demes graph to discoal demographics with name mapping */
int demes_convert_graph_with_names(SimulationParams *params, struct demes_graph *graph,
                                   char ***deme_name_map_out, int *num_demes_out) {
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
        
        /* Check for unsupported features in epochs */
        for (size_t e = 0; e < deme->n_epochs; e++) {
            struct demes_epoch *epoch = &deme->epochs[e];
            
            /* Warn about selfing rate */
            if (epoch->selfing_rate > 0) {
                fprintf(stderr, "Warning: Selfing rate %.3f in deme '%s' epoch %zu is not supported by discoal and will be ignored\n",
                        epoch->selfing_rate, deme->name, e);
            }
            
            /* Warn about cloning rate */
            if (epoch->cloning_rate > 0) {
                fprintf(stderr, "Warning: Cloning rate %.3f in deme '%s' epoch %zu is not supported by discoal and will be ignored\n",
                        epoch->cloning_rate, deme->name, e);
            }
            
            /* Warn about exponential growth */
            if (epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL) {
                fprintf(stderr, "Warning: Exponential growth in deme '%s' epoch %zu (%.0f to %.0f) will be approximated with discrete size changes\n",
                        deme->name, e, epoch->start_size, epoch->end_size);
            }
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
    
    /* Migration events - need to check if they form symmetric pattern first */
    bool count_symmetric = false;
    if (graph->n_migrations > 0 && params->core.num_populations > 1) {
        /* Quick check: if we have exactly n*(n-1) migrations, might be symmetric */
        int expected_symmetric = params->core.num_populations * (params->core.num_populations - 1);
        if (graph->n_migrations == expected_symmetric) {
            /* Assume symmetric for counting purposes - we'll verify later */
            count_symmetric = true;
            
            /* For symmetric migration starting at time > 0, we need start events */
            if (graph->migrations[0].start_time > 0) {
                num_events += expected_symmetric;
            }
            
            /* For time-bounded symmetric migration, we need stop events */
            if (graph->migrations[0].end_time > 0 && 
                graph->migrations[0].end_time != graph->migrations[0].start_time) {
                num_events += expected_symmetric;
            }
        }
    }
    
    if (!count_symmetric) {
        /* Count individual migration events */
        for (size_t i = 0; i < graph->n_migrations; i++) {
            num_events++; /* Start event */
            /* Check if we need an end event */
            if (graph->migrations[i].end_time > 0 && 
                graph->migrations[i].end_time != graph->migrations[i].start_time) {
                num_events++; /* End event */
            }
        }
    }
    
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
        int pop_idx = find_population_index((char*)deme->name,
                                          deme_name_map, next_pop_idx);
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
    
    /* Convert ancestry relationships to split/join events */
    for (size_t i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        
        if (deme->n_ancestors > 0) {
            int daughter_idx = find_population_index((char*)deme->name,
                                                   deme_name_map, next_pop_idx);
            if (daughter_idx < 0) continue;
            
            /* In discoal, splits are represented as joins going backward in time */
            /* We use -ed (EVENT_JOIN) to merge daughter into parent at start_time */
            
            if (deme->n_ancestors == 1) {
                /* Simple split: one parent, one daughter */
                int parent_idx = find_population_index((char*)deme->ancestors[0]->name,
                                                     deme_name_map, next_pop_idx);
                if (parent_idx < 0) continue;
                
                DemographicEvent *event = &params->demographics.events[event_idx++];
                event->time = demes_time_to_coalescent(deme->start_time, graph);
                event->type = EVENT_JOIN;
                event->params.join.pop1 = daughter_idx;
                event->params.join.pop2 = parent_idx;
                event->params.join.dest = parent_idx;
            } else if (deme->n_ancestors == 2) {
                /* Admixture event: two parents contribute to daughter */
                /* discoal uses -ea for admixture events */
                int parent1_idx = find_population_index((char*)deme->ancestors[0]->name,
                                                      deme_name_map, next_pop_idx);
                int parent2_idx = find_population_index((char*)deme->ancestors[1]->name,
                                                      deme_name_map, next_pop_idx);
                if (parent1_idx < 0 || parent2_idx < 0) continue;
                
                DemographicEvent *event = &params->demographics.events[event_idx++];
                event->time = demes_time_to_coalescent(deme->start_time, graph);
                event->type = EVENT_ADMIX;
                event->params.admix.from = daughter_idx;
                event->params.admix.to = parent1_idx;
                /* Proportion from parent1 */
                event->params.admix.prop = deme->proportions[0];
                
                /* Note: discoal's -ea only supports two-way admixture with one parent
                 * receiving all lineages. The second parent contributes via migration.
                 * This is a limitation compared to Demes' multi-way admixture. */
            } else {
                /* Multi-way admixture not supported by discoal */
                fprintf(stderr, "Warning: Multi-way admixture with %zu ancestors "
                       "not supported for deme '%s'\n", deme->n_ancestors, deme->name);
            }
        }
    }
    
    /* Warn about pulse events */
    if (graph->n_pulses > 0) {
        fprintf(stderr, "Warning: %zu pulse event(s) found in Demes file. Pulse events are not supported by discoal and will be ignored.\n", 
                graph->n_pulses);
        for (size_t i = 0; i < graph->n_pulses; i++) {
            struct demes_pulse *pulse = &graph->pulses[i];
            fprintf(stderr, "  - Pulse at time %.4f from ", pulse->time);
            for (size_t s = 0; s < pulse->n_sources; s++) {
                fprintf(stderr, "%s%.3f", 
                        s > 0 ? "+" : "",
                        pulse->proportions[s]);
                fprintf(stderr, "*%s", pulse->sources[s]->name);
            }
            fprintf(stderr, " into %s\n", pulse->dest->name);
        }
    }
    
    /* Check if all migrations form a symmetric pattern */
    bool all_symmetric = true;
    double symmetric_rate = -1;
    double symmetric_start = -1;
    double symmetric_end = -1;
    
    if (graph->n_migrations > 0 && params->core.num_populations > 1) {
        /* First pass: check if migrations form symmetric pairs */
        for (size_t i = 0; i < graph->n_migrations; i++) {
            struct demes_migration *mig1 = &graph->migrations[i];
            bool found_pair = false;
            
            /* Look for reciprocal migration */
            for (size_t j = 0; j < graph->n_migrations; j++) {
                if (i == j) continue;
                struct demes_migration *mig2 = &graph->migrations[j];
                
                /* Check if this is the reciprocal migration */
                if (mig1->source == mig2->dest && mig1->dest == mig2->source &&
                    mig1->rate == mig2->rate && mig1->start_time == mig2->start_time &&
                    mig1->end_time == mig2->end_time) {
                    found_pair = true;
                    
                    /* Check if all pairs have the same rate and timing */
                    if (symmetric_rate < 0) {
                        symmetric_rate = mig1->rate;
                        symmetric_start = mig1->start_time;
                        symmetric_end = mig1->end_time;
                    } else if (mig1->rate != symmetric_rate || 
                              mig1->start_time != symmetric_start ||
                              mig1->end_time != symmetric_end) {
                        all_symmetric = false;
                        break;
                    }
                    break;
                }
            }
            
            if (!found_pair) {
                all_symmetric = false;
                break;
            }
        }
        
        /* Check if we have migrations between all population pairs */
        if (all_symmetric) {
            int expected_migrations = params->core.num_populations * (params->core.num_populations - 1);
            if (graph->n_migrations != expected_migrations) {
                all_symmetric = false;
            }
        }
    } else {
        all_symmetric = false;
    }
    
    /* If all migrations are symmetric, use discoal's symmetric migration */
    if (all_symmetric && symmetric_rate > 0) {
        fprintf(stderr, "Note: Detected symmetric migration pattern (rate=%.6f) - using efficient symmetric migration\n", 
                symmetric_rate);
        
        /* Check if migrations start at present (time 0 or infinity in Demes) */
        bool starts_at_present = (symmetric_start == 0 || isinf(symmetric_start));
        
        /* Set up symmetric migration at time 0 */
        if (starts_at_present) {
            params->forces.symmetric_migration = true;
            params->forces.migration_rate = symmetric_rate;
            params_set_symmetric_migration(params, symmetric_rate);
            
            /* Skip creating individual migration events for present-day symmetric migration */
        } else {
            /* Create migration events to start symmetric migration at specified time */
            for (int i = 0; i < params->core.num_populations; i++) {
                for (int j = 0; j < params->core.num_populations; j++) {
                    if (i != j) {
                        DemographicEvent *event = &params->demographics.events[event_idx++];
                        event->time = demes_time_to_coalescent(symmetric_start, graph);
                        event->type = EVENT_MIGRATION_CHANGE;
                        event->params.migration.from = i;
                        event->params.migration.to = j;
                        event->params.migration.rate = symmetric_rate;
                    }
                }
            }
        }
        
        /* If time-bounded, create stop events */
        if (symmetric_end > 0 && symmetric_end != symmetric_start && !isinf(symmetric_end)) {
            fprintf(stderr, "Note: Time-bounded symmetric migration (%.4f-%.4f) will be stopped at end time\n",
                    symmetric_start, symmetric_end);
            for (int i = 0; i < params->core.num_populations; i++) {
                for (int j = 0; j < params->core.num_populations; j++) {
                    if (i != j) {
                        DemographicEvent *event = &params->demographics.events[event_idx++];
                        event->time = demes_time_to_coalescent(symmetric_end, graph);
                        event->type = EVENT_MIGRATION_CHANGE;
                        event->params.migration.from = i;
                        event->params.migration.to = j;
                        event->params.migration.rate = 0.0;
                    }
                }
            }
        }
        
        /* Skip processing individual migrations if symmetric pattern detected */
        goto skip_individual_migrations;
    }
    
    /* Convert migration events individually */
    for (size_t i = 0; i < graph->n_migrations; i++) {
        struct demes_migration *mig = &graph->migrations[i];
        
        int source_idx = find_population_index((char*)mig->source->name,
                                             deme_name_map, next_pop_idx);
        int dest_idx = find_population_index((char*)mig->dest->name,
                                           deme_name_map, next_pop_idx);
        
        if (source_idx < 0 || dest_idx < 0) {
            fprintf(stderr, "Warning: Skipping migration between %s and %s (population not found)\n",
                    mig->source->name, mig->dest->name);
            continue;
        }
        
        /* discoal uses -m for migration rate changes */
        /* Note: discoal migration is backward in time (dest to source) */
        DemographicEvent *event = &params->demographics.events[event_idx++];
        
        /* Use start_time as the event time (when migration begins) */
        event->time = demes_time_to_coalescent(mig->start_time, graph);
        event->type = EVENT_MIGRATION_CHANGE;
        /* In discoal, migration is from dest to source (backward in time) */
        event->params.migration.from = dest_idx;
        event->params.migration.to = source_idx;
        event->params.migration.rate = mig->rate;
        
        /* If migration has an end time, create another event to stop it */
        if (mig->end_time > 0 && mig->end_time != mig->start_time) {
            fprintf(stderr, "Note: Time-bounded migration from %s to %s (%.4f-%.4f) approximated with start/stop events\n",
                    mig->source->name, mig->dest->name, mig->start_time, mig->end_time);
            DemographicEvent *end_event = &params->demographics.events[event_idx++];
            end_event->time = demes_time_to_coalescent(mig->end_time, graph);
            end_event->type = EVENT_MIGRATION_CHANGE;
            end_event->params.migration.from = dest_idx;
            end_event->params.migration.to = source_idx;
            end_event->params.migration.rate = 0.0;  /* Stop migration */
        }
    }
    
skip_individual_migrations:
    /* TODO: Implement remaining conversions:
     * - Exponential growth to discrete size changes (discoal has no -eg)
     * - Pulses to admixture events
     * 
     * CRITICAL LIMITATION: discoal does not support exponential growth!
     * Only discrete population size changes via -en events are supported.
     * Exponential growth in Demes models must be approximated.
     */
    
    /* Sort events by time (required by discoal) */
    if (event_idx > 1) {
        qsort(params->demographics.events, event_idx, sizeof(DemographicEvent), 
              compare_demographic_events);
    }
    
    params->demographics.num_events = event_idx;
    
    /* Note: We don't set up symmetric migration here because:
     * 1. Demes specifies individual migration events, not symmetric rates
     * 2. We've already converted migrations to EVENT_MIGRATION_CHANGE events
     * 3. If users want symmetric migration, they should use discoal's -M flag
     */
    
    /* Return name map if requested */
    if (deme_name_map_out && num_demes_out) {
        *deme_name_map_out = calloc(next_pop_idx, sizeof(char*));
        if (*deme_name_map_out) {
            for (int i = 0; i < next_pop_idx; i++) {
                (*deme_name_map_out)[i] = deme_name_map[i];
                deme_name_map[i] = NULL; /* Transfer ownership */
            }
            *num_demes_out = next_pop_idx;
        }
    }
    
    /* Clean up (only if not transferred) */
    for (int i = 0; i < next_pop_idx; i++) {
        if (deme_name_map[i]) {
            free(deme_name_map[i]);
        }
    }
    
    return 0;

cleanup:
    for (int i = 0; i < next_pop_idx; i++) {
        if (deme_name_map[i]) {
            free(deme_name_map[i]);
        }
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
    char **deme_name_map = NULL;
    int num_demes = 0;
    
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
                        
                        /* Load demographics from Demes file with name mapping */
                        const char *demes_file = (char*)demo_value->data.scalar.value;
                        struct demes_graph *graph = NULL;
                        
                        /* Load the Demes file */
                        ret = demes_graph_load(demes_file, &graph);
                        if (ret != 0) {
                            fprintf(stderr, "Error loading Demes file '%s': error code %d\n", 
                                    demes_file, ret);
                            break;
                        }
                        
                        /* Convert to discoal format with name mapping */
                        ret = demes_convert_graph_with_names(params, graph, &deme_name_map, &num_demes);
                        
                        /* Free the Demes graph */
                        demes_graph_free(graph);
                        
                        if (ret != 0) {
                            fprintf(stderr, "Error converting demographics from '%s'\n", 
                                    demes_file);
                        }
                        break;
                    }
                }
                break;
            }
        }
    }
    
    /* Apply sample specifications if we have a Demes name map */
    if (deme_name_map && num_demes > 0 && params->demographics.sample_specs) {
        ret = demes_apply_sample_specs(params, deme_name_map, num_demes);
        if (ret != 0) {
            fprintf(stderr, "Error applying sample specifications\n");
        }
        
        /* Free the name map */
        for (int i = 0; i < num_demes; i++) {
            if (deme_name_map[i]) {
                free(deme_name_map[i]);
            }
        }
        free(deme_name_map);
    }
    
    yaml_document_delete(&document);
    yaml_parser_delete(&parser);
    fclose(file);
    
    return ret;
}

/* Apply sample specifications to the simulation parameters */
int demes_apply_sample_specs(SimulationParams *params, char **deme_name_map, int num_demes) {
    if (!params->demographics.sample_specs || params->demographics.num_sample_specs == 0) {
        /* No sample specs provided */
        return 0;
    }
    
    /* Reset all sample sizes to 0 */
    params->core.total_samples = 0;
    for (int i = 0; i < params->core.num_populations; i++) {
        params->core.sample_sizes[i] = 0;
    }
    
    /* Apply each sample specification */
    for (int i = 0; i < params->demographics.num_sample_specs; i++) {
        SampleSpec *spec = &params->demographics.sample_specs[i];
        
        /* Find population index for the named population */
        int pop_idx = -1;
        for (int j = 0; j < num_demes; j++) {
            if (deme_name_map[j] && strcmp(deme_name_map[j], spec->population) == 0) {
                pop_idx = j;
                break;
            }
        }
        
        if (pop_idx < 0) {
            fprintf(stderr, "Error: Sample specification references unknown population '%s'\n",
                    spec->population);
            return -1;
        }
        
        /* Apply the sample specification */
        if (spec->time == 0.0) {
            /* Present-day samples */
            params->core.sample_sizes[pop_idx] += spec->size;
            params->core.total_samples += spec->size;
        } else {
            /* Ancient samples - need to create an event */
            /* TODO: Implement ancient sample support */
            fprintf(stderr, "Warning: Ancient samples (time > 0) not yet implemented for population '%s'\n",
                    spec->population);
        }
    }
    
    return 0;
}

/* Get population index for a deme name */
int demes_get_population_index(SimulationParams *params, const char *deme_name) {
    /* This would need to maintain a mapping between deme names and indices */
    /* For now, return -1 (not implemented) */
    return -1;
}
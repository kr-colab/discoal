/**
 * yaml_loader.c - YAML configuration loading using libyaml
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "yaml_loader.h"
#include "params.h"

/* Helper to get scalar value from YAML node */
static const char* get_scalar_value(yaml_node_t *node) {
    if (node->type == YAML_SCALAR_NODE) {
        return (const char*)node->data.scalar.value;
    }
    return NULL;
}

/* Helper to parse a double from YAML node */
static int parse_double(yaml_node_t *node, double *value) {
    const char *str = get_scalar_value(node);
    if (!str) return -1;
    
    char *endptr;
    *value = strtod(str, &endptr);
    return (*endptr == '\0') ? 0 : -1;
}

/* Helper to parse an int from YAML node */
static int parse_int(yaml_node_t *node, int *value) {
    const char *str = get_scalar_value(node);
    if (!str) return -1;
    
    char *endptr;
    long val = strtol(str, &endptr, 10);
    if (*endptr != '\0' || val > INT_MAX || val < INT_MIN) {
        return -1;
    }
    *value = (int)val;
    return 0;
}

/* Helper to parse a boolean from YAML node */
static int parse_bool(yaml_node_t *node, bool *value) {
    const char *str = get_scalar_value(node);
    if (!str) return -1;
    
    if (strcasecmp(str, "true") == 0 || strcasecmp(str, "yes") == 0 || 
        strcasecmp(str, "on") == 0 || strcmp(str, "1") == 0) {
        *value = true;
    } else {
        *value = false;
    }
    return 0;
}

/* Parse simulation section */
static int parse_simulation_section(yaml_document_t *document, yaml_node_t *node, SimulationParams *params) {
    if (node->type != YAML_MAPPING_NODE) return -1;
    
    yaml_node_pair_t *pair;
    for (pair = node->data.mapping.pairs.start; pair < node->data.mapping.pairs.top; pair++) {
        yaml_node_t *key = yaml_document_get_node(document, pair->key);
        yaml_node_t *value = yaml_document_get_node(document, pair->value);
        
        const char *key_str = get_scalar_value(key);
        if (!key_str) continue;
        
        if (strcmp(key_str, "samples") == 0) {
            parse_int(value, &params->core.total_samples);
            /* For single population, set sample size */
            if (params->core.num_populations == 1) {
                params->core.sample_sizes[0] = params->core.total_samples;
            }
        } else if (strcmp(key_str, "replicates") == 0) {
            parse_int(value, &params->core.num_replicates);
        } else if (strcmp(key_str, "sites") == 0) {
            parse_int(value, &params->core.num_sites);
        } else if (strcmp(key_str, "populations") == 0) {
            parse_int(value, &params->core.num_populations);
        } else if (strcmp(key_str, "effective_size") == 0 || strcmp(key_str, "Ne") == 0) {
            parse_double(value, &params->core.Ne);
        } else if (strcmp(key_str, "sample_distribution") == 0 && value->type == YAML_SEQUENCE_NODE) {
            /* Parse array of sample sizes */
            int idx = 0;
            yaml_node_item_t *item;
            for (item = value->data.sequence.items.start; 
                 item < value->data.sequence.items.top && idx < MAXPOPS; item++) {
                yaml_node_t *item_node = yaml_document_get_node(document, *item);
                parse_int(item_node, &params->core.sample_sizes[idx++]);
            }
        }
    }
    
    return 0;
}

/* Parse evolution section */
static int parse_evolution_section(yaml_document_t *document, yaml_node_t *node, SimulationParams *params) {
    if (node->type != YAML_MAPPING_NODE) return -1;
    
    yaml_node_pair_t *pair;
    for (pair = node->data.mapping.pairs.start; pair < node->data.mapping.pairs.top; pair++) {
        yaml_node_t *key = yaml_document_get_node(document, pair->key);
        yaml_node_t *value = yaml_document_get_node(document, pair->value);
        
        const char *key_str = get_scalar_value(key);
        if (!key_str) continue;
        
        if (strcmp(key_str, "mutation") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse mutation subsection */
            yaml_node_pair_t *mut_pair;
            for (mut_pair = value->data.mapping.pairs.start; mut_pair < value->data.mapping.pairs.top; mut_pair++) {
                yaml_node_t *mut_key = yaml_document_get_node(document, mut_pair->key);
                yaml_node_t *mut_value = yaml_document_get_node(document, mut_pair->value);
                const char *mut_key_str = get_scalar_value(mut_key);
                
                if (mut_key_str && strcmp(mut_key_str, "theta") == 0) {
                    parse_double(mut_value, &params->forces.theta);
                }
            }
        } else if (strcmp(key_str, "recombination") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse recombination subsection */
            yaml_node_pair_t *rec_pair;
            for (rec_pair = value->data.mapping.pairs.start; rec_pair < value->data.mapping.pairs.top; rec_pair++) {
                yaml_node_t *rec_key = yaml_document_get_node(document, rec_pair->key);
                yaml_node_t *rec_value = yaml_document_get_node(document, rec_pair->value);
                const char *rec_key_str = get_scalar_value(rec_key);
                
                if (rec_key_str && strcmp(rec_key_str, "rho") == 0) {
                    parse_double(rec_value, &params->forces.rho);
                } else if (rec_key_str && strcmp(rec_key_str, "gene_conversion") == 0 && 
                          rec_value->type == YAML_MAPPING_NODE) {
                    /* Parse gene conversion sub-subsection */
                    yaml_node_pair_t *gc_pair;
                    for (gc_pair = rec_value->data.mapping.pairs.start; 
                         gc_pair < rec_value->data.mapping.pairs.top; gc_pair++) {
                        yaml_node_t *gc_key = yaml_document_get_node(document, gc_pair->key);
                        yaml_node_t *gc_value = yaml_document_get_node(document, gc_pair->value);
                        const char *gc_key_str = get_scalar_value(gc_key);
                        
                        if (gc_key_str && strcmp(gc_key_str, "ratio") == 0) {
                            parse_double(gc_value, &params->forces.gene_conversion_ratio);
                        } else if (gc_key_str && strcmp(gc_key_str, "tract_length") == 0) {
                            parse_int(gc_value, &params->forces.gc_tract_mean);
                        }
                    }
                }
            }
        } else if (strcmp(key_str, "migration") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse migration subsection */
            yaml_node_pair_t *mig_pair;
            for (mig_pair = value->data.mapping.pairs.start; mig_pair < value->data.mapping.pairs.top; mig_pair++) {
                yaml_node_t *mig_key = yaml_document_get_node(document, mig_pair->key);
                yaml_node_t *mig_value = yaml_document_get_node(document, mig_pair->value);
                const char *mig_key_str = get_scalar_value(mig_key);
                
                if (mig_key_str && strcmp(mig_key_str, "symmetric") == 0) {
                    parse_bool(mig_value, &params->forces.symmetric_migration);
                } else if (mig_key_str && strcmp(mig_key_str, "rate") == 0) {
                    double rate;
                    if (parse_double(mig_value, &rate) == 0 && params->forces.symmetric_migration) {
                        params_set_symmetric_migration(params, rate);
                    }
                }
            }
        }
    }
    
    return 0;
}

/* Parse selection section */
static int parse_selection_section(yaml_document_t *document, yaml_node_t *node, SimulationParams *params) {
    if (node->type != YAML_MAPPING_NODE) return -1;
    
    yaml_node_pair_t *pair;
    for (pair = node->data.mapping.pairs.start; pair < node->data.mapping.pairs.top; pair++) {
        yaml_node_t *key = yaml_document_get_node(document, pair->key);
        yaml_node_t *value = yaml_document_get_node(document, pair->value);
        
        const char *key_str = get_scalar_value(key);
        if (!key_str) continue;
        
        if (strcmp(key_str, "coefficient") == 0 || strcmp(key_str, "alpha") == 0) {
            parse_double(value, &params->selection.alpha);
        } else if (strcmp(key_str, "position") == 0) {
            parse_double(value, &params->selection.sweep_position);
        } else if (strcmp(key_str, "mode") == 0) {
            const char *mode = get_scalar_value(value);
            if (mode) {
                if (strcmp(mode, "deterministic") == 0) {
                    params->selection.sweep_mode = SWEEP_DETERMINISTIC;
                } else if (strcmp(mode, "stochastic") == 0) {
                    params->selection.sweep_mode = SWEEP_STOCHASTIC;
                } else if (strcmp(mode, "neutral") == 0) {
                    params->selection.sweep_mode = SWEEP_NEUTRAL;
                }
            }
        } else if (strcmp(key_str, "time") == 0 || strcmp(key_str, "tau") == 0) {
            double time;
            if (parse_double(value, &time) == 0) {
                params->selection.tau = time * 2.0;  /* Convert to coalescent units */
            }
        } else if (strcmp(key_str, "partial_sweep") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse partial sweep subsection */
            yaml_node_pair_t *ps_pair;
            for (ps_pair = value->data.mapping.pairs.start; ps_pair < value->data.mapping.pairs.top; ps_pair++) {
                yaml_node_t *ps_key = yaml_document_get_node(document, ps_pair->key);
                yaml_node_t *ps_value = yaml_document_get_node(document, ps_pair->value);
                const char *ps_key_str = get_scalar_value(ps_key);
                
                if (ps_key_str && strcmp(ps_key_str, "enabled") == 0) {
                    parse_bool(ps_value, &params->selection.partial_sweep);
                } else if (ps_key_str && strcmp(ps_key_str, "final_frequency") == 0) {
                    parse_double(ps_value, &params->selection.partial_sweep_final_freq);
                }
            }
        }
    }
    
    return 0;
}

/* Parse output section */
static int parse_output_section(yaml_document_t *document, yaml_node_t *node, SimulationParams *params) {
    if (node->type != YAML_MAPPING_NODE) return -1;
    
    yaml_node_pair_t *pair;
    for (pair = node->data.mapping.pairs.start; pair < node->data.mapping.pairs.top; pair++) {
        yaml_node_t *key = yaml_document_get_node(document, pair->key);
        yaml_node_t *value = yaml_document_get_node(document, pair->value);
        
        const char *key_str = get_scalar_value(key);
        if (!key_str) continue;
        
        if (strcmp(key_str, "format") == 0) {
            const char *format = get_scalar_value(value);
            if (format) {
                if (strcmp(format, "ms") == 0) {
                    params->output.format = OUTPUT_MS;
                } else if (strcmp(format, "fasta") == 0) {
                    params->output.format = OUTPUT_FASTA;
                } else if (strcmp(format, "nexus") == 0) {
                    params->output.format = OUTPUT_NEXUS;
                }
            }
        } else if (strcmp(key_str, "tree_sequence") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse tree sequence subsection */
            yaml_node_pair_t *ts_pair;
            for (ts_pair = value->data.mapping.pairs.start; ts_pair < value->data.mapping.pairs.top; ts_pair++) {
                yaml_node_t *ts_key = yaml_document_get_node(document, ts_pair->key);
                yaml_node_t *ts_value = yaml_document_get_node(document, ts_pair->value);
                const char *ts_key_str = get_scalar_value(ts_key);
                
                if (ts_key_str && strcmp(ts_key_str, "enabled") == 0) {
                    parse_bool(ts_value, &params->output.tree_sequence_output);
                } else if (ts_key_str && strcmp(ts_key_str, "filename") == 0) {
                    const char *filename = get_scalar_value(ts_value);
                    if (filename) {
                        strncpy(params->output.output_filename, filename, PATH_MAX - 1);
                        params->output.output_filename[PATH_MAX - 1] = '\0';
                    }
                } else if (ts_key_str && strcmp(ts_key_str, "minimal") == 0) {
                    parse_bool(ts_value, &params->output.minimal_tree_seq);
                }
            }
        }
    }
    
    return 0;
}

/* Parse a YAML document into simulation parameters */
int yaml_parse_document(yaml_parser_t *parser, SimulationParams *params) {
    yaml_document_t document;
    
    if (!yaml_parser_load(parser, &document)) {
        fprintf(stderr, "Error loading YAML document\n");
        return -1;
    }
    
    /* Get root node */
    yaml_node_t *root = yaml_document_get_root_node(&document);
    if (!root || root->type != YAML_MAPPING_NODE) {
        fprintf(stderr, "Error: YAML root must be a mapping\n");
        yaml_document_delete(&document);
        return -1;
    }
    
    /* Parse top-level sections */
    yaml_node_pair_t *pair;
    for (pair = root->data.mapping.pairs.start; pair < root->data.mapping.pairs.top; pair++) {
        yaml_node_t *key = yaml_document_get_node(&document, pair->key);
        yaml_node_t *value = yaml_document_get_node(&document, pair->value);
        
        const char *section = get_scalar_value(key);
        if (!section) continue;
        
        if (strcmp(section, "simulation") == 0) {
            parse_simulation_section(&document, value, params);
        } else if (strcmp(section, "evolution") == 0) {
            parse_evolution_section(&document, value, params);
        } else if (strcmp(section, "selection") == 0) {
            parse_selection_section(&document, value, params);
        } else if (strcmp(section, "output") == 0) {
            parse_output_section(&document, value, params);
        } else if (strcmp(section, "random") == 0 && value->type == YAML_MAPPING_NODE) {
            /* Parse random section */
            yaml_node_pair_t *rand_pair;
            for (rand_pair = value->data.mapping.pairs.start; 
                 rand_pair < value->data.mapping.pairs.top; rand_pair++) {
                yaml_node_t *rand_key = yaml_document_get_node(&document, rand_pair->key);
                yaml_node_t *rand_value = yaml_document_get_node(&document, rand_pair->value);
                const char *rand_key_str = get_scalar_value(rand_key);
                
                if (rand_key_str && strcmp(rand_key_str, "seed1") == 0) {
                    int seed;
                    if (parse_int(rand_value, &seed) == 0) {
                        params->random.seed1 = (unsigned long)seed;
                    }
                } else if (rand_key_str && strcmp(rand_key_str, "seed2") == 0) {
                    int seed;
                    if (parse_int(rand_value, &seed) == 0) {
                        params->random.seed2 = (unsigned long)seed;
                    }
                }
            }
        }
        /* TODO: Parse demographics section (will integrate with Demes) */
    }
    
    yaml_document_delete(&document);
    return 0;
}

/* Load simulation parameters from a YAML file */
int yaml_load_params(SimulationParams *params, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open YAML file '%s': %s\n", filename, strerror(errno));
        return -1;
    }
    
    yaml_parser_t parser;
    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Error: Failed to initialize YAML parser\n");
        fclose(file);
        return -1;
    }
    
    yaml_parser_set_input_file(&parser, file);
    
    int result = yaml_parse_document(&parser, params);
    
    yaml_parser_delete(&parser);
    fclose(file);
    
    if (result == 0) {
        /* Validate parameters */
        char error_msg[256];
        if (params_validate(params, error_msg, sizeof(error_msg)) != 0) {
            fprintf(stderr, "Parameter validation error: %s\n", error_msg);
            return -1;
        }
    }
    
    return result;
}

/* Save simulation parameters to a YAML file */
int yaml_save_params(const SimulationParams *params, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Cannot create YAML file '%s': %s\n", filename, strerror(errno));
        return -1;
    }
    
    yaml_emitter_t emitter;
    yaml_event_t event;
    
    if (!yaml_emitter_initialize(&emitter)) {
        fprintf(stderr, "Error: Failed to initialize YAML emitter\n");
        fclose(file);
        return -1;
    }
    
    yaml_emitter_set_output_file(&emitter, file);
    
    /* Start document */
    yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
    yaml_emitter_emit(&emitter, &event);
    
    yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 0);
    yaml_emitter_emit(&emitter, &event);
    
    /* Root mapping */
    yaml_mapping_start_event_initialize(&event, NULL, NULL, 1, YAML_BLOCK_MAPPING_STYLE);
    yaml_emitter_emit(&emitter, &event);
    
    /* TODO: Implement full YAML output using emitter API */
    /* For now, use simple fprintf approach */
    yaml_mapping_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    
    yaml_document_end_event_initialize(&event, 0);
    yaml_emitter_emit(&emitter, &event);
    
    yaml_stream_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    
    yaml_emitter_delete(&emitter);
    fclose(file);
    
    /* For now, fall back to simple text output */
    FILE *out = fopen(filename, "w");
    if (!out) return -1;
    
    fprintf(out, "# discoal simulation parameters\n\n");
    fprintf(out, "simulation:\n");
    fprintf(out, "  samples: %d\n", params->core.total_samples);
    fprintf(out, "  replicates: %d\n", params->core.num_replicates);
    fprintf(out, "  sites: %d\n", params->core.num_sites);
    if (params->core.num_populations > 1) {
        fprintf(out, "  populations: %d\n", params->core.num_populations);
    }
    
    fprintf(out, "\nevolution:\n");
    if (params->forces.theta > 0) {
        fprintf(out, "  mutation:\n");
        fprintf(out, "    theta: %.4f\n", params->forces.theta);
    }
    if (params->forces.rho > 0) {
        fprintf(out, "  recombination:\n");
        fprintf(out, "    rho: %.4f\n", params->forces.rho);
    }
    
    fclose(out);
    return 0;
}
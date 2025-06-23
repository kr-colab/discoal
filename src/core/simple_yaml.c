/**
 * simple_yaml.c - Simple YAML-like configuration parser for discoal
 * 
 * This implements a subset of YAML that's sufficient for our needs
 * without requiring an external dependency.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "params.h"
#include "version.h"

#define MAX_LINE_LENGTH 1024
#define MAX_KEY_LENGTH 256

/* Trim whitespace from both ends of a string */
static char* trim(char *str) {
    char *end;
    
    /* Trim leading whitespace */
    while (isspace((unsigned char)*str)) str++;
    
    if (*str == 0) return str;
    
    /* Trim trailing whitespace */
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    
    return str;
}

/* Parse a key: value line */
static int parse_key_value(char *line, char *key, char *value) {
    char *colon = strchr(line, ':');
    if (!colon) return -1;
    
    /* Split at colon */
    *colon = '\0';
    
    /* Copy and trim key */
    strncpy(key, trim(line), MAX_KEY_LENGTH - 1);
    key[MAX_KEY_LENGTH - 1] = '\0';
    
    /* Copy and trim value */
    strncpy(value, trim(colon + 1), MAX_LINE_LENGTH - 1);
    value[MAX_LINE_LENGTH - 1] = '\0';
    
    return 0;
}

/* Get the indentation level of a line (number of leading spaces) */
static int get_indent_level(const char *line) {
    int count = 0;
    while (*line == ' ') {
        count++;
        line++;
    }
    return count;
}

/* Parse a boolean value */
static int parse_bool(const char *value) {
    if (strcasecmp(value, "true") == 0 || 
        strcasecmp(value, "yes") == 0 ||
        strcasecmp(value, "on") == 0 ||
        strcmp(value, "1") == 0) {
        return 1;
    }
    return 0;
}

/* State for parsing nested structures */
typedef struct {
    FILE *file;
    int line_number;
    char current_line[MAX_LINE_LENGTH];
    int current_indent;
    SimulationParams *params;
} ParseState;

/* Read next non-empty, non-comment line */
static int read_next_line(ParseState *state) {
    while (fgets(state->current_line, MAX_LINE_LENGTH, state->file)) {
        state->line_number++;
        
        /* Remove newline */
        state->current_line[strcspn(state->current_line, "\n")] = '\0';
        
        /* Skip empty lines and comments */
        char *trimmed = trim(state->current_line);
        if (*trimmed == '\0' || *trimmed == '#') continue;
        
        state->current_indent = get_indent_level(state->current_line);
        return 1;
    }
    return 0;
}

/* Parse simulation core section */
static int parse_simulation_section(ParseState *state, int base_indent) {
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    while (read_next_line(state)) {
        if (state->current_indent <= base_indent) {
            /* Back up one line for parent parser */
            fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
            state->line_number--;
            break;
        }
        
        if (parse_key_value(state->current_line, key, value) < 0) continue;
        
        if (strcmp(key, "samples") == 0) {
            state->params->core.total_samples = atoi(value);
            /* If single population, set sample size */
            if (state->params->core.num_populations == 1) {
                state->params->core.sample_sizes[0] = state->params->core.total_samples;
            }
        } else if (strcmp(key, "replicates") == 0) {
            state->params->core.num_replicates = atoi(value);
        } else if (strcmp(key, "sites") == 0) {
            state->params->core.num_sites = atoi(value);
        } else if (strcmp(key, "populations") == 0) {
            state->params->core.num_populations = atoi(value);
        } else if (strcmp(key, "effective_size") == 0 || strcmp(key, "Ne") == 0) {
            state->params->core.Ne = atof(value);
        } else if (strcmp(key, "sample_distribution") == 0) {
            /* Parse array like [50, 50] */
            char *p = strchr(value, '[');
            if (p) {
                p++; /* Skip [ */
                int i = 0;
                char *token = strtok(p, ",]");
                while (token && i < MAXPOPS) {
                    state->params->core.sample_sizes[i++] = atoi(trim(token));
                    token = strtok(NULL, ",]");
                }
            }
        }
    }
    
    return 0;
}

/* Parse evolution section */
static int parse_evolution_section(ParseState *state, int base_indent) {
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    while (read_next_line(state)) {
        if (state->current_indent <= base_indent) {
            fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
            state->line_number--;
            break;
        }
        
        if (parse_key_value(state->current_line, key, value) < 0) continue;
        
        /* Check for subsections */
        if (strcmp(key, "mutation") == 0) {
            /* Handle mutation subsection */
            while (read_next_line(state) && state->current_indent > base_indent + 2) {
                if (parse_key_value(state->current_line, key, value) < 0) continue;
                if (strcmp(key, "theta") == 0) {
                    state->params->forces.theta = atof(value);
                }
            }
            /* Back up for parent */
            if (state->current_indent <= base_indent + 2) {
                fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                state->line_number--;
            }
        } else if (strcmp(key, "recombination") == 0) {
            /* Handle recombination subsection */
            while (read_next_line(state) && state->current_indent > base_indent + 2) {
                if (parse_key_value(state->current_line, key, value) < 0) continue;
                if (strcmp(key, "rho") == 0) {
                    state->params->forces.rho = atof(value);
                } else if (strcmp(key, "gene_conversion") == 0) {
                    /* Handle gene conversion sub-subsection */
                    while (read_next_line(state) && state->current_indent > base_indent + 4) {
                        if (parse_key_value(state->current_line, key, value) < 0) continue;
                        if (strcmp(key, "enabled") == 0) {
                            /* If enabled, we'll set ratio when we see it */
                        } else if (strcmp(key, "ratio") == 0) {
                            state->params->forces.gene_conversion_ratio = atof(value);
                        } else if (strcmp(key, "tract_length") == 0) {
                            state->params->forces.gc_tract_mean = atoi(value);
                        }
                    }
                    if (state->current_indent <= base_indent + 4) {
                        fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                        state->line_number--;
                    }
                }
            }
            if (state->current_indent <= base_indent + 2) {
                fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                state->line_number--;
            }
        } else if (strcmp(key, "migration") == 0) {
            /* Handle migration subsection */
            while (read_next_line(state) && state->current_indent > base_indent + 2) {
                if (parse_key_value(state->current_line, key, value) < 0) continue;
                if (strcmp(key, "symmetric") == 0) {
                    state->params->forces.symmetric_migration = parse_bool(value);
                } else if (strcmp(key, "rate") == 0) {
                    state->params->forces.migration_rate = atof(value);
                    if (state->params->forces.symmetric_migration) {
                        params_set_symmetric_migration(state->params, state->params->forces.migration_rate);
                    }
                }
            }
            if (state->current_indent <= base_indent + 2) {
                fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                state->line_number--;
            }
        }
    }
    
    return 0;
}

/* Parse selection section */
static int parse_selection_section(ParseState *state, int base_indent) {
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    while (read_next_line(state)) {
        if (state->current_indent <= base_indent) {
            fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
            state->line_number--;
            break;
        }
        
        if (parse_key_value(state->current_line, key, value) < 0) continue;
        
        if (strcmp(key, "enabled") == 0) {
            /* Just for clarity in config */
        } else if (strcmp(key, "coefficient") == 0 || strcmp(key, "alpha") == 0) {
            state->params->selection.alpha = atof(value);
        } else if (strcmp(key, "position") == 0) {
            state->params->selection.sweep_position = atof(value);
        } else if (strcmp(key, "mode") == 0) {
            if (strcmp(value, "deterministic") == 0) {
                state->params->selection.sweep_mode = SWEEP_DETERMINISTIC;
            } else if (strcmp(value, "stochastic") == 0) {
                state->params->selection.sweep_mode = SWEEP_STOCHASTIC;
            } else if (strcmp(value, "neutral") == 0) {
                state->params->selection.sweep_mode = SWEEP_NEUTRAL;
            }
        } else if (strcmp(key, "tau") == 0 || strcmp(key, "time") == 0) {
            state->params->selection.tau = atof(value) * 2.0;  /* Convert to coalescent units */
        } else if (strcmp(key, "partial_sweep") == 0) {
            /* Handle partial sweep subsection */
            while (read_next_line(state) && state->current_indent > base_indent + 2) {
                if (parse_key_value(state->current_line, key, value) < 0) continue;
                if (strcmp(key, "enabled") == 0) {
                    state->params->selection.partial_sweep = parse_bool(value);
                } else if (strcmp(key, "final_frequency") == 0) {
                    state->params->selection.partial_sweep_final_freq = atof(value);
                }
            }
            if (state->current_indent <= base_indent + 2) {
                fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                state->line_number--;
            }
        }
    }
    
    return 0;
}

/* Parse output section */
static int parse_output_section(ParseState *state, int base_indent) {
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    while (read_next_line(state)) {
        if (state->current_indent <= base_indent) {
            fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
            state->line_number--;
            break;
        }
        
        if (parse_key_value(state->current_line, key, value) < 0) continue;
        
        if (strcmp(key, "format") == 0) {
            if (strcmp(value, "ms") == 0) {
                state->params->output.format = OUTPUT_MS;
            } else if (strcmp(value, "fasta") == 0) {
                state->params->output.format = OUTPUT_FASTA;
            } else if (strcmp(value, "nexus") == 0) {
                state->params->output.format = OUTPUT_NEXUS;
            }
        } else if (strcmp(key, "tree_sequence") == 0) {
            /* Handle tree sequence subsection */
            while (read_next_line(state) && state->current_indent > base_indent + 2) {
                if (parse_key_value(state->current_line, key, value) < 0) continue;
                if (strcmp(key, "enabled") == 0) {
                    state->params->output.tree_sequence_output = parse_bool(value);
                } else if (strcmp(key, "filename") == 0) {
                    /* Remove quotes if present */
                    if (value[0] == '"' || value[0] == '\'') {
                        int len = strlen(value);
                        if (len > 2 && value[len-1] == value[0]) {
                            value[len-1] = '\0';
                            memmove(value, value + 1, len - 1);
                        }
                    }
                    strncpy(state->params->output.output_filename, value, PATH_MAX - 1);
                } else if (strcmp(key, "minimal") == 0) {
                    state->params->output.minimal_tree_seq = parse_bool(value);
                }
            }
            if (state->current_indent <= base_indent + 2) {
                fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                state->line_number--;
            }
        }
    }
    
    return 0;
}

/* Parse demographic events array */
static int parse_demographics_section(ParseState *state, int base_indent) {
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    while (read_next_line(state)) {
        if (state->current_indent <= base_indent) {
            fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
            state->line_number--;
            break;
        }
        
        if (parse_key_value(state->current_line, key, value) < 0) continue;
        
        if (strcmp(key, "events") == 0) {
            /* Parse events array */
            while (read_next_line(state)) {
                /* Check if we've left the events section */
                if (state->current_indent <= base_indent) {
                    fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                    state->line_number--;
                    break;
                }
                
                /* Look for array item marker (- ) */
                char *trimmed = trim(state->current_line);
                if (trimmed[0] == '-') {
                    /* New event */
                    DemographicEvent event;
                    memset(&event, 0, sizeof(event));
                    
                    /* Read event properties */
                    while (read_next_line(state) && state->current_indent > base_indent + 2) {
                        if (parse_key_value(state->current_line, key, value) < 0) continue;
                        
                        if (strcmp(key, "time") == 0) {
                            event.time = atof(value) * 2.0;  /* Convert to coalescent units */
                        } else if (strcmp(key, "type") == 0) {
                            if (strcmp(value, "size_change") == 0) {
                                event.type = EVENT_SIZE_CHANGE;
                            } else if (strcmp(value, "split") == 0) {
                                event.type = EVENT_SPLIT;
                            } else if (strcmp(value, "join") == 0) {
                                event.type = EVENT_JOIN;
                            } else if (strcmp(value, "admix") == 0) {
                                event.type = EVENT_ADMIX;
                            }
                        } else if (strcmp(key, "population") == 0) {
                            event.params.size_change.pop = atoi(value);
                        } else if (strcmp(key, "size") == 0) {
                            event.params.size_change.size = atof(value);
                        } else if (strcmp(key, "from") == 0) {
                            if (event.type == EVENT_SPLIT) {
                                event.params.split.from = atoi(value);
                            } else if (event.type == EVENT_ADMIX) {
                                event.params.admix.from = atoi(value);
                            }
                        } else if (strcmp(key, "to") == 0) {
                            if (event.type == EVENT_ADMIX) {
                                event.params.admix.to = atoi(value);
                            } else if (event.type == EVENT_SPLIT) {
                                /* Parse array like [1, 2] */
                                char *p = strchr(value, '[');
                                if (p) {
                                    p++;
                                    char *token = strtok(p, ",]");
                                    if (token) event.params.split.to1 = atoi(trim(token));
                                    token = strtok(NULL, ",]");
                                    if (token) event.params.split.to2 = atoi(trim(token));
                                }
                            }
                        } else if (strcmp(key, "proportions") == 0) {
                            /* Parse array like [0.3, 0.7] */
                            char *p = strchr(value, '[');
                            if (p) {
                                p++;
                                char *token = strtok(p, ",]");
                                if (token) event.params.split.prop = atof(trim(token));
                            }
                        } else if (strcmp(key, "populations") == 0 && event.type == EVENT_JOIN) {
                            /* Parse array like [1, 2] for join event */
                            char *p = strchr(value, '[');
                            if (p) {
                                p++;
                                char *token = strtok(p, ",]");
                                if (token) event.params.join.pop1 = atoi(trim(token));
                                token = strtok(NULL, ",]");
                                if (token) event.params.join.pop2 = atoi(trim(token));
                            }
                        } else if (strcmp(key, "destination") == 0) {
                            event.params.join.dest = atoi(value);
                        }
                    }
                    
                    /* Add the event */
                    params_add_demographic_event(state->params, &event);
                    
                    /* Continue to next iteration to look for more events */
                } else {
                    /* Not an array item, back up */
                    fseek(state->file, -strlen(state->current_line) - 1, SEEK_CUR);
                    state->line_number--;
                    break;
                }
            }
        }
    }
    
    return 0;
}

/* Main YAML loading function */
int params_load_from_yaml(SimulationParams *params, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open YAML file '%s': %s\n", filename, strerror(errno));
        return -1;
    }
    
    ParseState state = {
        .file = file,
        .line_number = 0,
        .params = params
    };
    
    char key[MAX_KEY_LENGTH], value[MAX_LINE_LENGTH];
    
    /* Parse top-level sections */
    while (read_next_line(&state)) {
        if (state.current_indent > 0) continue;  /* Skip indented lines at top level */
        
        if (parse_key_value(state.current_line, key, value) < 0) continue;
        
        /* Handle top-level sections */
        if (strcmp(key, "simulation") == 0) {
            parse_simulation_section(&state, 0);
        } else if (strcmp(key, "evolution") == 0) {
            parse_evolution_section(&state, 0);
        } else if (strcmp(key, "selection") == 0) {
            parse_selection_section(&state, 0);
        } else if (strcmp(key, "demographics") == 0) {
            parse_demographics_section(&state, 0);
        } else if (strcmp(key, "output") == 0) {
            parse_output_section(&state, 0);
        } else if (strcmp(key, "random") == 0) {
            /* Parse random section */
            while (read_next_line(&state) && state.current_indent > 0) {
                if (parse_key_value(state.current_line, key, value) < 0) continue;
                if (strcmp(key, "seed1") == 0) {
                    params->random.seed1 = atol(value);
                } else if (strcmp(key, "seed2") == 0) {
                    params->random.seed2 = atol(value);
                } else if (strcmp(key, "engine") == 0) {
                    if (strstr(value, "xoshiro") != NULL) {
                        params->random.use_xoshiro = true;
                    }
                }
            }
            if (state.current_indent <= 0) {
                fseek(state.file, -strlen(state.current_line) - 1, SEEK_CUR);
                state.line_number--;
            }
        }
    }
    
    fclose(file);
    
    /* Validate the loaded parameters */
    char error_msg[256];
    if (params_validate(params, error_msg, sizeof(error_msg)) != 0) {
        fprintf(stderr, "YAML validation error: %s\n", error_msg);
        return -1;
    }
    
    return 0;
}

/* Save parameters to YAML file */
int params_save_to_yaml(const SimulationParams *params, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Cannot create YAML file '%s': %s\n", filename, strerror(errno));
        return -1;
    }
    
    /* Write header */
    fprintf(file, "# discoal simulation parameters\n");
    fprintf(file, "# Generated by discoal version " DISCOAL_VERSION "\n\n");
    
    /* Simulation section */
    fprintf(file, "simulation:\n");
    fprintf(file, "  samples: %d\n", params->core.total_samples);
    fprintf(file, "  replicates: %d\n", params->core.num_replicates);
    fprintf(file, "  sites: %d\n", params->core.num_sites);
    if (params->core.num_populations > 1) {
        fprintf(file, "  populations: %d\n", params->core.num_populations);
        fprintf(file, "  sample_distribution: [");
        for (int i = 0; i < params->core.num_populations; i++) {
            if (i > 0) fprintf(file, ", ");
            fprintf(file, "%d", params->core.sample_sizes[i]);
        }
        fprintf(file, "]\n");
    }
    if (params->core.Ne != 1.0) {
        fprintf(file, "  effective_size: %.0f\n", params->core.Ne);
    }
    
    /* Evolution section */
    fprintf(file, "\nevolution:\n");
    
    if (params->forces.theta > 0) {
        fprintf(file, "  mutation:\n");
        fprintf(file, "    theta: %.4f\n", params->forces.theta);
    }
    
    if (params->forces.rho > 0) {
        fprintf(file, "  recombination:\n");
        fprintf(file, "    rho: %.4f\n", params->forces.rho);
        if (params->forces.gene_conversion_ratio > 0) {
            fprintf(file, "    gene_conversion:\n");
            fprintf(file, "      enabled: true\n");
            fprintf(file, "      ratio: %.4f\n", params->forces.gene_conversion_ratio);
            fprintf(file, "      tract_length: %d\n", params->forces.gc_tract_mean);
        }
    }
    
    if (params->forces.symmetric_migration && params->forces.migration_rate > 0) {
        fprintf(file, "  migration:\n");
        fprintf(file, "    symmetric: true\n");
        fprintf(file, "    rate: %.4f\n", params->forces.migration_rate);
    }
    
    /* Selection section */
    if (params->selection.alpha != 0) {
        fprintf(file, "\nselection:\n");
        fprintf(file, "  enabled: true\n");
        fprintf(file, "  coefficient: %.4f\n", params->selection.alpha);
        fprintf(file, "  position: %.4f\n", params->selection.sweep_position);
        
        const char *mode_str = "neutral";
        switch (params->selection.sweep_mode) {
            case SWEEP_DETERMINISTIC: mode_str = "deterministic"; break;
            case SWEEP_STOCHASTIC: mode_str = "stochastic"; break;
            case SWEEP_NEUTRAL: mode_str = "neutral"; break;
        }
        fprintf(file, "  mode: %s\n", mode_str);
        
        if (params->selection.tau > 0) {
            fprintf(file, "  time: %.4f\n", params->selection.tau / 2.0);  /* Convert from coalescent units */
        }
        
        if (params->selection.partial_sweep) {
            fprintf(file, "  partial_sweep:\n");
            fprintf(file, "    enabled: true\n");
            fprintf(file, "    final_frequency: %.4f\n", params->selection.partial_sweep_final_freq);
        }
    }
    
    /* Demographics section */
    if (params->demographics.num_events > 0) {
        fprintf(file, "\ndemographics:\n");
        fprintf(file, "  events:\n");
        
        for (int i = 0; i < params->demographics.num_events; i++) {
            const DemographicEvent *event = &params->demographics.events[i];
            fprintf(file, "    - time: %.6f\n", event->time / 2.0);  /* Convert from coalescent units */
            
            switch (event->type) {
                case EVENT_SIZE_CHANGE:
                    fprintf(file, "      type: size_change\n");
                    fprintf(file, "      population: %d\n", event->params.size_change.pop);
                    fprintf(file, "      size: %.4f\n", event->params.size_change.size);
                    break;
                case EVENT_SPLIT:
                    fprintf(file, "      type: split\n");
                    fprintf(file, "      from: %d\n", event->params.split.from);
                    fprintf(file, "      to: [%d, %d]\n", event->params.split.to1, event->params.split.to2);
                    fprintf(file, "      proportions: [%.4f, %.4f]\n", 
                            event->params.split.prop, 1.0 - event->params.split.prop);
                    break;
                case EVENT_JOIN:
                    fprintf(file, "      type: join\n");
                    fprintf(file, "      populations: [%d, %d]\n", 
                            event->params.join.pop1, event->params.join.pop2);
                    fprintf(file, "      destination: %d\n", event->params.join.dest);
                    break;
                case EVENT_ADMIX:
                    fprintf(file, "      type: admix\n");
                    fprintf(file, "      to: %d\n", event->params.admix.to);
                    fprintf(file, "      from: %d\n", event->params.admix.from);
                    fprintf(file, "      proportion: %.4f\n", event->params.admix.prop);
                    break;
                default:
                    break;
            }
        }
    }
    
    /* Output section */
    fprintf(file, "\noutput:\n");
    const char *format_names[] = {"ms", "fasta", "nexus"};
    fprintf(file, "  format: %s\n", format_names[params->output.format]);
    
    if (params->output.tree_sequence_output) {
        fprintf(file, "  tree_sequence:\n");
        fprintf(file, "    enabled: true\n");
        fprintf(file, "    filename: \"%s\"\n", params->output.output_filename);
        fprintf(file, "    minimal: %s\n", params->output.minimal_tree_seq ? "true" : "false");
    }
    
    /* Random section */
    if (params->random.seed1 != 0 || params->random.seed2 != 0) {
        fprintf(file, "\nrandom:\n");
        if (params->random.seed1 != 0) {
            fprintf(file, "  seed1: %lu\n", params->random.seed1);
        }
        if (params->random.seed2 != 0) {
            fprintf(file, "  seed2: %lu\n", params->random.seed2);
        }
        if (params->random.use_xoshiro) {
            fprintf(file, "  engine: xoshiro256++\n");
        }
    }
    
    fclose(file);
    return 0;
}
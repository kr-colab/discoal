#include "configInterface.h"
#include "demesInterface.h"
#include "discoalFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>

extern int eventNumber;
extern int eventsCapacity;
extern struct event *events;
extern double *currentSize;

// Initialize configuration with default values
void initializeDefaultConfig(SimulationConfig *config) {
    memset(config, 0, sizeof(SimulationConfig));
    
    // Set default values that match discoal defaults
    config->sample_size = 0;
    config->num_replicates = 0;
    config->num_sites = 0;
    config->seed1 = 0;
    config->seed2 = 0;
    config->has_seed = 0;
    
    config->npops = 1;
    config->effective_popn_size = 1000000;  // EFFECTIVE_POPN_SIZE default
    config->has_populations = 0;
    
    config->theta = 0.0;
    config->rho = 0.0;
    config->gamma = 0.0;
    config->gc_mean = 0;
    config->gamma_co_ratio = 0.0;
    config->gamma_co_ratio_mode = 0;
    config->has_genetics = 0;
    
    config->alpha = 0.0;
    config->sweep_site = 0.5;
    config->sweep_mode = 's';  // default to stochastic
    config->tau = 0.0;
    config->f0 = 0.0;
    config->ua = 0.0;
    config->partial_sweep_final_freq = 0.0;
    config->recur_sweep_mode = 0;
    config->recur_sweep_rate = 0.0;
    config->partial_sweep_mode = 0;
    config->soft_sweep_mode = 0;
    config->has_selection = 0;
    
    config->output_style = 'h';
    config->finite_output_flag = 0;
    config->hide_partial_snp = 0;
    config->tskit_output_mode = 0;
    config->tskit_output_filename[0] = '\0';
    config->minimal_tree_seq = 1;
    config->has_output = 0;
    
    config->mig_flag = 0;
    for (int i = 0; i < MAXPOPS; i++) {
        for (int j = 0; j < MAXPOPS; j++) {
            config->mig_mat_const[i][j] = 0.0;
        }
    }
    
    config->demes_file[0] = '\0';
    config->use_demes = 0;
    config->num_explicit_events = 0;
    
    config->anc_sample_flag = 0;
    config->anc_sample_size = 0;
}

// Helper function to parse YAML scalar values
static int parseScalarValue(yaml_event_t *event, char *buffer, size_t buffer_size) {
    if (event->type != YAML_SCALAR_EVENT) return 0;
    
    size_t length = event->data.scalar.length;
    if (length >= buffer_size) length = buffer_size - 1;
    
    strncpy(buffer, (char*)event->data.scalar.value, length);
    buffer[length] = '\0';
    
    return 1;
}

// Parse simulation section
static int parseSimulationSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    
    // fprintf(stderr, "DEBUG: Entering parseSimulationSection\n");
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in simulation section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (!parseScalarValue(&event, value, sizeof(value))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Process key-value pairs
                // fprintf(stderr, "DEBUG: Processing key='%s', value='%s'\n", key, value);
                if (strcmp(key, "sample_size") == 0) {
                    config->sample_size = atoi(value);
                } else if (strcmp(key, "num_replicates") == 0) {
                    config->num_replicates = atoi(value);
                } else if (strcmp(key, "num_sites") == 0) {
                    config->num_sites = atoi(value);
                } else if (strcmp(key, "effective_popn_size") == 0) {
                    config->effective_popn_size = atof(value);
                } else if (strcmp(key, "seed") == 0) {
                    // Handle seed as array - for now just parse first value
                    // fprintf(stderr, "DEBUG: Parsing scalar seed: %s\n", value);
                    config->seed1 = atoi(value);
                    config->has_seed = 1;
                } else {
                    // Check if this might be seed array elements (numeric keys)
                    char *endptr;
                    long first_seed = strtol(key, &endptr, 10);
                    if (*endptr == '\0') {  // key is a valid number
                        long second_seed = strtol(value, &endptr, 10);
                        if (*endptr == '\0') {  // value is also a valid number
                            // fprintf(stderr, "DEBUG: Detected seed array elements: %ld, %ld\n", first_seed, second_seed);
                            config->seed1 = (int)first_seed;
                            config->seed2 = (int)second_seed;
                            config->has_seed = 1;
                        }
                    }
                }
                break;
                
            case YAML_SEQUENCE_START_EVENT:
                // Handle sequences (like seed array)
                if (strcmp(key, "seed") == 0) {
                    // fprintf(stderr, "DEBUG: Parsing seed sequence\n");
                    yaml_event_delete(&event);
                    
                    // Parse first seed value
                    if (yaml_parser_parse(parser, &event) && event.type == YAML_SCALAR_EVENT) {
                        parseScalarValue(&event, value, sizeof(value));
                        config->seed1 = atoi(value);
                        // fprintf(stderr, "DEBUG: First seed: %s -> %d\n", value, config->seed1);
                        yaml_event_delete(&event);
                        
                        // Parse second seed value
                        if (yaml_parser_parse(parser, &event) && event.type == YAML_SCALAR_EVENT) {
                            parseScalarValue(&event, value, sizeof(value));
                            config->seed2 = atoi(value);
                            config->has_seed = 1;
                            // fprintf(stderr, "DEBUG: Second seed: %s -> %d, has_seed set to 1\n", value, config->seed2);
                        }
                    }
                    
                    // Skip to end of sequence
                    while (yaml_parser_parse(parser, &event) && event.type != YAML_SEQUENCE_END_EVENT) {
                        yaml_event_delete(&event);
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Parse genetics section
static int parseGeneticsSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    
    config->has_genetics = 1;
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in genetics section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&event, value, sizeof(value));
                    
                    // Process key-value pairs
                    if (strcmp(key, "mutation_rate") == 0) {
                        config->theta = atof(value);
                    } else if (strcmp(key, "recombination_rate") == 0) {
                        config->rho = atof(value);
                    }
                } else if (event.type == YAML_MAPPING_START_EVENT) {
                    // Handle nested mappings like gene_conversion
                    if (strcmp(key, "gene_conversion") == 0) {
                        // Parse gene conversion sub-section
                        yaml_event_t gc_event;
                        char gc_key[256], gc_value[256];
                        
                        while (yaml_parser_parse(parser, &gc_event)) {
                            if (gc_event.type == YAML_MAPPING_END_EVENT) {
                                yaml_event_delete(&gc_event);
                                break;
                            }
                            
                            if (gc_event.type == YAML_SCALAR_EVENT) {
                                parseScalarValue(&gc_event, gc_key, sizeof(gc_key));
                                yaml_event_delete(&gc_event);
                                
                                if (yaml_parser_parse(parser, &gc_event) && gc_event.type == YAML_SCALAR_EVENT) {
                                    parseScalarValue(&gc_event, gc_value, sizeof(gc_value));
                                    
                                    if (strcmp(gc_key, "rate") == 0) {
                                        config->gamma = atof(gc_value);
                                    } else if (strcmp(gc_key, "tract_length") == 0) {
                                        config->gc_mean = atoi(gc_value);
                                    } else if (strcmp(gc_key, "crossover_ratio") == 0) {
                                        config->gamma_co_ratio = atof(gc_value);
                                        config->gamma_co_ratio_mode = 1;
                                    }
                                }
                            }
                            yaml_event_delete(&gc_event);
                        }
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Parse populations section
static int parsePopulationsSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    
    config->has_populations = 1;
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in populations section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&event, value, sizeof(value));
                    
                    if (strcmp(key, "count") == 0) {
                        config->npops = atoi(value);
                    }
                } else if (event.type == YAML_SEQUENCE_START_EVENT) {
                    if (strcmp(key, "sample_sizes") == 0) {
                        int pop_index = 0;
                        yaml_event_t seq_event;
                        
                        while (yaml_parser_parse(parser, &seq_event) && seq_event.type != YAML_SEQUENCE_END_EVENT) {
                            if (seq_event.type == YAML_SCALAR_EVENT && pop_index < MAXPOPS) {
                                parseScalarValue(&seq_event, value, sizeof(value));
                                config->sample_sizes[pop_index] = atoi(value);
                                pop_index++;
                            }
                            yaml_event_delete(&seq_event);
                        }
                        yaml_event_delete(&seq_event);
                        
                        // Set npops based on number of sample sizes
                        config->npops = pop_index;
                        
                        // Calculate total sample size
                        config->sample_size = 0;
                        for (int i = 0; i < pop_index; i++) {
                            config->sample_size += config->sample_sizes[i];
                        }
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Parse events section
static int parseEventsSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    int in_sequence = 0;
    
    // Check if we're already in a sequence (called from main parser after YAML_SEQUENCE_START_EVENT)
    if (!yaml_parser_parse(parser, &event)) {
        fprintf(stderr, "Error parsing YAML in events section\n");
        return -1;
    }
    
    // If we got a mapping start, we're directly in the event sequence
    if (event.type == YAML_MAPPING_START_EVENT) {
        in_sequence = 1;
        yaml_event_delete(&event);
        
        // Parse events directly
        yaml_event_t ev_event;
        do {
            // Parse individual event
            char item_key[256], item_value[256];
            char event_type[256] = "";
            double event_time = 0.0;
            int pop_id = 0, pop_id2 = 0;
            double size_or_rate = 0.0;
            
            while (yaml_parser_parse(parser, &ev_event)) {
                if (ev_event.type == YAML_MAPPING_END_EVENT) {
                    yaml_event_delete(&ev_event);
                    break;
                }
                
                if (ev_event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&ev_event, item_key, sizeof(item_key));
                    yaml_event_delete(&ev_event);
                    
                    if (yaml_parser_parse(parser, &ev_event) && ev_event.type == YAML_SCALAR_EVENT) {
                        parseScalarValue(&ev_event, item_value, sizeof(item_value));
                        
                        if (strcmp(item_key, "type") == 0) {
                            strncpy(event_type, item_value, sizeof(event_type) - 1);
                        } else if (strcmp(item_key, "time") == 0) {
                            event_time = atof(item_value);
                        } else if (strcmp(item_key, "population") == 0) {
                            pop_id = atoi(item_value);
                        } else if (strcmp(item_key, "source") == 0) {
                            pop_id = atoi(item_value);
                        } else if (strcmp(item_key, "destination") == 0) {
                            pop_id2 = atoi(item_value);
                        } else if (strcmp(item_key, "derived") == 0) {
                            pop_id = atoi(item_value);
                        } else if (strcmp(item_key, "ancestral") == 0) {
                            pop_id2 = atoi(item_value);
                        } else if (strcmp(item_key, "new_size") == 0) {
                            size_or_rate = atof(item_value);
                        } else if (strcmp(item_key, "rate") == 0) {
                            size_or_rate = atof(item_value);
                        }
                    }
                }
                yaml_event_delete(&ev_event);
            }
            
            // Add event to config
            if (config->num_explicit_events < 1000) {
                struct event *ev = &config->explicit_events[config->num_explicit_events];
                ev->time = event_time;
                ev->popID = pop_id;
                ev->popID2 = pop_id2;
                
                if (strcmp(event_type, "population_size_change") == 0) {
                    ev->type = 'g';
                    ev->popnSize = size_or_rate;
                } else if (strcmp(event_type, "population_split") == 0) {
                    ev->type = 'p';
                } else if (strcmp(event_type, "migration_rate_change") == 0) {
                    ev->type = 'm';
                    ev->popnSize = size_or_rate; // Migration rate stored in popnSize
                }
                
                config->num_explicit_events++;
            }
            
            // Check for next event
            if (!yaml_parser_parse(parser, &ev_event)) {
                return -1;
            }
            
            if (ev_event.type == YAML_SEQUENCE_END_EVENT) {
                yaml_event_delete(&ev_event);
                return 0;
            } else if (ev_event.type == YAML_MAPPING_START_EVENT) {
                yaml_event_delete(&ev_event);
                // Continue to parse next event
            } else {
                yaml_event_delete(&ev_event);
                return -1;
            }
        } while (1);
        
        return 0;
    }
    
    // Handle original parsing logic (sequence might be nested under a key)
    yaml_event_delete(&event);
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in events section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SEQUENCE_START_EVENT:
                // Handle direct event sequence under "events:" key
                {
                    yaml_event_t ev_event;
                    while (yaml_parser_parse(parser, &ev_event)) {
                        if (ev_event.type == YAML_SEQUENCE_END_EVENT) {
                            yaml_event_delete(&ev_event);
                            break;
                        }
                        
                        if (ev_event.type == YAML_MAPPING_START_EVENT) {
                            // Parse individual event
                            yaml_event_t item_event;
                            char item_key[256], item_value[256];
                            char event_type[256] = "";
                            double event_time = 0.0;
                            int pop_id = 0, pop_id2 = 0;
                            double size_or_rate = 0.0;
                            
                            while (yaml_parser_parse(parser, &item_event)) {
                                if (item_event.type == YAML_MAPPING_END_EVENT) {
                                    yaml_event_delete(&item_event);
                                    break;
                                }
                                
                                if (item_event.type == YAML_SCALAR_EVENT) {
                                    parseScalarValue(&item_event, item_key, sizeof(item_key));
                                    yaml_event_delete(&item_event);
                                    
                                    if (yaml_parser_parse(parser, &item_event) && item_event.type == YAML_SCALAR_EVENT) {
                                        parseScalarValue(&item_event, item_value, sizeof(item_value));
                                        
                                        if (strcmp(item_key, "type") == 0) {
                                            strncpy(event_type, item_value, sizeof(event_type) - 1);
                                        } else if (strcmp(item_key, "time") == 0) {
                                            event_time = atof(item_value);
                                        } else if (strcmp(item_key, "population") == 0) {
                                            pop_id = atoi(item_value);
                                        } else if (strcmp(item_key, "source") == 0) {
                                            pop_id = atoi(item_value);
                                        } else if (strcmp(item_key, "destination") == 0) {
                                            pop_id2 = atoi(item_value);
                                        } else if (strcmp(item_key, "derived") == 0) {
                                            pop_id = atoi(item_value);
                                        } else if (strcmp(item_key, "ancestral") == 0) {
                                            pop_id2 = atoi(item_value);
                                        } else if (strcmp(item_key, "new_size") == 0) {
                                            size_or_rate = atof(item_value);
                                        } else if (strcmp(item_key, "rate") == 0) {
                                            size_or_rate = atof(item_value);
                                        }
                                    }
                                }
                                yaml_event_delete(&item_event);
                            }
                            
                            // Add event to config
                            if (config->num_explicit_events < 1000) {
                                struct event *ev = &config->explicit_events[config->num_explicit_events];
                                ev->time = event_time;
                                ev->popID = pop_id;
                                ev->popID2 = pop_id2;
                                
                                if (strcmp(event_type, "population_size_change") == 0) {
                                    ev->type = 'g';
                                    ev->popnSize = size_or_rate;
                                } else if (strcmp(event_type, "population_split") == 0) {
                                    ev->type = 'p';
                                } else if (strcmp(event_type, "migration_rate_change") == 0) {
                                    ev->type = 'm';
                                    ev->popnSize = size_or_rate; // Migration rate stored in popnSize
                                }
                                
                                config->num_explicit_events++;
                            }
                        }
                        yaml_event_delete(&ev_event);
                    }
                }
                break;
                
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&event, value, sizeof(value));
                    
                    if (strcmp(key, "demes_file") == 0) {
                        strncpy(config->demes_file, value, sizeof(config->demes_file) - 1);
                        config->demes_file[sizeof(config->demes_file) - 1] = '\0';
                        config->use_demes = 1;
                    }
                } else if (event.type == YAML_SEQUENCE_START_EVENT) {
                    // Handle demographic or selection events sequences
                    if (strcmp(key, "selection") == 0) {
                        config->has_selection = 1;
                        
                        // Parse selection events sequence
                        yaml_event_t sel_event;
                        while (yaml_parser_parse(parser, &sel_event)) {
                            if (sel_event.type == YAML_SEQUENCE_END_EVENT) {
                                yaml_event_delete(&sel_event);
                                break;
                            }
                            
                            if (sel_event.type == YAML_MAPPING_START_EVENT) {
                                // Parse individual selection event
                                yaml_event_t sweep_event;
                                char sweep_key[256], sweep_value[256];
                                char sweep_type[256] = "";
                                char sweep_mode[256] = "stochastic";
                                double sweep_time = 0.0;
                                double selection_coeff = 0.0;
                                double sweep_position = 0.5;
                                double initial_freq = 0.0;
                                
                                while (yaml_parser_parse(parser, &sweep_event)) {
                                    if (sweep_event.type == YAML_MAPPING_END_EVENT) {
                                        yaml_event_delete(&sweep_event);
                                        break;
                                    }
                                    
                                    if (sweep_event.type == YAML_SCALAR_EVENT) {
                                        parseScalarValue(&sweep_event, sweep_key, sizeof(sweep_key));
                                        yaml_event_delete(&sweep_event);
                                        
                                        if (yaml_parser_parse(parser, &sweep_event) && sweep_event.type == YAML_SCALAR_EVENT) {
                                            parseScalarValue(&sweep_event, sweep_value, sizeof(sweep_value));
                                            
                                            if (strcmp(sweep_key, "type") == 0) {
                                                strncpy(sweep_type, sweep_value, sizeof(sweep_type) - 1);
                                            } else if (strcmp(sweep_key, "mode") == 0) {
                                                strncpy(sweep_mode, sweep_value, sizeof(sweep_mode) - 1);
                                            } else if (strcmp(sweep_key, "time") == 0) {
                                                sweep_time = atof(sweep_value);
                                            } else if (strcmp(sweep_key, "selection_coeff") == 0) {
                                                selection_coeff = atof(sweep_value);
                                            } else if (strcmp(sweep_key, "position") == 0) {
                                                sweep_position = atof(sweep_value);
                                            } else if (strcmp(sweep_key, "initial_freq") == 0) {
                                                initial_freq = atof(sweep_value);
                                            }
                                        }
                                    }
                                    yaml_event_delete(&sweep_event);
                                }
                                
                                // Apply parsed selection event to config
                                if (strcmp(sweep_type, "sweep") == 0) {
                                    config->alpha = selection_coeff;
                                    config->tau = sweep_time * 2.0;  // Convert to discoal time units
                                    config->sweep_site = sweep_position;
                                    
                                    if (strcmp(sweep_mode, "deterministic") == 0) {
                                        config->sweep_mode = 'd';
                                    } else if (strcmp(sweep_mode, "neutral") == 0) {
                                        config->sweep_mode = 'N';
                                    } else {
                                        config->sweep_mode = 's';  // default to stochastic
                                    }
                                    
                                    if (initial_freq > 0.0) {
                                        config->f0 = initial_freq;
                                        config->soft_sweep_mode = 1;
                                    }
                                }
                            }
                            yaml_event_delete(&sel_event);
                        }
                    } else if (strcmp(key, "demographic") == 0) {
                        // Parse demographic events sequence (only when not using demes)
                        if (!config->use_demes) {
                            yaml_event_t demo_event;
                            while (yaml_parser_parse(parser, &demo_event)) {
                                if (demo_event.type == YAML_SEQUENCE_END_EVENT) {
                                    yaml_event_delete(&demo_event);
                                    break;
                                }
                                
                                if (demo_event.type == YAML_MAPPING_START_EVENT) {
                                    // Parse individual demographic event
                                    yaml_event_t dem_event;
                                    char dem_key[256], dem_value[256];
                                    char event_type[256] = "";
                                    double event_time = 0.0;
                                    int pop_id = 0, pop_id2 = 0;
                                    double size_or_rate = 0.0;
                                    
                                    while (yaml_parser_parse(parser, &dem_event)) {
                                        if (dem_event.type == YAML_MAPPING_END_EVENT) {
                                            yaml_event_delete(&dem_event);
                                            break;
                                        }
                                        
                                        if (dem_event.type == YAML_SCALAR_EVENT) {
                                            parseScalarValue(&dem_event, dem_key, sizeof(dem_key));
                                            yaml_event_delete(&dem_event);
                                            
                                            if (yaml_parser_parse(parser, &dem_event) && dem_event.type == YAML_SCALAR_EVENT) {
                                                parseScalarValue(&dem_event, dem_value, sizeof(dem_value));
                                                
                                                if (strcmp(dem_key, "type") == 0) {
                                                    strncpy(event_type, dem_value, sizeof(event_type) - 1);
                                                } else if (strcmp(dem_key, "time") == 0) {
                                                    event_time = atof(dem_value);
                                                } else if (strcmp(dem_key, "population") == 0) {
                                                    pop_id = atoi(dem_value);
                                                } else if (strcmp(dem_key, "source_pop") == 0) {
                                                    pop_id = atoi(dem_value);
                                                } else if (strcmp(dem_key, "dest_pop") == 0) {
                                                    pop_id2 = atoi(dem_value);
                                                } else if (strcmp(dem_key, "size") == 0) {
                                                    size_or_rate = atof(dem_value);
                                                } else if (strcmp(dem_key, "rate") == 0) {
                                                    size_or_rate = atof(dem_value);
                                                }
                                            }
                                        }
                                        yaml_event_delete(&dem_event);
                                    }
                                    
                                    // Add demographic event to config
                                    if (config->num_explicit_events < 1000) {
                                        struct event *ev = &config->explicit_events[config->num_explicit_events];
                                        ev->time = event_time * 2.0;  // Convert to discoal time units
                                        ev->popID = pop_id;
                                        ev->popID2 = pop_id2;
                                        
                                        if (strcmp(event_type, "size_change") == 0) {
                                            ev->type = 'n';
                                            ev->popnSize = size_or_rate;
                                        } else if (strcmp(event_type, "population_split") == 0) {
                                            ev->type = 'p';
                                        } else if (strcmp(event_type, "migration_change") == 0) {
                                            ev->type = 'M';
                                            // Migration rate handling would be more complex
                                        }
                                        
                                        config->num_explicit_events++;
                                    }
                                }
                                yaml_event_delete(&demo_event);
                            }
                        } else {
                            // Skip demographic events if using demes file
                            fprintf(stderr, "Warning: Ignoring demographic events in YAML - using demes file instead\n");
                            int sequence_depth = 1;
                            while (sequence_depth > 0 && yaml_parser_parse(parser, &event)) {
                                if (event.type == YAML_SEQUENCE_START_EVENT) sequence_depth++;
                                else if (event.type == YAML_SEQUENCE_END_EVENT) sequence_depth--;
                                yaml_event_delete(&event);
                            }
                        }
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Parse selection section
static int parseSelectionSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    
    config->has_selection = 1;
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in selection section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&event, value, sizeof(value));
                    
                    // Process key-value pairs
                    if (strcmp(key, "alpha") == 0) {
                        config->alpha = atof(value);
                    } else if (strcmp(key, "sweep_site") == 0) {
                        config->sweep_site = atof(value);
                    } else if (strcmp(key, "sweep_mode") == 0) {
                        if (strcmp(value, "deterministic") == 0) {
                            config->sweep_mode = 'd';
                        } else if (strcmp(value, "neutral") == 0) {
                            config->sweep_mode = 'N';
                        } else {
                            config->sweep_mode = 's';  // default to stochastic
                        }
                    } else if (strcmp(key, "tau") == 0) {
                        config->tau = atof(value);
                    }
                } else if (event.type == YAML_MAPPING_START_EVENT) {
                    // Handle nested mappings like soft_sweep
                    if (strcmp(key, "soft_sweep") == 0) {
                        // Parse soft sweep sub-section
                        yaml_event_t ss_event;
                        char ss_key[256], ss_value[256];
                        
                        while (yaml_parser_parse(parser, &ss_event)) {
                            if (ss_event.type == YAML_MAPPING_END_EVENT) {
                                yaml_event_delete(&ss_event);
                                break;
                            }
                            
                            if (ss_event.type == YAML_SCALAR_EVENT) {
                                parseScalarValue(&ss_event, ss_key, sizeof(ss_key));
                                yaml_event_delete(&ss_event);
                                
                                if (yaml_parser_parse(parser, &ss_event) && ss_event.type == YAML_SCALAR_EVENT) {
                                    parseScalarValue(&ss_event, ss_value, sizeof(ss_value));
                                    
                                    if (strcmp(ss_key, "initial_frequency") == 0) {
                                        config->f0 = atof(ss_value);
                                        config->soft_sweep_mode = 1;
                                    }
                                }
                            }
                            yaml_event_delete(&ss_event);
                        }
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Parse output section
static int parseOutputSection(yaml_parser_t *parser, SimulationConfig *config) {
    yaml_event_t event;
    char key[256], value[256];
    int done = 0;
    
    config->has_output = 1;
    
    while (!done) {
        if (!yaml_parser_parse(parser, &event)) {
            fprintf(stderr, "Error parsing YAML in output section\n");
            return -1;
        }
        
        switch (event.type) {
            case YAML_SCALAR_EVENT:
                if (!parseScalarValue(&event, key, sizeof(key))) {
                    yaml_event_delete(&event);
                    continue;
                }
                
                // Get the value
                yaml_event_delete(&event);
                if (!yaml_parser_parse(parser, &event)) {
                    fprintf(stderr, "Error parsing YAML value for key: %s\n", key);
                    return -1;
                }
                
                if (event.type == YAML_SCALAR_EVENT) {
                    parseScalarValue(&event, value, sizeof(value));
                    
                    // Process key-value pairs
                    if (strcmp(key, "style") == 0) {
                        if (strcmp(value, "snp_matrix") == 0) {
                            config->output_style = 's';
                        } else if (strcmp(value, "haplotype") == 0) {
                            config->output_style = 'h';
                        }
                    } else if (strcmp(key, "finite_output") == 0) {
                        config->finite_output_flag = (strcmp(value, "true") == 0) ? 1 : 0;
                    } else if (strcmp(key, "hide_partial_snp") == 0) {
                        config->hide_partial_snp = (strcmp(value, "true") == 0) ? 1 : 0;
                    }
                } else if (event.type == YAML_MAPPING_START_EVENT) {
                    // Handle nested mappings like tskit
                    if (strcmp(key, "tskit") == 0) {
                        // Parse tskit sub-section
                        yaml_event_t ts_event;
                        char ts_key[256], ts_value[256];
                        
                        while (yaml_parser_parse(parser, &ts_event)) {
                            if (ts_event.type == YAML_MAPPING_END_EVENT) {
                                yaml_event_delete(&ts_event);
                                break;
                            }
                            
                            if (ts_event.type == YAML_SCALAR_EVENT) {
                                parseScalarValue(&ts_event, ts_key, sizeof(ts_key));
                                yaml_event_delete(&ts_event);
                                
                                if (yaml_parser_parse(parser, &ts_event) && ts_event.type == YAML_SCALAR_EVENT) {
                                    parseScalarValue(&ts_event, ts_value, sizeof(ts_value));
                                    
                                    if (strcmp(ts_key, "enabled") == 0) {
                                        config->tskit_output_mode = (strcmp(ts_value, "true") == 0) ? 1 : 0;
                                    } else if (strcmp(ts_key, "filename") == 0) {
                                        strncpy(config->tskit_output_filename, ts_value, sizeof(config->tskit_output_filename) - 1);
                                        config->tskit_output_filename[sizeof(config->tskit_output_filename) - 1] = '\0';
                                    } else if (strcmp(ts_key, "minimal") == 0) {
                                        config->minimal_tree_seq = (strcmp(ts_value, "true") == 0) ? 1 : 0;
                                    }
                                }
                            }
                            yaml_event_delete(&ts_event);
                        }
                    }
                }
                break;
                
            case YAML_MAPPING_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    return 0;
}

// Main function to load configuration from YAML file
int loadConfigFile(const char *filename, SimulationConfig *config) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open configuration file '%s'\n", filename);
        return -1;
    }
    
    yaml_parser_t parser;
    yaml_event_t event;
    
    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Error: Failed to initialize YAML parser\n");
        fclose(file);
        return -1;
    }
    
    yaml_parser_set_input_file(&parser, file);
    
    // Initialize config with defaults
    initializeDefaultConfig(config);
    
    // Parse YAML document
    int done = 0;
    char current_section[256] = "";
    
    while (!done) {
        if (!yaml_parser_parse(&parser, &event)) {
            fprintf(stderr, "Error: YAML parsing failed\n");
            yaml_parser_delete(&parser);
            fclose(file);
            return -1;
        }
        
        switch (event.type) {
            case YAML_STREAM_START_EVENT:
            case YAML_DOCUMENT_START_EVENT:
            case YAML_MAPPING_START_EVENT:
                break;
                
            case YAML_SCALAR_EVENT:
                {
                    char key[256];
                    if (parseScalarValue(&event, key, sizeof(key))) {
                        strcpy(current_section, key);
                        
                        yaml_event_delete(&event);
                        if (!yaml_parser_parse(&parser, &event)) {
                            fprintf(stderr, "Error parsing section: %s\n", key);
                            break;
                        }
                        
                        if (event.type == YAML_MAPPING_START_EVENT) {
                            if (strcmp(key, "simulation") == 0) {
                                yaml_event_delete(&event);
                                if (parseSimulationSection(&parser, config) != 0) {
                                    done = 1;
                                }
                                continue;
                            } else if (strcmp(key, "genetics") == 0) {
                                yaml_event_delete(&event);
                                if (parseGeneticsSection(&parser, config) != 0) {
                                    done = 1;
                                }
                                continue;
                            } else if (strcmp(key, "populations") == 0) {
                                yaml_event_delete(&event);
                                if (parsePopulationsSection(&parser, config) != 0) {
                                    done = 1;
                                }
                                continue;
                            } else if (strcmp(key, "selection") == 0) {
                                yaml_event_delete(&event);
                                if (parseSelectionSection(&parser, config) != 0) {
                                    done = 1;
                                }
                                continue;
                            } else if (strcmp(key, "output") == 0) {
                                yaml_event_delete(&event);
                                if (parseOutputSection(&parser, config) != 0) {
                                    done = 1;
                                }
                                continue;
                            }
                        } else if (event.type == YAML_SEQUENCE_START_EVENT && strcmp(key, "events") == 0) {
                            // Handle events as a direct sequence
                            yaml_event_delete(&event);
                            if (parseEventsSection(&parser, config) != 0) {
                                done = 1;
                            }
                            continue;
                        } else if (event.type == YAML_SCALAR_EVENT && strcmp(key, "demes") == 0) {
                            // Handle demes: filename directly
                            char value[256];
                            if (parseScalarValue(&event, value, sizeof(value))) {
                                strncpy(config->demes_file, value, sizeof(config->demes_file) - 1);
                                config->demes_file[sizeof(config->demes_file) - 1] = '\0';
                                config->use_demes = 1;
                            }
                        }
                    }
                }
                break;
                
            case YAML_DOCUMENT_END_EVENT:
            case YAML_STREAM_END_EVENT:
                done = 1;
                break;
                
            default:
                break;
        }
        
        yaml_event_delete(&event);
    }
    
    yaml_parser_delete(&parser);
    fclose(file);
    
    return 0;
}

// Apply configuration to global variables (called after command line parsing for backwards compatibility)
int applyConfiguration(const SimulationConfig *config) {
    extern int sampleSize, sampleNumber, nSites;
    extern long seed1, seed2;
    extern int npops;
    extern int sampleSizes[MAXPOPS];
    extern int EFFECTIVE_POPN_SIZE;
    extern double theta, rho, my_gamma, gammaCoRatio;
    extern int gcMean;
    extern double gammaCoRatioMode;
    extern double alpha, sweepSite, tau, f0, uA, partialSweepFinalFreq, recurSweepRate;
    extern char sweepMode;
    extern int recurSweepMode, partialSweepMode, softSweepMode;
    extern int tskitOutputMode, minimalTreeSeq, hidePartialSNP;
    extern char tskitOutputFilename[1024];
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern int migFlag;
    
    // Apply simulation settings from YAML (these are defaults that can be overridden by command line)
    if (config->sample_size > 0) {
        sampleSize = config->sample_size;
    }
    if (config->num_replicates > 0) {
        sampleNumber = config->num_replicates;
    }
    if (config->num_sites > 0) {
        nSites = config->num_sites;
    }
    if (config->has_seed) {
        // fprintf(stderr, "DEBUG: Applying seeds from YAML: %d, %d\n", config->seed1, config->seed2);
        seed1 = config->seed1;
        seed2 = config->seed2;
        // fprintf(stderr, "DEBUG: Seeds after application: %ld, %ld\n", seed1, seed2);
    }
    if (config->effective_popn_size > 0) {
        EFFECTIVE_POPN_SIZE = config->effective_popn_size;
    }
    
    // Apply genetics parameters from YAML
    if (config->has_genetics) {
        if (config->theta > 0) {
            // fprintf(stderr, "DEBUG: Applying theta from YAML: %f\n", config->theta);
            theta = config->theta;
        }
        if (config->rho > 0) {
            rho = config->rho;
        }
        if (config->gamma_co_ratio_mode) {
            gammaCoRatio = config->gamma_co_ratio;
            gammaCoRatioMode = 1;
            gcMean = config->gc_mean;
            // Gene conversion rate will be calculated as rho * gammaCoRatio later
        } else if (config->gamma > 0) {
            my_gamma = config->gamma;
            gcMean = config->gc_mean;
        }
    }
    
    // Apply population configuration from YAML
    if (config->has_populations) {
        npops = config->npops;
        for (int i = 0; i < npops && i < MAXPOPS; i++) {
            sampleSizes[i] = config->sample_sizes[i];
        }
    }
    
    // Apply selection parameters from YAML
    if (config->has_selection) {
        if (config->alpha > 0) {
            alpha = config->alpha;
        }
        if (config->sweep_site >= 0) {
            sweepSite = config->sweep_site;
        }
        if (config->tau > 0) {
            tau = config->tau;
        }
        sweepMode = config->sweep_mode;
        if (config->soft_sweep_mode) {
            softSweepMode = config->soft_sweep_mode;
            if (config->f0 > 0) {
                f0 = config->f0;
            }
        }
        
        // Add sweep event to events array
        ensureEventsCapacity();
        events[eventNumber].time = tau;
        events[eventNumber].type = 's'; // sweep event
        eventNumber++;
    }
    
    // Apply explicit demographic events (only if not using demes)
    if (!config->use_demes && config->num_explicit_events > 0) {
        for (int i = 0; i < config->num_explicit_events; i++) {
            ensureEventsCapacity();
            events[eventNumber] = config->explicit_events[i];
            eventNumber++;
        }
    }
    
    // Handle demes file if specified
    if (config->use_demes && strlen(config->demes_file) > 0) {
        int ret = loadDemesFile(config->demes_file, &events, &eventNumber, &eventsCapacity, 
                               currentSize, &npops, sampleSizes, EFFECTIVE_POPN_SIZE);
        if (ret != 0) {
            fprintf(stderr, "Error: Failed to load demes file '%s' from YAML config\n", config->demes_file);
            return -1;
        }
        
        fprintf(stderr, "Loaded %d populations and %d events from demes file '%s' (via YAML config)\n", 
                npops, eventNumber - 1, config->demes_file);
    }
    
    return 0;
}
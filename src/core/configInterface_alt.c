#include "configInterface_alt.h"
#include <assert.h>
#include <string.h>

void ensureEventsCapacity();

/* for string-valued options, use enums so that CYAML can automatically check
 * for invalid values */
static const cyaml_strval_t output_types_strings[] = {
    {"haplotype",     OUTPUT_HAPLOTYPE},
    {"snp_array",     OUTPUT_SNP_ARRAY},
    {"tree_sequence", OUTPUT_TREE_SEQN},
};

static const cyaml_strval_t sweep_modes_strings[] = {
    {"stochastic",    SWEEP_STOCHASTIC},
    {"deterministic", SWEEP_DETERMINISTIC},
    {"neutral",       SWEEP_NEUTRAL},
};


/* for sequence valued options, define a schema for the data type used
 * held by the sequence */
static const cyaml_schema_value_t int_array_schema = {
    CYAML_VALUE_INT(CYAML_FLAG_DEFAULT, int),
};

static const cyaml_schema_value_t float_array_schema = {
    CYAML_VALUE_FLOAT(CYAML_FLAG_DEFAULT, double),
};


/* schema for fields within simulation block */
static const cyaml_schema_field_t simulation_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_INT("sample_size", CYAML_FLAG_DEFAULT, 
        struct simulation_config, sample_size),
    CYAML_FIELD_INT("num_replicates", CYAML_FLAG_DEFAULT, 
        struct simulation_config, num_replicates),
    CYAML_FIELD_INT("num_sites", CYAML_FLAG_DEFAULT, 
        struct simulation_config, num_sites),
    /* optional arguments */
    CYAML_FIELD_SEQUENCE_FIXED("seed", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL, 
        struct simulation_config, seed, &int_array_schema, 2),
    CYAML_FIELD_END
};


/* schema for fields within genetics block */
static const cyaml_schema_field_t genetics_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("mutation_rate", CYAML_FLAG_DEFAULT, 
        struct genetics_config, mutation_rate),
    CYAML_FIELD_FLOAT("recombination_rate", CYAML_FLAG_DEFAULT, 
        struct genetics_config, recombination_rate),
    /* optional arguments */
    CYAML_FIELD_FLOAT("gene_conversion_rate", CYAML_FLAG_OPTIONAL, 
        struct genetics_config, gene_conversion_rate),
    CYAML_FIELD_INT("gene_conversion_tract_length", CYAML_FLAG_OPTIONAL, 
        struct genetics_config, gene_conversion_tract_length),
    CYAML_FIELD_FLOAT("crossover_ratio", CYAML_FLAG_OPTIONAL, 
        struct genetics_config, crossover_ratio),
    CYAML_FIELD_END
};


/* schema for specific demographic event types and event sub-block */
static const cyaml_schema_field_t population_size_change_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct population_size_change, time),
    CYAML_FIELD_FLOAT("size", CYAML_FLAG_DEFAULT, struct population_size_change, size),
    CYAML_FIELD_INT("population", CYAML_FLAG_DEFAULT, struct population_size_change, population),
    /* optional arguments */
    CYAML_FIELD_END
};
static const cyaml_schema_value_t population_size_change_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct population_size_change, 
        population_size_change_fields_schema),
};

static const cyaml_schema_field_t migration_rate_change_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct migration_rate_change, time),
    CYAML_FIELD_FLOAT("rate", CYAML_FLAG_DEFAULT, struct migration_rate_change, rate),
    CYAML_FIELD_INT("source", CYAML_FLAG_DEFAULT, struct migration_rate_change, source),
    CYAML_FIELD_INT("destination", CYAML_FLAG_DEFAULT, struct migration_rate_change, destination),
    /* optional arguments */
    CYAML_FIELD_END
};
static const cyaml_schema_value_t migration_rate_change_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct migration_rate_change, 
        migration_rate_change_fields_schema),
};

static const cyaml_schema_field_t population_split_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct population_split, time),
    CYAML_FIELD_INT("ancestral", CYAML_FLAG_DEFAULT, struct population_split, ancestral),
    CYAML_FIELD_INT("derived", CYAML_FLAG_DEFAULT, struct population_split, derived),
    /* optional arguments */
    CYAML_FIELD_END
};
static const cyaml_schema_value_t population_split_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct population_split, 
        population_split_fields_schema),
};

static const cyaml_schema_field_t ancient_sample_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct ancient_sample, time),
    CYAML_FIELD_INT("population", CYAML_FLAG_DEFAULT, struct ancient_sample, population),
    CYAML_FIELD_INT("sample_size", CYAML_FLAG_DEFAULT, struct ancient_sample, sample_size),
    /* optional arguments */
    CYAML_FIELD_END
};
static const cyaml_schema_value_t ancient_sample_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct ancient_sample, 
        ancient_sample_fields_schema),
};

static const cyaml_schema_field_t demographic_events_fields_schema[] = { /* sub-block */
    /* optional arguments */
    CYAML_FIELD_SEQUENCE_COUNT("population_size_changes", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demographic_events, population_size_changes, num_population_size_changes,
        &population_size_change_schema, 1, CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE_COUNT("migration_rate_changes", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demographic_events, migration_rate_changes, num_migration_rate_changes,
        &migration_rate_change_schema, 1, CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE_COUNT("population_splits", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demographic_events, population_splits, num_population_splits,
        &population_split_schema, 1, CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE_COUNT("ancient_samples", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demographic_events, ancient_samples, num_ancient_samples,
        &ancient_sample_schema, 1, CYAML_UNLIMITED),
    CYAML_FIELD_END
};


/* schema for migration matrix */
static const cyaml_schema_field_t migration_matrix_row_fields_schema[] = { 
    CYAML_FIELD_SEQUENCE_COUNT("row", CYAML_FLAG_POINTER,
        struct migration_matrix_row, rates, num_cols,
        &float_array_schema, 0, CYAML_UNLIMITED),
    CYAML_FIELD_END
};
static const cyaml_schema_value_t migration_matrix_row_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct migration_matrix_row, 
        migration_matrix_row_fields_schema),
};


/* schema for fields within demography block */
static const cyaml_schema_field_t demography_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_SEQUENCE_COUNT("deme_sample_size", CYAML_FLAG_POINTER,
        struct demography_config, deme_sample_size, num_demes, &int_array_schema, 
        1, MAXPOPS),
    /* optional arguments */
    CYAML_FIELD_FLOAT("effective_population_size", CYAML_FLAG_OPTIONAL, 
        struct demography_config, effective_population_size),
    CYAML_FIELD_MAPPING_PTR("demographic_events", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demography_config, demographic_events, demographic_events_fields_schema),
    CYAML_FIELD_SEQUENCE_COUNT("migration_matrix", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demography_config, migration_matrix, num_migration_matrix_rows,
        &migration_matrix_row_schema, 1, MAXPOPS),
    CYAML_FIELD_STRING_PTR("demes_filename", CYAML_FLAG_OPTIONAL,
        struct demography_config, demes_filename, 0, PATH_MAX),
    CYAML_FIELD_END
};


/* schema for fields within selection block */
static const cyaml_schema_field_t selection_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_ENUM("sweep_mode", CYAML_FLAG_DEFAULT, 
        struct selection_config, sweep_mode, 
        sweep_modes_strings, CYAML_ARRAY_LEN(sweep_modes_strings)),
    CYAML_FIELD_FLOAT("selection_coefficient", CYAML_FLAG_DEFAULT,
        struct selection_config, selection_coefficient),
    CYAML_FIELD_FLOAT("sweep_position", CYAML_FLAG_DEFAULT,
        struct selection_config, sweep_position),
    /* optional arguments */
    CYAML_FIELD_FLOAT("fixation_time_ago", CYAML_FLAG_OPTIONAL,
        struct selection_config, fixation_time_ago),
    CYAML_FIELD_FLOAT("initial_frequency", CYAML_FLAG_OPTIONAL,
        struct selection_config, initial_frequency),
    CYAML_FIELD_FLOAT("final_frequency", CYAML_FLAG_OPTIONAL,
        struct selection_config, final_frequency),
    CYAML_FIELD_FLOAT("beneficial_mutation_rate", CYAML_FLAG_OPTIONAL,
        struct selection_config, beneficial_mutation_rate),
    CYAML_FIELD_FLOAT("recurrent_sweep_rate", CYAML_FLAG_OPTIONAL,
        struct selection_config, recurrent_sweep_rate),
    CYAML_FIELD_END
};


/* schema for fields within output block */
static const cyaml_schema_field_t output_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_ENUM("output_type", CYAML_FLAG_DEFAULT, 
        struct output_config, output_type, 
        output_types_strings, CYAML_ARRAY_LEN(output_types_strings)),
    /* optional arguments */
    CYAML_FIELD_BOOL("finite_output", CYAML_FLAG_OPTIONAL,
        struct output_config, finite_output),
    CYAML_FIELD_BOOL("hide_partial_snp", CYAML_FLAG_OPTIONAL,
        struct output_config, hide_partial_snp),
    CYAML_FIELD_BOOL("unsimplified_tree_sequence", CYAML_FLAG_OPTIONAL,
        struct output_config, unsimplified_tree_sequence),
    CYAML_FIELD_STRING_PTR("tree_sequence_filename", CYAML_FLAG_OPTIONAL,
        struct output_config, tree_sequence_filename, 0, PATH_MAX),
    CYAML_FIELD_END
};


/* top level schema */
static const cyaml_schema_field_t discoal_config_fields_schema[] = {
    /* required blocks */
    CYAML_FIELD_MAPPING_PTR("simulation", CYAML_FLAG_POINTER,
        struct discoal_config, simulation, simulation_fields_schema),
    CYAML_FIELD_MAPPING_PTR("genetics", CYAML_FLAG_POINTER,
        struct discoal_config, genetics, genetics_fields_schema),
    /* optional blocks */
    CYAML_FIELD_MAPPING_PTR("demography", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct discoal_config, demography, demography_fields_schema),
    CYAML_FIELD_MAPPING_PTR("selection", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct discoal_config, selection, selection_fields_schema),
    CYAML_FIELD_MAPPING_PTR("output", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct discoal_config, output, output_fields_schema),
    CYAML_FIELD_END
};
static const cyaml_schema_value_t discoal_config_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER,
        struct discoal_config, discoal_config_fields_schema),
};


/* options for cyaml parsing */
static const cyaml_config_t cyaml_config = {
    .log_fn = cyaml_log,            
    .mem_fn = cyaml_mem,
    .log_level = CYAML_LOG_WARNING,
};


/* FIXME: initialization should be explicit, and we probably don't want to
 * constrain yaml to positive values
 *
 * libcyaml is guaranteed to initialize primitive types to zero and pointers to NULL.
 * Hence, we do not worry about initializing the config structs explicitly, and instead
 * only copy into extern fields if the config members are positive.
 *
 * This might need more careful handling if the defaults are non-positive, and a
 * particular member of the config is explicitly set to zero or less in the yaml.
 *
 * For booleans, we need to ensure that the option is always false by default.
 */

int parse_simulation_block(struct simulation_config *cfg) 
{
    extern int sampleSize, sampleNumber, nSites;
    extern long seed1, seed2;
    if (cfg != NULL) {
        if (cfg->sample_size > 0) {
            sampleSize = cfg->sample_size;
        }
        if (cfg->num_replicates > 0) {
            sampleNumber = cfg->num_replicates;
        }
        if (cfg->num_sites > 0) {
            nSites = cfg->num_sites;
        }
        if (cfg->seed != NULL) {
            seed1 = cfg->seed[0];
            seed2 = cfg->seed[1];
        }
    }
    return EXIT_SUCCESS;
}

int parse_genetics_block(struct genetics_config *cfg) 
{
    extern double theta, rho;
    extern double gammaCoRatio, my_gamma, gammaCoRatioMode;
    extern int gcMean;
    if (cfg != NULL) {
        if (cfg->mutation_rate > 0) {
            theta = cfg->mutation_rate;
        }
        if (cfg->recombination_rate > 0) {
            rho = cfg->recombination_rate;
        }
        if (cfg->crossover_ratio > 0) {
            if (cfg->gene_conversion_rate > 0) {
                fprintf(stderr, 
                  "Error parsing config: `gene_conversion_rate` "
                  "and `crossover_ratio` cannot both be set\n"
                );
                return EXIT_FAILURE;
            }
            gammaCoRatioMode = 1;
            gammaCoRatio = cfg->crossover_ratio;
        }
        if (cfg->gene_conversion_rate > 0) {
            gammaCoRatioMode = 0; // FIXME: is this correct?
            my_gamma = cfg->gene_conversion_rate;
        }
        if (cfg->gene_conversion_tract_length > 0) {
            // FIXME: check that is positive, if one of the gc rates is set?
            gcMean = cfg->gene_conversion_tract_length;
        }
    }
    return EXIT_SUCCESS;
}

int parse_demography_block(struct demography_config *cfg) 
{
    extern int sampleSizes[MAXPOPS];
    extern int npops;
    extern int migFlag;
    extern int EFFECTIVE_POPN_SIZE;
    extern int eventNumber, eventsCapacity;
    extern struct event *events;
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern double *currentSize;
    extern double tDiv;
    if (cfg != NULL) {
        npops = cfg->num_demes;
        assert(cfg->num_demes > 0);
        for (int i = 0; i < cfg->num_demes; ++i) {
            sampleSizes[i] = cfg->deme_sample_size[i];
            currentSize[i] = 1.0;
            // FIXME: do we need to check that these sum to sampleSize or is that done elsewhere?
        }
        if (cfg->effective_population_size > 0) {
            EFFECTIVE_POPN_SIZE = (int)(cfg->effective_population_size);
            // FIXME: should this really be an int?
        }
        /* parse demographic events; these will be sorted into time order later */
        if (cfg->demographic_events != NULL) {
            if (cfg->demes_filename != NULL) {
                fprintf(stderr, "Cannot use demographic_events if demes_filename is provided\n");
                return EXIT_FAILURE;
            }
            struct demographic_events *dmo = cfg->demographic_events;
            if (dmo->num_ancient_samples > 0) {
                /* FIXME: not sure how to implement parsing */
                /* FIXME: are these counted in deme_sample_sizes? */
                fprintf(stderr, "Ancient sample events not yet implemented\n");
                return EXIT_FAILURE;
            }
            /* FIXME: do we need to check population indices < npops, times are positive etc */
            /* Time-varying migration rate changes are not implemented in the
             * main event loop yet (see discoalFunctions.c); reject rather than
             * silently dropping them. */
            if (dmo->num_migration_rate_changes > 0) {
                fprintf(stderr,
                    "Error parsing config: migration_rate_changes are not yet "
                    "implemented (no handler in the main event loop)\n");
                return EXIT_FAILURE;
            }
            /* Times in YAML follow the same convention as the command line:
             * the user supplies them in 2N units and the parser scales by 2 to
             * convert to discoal's internal 4N units. Keep this in sync with
             * `-en`, `-ed`, and `-ws` time scaling in getParameters(). */
            for (int i = 0; i < dmo->num_population_size_changes; ++i) {
                ensureEventsCapacity();
                events[eventNumber].type = 'n';
                events[eventNumber].time = dmo->population_size_changes[i].time * 2.0;
                events[eventNumber].popID = dmo->population_size_changes[i].population;
                events[eventNumber].popnSize = dmo->population_size_changes[i].size;
                eventNumber++;
            }
            for (int i = 0; i < dmo->num_population_splits; ++i) {
                ensureEventsCapacity();
                events[eventNumber].type = 'p';
                events[eventNumber].time = dmo->population_splits[i].time * 2.0;
                events[eventNumber].popID = dmo->population_splits[i].derived;
                events[eventNumber].popID2 = dmo->population_splits[i].ancestral;
                eventNumber++;
                tDiv = dmo->population_splits[i].time;  /* mark merger model active */
            }
        }
        /* parse migration matrix */
        if (cfg->migration_matrix != NULL) {
            if (cfg->demes_filename != NULL) {
                fprintf(stderr, "Cannot use migration_matrix if demes_filename "
                    "is provided\n");
                return EXIT_FAILURE;
            }
            if (cfg->num_migration_matrix_rows != cfg->num_demes) {
                fprintf(stderr, "Number of migration matrix rows does not "
                    "match number of demes\n");
                return EXIT_FAILURE;
            }
            for (int i = 0; i < cfg->num_migration_matrix_rows; ++i) {
                struct migration_matrix_row *row = &cfg->migration_matrix[i];
                if (row->num_cols != cfg->num_demes) {
                    fprintf(stderr, "Number of elements in migration matrix row "
                        "does not match number of demes\n");
                    return EXIT_FAILURE;
                }
                for (int j = 0; j < row->num_cols; ++j) {
                    /* FIXME: assert rate > 0? */
                    migMatConst[i][j] = row->rates[j];
                }
            }
            migFlag = 1;  /* set migration mode */
        }
        /* parse demes YAML into events */
        if (cfg->demes_filename != NULL) {
            int ret = loadDemesFile(cfg->demes_filename, &events, &eventNumber, 
                &eventsCapacity, currentSize, &npops, sampleSizes, EFFECTIVE_POPN_SIZE);
            if (ret != 0) {
                fprintf(stderr, "Error: Failed to load demes file '%s' from YAML config\n", 
                    cfg->demes_filename);
                return EXIT_FAILURE;
            }
            fprintf(stderr, 
                "Loaded %d populations and %d events from demes file '%s' "
                "(via YAML config)\n", npops, eventNumber - 1, 
                cfg->demes_filename);
        }
    }
    return EXIT_SUCCESS;
}

int parse_selection_block(struct selection_config *cfg) 
{
    extern double alpha, sweepSite, tau, f0, uA;
    extern double partialSweepFinalFreq, recurSweepRate;
    extern int recurSweepMode, partialSweepMode, softSweepMode;
    extern char sweepMode;
    extern int eventNumber, eventsCapacity;
    extern struct event *events;
    if (cfg != NULL) {
        switch (cfg->sweep_mode) {
            case SWEEP_STOCHASTIC:
                sweepMode = 's';
                break;
            case SWEEP_DETERMINISTIC:
                sweepMode = 'd';
                break;
            case SWEEP_NEUTRAL:
                sweepMode = 'N';
                break;
            default:
                break;
        } /* FIXME: need to add recurrent sweep modes */
        if (cfg->selection_coefficient > 0) {
            alpha = cfg->selection_coefficient;
        }
        if (cfg->sweep_position >= 0) { /* FIXME: need initialization as default is 0.5 */
            /* FIXME: check 0 <= x <= 1 here? */
            sweepSite = cfg->sweep_position;
        }
        if (cfg->fixation_time_ago > 0) {
            /* User supplies tau in 2N units (matching `-ws`); scale to 4N. */
            tau = cfg->fixation_time_ago * 2.0;
        }
        if (cfg->initial_frequency > 0) {
            /* FIXME: check 0 <= f <= 1 here? */
            f0 = cfg->initial_frequency;
            softSweepMode = 1;
        }
        if (cfg->final_frequency > 0) {
            /* FIXME: check 0 <= f <= 1 here? */
            partialSweepFinalFreq = cfg->final_frequency;
            partialSweepMode = 1;
        }
        if (cfg->beneficial_mutation_rate > 0) {
            uA = cfg->beneficial_mutation_rate;
        }
        if (cfg->recurrent_sweep_rate > 0) {
            recurSweepRate = cfg->recurrent_sweep_rate;
            recurSweepMode = 1;
        }

        if (cfg->fixation_time_ago > 0 && cfg->recurrent_sweep_rate > 0) {
            fprintf(stderr,
                "Error parsing config: fixation_time_ago and "
                "recurrent_sweep_rate cannot both be set\n");
            return EXIT_FAILURE;
        }

        /* Single sweep (matching `-ws`/`-wd`/`-wn`) requires fixation_time_ago
         * and emits an `'s'` event at tau. Recurrent sweeps (matching `-R`)
         * are driven by recurSweepMode + recurSweepRate alone and do not
         * produce an event. */
        if (cfg->recurrent_sweep_rate <= 0) {
            if (cfg->fixation_time_ago <= 0) {
                fprintf(stderr,
                    "Error parsing config: selection block requires "
                    "fixation_time_ago > 0 unless recurrent_sweep_rate "
                    "is set\n");
                return EXIT_FAILURE;
            }
            ensureEventsCapacity();
            events[eventNumber].time = tau;
            events[eventNumber].type = 's';
            eventNumber++;
        }
    }
    return EXIT_SUCCESS;
}

int parse_output_block(struct output_config *cfg) 
{
    extern int tskitOutputMode, minimalTreeSeq, hidePartialSNP;
    extern int finiteOutputFlag;
    extern char tskitOutputFilename[1024];
    if (cfg != NULL) {
        switch (cfg->output_type) {
            case OUTPUT_HAPLOTYPE:
                // FIXME: it is not clear what should be done here
                break;
            case OUTPUT_SNP_ARRAY:
                // FIXME: it is not clear what should be done here
                break;
            case OUTPUT_TREE_SEQN:
                if (cfg->tree_sequence_filename == NULL) {
                    fprintf(stderr, 
                        "Must provide tree_sequence_filename if using "
                        "output mode tree_sequence\n");
                    return EXIT_FAILURE;
                }
                tskitOutputMode = 1;
                break;
            default:
                break;
        }
        if (cfg->finite_output) {
            finiteOutputFlag = 1;
        }
        if (cfg->hide_partial_snp) {
            hidePartialSNP = 1;
        }
        if (cfg->unsimplified_tree_sequence) {
            minimalTreeSeq = 0;
        }
        if (cfg->tree_sequence_filename != NULL) {
            if (cfg->output_type != OUTPUT_TREE_SEQN) {
                fprintf(stderr, 
                    "Can only provide tree_sequence_filename if using "
                    "output mode tree_sequence\n");
                return EXIT_FAILURE;
            }
            // FIXME: safer way to do this?
            strncpy(tskitOutputFilename, cfg->tree_sequence_filename, 
                sizeof(tskitOutputFilename) - 1);
            tskitOutputFilename[sizeof(tskitOutputFilename) - 1] = '\0';
        }
    }
    return EXIT_SUCCESS;
}

int apply_yaml_config(struct discoal_config *config)
{
    int ret = EXIT_SUCCESS;
    assert(config != NULL);
    ret = parse_simulation_block(config->simulation);
    if (ret != EXIT_SUCCESS) { goto out; }
    ret = parse_genetics_block(config->genetics);
    if (ret != EXIT_SUCCESS) { goto out; }
    ret = parse_demography_block(config->demography);
    if (ret != EXIT_SUCCESS) { goto out; }
    ret = parse_selection_block(config->selection);
    if (ret != EXIT_SUCCESS) { goto out; }
    ret = parse_output_block(config->output);
    if (ret != EXIT_SUCCESS) { goto out; }
out:
    cyaml_free(&cyaml_config, &discoal_config_schema, config, 0);
    return ret;
}

int load_yaml_config(const char *yaml_path, struct discoal_config **config)
{
    int err;
    assert(*config == NULL);
    err = cyaml_load_file(yaml_path, &cyaml_config,
        &discoal_config_schema, (void **) config, NULL);
    if (err != CYAML_OK) {
        fprintf(stderr, "ERROR: %s\n", cyaml_strerror(err));
        cyaml_free(&cyaml_config, &discoal_config_schema, config, 0);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

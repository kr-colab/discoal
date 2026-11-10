#include "configInterface.h"
#include "discoalFunctions.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>

void ensureEventsCapacity();


/* input validation */
static int check_pop_index(const char *arr, const char *field, int idx,
    int pop, unsigned num_demes)
{
    if (pop < 0 || pop >= (int)num_demes) {
        fprintf(stderr,
            "Error parsing config: %s[%d].%s (%d) must be in [0, %u)\n",
            arr, idx, field, pop, num_demes);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int check_positive_double(const char *arr, const char *field, int idx,
    double val)
{
    if (val <= 0) {
        fprintf(stderr,
            "Error parsing config: %s[%d].%s (%g) must be > 0\n",
            arr, idx, field, val);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int check_nonneg_double(const char *arr, const char *field, int idx,
    double val)
{
    if (val < 0) {
        fprintf(stderr,
            "Error parsing config: %s[%d].%s (%g) must be >= 0\n",
            arr, idx, field, val);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/* A population no longer exists from the time a population_splits entry
 * merges it (as `derived`) into another one, or an admixtures entry moves
 * all its lineages into two others (as `population`). Reject an entry of
 * `arr` that involves population `pop` at or after that time.
 * skip_admixture is the index of an admixtures entry to ignore (the entry
 * being checked, which ends `pop` itself), or -1. */
static int check_not_merged(const char *arr, int idx, int pop, double time,
                            const struct demographic_events *dmo,
                            int skip_admixture) {
    for (int k = 0; k < dmo->num_population_splits; ++k) {
        const struct population_split *sp = &dmo->population_splits[k];
        if (pop == sp->derived && time >= sp->time) {
            fprintf(stderr,
                "Error parsing config: %s[%d] at time %g involves population "
                "%d, which merges into population %d at time %g "
                "(population_splits[%d])\n",
                arr, idx, time, pop, sp->ancestral, sp->time, k);
            return EXIT_FAILURE;
        }
    }
    for (int k = 0; k < dmo->num_admixtures; ++k) {
        const struct admixture *ad = &dmo->admixtures[k];
        if (k == skip_admixture) continue;
        if (pop == ad->population && time >= ad->time) {
            fprintf(stderr,
                "Error parsing config: %s[%d] at time %g involves population "
                "%d, which is admixed from populations %d and %d at time %g "
                "(admixtures[%d])\n",
                arr, idx, time, pop, ad->source_1, ad->source_2, ad->time, k);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

/* Unlike -ed and -ea, a YAML split or admixture ends a population, as in
 * demes: set all migration into and out of it to zero at that time.
 * Otherwise lineages can migrate into the ended population and be
 * stranded there. 'm' events set migShape[popID2][popID] (lineages move
 * popID2 -> popID). */
static void end_population_migration(int pop, double time, unsigned num_demes) {
    for (int k = 0; k < (int)num_demes; ++k) {
        if (k == pop) continue;
        for (int dir = 0; dir < 2; ++dir) {
            ensureEventsCapacity();
            events[eventNumber].type = 'm';
            events[eventNumber].time = userTimeToInternal(time);
            events[eventNumber].popID2 = dir ? pop : k;
            events[eventNumber].popID = dir ? k : pop;
            events[eventNumber].popnSize = 0.0;
            eventNumber++;
        }
    }
}

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
    {"trajectory_file", SWEEP_TRAJECTORY_FILE},
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

static const cyaml_schema_field_t admixture_fields_schema[] = {
    /* required arguments */
    CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct admixture, time),
    CYAML_FIELD_INT("population", CYAML_FLAG_DEFAULT, struct admixture, population),
    CYAML_FIELD_INT("source_1", CYAML_FLAG_DEFAULT, struct admixture, source_1),
    CYAML_FIELD_INT("source_2", CYAML_FLAG_DEFAULT, struct admixture, source_2),
    CYAML_FIELD_FLOAT("proportion", CYAML_FLAG_DEFAULT, struct admixture, proportion),
    CYAML_FIELD_END
};
static const cyaml_schema_value_t admixture_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct admixture,
        admixture_fields_schema),
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
    CYAML_FIELD_SEQUENCE_COUNT("admixtures", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct demographic_events, admixtures, num_admixtures,
        &admixture_schema, 1, CYAML_UNLIMITED),
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
    /* optional arguments; selection_coefficient is required unless
     * sweep_mode is trajectory_file */
    CYAML_FIELD_FLOAT_PTR("selection_coefficient",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct selection_config, selection_coefficient),
    CYAML_FIELD_STRING_PTR("trajectory_file", CYAML_FLAG_OPTIONAL,
        struct selection_config, trajectory_file, 1, CYAML_UNLIMITED),
    CYAML_FIELD_FLOAT_PTR("sweep_position",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct selection_config, sweep_position),
    CYAML_FIELD_FLOAT_PTR("fixation_time_ago",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct selection_config, fixation_time_ago),
    CYAML_FIELD_FLOAT("initial_frequency", CYAML_FLAG_OPTIONAL,
        struct selection_config, initial_frequency),
    CYAML_FIELD_FLOAT("final_frequency", CYAML_FLAG_OPTIONAL,
        struct selection_config, final_frequency),
    CYAML_FIELD_FLOAT("beneficial_mutation_rate", CYAML_FLAG_OPTIONAL,
        struct selection_config, beneficial_mutation_rate),
    CYAML_FIELD_FLOAT("recurrent_sweep_rate", CYAML_FLAG_OPTIONAL,
        struct selection_config, recurrent_sweep_rate),
    CYAML_FIELD_FLOAT_PTR("left_rho", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct selection_config, left_rho),
    CYAML_FIELD_BOOL("recurrent_sweeps_left_of_locus", CYAML_FLAG_OPTIONAL,
        struct selection_config, recurrent_sweeps_left_of_locus),
    CYAML_FIELD_FLOAT_PTR("sweep_time_step",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct selection_config, sweep_time_step),
    CYAML_FIELD_END
};


/* schema for fields within output block */
static const cyaml_schema_field_t output_fields_schema[] = {
    /* required arguments */
    /* optional arguments; output_type defaults to haplotype (ms output) */
    CYAML_FIELD_ENUM("output_type", CYAML_FLAG_OPTIONAL,
        struct output_config, output_type,
        output_types_strings, CYAML_ARRAY_LEN(output_types_strings)),
    CYAML_FIELD_BOOL("hide_partial_snp", CYAML_FLAG_OPTIONAL,
        struct output_config, hide_partial_snp),
    CYAML_FIELD_BOOL("unsimplified_tree_sequence", CYAML_FLAG_OPTIONAL,
        struct output_config, unsimplified_tree_sequence),
    CYAML_FIELD_STRING_PTR("tree_sequence_filename", CYAML_FLAG_OPTIONAL,
        struct output_config, tree_sequence_filename, 0, PATH_MAX),
    CYAML_FIELD_END
};


/* top level schema */
static const cyaml_schema_field_t uniform_prior_fields_schema[] = {
    CYAML_FIELD_FLOAT("low", CYAML_FLAG_DEFAULT, struct uniform_prior, low),
    CYAML_FIELD_FLOAT("high", CYAML_FLAG_DEFAULT, struct uniform_prior, high),
    CYAML_FIELD_END
};
static const cyaml_schema_field_t exponential_prior_fields_schema[] = {
    CYAML_FIELD_FLOAT("mean", CYAML_FLAG_DEFAULT, struct exponential_prior, mean),
    CYAML_FIELD_FLOAT("upper", CYAML_FLAG_DEFAULT, struct exponential_prior, upper),
    CYAML_FIELD_END
};
static const cyaml_schema_field_t prior_size_change_fields_schema[] = {
    CYAML_FIELD_FLOAT("time_low", CYAML_FLAG_DEFAULT,
        struct prior_size_change, time_low),
    CYAML_FIELD_FLOAT("time_high", CYAML_FLAG_DEFAULT,
        struct prior_size_change, time_high),
    CYAML_FIELD_FLOAT("size_low", CYAML_FLAG_DEFAULT,
        struct prior_size_change, size_low),
    CYAML_FIELD_FLOAT("size_high", CYAML_FLAG_DEFAULT,
        struct prior_size_change, size_high),
    CYAML_FIELD_INT("population", CYAML_FLAG_OPTIONAL,
        struct prior_size_change, population),
    CYAML_FIELD_END
};
static const cyaml_schema_value_t prior_size_change_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, struct prior_size_change,
        prior_size_change_fields_schema),
};
#define UNIFORM_PRIOR_FIELD(key, member) \
    CYAML_FIELD_MAPPING_PTR(key, CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL, \
        struct priors_config, member, uniform_prior_fields_schema)
static const cyaml_schema_field_t priors_fields_schema[] = {
    UNIFORM_PRIOR_FIELD("mutation_rate", mutation_rate),
    UNIFORM_PRIOR_FIELD("recombination_rate", recombination_rate),
    CYAML_FIELD_MAPPING_PTR("recombination_rate_exponential",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL, struct priors_config,
        recombination_rate_exponential, exponential_prior_fields_schema),
    UNIFORM_PRIOR_FIELD("selection_coefficient", selection_coefficient),
    UNIFORM_PRIOR_FIELD("fixation_time_ago", fixation_time_ago),
    UNIFORM_PRIOR_FIELD("beneficial_mutation_rate", beneficial_mutation_rate),
    UNIFORM_PRIOR_FIELD("sweep_position", sweep_position),
    UNIFORM_PRIOR_FIELD("initial_frequency", initial_frequency),
    UNIFORM_PRIOR_FIELD("final_frequency", final_frequency),
    CYAML_FIELD_SEQUENCE_COUNT("population_size_changes",
        CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL, struct priors_config,
        population_size_changes, num_population_size_changes,
        &prior_size_change_schema, 1, CYAML_UNLIMITED),
    CYAML_FIELD_END
};
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
    CYAML_FIELD_MAPPING_PTR("priors", CYAML_FLAG_POINTER | CYAML_FLAG_OPTIONAL,
        struct discoal_config, priors, priors_fields_schema),
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


/* libcyaml zeros every primitive struct field before populating it from
 * the YAML, so a user-supplied 0 looks indistinguishable from a missing
 * field. Each block parser below handles three flavours of input:
 *
 *   - Required fields with no semantic "unset" value (sample_size,
 *     num_sites, mutation_rate, sweep_position, ...): cyaml guarantees
 *     the field is set, so validate the range and assign unconditionally.
 *
 *   - Optional fields where 0 is semantically equivalent to "not
 *     configured" (gene_conversion_rate, initial_frequency, ...): keep
 *     the `if (x > 0)` copy guard, but reject explicit out-of-range
 *     values like negatives or >= 1 for frequencies.
 *
 *   - Optional fields where 0 is a meaningful user choice distinct from
 *     "absent" (none today): would need cyaml's pointer-typed optional
 *     machinery so absence -> NULL. Left for follow-up.
 *
 * Booleans are assigned unconditionally; cyaml zeroing matches our
 * default-false convention.
 */

int parse_simulation_block(struct simulation_config *cfg)
{
    extern int sampleSize, sampleNumber, nSites;
    extern int sampleSizes[MAXPOPS];
    extern int popnSizes[MAXPOPS];
    extern int npops;
    extern double *currentSize;
    extern long seed1, seed2;
    if (cfg != NULL) {
        if (cfg->sample_size <= 0) {
            fprintf(stderr,
                "Error parsing config: sample_size (%d) must be > 0\n",
                cfg->sample_size);
            return EXIT_FAILURE;
        }
        sampleSize = cfg->sample_size;
        if (cfg->num_replicates <= 0) {
            fprintf(stderr,
                "Error parsing config: num_replicates (%d) must be > 0\n",
                cfg->num_replicates);
            return EXIT_FAILURE;
        }
        sampleNumber = cfg->num_replicates;
        if (cfg->num_sites <= 0) {
            fprintf(stderr,
                "Error parsing config: num_sites (%d) must be > 0\n",
                cfg->num_sites);
            return EXIT_FAILURE;
        }
        nSites = cfg->num_sites;
        if (cfg->seed != NULL) {
            seed1 = cfg->seed[0];
            seed2 = cfg->seed[1];
        }
        /* Establish the single-population defaults implied by sample_size:
         * one population containing every sample. parse_demography_block
         * overrides these when an explicit demography block is present.
         * Without this, YAMLs that omit the optional demography block left
         * sampleSizes[]/currentSize[]/npops at zero, and the simulation
         * crashed with no sample nodes. */
        sampleSizes[0] = sampleSize;
        popnSizes[0]   = sampleSize;
        currentSize[0] = 1.0;
        npops          = 1;
    }
    return EXIT_SUCCESS;
}

int parse_genetics_block(struct genetics_config *cfg)
{
    extern double theta, rho;
    extern double gammaCoRatio, my_gamma, gammaCoRatioMode;
    extern int gcMean;
    if (cfg != NULL) {
        if (cfg->mutation_rate < 0) {
            fprintf(stderr,
                "Error parsing config: mutation_rate (%g) must be >= 0\n",
                cfg->mutation_rate);
            return EXIT_FAILURE;
        }
        theta = cfg->mutation_rate;
        if (cfg->recombination_rate < 0) {
            fprintf(stderr,
                "Error parsing config: recombination_rate (%g) must be >= 0\n",
                cfg->recombination_rate);
            return EXIT_FAILURE;
        }
        rho = cfg->recombination_rate;
        if (cfg->crossover_ratio < 0) {
            fprintf(stderr,
                "Error parsing config: crossover_ratio (%g) must be >= 0\n",
                cfg->crossover_ratio);
            return EXIT_FAILURE;
        }
        if (cfg->gene_conversion_rate < 0) {
            fprintf(stderr,
                "Error parsing config: gene_conversion_rate (%g) must be >= 0\n",
                cfg->gene_conversion_rate);
            return EXIT_FAILURE;
        }
        if (cfg->gene_conversion_tract_length < 0) {
            fprintf(stderr,
                "Error parsing config: gene_conversion_tract_length (%d) "
                "must be >= 0\n", cfg->gene_conversion_tract_length);
            return EXIT_FAILURE;
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
            my_gamma = cfg->gene_conversion_rate;
        }
        if (cfg->gene_conversion_tract_length > 0) {
            gcMean = cfg->gene_conversion_tract_length;
        }
    }
    return EXIT_SUCCESS;
}

int parse_demography_block(struct demography_config *cfg)
{
    extern int sampleSize;
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
        int deme_sum = 0;
        for (int i = 0; i < cfg->num_demes; ++i) {
            int n = cfg->deme_sample_size[i];
            if (n < 0) {
                fprintf(stderr,
                    "Error parsing config: deme_sample_size[%d] (%d) "
                    "must be >= 0\n", i, n);
                return EXIT_FAILURE;
            }
            sampleSizes[i] = n;
            currentSize[i] = 1.0;
            deme_sum += n;
        }
        /* initialize() creates sum(sampleSizes) sample nodes but sets
         * alleleNumber to sampleSize; a mismatch leaves popLists[] and
         * nodes[] out of sync and later corrupts coalescence. */
        if (deme_sum != sampleSize) {
            fprintf(stderr,
                "Error parsing config: sum of deme_sample_size (%d) does "
                "not match sample_size (%d)\n", deme_sum, sampleSize);
            return EXIT_FAILURE;
        }
        /* effective_population_size is optional; cyaml zero-fills, so 0
         * means "not set". Reject any explicitly invalid value (negative,
         * NaN) before falling through to the "is set" guard. */
        if (cfg->effective_population_size < 0 ||
            isnan(cfg->effective_population_size)) {
            fprintf(stderr,
                "Error parsing config: effective_population_size (%g) "
                "must be >= 0\n", cfg->effective_population_size);
            return EXIT_FAILURE;
        }
        if (cfg->effective_population_size > 0) {
            /* EFFECTIVE_POPN_SIZE is int; cmdline `-N` uses strtol with a
             * range check. Mirror that here so a fractional or out-of-range
             * YAML value is rejected loudly rather than silently truncated. */
            double ne = cfg->effective_population_size;
            if (ne > (double)INT_MAX) {
                fprintf(stderr,
                    "Error parsing config: effective_population_size (%g) "
                    "exceeds INT_MAX (%d)\n", ne, INT_MAX);
                return EXIT_FAILURE;
            }
            if (ne != (double)(long)ne) {
                fprintf(stderr,
                    "Error parsing config: effective_population_size (%g) "
                    "must be an integer\n", ne);
                return EXIT_FAILURE;
            }
            EFFECTIVE_POPN_SIZE = (int)ne;
        }
        /* parse demographic events; these will be sorted into time order later */
        if (cfg->demographic_events != NULL) {
            struct demographic_events *dmo = cfg->demographic_events;
            /* A demes file gives the whole demography, so with one only
             * ancient_samples may be set; the demes importer checks that
             * each sample's deme exists at its sampling time. */
            if (cfg->demes_filename != NULL &&
                (dmo->num_population_size_changes > 0 ||
                 dmo->num_population_splits > 0 ||
                 dmo->num_migration_rate_changes > 0 ||
                 dmo->num_admixtures > 0)) {
                fprintf(stderr, "Error parsing config: with demes_filename, "
                    "demographic_events may contain only ancient_samples; the "
                    "demes file gives the demography\n");
                return EXIT_FAILURE;
            }
            /* Ancient samples, as -A sample_size population time: that many
             * of the population's deme_sample_size samples are taken at
             * `time` instead of at time 0. The count per population is
             * checked against deme_sample_size after parsing, for the CLI
             * and YAML alike. A sample from a population at or after the
             * time it merges away is rejected, as for migration changes. */
            for (int i = 0; i < dmo->num_ancient_samples; ++i) {
                struct ancient_sample *as = &dmo->ancient_samples[i];
                if (check_pop_index("ancient_samples", "population",
                        i, as->population, cfg->num_demes) != EXIT_SUCCESS ||
                    check_nonneg_double("ancient_samples", "time",
                        i, as->time) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (as->sample_size < 1) {
                    fprintf(stderr,
                        "Error parsing config: ancient_samples[%d] "
                        "sample_size (%d) must be >= 1\n",
                        i, as->sample_size);
                    return EXIT_FAILURE;
                }
                if (check_not_merged("ancient_samples", i, as->population,
                        as->time, dmo, -1) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                ensureEventsCapacity();
                events[eventNumber].type = 'A';
                events[eventNumber].time = userTimeToInternal(as->time);
                events[eventNumber].popID = as->population;
                events[eventNumber].lineageNumber = as->sample_size;
                eventNumber++;
            }
            /* Migration rate changes, as -em time source destination rate:
             * from `time` back into the past, lineages in `source` move to
             * `destination` at `rate` (4 N0 m). A change that involves a
             * population at or after the time it merges away (a
             * population_splits entry with it as `derived`) is rejected: the
             * population no longer exists then, and the split already stops
             * its migration. */
            for (int i = 0; i < dmo->num_migration_rate_changes; ++i) {
                struct migration_rate_change *mrc = &dmo->migration_rate_changes[i];
                if (check_pop_index("migration_rate_changes", "source",
                        i, mrc->source, cfg->num_demes) != EXIT_SUCCESS ||
                    check_pop_index("migration_rate_changes", "destination",
                        i, mrc->destination, cfg->num_demes) != EXIT_SUCCESS ||
                    check_nonneg_double("migration_rate_changes", "time",
                        i, mrc->time) != EXIT_SUCCESS ||
                    check_nonneg_double("migration_rate_changes", "rate",
                        i, mrc->rate) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (mrc->source == mrc->destination) {
                    fprintf(stderr,
                        "Error parsing config: migration_rate_changes[%d] "
                        "source and destination (%d) must differ\n",
                        i, mrc->source);
                    return EXIT_FAILURE;
                }
                /* The split already stops a merged population's migration. */
                if (check_not_merged("migration_rate_changes", i, mrc->source,
                        mrc->time, dmo, -1) != EXIT_SUCCESS ||
                    check_not_merged("migration_rate_changes", i,
                        mrc->destination, mrc->time, dmo, -1) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                ensureEventsCapacity();
                events[eventNumber].type = 'm';
                events[eventNumber].time = userTimeToInternal(mrc->time);
                events[eventNumber].popID2 = mrc->source;
                events[eventNumber].popID = mrc->destination;
                events[eventNumber].popnSize = mrc->rate;
                eventNumber++;
                /* Migration that starts later still counts as migration for
                 * initialize()'s check on multi-population models. */
                if (mrc->rate > 0) migFlag = 1;
            }
            /* YAML times use the command-line units; see userTimeToInternal. */
            for (int i = 0; i < dmo->num_population_size_changes; ++i) {
                struct population_size_change *psc = &dmo->population_size_changes[i];
                if (check_pop_index("population_size_changes", "population",
                        i, psc->population, cfg->num_demes) != EXIT_SUCCESS ||
                    check_nonneg_double("population_size_changes", "time",
                        i, psc->time) != EXIT_SUCCESS ||
                    check_positive_double("population_size_changes", "size",
                        i, psc->size) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                ensureEventsCapacity();
                events[eventNumber].type = 'n';
                events[eventNumber].time = userTimeToInternal(psc->time);
                events[eventNumber].popID = psc->population;
                events[eventNumber].popnSize = psc->size;
                eventNumber++;
            }
            for (int i = 0; i < dmo->num_population_splits; ++i) {
                struct population_split *split = &dmo->population_splits[i];
                if (check_pop_index("population_splits", "derived",
                        i, split->derived, cfg->num_demes) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (check_pop_index("population_splits", "ancestral",
                        i, split->ancestral, cfg->num_demes) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (split->derived == split->ancestral) {
                    fprintf(stderr,
                        "Error parsing config: population_splits[%d] derived "
                        "(%d) and ancestral (%d) must be distinct\n",
                        i, split->derived, split->ancestral);
                    return EXIT_FAILURE;
                }
                if (check_positive_double("population_splits", "time",
                        i, split->time) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                ensureEventsCapacity();
                events[eventNumber].type = 'p';
                events[eventNumber].time = userTimeToInternal(split->time);
                events[eventNumber].popID = split->derived;
                events[eventNumber].popID2 = split->ancestral;
                eventNumber++;
                tDiv = split->time;  /* mark merger model active */
                end_population_migration(split->derived, split->time,
                    cfg->num_demes);
            }
            for (int i = 0; i < dmo->num_admixtures; ++i) {
                struct admixture *ad = &dmo->admixtures[i];
                if (check_pop_index("admixtures", "population",
                        i, ad->population, cfg->num_demes) != EXIT_SUCCESS ||
                    check_pop_index("admixtures", "source_1",
                        i, ad->source_1, cfg->num_demes) != EXIT_SUCCESS ||
                    check_pop_index("admixtures", "source_2",
                        i, ad->source_2, cfg->num_demes) != EXIT_SUCCESS ||
                    check_positive_double("admixtures", "time",
                        i, ad->time) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (ad->population == ad->source_1 ||
                    ad->population == ad->source_2 ||
                    ad->source_1 == ad->source_2) {
                    fprintf(stderr,
                        "Error parsing config: admixtures[%d] population (%d), "
                        "source_1 (%d) and source_2 (%d) must be distinct\n",
                        i, ad->population, ad->source_1, ad->source_2);
                    return EXIT_FAILURE;
                }
                if (ad->proportion < 0.0 || ad->proportion > 1.0) {
                    fprintf(stderr,
                        "Error parsing config: admixtures[%d].proportion (%g) "
                        "must be in [0, 1]\n", i, ad->proportion);
                    return EXIT_FAILURE;
                }
                /* The sources and the admixed population must still exist:
                 * not merged away by a split, and not admixed by another
                 * entry (this entry ends the admixed population itself). */
                if (check_not_merged("admixtures", i, ad->source_1, ad->time,
                        dmo, -1) != EXIT_SUCCESS ||
                    check_not_merged("admixtures", i, ad->source_2, ad->time,
                        dmo, -1) != EXIT_SUCCESS ||
                    check_not_merged("admixtures", i, ad->population,
                        ad->time, dmo, i) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                ensureEventsCapacity();
                events[eventNumber].type = 'a';
                events[eventNumber].time = userTimeToInternal(ad->time);
                events[eventNumber].popID = ad->population;
                events[eventNumber].popID2 = ad->source_1;
                events[eventNumber].popID3 = ad->source_2;
                events[eventNumber].admixProp = ad->proportion;
                eventNumber++;
                tDiv = ad->time;  /* populations are joined, as by a split */
                end_population_migration(ad->population, ad->time,
                    cfg->num_demes);
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
                /* Self-migration is not meaningful in a coalescent; cmdline
                 * `-M` forces the diagonal to zero. Reject rather than
                 * silently massage the user's input. */
                for (int j = 0; j < row->num_cols; ++j) {
                    if (j == i) {
                        if (row->rates[j] != 0.0) {
                            fprintf(stderr,
                                "Error parsing config: migration_matrix diagonal "
                                "must be zero (row %d column %d = %g)\n",
                                i, j, row->rates[j]);
                            return EXIT_FAILURE;
                        }
                    } else if (row->rates[j] < 0.0) {
                        fprintf(stderr,
                            "Error parsing config: migration_matrix rate "
                            "(row %d column %d = %g) must be >= 0\n",
                            i, j, row->rates[j]);
                        return EXIT_FAILURE;
                    }
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
            VERBOSE_LOG(
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
        int recurrent = cfg->recurrent_sweep_rate > 0;
        /* left_rho: 0 is a sweep right at the left end of the locus,
         * as -ls tau 0. */
        int linked = cfg->left_rho != NULL;
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
        }
        /* sweep_mode trajectory_file is -S: the file gives the sweep's
         * time and frequencies, so the fields that set these are not
         * allowed. The other combinations (priors, for example) are
         * checked with the command-line options. */
        const int from_file = cfg->sweep_mode == SWEEP_TRAJECTORY_FILE;
        if (from_file != (cfg->trajectory_file != NULL)) {
            fprintf(stderr,
                "Error parsing config: trajectory_file goes with sweep_mode "
                "'trajectory_file', and that mode needs it\n");
            return EXIT_FAILURE;
        }
        if (from_file && (cfg->selection_coefficient != NULL ||
                          cfg->fixation_time_ago != NULL ||
                          cfg->initial_frequency != 0 ||
                          cfg->final_frequency != 0 ||
                          cfg->left_rho != NULL ||
                          cfg->recurrent_sweep_rate != 0 ||
                          cfg->recurrent_sweeps_left_of_locus)) {
            fprintf(stderr,
                "Error parsing config: with sweep_mode 'trajectory_file', "
                "the file gives the sweep's time and frequencies; remove "
                "selection_coefficient, fixation_time_ago, "
                "initial_frequency, final_frequency, left_rho and the "
                "recurrent sweep fields\n");
            return EXIT_FAILURE;
        }
        if (!from_file) {
            if (cfg->selection_coefficient == NULL) {
                fprintf(stderr,
                    "Error parsing config: selection block requires "
                    "selection_coefficient\n");
                return EXIT_FAILURE;
            }
            if (*cfg->selection_coefficient < 0) {
                fprintf(stderr,
                    "Error parsing config: selection_coefficient (%g) "
                    "must be >= 0\n", *cfg->selection_coefficient);
                return EXIT_FAILURE;
            }
            if (cfg->sweep_mode == SWEEP_NEUTRAL && *cfg->selection_coefficient > 0) {
                fprintf(stderr,
                    "Error parsing config: selection_coefficient (%g) must be 0 "
                    "when sweep_mode is 'neutral'\n", *cfg->selection_coefficient);
                return EXIT_FAILURE;
            }
            alpha = *cfg->selection_coefficient;
        }
        if (cfg->sweep_position != NULL) {
            if (*cfg->sweep_position < 0.0 || *cfg->sweep_position > 1.0) {
                fprintf(stderr,
                    "Error parsing config: sweep_position (%g) must be in "
                    "[0, 1]\n", *cfg->sweep_position);
                return EXIT_FAILURE;
            }
            sweepSite = *cfg->sweep_position;
        }
        if (cfg->fixation_time_ago != NULL) {
            if (*cfg->fixation_time_ago < 0) {
                fprintf(stderr,
                    "Error parsing config: fixation_time_ago (%g) must be "
                    ">= 0\n", *cfg->fixation_time_ago);
                return EXIT_FAILURE;
            }
            tau = userTimeToInternal(*cfg->fixation_time_ago);
        }
        if (cfg->initial_frequency < 0 || cfg->initial_frequency >= 1.0) {
            fprintf(stderr,
                "Error parsing config: initial_frequency (%g) "
                "must be in (0, 1)\n", cfg->initial_frequency);
            return EXIT_FAILURE;
        }
        if (cfg->initial_frequency > 0) {
            f0 = cfg->initial_frequency;
            softSweepMode = 1;
        }
        if (cfg->final_frequency < 0 || cfg->final_frequency >= 1.0) {
            fprintf(stderr,
                "Error parsing config: final_frequency (%g) "
                "must be in (0, 1)\n", cfg->final_frequency);
            return EXIT_FAILURE;
        }
        if (cfg->final_frequency > 0) {
            partialSweepFinalFreq = cfg->final_frequency;
            partialSweepMode = 1;
        }
        /* Fixed values only: with a prior on final_frequency (parsed
         * before this block), final_frequency here is 0; getParameters()
         * checks the priors' bounds. */
        if (softSweepMode && cfg->final_frequency > 0 &&
            cfg->initial_frequency >= cfg->final_frequency) {
            fprintf(stderr,
                "Error parsing config: initial_frequency (%g) must be less "
                "than final_frequency (%g) when both are set\n",
                cfg->initial_frequency, cfg->final_frequency);
            return EXIT_FAILURE;
        }
        if (cfg->beneficial_mutation_rate < 0) {
            fprintf(stderr,
                "Error parsing config: beneficial_mutation_rate (%g) "
                "must be >= 0\n", cfg->beneficial_mutation_rate);
            return EXIT_FAILURE;
        }
        if (cfg->beneficial_mutation_rate > 0) {
            uA = cfg->beneficial_mutation_rate;
        }
        if (cfg->recurrent_sweep_rate < 0) {
            fprintf(stderr,
                "Error parsing config: recurrent_sweep_rate (%g) "
                "must be >= 0\n", cfg->recurrent_sweep_rate);
            return EXIT_FAILURE;
        }
        if (linked && *cfg->left_rho < 0) {
            fprintf(stderr,
                "Error parsing config: left_rho (%g) must be >= 0\n",
                *cfg->left_rho);
            return EXIT_FAILURE;
        }
        if (cfg->sweep_time_step != NULL && *cfg->sweep_time_step <= 0) {
            fprintf(stderr,
                "Error parsing config: sweep_time_step (%g) must be > 0\n",
                *cfg->sweep_time_step);
            return EXIT_FAILURE;
        }
        if (cfg->sweep_time_step != NULL) {
            deltaTMod = *cfg->sweep_time_step;
        }
        /* A sweep outside the locus has no position inside it: -ls, -ld,
         * -ln and -L set the sweep site to -1. */
        if ((linked || cfg->recurrent_sweeps_left_of_locus) &&
            cfg->sweep_position != NULL) {
            fprintf(stderr,
                "Error parsing config: sweep_position cannot be set with "
                "left_rho or recurrent_sweeps_left_of_locus; the sweep is "
                "outside the locus\n");
            return EXIT_FAILURE;
        }
        if (cfg->recurrent_sweeps_left_of_locus && !recurrent) {
            fprintf(stderr,
                "Error parsing config: recurrent_sweeps_left_of_locus needs "
                "recurrent_sweep_rate > 0\n");
            return EXIT_FAILURE;
        }

        if (recurrent) {
            /* -R and -L: recurrent sweeps, which are always stochastic. */
            if (cfg->fixation_time_ago != NULL) {
                fprintf(stderr,
                    "Error parsing config: fixation_time_ago and "
                    "recurrent_sweep_rate cannot both be set\n");
                return EXIT_FAILURE;
            }
            if (linked) {
                fprintf(stderr,
                    "Error parsing config: left_rho is for a single sweep; "
                    "with recurrent_sweep_rate use "
                    "recurrent_sweeps_left_of_locus\n");
                return EXIT_FAILURE;
            }
            if (cfg->sweep_mode != SWEEP_STOCHASTIC) {
                fprintf(stderr,
                    "Error parsing config: recurrent sweeps are stochastic; "
                    "set sweep_mode to 'stochastic' with "
                    "recurrent_sweep_rate\n");
                return EXIT_FAILURE;
            }
            recurSweepRate = cfg->recurrent_sweep_rate;
            recurSweepMode = 1;
            if (cfg->recurrent_sweeps_left_of_locus) {
                sweepSite = -1.0;
                recurrentLeftOfLocus = 1;
            }
            return EXIT_SUCCESS;
        }

        if (from_file) {
            return setSweepTrajectoryFile(cfg->trajectory_file) == 0
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        /* A single sweep (-ws/-wd/-wn, or -ls/-ld/-ln with left_rho) is an
         * 's' event at tau. With a prior on fixation_time_ago, the event's
         * time is drawn for each replicate; until then it is the lower
         * bound, set by the priors block. */
        if (cfg->fixation_time_ago == NULL && !priorTau) {
            fprintf(stderr,
                "Error parsing config: selection block requires "
                "fixation_time_ago (or a prior on it) unless "
                "recurrent_sweep_rate is set\n");
            return EXIT_FAILURE;
        }
        if (cfg->fixation_time_ago == NULL) {
            /* Only a prior: the event holds the lower bound until the
             * first draw. */
            tau = pTauLow;
        }
        if (linked) {
            sweepSite = -1.0;
            leftRho = *cfg->left_rho * 2.0;  /* stored as -ls stores it */
            leftRhoFlag = 1;
        }
        ensureEventsCapacity();
        events[eventNumber].time = tau;
        events[eventNumber].type = 's';
        eventNumber++;
    }
    return EXIT_SUCCESS;
}

int parse_output_block(struct output_config *cfg) 
{
    extern int tskitOutputMode, minimalTreeSeq, hidePartialSNP;
    extern char tskitOutputFilename[1024];
    if (cfg != NULL) {
        switch (cfg->output_type) {
            case OUTPUT_DEFAULT:
            case OUTPUT_HAPLOTYPE:
                /* ms-format haplotypes on stdout, discoal's normal output */
                break;
            case OUTPUT_SNP_ARRAY:
                /* Kept in the value table only to give a clear error. */
                fprintf(stderr,
                    "Error parsing config: output_type \"snp_array\" is not "
                    "supported; use \"haplotype\" (ms output, the default) "
                    "or \"tree_sequence\"\n");
                return EXIT_FAILURE;
            case OUTPUT_TREE_SEQN:
                if (cfg->tree_sequence_filename == NULL ||
                    cfg->tree_sequence_filename[0] == '\0') {
                    fprintf(stderr,
                        "Must provide non-empty tree_sequence_filename if "
                        "using output mode tree_sequence\n");
                    return EXIT_FAILURE;
                }
                tskitOutputMode = 1;
                break;
            default:
                break;
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

/* One uniform prior of the priors block and the globals its -P flag sets.
 * The bounds must satisfy lo <= low <= high <= hi, with strict
 * inequalities where open_lo / open_hi is set. */
/* The range of each uniform prior: its YAML field, its command-line flag,
 * and the globals it sets. low <= high must lie inside [lo, hi], with lo or
 * hi excluded if open_lo or open_hi. Both the YAML priors and the -P flags
 * are checked against this table. */
struct prior_range {
    const char *name, *flag;
    int needs_selection;
    double lo, hi;
    int open_lo, open_hi;
    int is_time;                         /* convert to internal units */
    int *on;
    double *low, *high;
};
static const struct prior_range prior_ranges[] = {
    {"mutation_rate", "-Pt", 0, 0, INFINITY, 0, 0, 0,
        &priorTheta, &pThetaLow, &pThetaUp},
    {"recombination_rate", "-Pr", 0, 0, INFINITY, 0, 0, 0,
        &priorRho, &pRhoLow, &pRhoUp},
    {"selection_coefficient", "-Pa", 1, 0, INFINITY, 0, 0, 0,
        &priorAlpha, &pAlphaLow, &pAlphaUp},
    {"fixation_time_ago", "-Pu", 1, 0, INFINITY, 0, 0, 1,
        &priorTau, &pTauLow, &pTauUp},
    {"beneficial_mutation_rate", "-PuA", 1, 0, INFINITY, 0, 0, 0,
        &priorUA, &pUALow, &pUAUp},
    {"sweep_position", "-Px", 1, 0, 1, 0, 0, 0,
        &priorX, &pXLow, &pXUp},
    {"initial_frequency", "-Pf", 1, 0, 1, 0, 1, 0,
        &priorF0, &pF0Low, &pF0Up},
    {"final_frequency", "-Pc", 1, 0, 1, 0, 1, 0,
        &priorC, &pCLow, &pCUp},
};
#define N_PRIOR_RANGES (sizeof(prior_ranges) / sizeof(prior_ranges[0]))

static int uniformBoundsOk(const struct prior_range *r, double low, double high) {
    int lo_ok = r->open_lo ? low > r->lo : low >= r->lo;
    int hi_ok = r->open_hi ? high < r->hi : high <= r->hi;
    return lo_ok && hi_ok && low <= high;
}

/* "lo <= low <= high <= hi", with < for an open end and no upper part for
 * an infinite hi. */
static void describeRange(const struct prior_range *r, char *out, size_t n) {
    int k = snprintf(out, n, "%g %s low <= high", r->lo, r->open_lo ? "<" : "<=");
    if (isfinite(r->hi) && k > 0 && (size_t)k < n)
        snprintf(out + k, n - k, " %s %g", r->open_hi ? "<" : "<=", r->hi);
}

static int apply_uniform_prior(const struct prior_range *u,
                               const struct uniform_prior *p) {
    if (!uniformBoundsOk(u, p->low, p->high)) {
        char range[64];
        describeRange(u, range, sizeof(range));
        fprintf(stderr,
            "Error parsing config: priors.%s must satisfy %s (got low %g, "
            "high %g)\n", u->name, range, p->low, p->high);
        return EXIT_FAILURE;
    }
    *u->on = 1;
    *u->low = u->is_time ? userTimeToInternal(p->low) : p->low;
    *u->high = u->is_time ? userTimeToInternal(p->high) : p->high;
    return EXIT_SUCCESS;
}

/* Check every uniform prior that is set, by a -P flag or by YAML, against
 * prior_ranges. Times are checked in 4N0 generations, as given. Returns 0,
 * or -1 after printing an error that names the flag. (-Pre sets priorRho to
 * 2 and is checked by the caller.) */
int checkUniformPriorBounds(void) {
    for (size_t k = 0; k < N_PRIOR_RANGES; ++k) {
        const struct prior_range *r = &prior_ranges[k];
        if (*r->on != 1) continue;
        double low = *r->low, high = *r->high;
        if (r->is_time) {
            low /= USER_TO_INTERNAL_TIME;
            high /= USER_TO_INTERNAL_TIME;
        }
        if (!uniformBoundsOk(r, low, high)) {
            char range[64];
            describeRange(r, range, sizeof(range));
            fprintf(stderr, "Error: %s %g %g: the bounds must satisfy %s\n",
                    r->flag, low, high, range);
            return -1;
        }
    }
    return 0;
}

/* The priors block, as the -P flags: each entry sets the same globals as
 * its flag, and drawPriors() draws the values for every replicate. Parsed
 * before the selection block, which emits the sweep event when
 * fixation_time_ago has only a prior. */
int parse_priors_block(struct priors_config *cfg, int has_selection,
                       unsigned num_demes, const struct demographic_events *dmo)
{
    if (cfg == NULL) {
        return EXIT_SUCCESS;
    }
    /* In the order of prior_ranges. */
    const struct uniform_prior *given[] = {
        cfg->mutation_rate, cfg->recombination_rate, cfg->selection_coefficient,
        cfg->fixation_time_ago, cfg->beneficial_mutation_rate,
        cfg->sweep_position, cfg->initial_frequency, cfg->final_frequency,
    };
    _Static_assert(sizeof(given) / sizeof(given[0]) == N_PRIOR_RANGES,
                   "given[] must match prior_ranges[]");
    for (size_t k = 0; k < N_PRIOR_RANGES; ++k) {
        if (given[k] == NULL) continue;
        if (prior_ranges[k].needs_selection && !has_selection) {
            fprintf(stderr,
                "Error parsing config: priors.%s needs a selection block\n",
                prior_ranges[k].name);
            return EXIT_FAILURE;
        }
        if (apply_uniform_prior(&prior_ranges[k], given[k]) != EXIT_SUCCESS)
            return EXIT_FAILURE;
    }
    /* A prior on the final frequency makes the sweep partial, as -Pc. */
    if (cfg->final_frequency) {
        partialSweepMode = 1;
    }
    if (cfg->recombination_rate_exponential) {
        struct exponential_prior *e = cfg->recombination_rate_exponential;
        if (cfg->recombination_rate) {
            fprintf(stderr,
                "Error parsing config: priors.recombination_rate and "
                "priors.recombination_rate_exponential cannot both be set\n");
            return EXIT_FAILURE;
        }
        if (!(e->mean > 0 && e->upper > 0)) {
            fprintf(stderr,
                "Error parsing config: priors.recombination_rate_exponential "
                "mean (%g) and upper (%g) must be > 0\n", e->mean, e->upper);
            return EXIT_FAILURE;
        }
        priorRho = 2;
        pRhoMean = e->mean;
        pRhoUp = e->upper;
    }
    for (unsigned i = 0; i < cfg->num_population_size_changes; ++i) {
        struct prior_size_change *p = &cfg->population_size_changes[i];
        const char *what = "priors.population_size_changes";
        if (check_pop_index(what, "population", (int)i, p->population,
                num_demes) != EXIT_SUCCESS)
            return EXIT_FAILURE;
        /* The drawn time is not known yet; reject a population that has
         * already ended at the lower bound. */
        if (dmo != NULL &&
            check_not_merged(what, (int)i, p->population, p->time_low, dmo,
                -1) != EXIT_SUCCESS)
            return EXIT_FAILURE;
        if (addPriorSizeChange(what, userTimeToInternal(p->time_low),
                userTimeToInternal(p->time_high), p->size_low, p->size_high,
                p->population) != 0)
            return EXIT_FAILURE;
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
    ret = parse_priors_block(config->priors, config->selection != NULL,
        config->demography ? config->demography->num_demes : 1,
        config->demography ? config->demography->demographic_events : NULL);
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
        cyaml_free(&cyaml_config, &discoal_config_schema, *config, 0);
        return EXIT_FAILURE;
    }
    /* cyaml returns OK with a NULL top-level when the file is empty or
     * contains only comments; downstream parsers assume non-NULL. */
    if (*config == NULL) {
        fprintf(stderr, "Error parsing config: YAML file '%s' is empty or "
            "contains no top-level mapping\n", yaml_path);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

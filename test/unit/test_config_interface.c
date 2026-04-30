/*
 * test_config_interface.c
 *
 * Unit tests for src/core/configInterface.[ch], the YAML config
 * parser invoked by discoal's -Y option. These tests check the
 * parse_* validation functions and apply_yaml_config; they do not
 * try to test libcyaml's own schema parser (missing required fields,
 * wrong enum strings, etc. are covered by the integration tests in
 * testing/yaml_validation_suite.sh).
 *
 * How apply_yaml_config disposes of its argument
 * ----------------------------------------------
 * apply_yaml_config calls cyaml_free on the struct discoal_config *
 * it receives, even on failure. After apply_yaml_config returns, that
 * pointer is no longer valid. Each test reads the discoal globals the
 * parsers write into (theta, sampleSize, events[], etc.) rather than
 * the parsed config struct.
 *
 * Resetting global state between tests
 * ------------------------------------
 * The parse_* functions write into roughly 25 globals declared across
 * discoal.h and test_globals.c. Each test starts from a known state
 * by calling reset_config_globals() (defined in test_globals.c) from
 * setUp(). If a parse_* function gains a new global write, add it to
 * reset_config_globals() too; otherwise the next test will see the
 * previous test's value.
 *
 * Working directory
 * -----------------
 * setUp() creates a fresh directory under /tmp and chdir's into it.
 * tearDown() chdir's back to the original directory and removes the
 * temporary one. Tests that need an external file (such as the
 * demes-import test) must either use an absolute path or copy the
 * file into the working directory; a relative path that points back
 * into the source tree will not resolve.
 */

/* TEST_ASSERT_EQUAL_DOUBLE requires UNITY_INCLUDE_DOUBLE to be defined
 * for both this translation unit and unity.c; the Makefile target adds
 * -DUNITY_INCLUDE_DOUBLE so both see it. */
#include "unity.h"
#include "configInterface.h"
#include "discoal.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* declared in test_globals.c */
extern void reset_config_globals(void);

/* discoal globals the tests inspect after apply_yaml_config */
extern int sampleSize, sampleNumber, nSites;
extern int sampleSizes[MAXPOPS];
extern int popnSizes[MAXPOPS];
extern int npops;
extern int migFlag;
extern int EFFECTIVE_POPN_SIZE;
extern int gcMean, finiteOutputFlag;
extern int recurSweepMode, partialSweepMode, softSweepMode;
extern int tskitOutputMode, minimalTreeSeq, hidePartialSNP;
extern double theta, rho;
extern double gammaCoRatio, my_gamma, gammaCoRatioMode;
extern double alpha, sweepSite, tau, f0, uA;
extern double partialSweepFinalFreq, recurSweepRate;
extern double tDiv;
extern double migMatConst[MAXPOPS][MAXPOPS];
extern char sweepMode;
extern char tskitOutputFilename[1024];
extern long seed1, seed2;
extern double *currentSize;
extern struct event *events;
extern int eventNumber, eventsCapacity;

/* per-test fixture state */
static char workdir_path[PATH_MAX];
static char yaml_path[PATH_MAX];
static char original_cwd[PATH_MAX];

void setUp(void) {
    if (getcwd(original_cwd, sizeof(original_cwd)) == NULL) {
        TEST_FAIL_MESSAGE("setUp: getcwd failed");
    }
    strcpy(workdir_path, "/tmp/test_config_XXXXXX");
    if (mkdtemp(workdir_path) == NULL) {
        TEST_FAIL_MESSAGE("setUp: mkdtemp failed");
    }
    if (chdir(workdir_path) != 0) {
        TEST_FAIL_MESSAGE("setUp: chdir to workdir failed");
    }
    snprintf(yaml_path, sizeof(yaml_path), "%s/test.yaml", workdir_path);

    if (currentSize != NULL) {
        free(currentSize);
        currentSize = NULL;
    }
    currentSize = (double *)calloc(MAXPOPS, sizeof(double));
    TEST_ASSERT_NOT_NULL_MESSAGE(currentSize, "setUp: currentSize calloc failed");

    if (events != NULL) {
        free(events);
        events = NULL;
    }
    eventNumber = 0;
    eventsCapacity = 0;

    reset_config_globals();
}

void tearDown(void) {
    unlink(yaml_path);
    if (chdir(original_cwd) != 0) {
        /* nothing useful we can do here; test already over */
    }
    rmdir(workdir_path);

    if (currentSize != NULL) {
        free(currentSize);
        currentSize = NULL;
    }
    if (events != NULL) {
        free(events);
        events = NULL;
    }
    eventNumber = 0;
    eventsCapacity = 0;
}

/*
 * Write `yaml` to the per-test temporary file, then run
 * load_yaml_config followed by apply_yaml_config. Returns
 * EXIT_SUCCESS only if both calls succeed; otherwise returns the
 * status from the first one that failed. The parsed config is
 * freed inside apply_yaml_config (see the file-level comment), so
 * the caller checks the discoal globals to see what was set.
 */
static int apply_from_yaml_string(const char *yaml) {
    FILE *fp = fopen(yaml_path, "w");
    if (fp == NULL) {
        return EXIT_FAILURE;
    }
    fputs(yaml, fp);
    fclose(fp);

    struct discoal_config *cfg = NULL;
    int rc = load_yaml_config(yaml_path, &cfg);
    if (rc != EXIT_SUCCESS) {
        return rc;
    }
    return apply_yaml_config(cfg);
}

/*
 * Copy a file from <repo_root>/test/unit/fixtures/<basename> into
 * the per-test working directory under the same basename. This
 * lets a relative demes_filename: in the YAML resolve through the
 * current working directory, which setUp has already chdir'd into.
 *
 * original_cwd is captured by setUp and holds the directory the
 * test binary was launched from. The Makefile's run_tests target
 * runs ./build/test_config_interface from the repo root, so
 * <original_cwd>/test/unit/fixtures/... is the correct source
 * path both there and when the binary is run directly.
 */
static void copy_fixture_to_workdir(const char *fixture_basename) {
    char src_path[PATH_MAX];
    char dst_path[PATH_MAX];
    int written;
    written = snprintf(src_path, sizeof(src_path),
        "%s/test/unit/fixtures/%s", original_cwd, fixture_basename);
    TEST_ASSERT_TRUE_MESSAGE(written > 0 && (size_t)written < sizeof(src_path),
        "fixture src path overflow");
    written = snprintf(dst_path, sizeof(dst_path),
        "%s/%s", workdir_path, fixture_basename);
    TEST_ASSERT_TRUE_MESSAGE(written > 0 && (size_t)written < sizeof(dst_path),
        "fixture dst path overflow");

    FILE *src = fopen(src_path, "rb");
    TEST_ASSERT_NOT_NULL_MESSAGE(src, "fixture source open failed");
    FILE *dst = fopen(dst_path, "wb");
    TEST_ASSERT_NOT_NULL_MESSAGE(dst, "fixture dest open failed");

    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), src)) > 0) {
        TEST_ASSERT_EQUAL_MESSAGE(n, fwrite(buf, 1, n, dst),
            "fixture copy write short");
    }
    fclose(src);
    fclose(dst);
}

/*
 * The simplest valid input: a YAML with only the two required blocks
 * (simulation and genetics) is loaded and applied without error.
 * Useful as a quick sanity check; if this fails, the problem is most
 * likely in setUp / tearDown rather than in any specific validation.
 */
void test_smoke_minimum_viable_config(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(10, sampleSize);
    TEST_ASSERT_EQUAL_INT(1, sampleNumber);
    TEST_ASSERT_EQUAL_INT(1000, nSites);
}

/* ----- Simulation block ----- */

/*
 * A simulation block with explicit sample_size, num_replicates, and
 * num_sites populates the corresponding globals.
 */
void test_load_simulation_section(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 20\n"
        "  num_replicates: 10\n"
        "  num_sites: 100000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(20, sampleSize);
    TEST_ASSERT_EQUAL_INT(10, sampleNumber);
    TEST_ASSERT_EQUAL_INT(100000, nSites);
}

/*
 * A seed: [seed1, seed2] sequence populates the seed1 and seed2
 * globals.
 */
void test_load_seed_array(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "  seed: [12345, 67890]\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(12345, seed1);
    TEST_ASSERT_EQUAL_INT(67890, seed2);
}

/*
 * parse_simulation_block requires sample_size > 0; both 0 and a
 * negative value cause apply_yaml_config to fail. The two values
 * are exercised in separate functions so that each case runs in a
 * fresh setUp / tearDown cycle, rather than resetting globals by
 * hand mid-test.
 */
void test_simulation_rejects_zero_sample_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 0\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

void test_simulation_rejects_negative_sample_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: -1\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_simulation_block requires num_replicates > 0.
 */
void test_simulation_rejects_nonpositive_num_replicates(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 0\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_simulation_block requires num_sites > 0.
 */
void test_simulation_rejects_nonpositive_num_sites(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 0\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * When no demography block is present, parse_simulation_block sets
 * up a single population containing every sample and currentSize[0]
 * = 1.0. This is what downstream code in initialize() assumes; see
 * the comment in configInterface.c near "without this, YAMLs
 * that omit the optional demography block left ... at zero".
 */
void test_simulation_sets_single_pop_defaults(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 7\n"
        "  num_replicates: 1\n"
        "  num_sites: 100\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, npops);
    TEST_ASSERT_EQUAL_INT(7, sampleSizes[0]);
    TEST_ASSERT_EQUAL_INT(7, popnSizes[0]);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, currentSize[0]);
}

/* ----- Genetics block ----- */

/*
 * A genetics block with mutation_rate and recombination_rate
 * populates theta and rho.
 */
void test_load_genetics_section(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(10.0, theta);
    TEST_ASSERT_EQUAL_DOUBLE(5.0, rho);
}

/*
 * gene_conversion_rate and gene_conversion_tract_length populate
 * my_gamma and gcMean respectively.
 */
void test_load_gene_conversion(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  gene_conversion_rate: 2.5\n"
        "  gene_conversion_tract_length: 500\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(2.5, my_gamma);
    TEST_ASSERT_EQUAL_INT(500, gcMean);
}

/*
 * Setting genetics.crossover_ratio activates ratio mode:
 * gammaCoRatioMode is set to 1 and gammaCoRatio receives the value.
 * gammaCoRatioMode is declared as a double in discoal.h, so it is
 * compared as one.
 */
void test_load_gene_conversion_crossover_ratio(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  crossover_ratio: 1.5\n"
        "  gene_conversion_tract_length: 500\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(1.5, gammaCoRatio);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, gammaCoRatioMode);
    TEST_ASSERT_EQUAL_INT(500, gcMean);
}

/*
 * parse_genetics_block rejects a negative mutation_rate. Zero is
 * allowed (it is the documented "no mutation" case); only negatives
 * are an error.
 */
void test_genetics_rejects_negative_mutation_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: -1.0\n"
        "  recombination_rate: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_genetics_block rejects a negative recombination_rate.
 */
void test_genetics_rejects_negative_recombination_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: -1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_genetics_block rejects a negative gene_conversion_rate.
 */
void test_genetics_rejects_negative_gene_conversion_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "  gene_conversion_rate: -1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_genetics_block rejects a negative crossover_ratio.
 */
void test_genetics_rejects_negative_crossover_ratio(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "  crossover_ratio: -1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * parse_genetics_block rejects a negative gene_conversion_tract_length.
 */
void test_genetics_rejects_negative_tract_length(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "  gene_conversion_tract_length: -1\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * gene_conversion_rate and crossover_ratio cannot both be set; the
 * two are different parameterisations of the same quantity.
 */
void test_genetics_rejects_rate_and_crossover_ratio_both_set(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  gene_conversion_rate: 2.5\n"
        "  crossover_ratio: 1.5\n"
        "  gene_conversion_tract_length: 500\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/* ----- Demography block ----- */

/*
 * A demography block with an explicit deme_sample_size sequence sets
 * npops, populates sampleSizes, and the per-deme values sum to
 * sample_size.
 */
void test_load_populations_section(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 23\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10, 5, 8]\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(3, npops);
    TEST_ASSERT_EQUAL_INT(10, sampleSizes[0]);
    TEST_ASSERT_EQUAL_INT(5, sampleSizes[1]);
    TEST_ASSERT_EQUAL_INT(8, sampleSizes[2]);
    TEST_ASSERT_EQUAL_INT(23, sampleSize);
}

/*
 * A population_size_change event is recorded with type 'n', popID,
 * size, and time scaled by 2.0. The user supplies time in 2N units;
 * the parser scales to discoal's internal 4N units, matching the
 * convention used by the -en CLI flag.
 */
void test_load_population_size_change_event(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.05\n"
        "        population: 0\n"
        "        size: 0.5\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, eventNumber);
    TEST_ASSERT_EQUAL_INT('n', events[0].type);
    TEST_ASSERT_EQUAL_DOUBLE(0.10, events[0].time);  /* 0.05 * 2.0 */
    TEST_ASSERT_EQUAL_INT(0, events[0].popID);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, events[0].popnSize);
}

/*
 * A population_split event is recorded with type 'p', popID =
 * derived, popID2 = ancestral, and time scaled by 2.0.  tDiv is set
 * to the unscaled split time; parse_demography_block uses it as a
 * flag for downstream code to switch on the merger model.
 */
void test_load_population_split_event(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    population_splits:\n"
        "      - time: 0.2\n"
        "        derived: 1\n"
        "        ancestral: 0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, eventNumber);
    TEST_ASSERT_EQUAL_INT('p', events[0].type);
    TEST_ASSERT_EQUAL_DOUBLE(0.4, events[0].time);  /* 0.2 * 2.0 */
    TEST_ASSERT_EQUAL_INT(1, events[0].popID);
    TEST_ASSERT_EQUAL_INT(0, events[0].popID2);
    TEST_ASSERT_EQUAL_DOUBLE(0.2, tDiv);
}

/*
 * A demography.demes_filename pointing at the minimal_two_demes.yaml
 * file under test/unit/fixtures loads through demesInterface,
 * populates events, and sets npops to match the file (2). The file
 * is copied into the per-test working directory so the relative
 * path resolves correctly.
 */
void test_demography_demes_filename_loads_events(void) {
    copy_fixture_to_workdir("minimal_two_demes.yaml");

    const char *yaml =
        "simulation:\n"
        "  sample_size: 4\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [2, 2]\n"
        "  demes_filename: \"minimal_two_demes.yaml\"\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(2, npops);
    /* sampleSizes are written by parse_demography_block before
     * loadDemesFile is called and must survive the demes load. This
     * is the regression covered by kr-colab/discoal#75. */
    TEST_ASSERT_EQUAL_INT(2, sampleSizes[0]);
    TEST_ASSERT_EQUAL_INT(2, sampleSizes[1]);
    /* The demes importer emits at least one event for the migration
     * band between A and B; the test does not check the exact count
     * so it does not become coupled to demes-c implementation
     * details. */
    TEST_ASSERT_GREATER_THAN_INT(0, eventNumber);
}

/*
 * parse_demography_block rejects negative deme_sample_size values
 * inside the per-deme loop, before the sample-size sum check.
 */
void test_demography_rejects_negative_deme_sample_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 9\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10, -1]\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * The sum of deme_sample_size must match simulation.sample_size,
 * otherwise initialize() creates the wrong number of sample nodes.
 */
void test_demography_rejects_deme_sum_mismatch(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 20\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10, 5]\n";  /* sum 15 != sample_size 20 */
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * effective_population_size must be >= 0; negatives are rejected.
 */
void test_demography_rejects_negative_effective_population_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  effective_population_size: -100.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * effective_population_size must convert to int without loss.
 * Fractional values are rejected, mirroring the strtol check used
 * by the -N CLI flag.
 */
void test_demography_rejects_nonintegral_effective_population_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  effective_population_size: 1000.5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * effective_population_size > INT_MAX is rejected. The INT_MAX
 * check runs before the integer-conversion check, so a value that
 * is integer-valued but too large is caught here.
 */
void test_demography_rejects_effective_population_size_over_intmax(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  effective_population_size: 1.0e10\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * The recorded event time for a population_size_change is exactly 2
 * times the YAML time, with no additional rounding. Uses a value
 * with no obvious round-number factor so a missing scale would
 * produce an obviously-wrong recorded time.
 */
void test_demography_population_size_change_scales_time_by_2(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.123\n"
        "        population: 0\n"
        "        size: 1.5\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, eventNumber);
    TEST_ASSERT_EQUAL_DOUBLE(0.246, events[0].time);  /* 0.123 * 2.0 */
}

/*
 * population_size_changes[i].population must be in [0, num_demes).
 */
void test_demography_population_size_change_rejects_bad_pop_index(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.05\n"
        "        population: 5\n"  /* num_demes is 2 */
        "        size: 0.5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * population_size_changes[i].time must be >= 0.
 */
void test_demography_population_size_change_rejects_negative_time(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: -0.1\n"
        "        population: 0\n"
        "        size: 0.5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * population_size_changes[i].size must be > 0 (a zero-size population
 * is degenerate; check_positive_double rejects).
 */
void test_demography_population_size_change_rejects_nonpositive_size(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.05\n"
        "        population: 0\n"
        "        size: 0.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * population_splits[i].derived and .ancestral must be distinct;
 * a deme cannot split from itself.
 */
void test_demography_population_split_rejects_same_derived_ancestral(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    population_splits:\n"
        "      - time: 0.2\n"
        "        derived: 0\n"
        "        ancestral: 0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * population_splits[i].derived (and .ancestral) must be in
 * [0, num_demes); out-of-range values are rejected.
 */
void test_demography_population_split_rejects_bad_pop_index(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    population_splits:\n"
        "      - time: 0.2\n"
        "        derived: 5\n"  /* num_demes is 2 */
        "        ancestral: 0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * population_splits[i].time must be > 0 (check_positive_double).
 */
void test_demography_population_split_rejects_nonpositive_time(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    population_splits:\n"
        "      - time: 0.0\n"
        "        derived: 1\n"
        "        ancestral: 0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * migration_rate_changes are parsed by libcyaml but rejected by
 * parse_demography_block, because the main event loop in
 * discoalFunctions.c does not yet implement the handler. This test
 * exists so the rejection is exercised; if someone later wires up
 * the handler, this test will start failing as a reminder to
 * remove the explicit error path.
 */
void test_demography_rejects_migration_rate_changes_not_implemented(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  demographic_events:\n"
        "    migration_rate_changes:\n"
        "      - time: 0.1\n"
        "        rate: 1.0\n"
        "        source: 0\n"
        "        destination: 1\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * ancient_samples are parsed by libcyaml but rejected by
 * parse_demography_block; the handler is not yet implemented.
 */
void test_demography_rejects_ancient_samples_not_implemented(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [10]\n"
        "  demographic_events:\n"
        "    ancient_samples:\n"
        "      - time: 0.1\n"
        "        population: 0\n"
        "        sample_size: 5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * migration_matrix must have exactly num_demes rows.
 */
void test_demography_migration_matrix_requires_square(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [3, 3, 4]\n"  /* 3 demes */
        "  migration_matrix:\n"
        "    - row: [0, 1, 0]\n"
        "    - row: [1, 0, 0]\n";  /* only 2 rows */
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * Each migration_matrix row must have exactly num_demes columns.
 */
void test_demography_migration_matrix_requires_square_rows(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"  /* 2 demes */
        "  migration_matrix:\n"
        "    - row: [0, 1, 1]\n"  /* 3 columns */
        "    - row: [1, 0, 0]\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * migration_matrix diagonal entries must be exactly zero. The -M
 * CLI flag forces the diagonal to zero silently; the YAML parser
 * rejects rather than silently massage the user's input.
 */
void test_demography_migration_matrix_rejects_nonzero_diagonal(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  migration_matrix:\n"
        "    - row: [0.5, 1.0]\n"  /* diag != 0 */
        "    - row: [1.0, 0.0]\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * migration_matrix off-diagonal rates must be >= 0.
 */
void test_demography_migration_matrix_rejects_negative_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  migration_matrix:\n"
        "    - row: [0.0, -1.0]\n"
        "    - row: [1.0,  0.0]\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * A valid migration_matrix populates migMatConst and sets migFlag = 1.
 */
void test_demography_migration_matrix_sets_mig_flag(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [5, 5]\n"
        "  migration_matrix:\n"
        "    - row: [0.0, 1.5]\n"
        "    - row: [2.5, 0.0]\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, migFlag);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, migMatConst[0][0]);
    TEST_ASSERT_EQUAL_DOUBLE(1.5, migMatConst[0][1]);
    TEST_ASSERT_EQUAL_DOUBLE(2.5, migMatConst[1][0]);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, migMatConst[1][1]);
}

/*
 * demographic_events and demes_filename are mutually exclusive: the
 * demes file owns event generation when present, so explicit events
 * would conflict.
 */
void test_demography_rejects_events_and_demes_filename_both_set(void) {
    /* parse_demography_block rejects the combination before it ever
     * tries to read the demes file, so loading the file is not
     * strictly necessary here. We copy it in anyway so that if the
     * code is ever reordered to read the file first, a missing-file
     * error would not be mistaken for the rejection we want to
     * test. */
    copy_fixture_to_workdir("minimal_two_demes.yaml");

    const char *yaml =
        "simulation:\n"
        "  sample_size: 4\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [2, 2]\n"
        "  demes_filename: \"minimal_two_demes.yaml\"\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.05\n"
        "        population: 0\n"
        "        size: 0.5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * migration_matrix and demes_filename are mutually exclusive for the
 * same reason: demes-derived migrations would conflict with an
 * explicit matrix.
 */
void test_demography_rejects_migration_matrix_and_demes_filename_both_set(void) {
    copy_fixture_to_workdir("minimal_two_demes.yaml");

    const char *yaml =
        "simulation:\n"
        "  sample_size: 4\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "demography:\n"
        "  deme_sample_size: [2, 2]\n"
        "  demes_filename: \"minimal_two_demes.yaml\"\n"
        "  migration_matrix:\n"
        "    - row: [0.0, 1.0]\n"
        "    - row: [1.0, 0.0]\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/* ----- Selection block ----- */

/*
 * A single deterministic sweep block populates sweepMode, alpha,
 * sweepSite, and tau (= fixation_time_ago * 2.0). The sweep event
 * itself is checked by test_selection_emits_sweep_event_at_tau.
 */
void test_load_selection_section(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: deterministic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.6\n"
        "  fixation_time_ago: 0.05\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT('d', sweepMode);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.6, sweepSite);
    TEST_ASSERT_EQUAL_DOUBLE(0.10, tau);  /* 0.05 * 2.0 */
}

/*
 * Setting selection.initial_frequency activates softSweepMode and
 * writes the value into f0.
 */
void test_load_soft_sweep(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  initial_frequency: 0.01\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, sweepSite);
    TEST_ASSERT_EQUAL_DOUBLE(0.01, f0);
    TEST_ASSERT_EQUAL_INT(1, softSweepMode);
}

/*
 * selection_coefficient must be >= 0.
 */
void test_selection_rejects_negative_selection_coefficient(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: -1.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * sweep_mode: neutral with selection_coefficient > 0 is
 * contradictory; a neutral sweep has zero selection by definition.
 */
void test_selection_rejects_nonzero_coefficient_for_neutral_mode(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: neutral\n"
        "  selection_coefficient: 1.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * sweep_position must be in [0, 1]. The above-one and below-zero
 * cases are exercised in separate functions so each runs in a
 * fresh setUp / tearDown cycle.
 */
void test_selection_rejects_sweep_position_above_one(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 1.1\n"
        "  fixation_time_ago: 0.05\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

void test_selection_rejects_sweep_position_below_zero(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: -0.1\n"
        "  fixation_time_ago: 0.05\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * initial_frequency must be in [0, 1) — 1.0 is rejected.
 */
void test_selection_rejects_initial_frequency_ge_1(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  initial_frequency: 1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * final_frequency must be in [0, 1) — 1.0 is rejected.
 */
void test_selection_rejects_final_frequency_ge_1(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  final_frequency: 1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * When both initial and final frequencies are set (i.e. a partial
 * soft sweep), initial must be strictly less than final.
 */
void test_selection_rejects_initial_ge_final_when_both_set(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  initial_frequency: 0.5\n"
        "  final_frequency: 0.4\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * beneficial_mutation_rate must be >= 0.
 */
void test_selection_rejects_negative_beneficial_mutation_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  beneficial_mutation_rate: -1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * recurrent_sweep_rate must be >= 0. The negative-rate check runs
 * before the "non-recurrent block requires fixation_time_ago"
 * check, so this YAML need not supply fixation_time_ago to reach
 * the rejection.
 */
void test_selection_rejects_negative_recurrent_sweep_rate(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  recurrent_sweep_rate: -1.0\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * fixation_time_ago and recurrent_sweep_rate are mutually
 * exclusive; one specifies a single sweep, the other a recurrent
 * one.
 */
void test_selection_rejects_fixation_and_recurrent_both_set(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.05\n"
        "  recurrent_sweep_rate: 0.1\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * A non-recurrent selection block requires fixation_time_ago > 0;
 * without it, parse_selection_block has nothing to anchor the sweep
 * event on.
 */
void test_selection_requires_fixation_time_ago_when_not_recurrent(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * A non-recurrent selection block emits exactly one event of type
 * 's' at tau (= fixation_time_ago * 2.0).
 */
void test_selection_emits_sweep_event_at_tau(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 50.0\n"
        "  sweep_position: 0.4\n"
        "  fixation_time_ago: 0.025\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(0.05, tau);  /* 0.025 * 2.0 */
    TEST_ASSERT_EQUAL_INT(1, eventNumber);
    TEST_ASSERT_EQUAL_INT('s', events[0].type);
    TEST_ASSERT_EQUAL_DOUBLE(0.05, events[0].time);
}

/*
 * A recurrent-sweep selection block sets recurSweepRate and
 * recurSweepMode but does not emit a sweep event; recurrent sweeps
 * are driven by the rate, not anchored at tau.
 */
void test_selection_recurrent_sets_flags_without_event(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  recurrent_sweep_rate: 0.5\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, recurSweepRate);
    TEST_ASSERT_EQUAL_INT(1, recurSweepMode);
    TEST_ASSERT_EQUAL_INT(0, eventNumber);
}

/* ----- Output block + load-level + end-to-end ----- */

/*
 * An output block with finite_output and hide_partial_snp set
 * populates finiteOutputFlag and hidePartialSNP. The output_type is
 * haplotype here only because the boolean flags are independent of
 * which output_type is chosen.
 */
void test_load_output_section(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "output:\n"
        "  output_type: haplotype\n"
        "  finite_output: true\n"
        "  hide_partial_snp: true\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, finiteOutputFlag);
    TEST_ASSERT_EQUAL_INT(1, hidePartialSNP);
}

/*
 * output_type: tree_sequence with a tree_sequence_filename activates
 * tskit output mode and writes the filename into tskitOutputFilename.
 * unsimplified_tree_sequence: true unsets the minimal flag.
 */
void test_load_tskit_output(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "output:\n"
        "  output_type: tree_sequence\n"
        "  tree_sequence_filename: test_output.trees\n"
        "  unsimplified_tree_sequence: true\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(1, tskitOutputMode);
    TEST_ASSERT_EQUAL_STRING("test_output.trees", tskitOutputFilename);
    TEST_ASSERT_EQUAL_INT(0, minimalTreeSeq);
}

/*
 * output_type: tree_sequence requires tree_sequence_filename to be
 * set and non-empty; otherwise downstream tskit code has no place
 * to write.
 */
void test_output_tree_sequence_requires_filename(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "output:\n"
        "  output_type: tree_sequence\n";  /* no filename */
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * tree_sequence_filename is only meaningful when output_type is
 * tree_sequence. Setting it for haplotype or snp_array output is a
 * configuration error rather than a silently-ignored field.
 */
void test_output_filename_without_tree_sequence_mode_rejected(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "output:\n"
        "  output_type: haplotype\n"
        "  tree_sequence_filename: should_not_be_used.trees\n";
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, apply_from_yaml_string(yaml));
}

/*
 * unsimplified_tree_sequence: true causes minimalTreeSeq to be set
 * to 0. In a real discoal run, getParameters() in discoal_multipop.c
 * sets minimalTreeSeq = 1 before reading any flag, so the field is
 * normally 1 by the time the YAML parser runs and the YAML setting
 * flips it down to 0. The unit-test setUp / reset_config_globals
 * zeroes everything instead, so the test sets minimalTreeSeq to 1
 * by hand here to mirror the real starting state and make the
 * change observable.
 */
void test_output_unsimplified_tree_sequence_unsets_minimal(void) {
    minimalTreeSeq = 1;  /* mimic getParameters() initialisation */

    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n"
        "output:\n"
        "  output_type: tree_sequence\n"
        "  tree_sequence_filename: out.trees\n"
        "  unsimplified_tree_sequence: true\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(0, minimalTreeSeq);
}

/*
 * load_yaml_config rejects empty input. libcyaml itself returns OK
 * with a NULL top-level for empty files; the wrapper in
 * configInterface.c (around line 819) turns that case into
 * EXIT_FAILURE explicitly. Globals must remain unchanged since no
 * parsing took place.
 */
void test_load_empty_yaml(void) {
    int rc = apply_from_yaml_string("");
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, rc);
    TEST_ASSERT_EQUAL_INT(0, sampleSize);
    TEST_ASSERT_EQUAL_INT(0, npops);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, theta);
}

/*
 * load_yaml_config returns failure when the file does not exist.
 * apply_from_yaml_string is deliberately not used so that the
 * failure comes from cyaml_load_file rather than the temp-file
 * write step.
 */
void test_load_nonexistent_file(void) {
    struct discoal_config *cfg = NULL;
    int rc = load_yaml_config(
        "/tmp/configInterface_nonexistent_xyzabc.yaml", &cfg);
    TEST_ASSERT_NOT_EQUAL(EXIT_SUCCESS, rc);
}

/*
 * A YAML exercising every block is loaded and applied; the
 * cumulative effect on the discoal globals matches the YAML
 * inputs. Event ordering matters: parse_demography_block runs
 * before parse_selection_block, so the population_size_change
 * event lands at events[0] and the sweep event at events[1].
 */
void test_apply_full_config(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 20\n"
        "  num_replicates: 5\n"
        "  num_sites: 10000\n"
        "  seed: [42, 99]\n"
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  gene_conversion_rate: 1.0\n"
        "  gene_conversion_tract_length: 100\n"
        "demography:\n"
        "  deme_sample_size: [10, 10]\n"
        "  effective_population_size: 1000000\n"
        "  demographic_events:\n"
        "    population_size_changes:\n"
        "      - time: 0.05\n"
        "        population: 0\n"
        "        size: 0.5\n"
        "  migration_matrix:\n"
        "    - row: [0.0, 0.5]\n"
        "    - row: [0.7, 0.0]\n"
        "selection:\n"
        "  sweep_mode: stochastic\n"
        "  selection_coefficient: 100.0\n"
        "  sweep_position: 0.5\n"
        "  fixation_time_ago: 0.025\n"
        "output:\n"
        "  output_type: haplotype\n"
        "  finite_output: true\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);

    /* simulation */
    TEST_ASSERT_EQUAL_INT(20, sampleSize);
    TEST_ASSERT_EQUAL_INT(5, sampleNumber);
    TEST_ASSERT_EQUAL_INT(10000, nSites);
    TEST_ASSERT_EQUAL_INT(42, seed1);
    TEST_ASSERT_EQUAL_INT(99, seed2);

    /* genetics */
    TEST_ASSERT_EQUAL_DOUBLE(10.0, theta);
    TEST_ASSERT_EQUAL_DOUBLE(5.0, rho);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, my_gamma);
    TEST_ASSERT_EQUAL_INT(100, gcMean);

    /* demography */
    TEST_ASSERT_EQUAL_INT(2, npops);
    TEST_ASSERT_EQUAL_INT(10, sampleSizes[0]);
    TEST_ASSERT_EQUAL_INT(10, sampleSizes[1]);
    TEST_ASSERT_EQUAL_INT(1000000, EFFECTIVE_POPN_SIZE);
    TEST_ASSERT_EQUAL_INT(1, migFlag);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, migMatConst[0][1]);
    TEST_ASSERT_EQUAL_DOUBLE(0.7, migMatConst[1][0]);

    /* selection */
    TEST_ASSERT_EQUAL_INT('s', sweepMode);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, sweepSite);
    TEST_ASSERT_EQUAL_DOUBLE(0.05, tau);  /* 0.025 * 2.0 */

    /* output */
    TEST_ASSERT_EQUAL_INT(1, finiteOutputFlag);

    /* events: size-change at events[0] (demography first), sweep at events[1] */
    TEST_ASSERT_EQUAL_INT(2, eventNumber);
    TEST_ASSERT_EQUAL_INT('n', events[0].type);
    TEST_ASSERT_EQUAL_DOUBLE(0.10, events[0].time);  /* 0.05 * 2.0 */
    TEST_ASSERT_EQUAL_INT(0, events[0].popID);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, events[0].popnSize);
    TEST_ASSERT_EQUAL_INT('s', events[1].type);
    TEST_ASSERT_EQUAL_DOUBLE(0.05, events[1].time);
}

/*
 * A YAML with only the required blocks (simulation and genetics)
 * leaves every optional-block global at its initial value. Guards
 * against a parser for an absent block accidentally writing
 * something anyway.
 */
void test_apply_zero_block_yaml(void) {
    const char *yaml =
        "simulation:\n"
        "  sample_size: 10\n"
        "  num_replicates: 1\n"
        "  num_sites: 1000\n"
        "genetics:\n"
        "  mutation_rate: 0.0\n"
        "  recombination_rate: 0.0\n";

    int rc = apply_from_yaml_string(yaml);
    TEST_ASSERT_EQUAL_INT(EXIT_SUCCESS, rc);

    /* No demography block: EFFECTIVE_POPN_SIZE and migFlag untouched. */
    TEST_ASSERT_EQUAL_INT(0, EFFECTIVE_POPN_SIZE);
    TEST_ASSERT_EQUAL_INT(0, migFlag);

    /* No selection block: every selection global unchanged. */
    TEST_ASSERT_EQUAL_INT('\0', sweepMode);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, sweepSite);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, tau);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, f0);
    TEST_ASSERT_EQUAL_INT(0, softSweepMode);
    TEST_ASSERT_EQUAL_INT(0, partialSweepMode);
    TEST_ASSERT_EQUAL_INT(0, recurSweepMode);

    /* No output block: every output global unchanged. */
    TEST_ASSERT_EQUAL_INT(0, tskitOutputMode);
    TEST_ASSERT_EQUAL_INT(0, finiteOutputFlag);
    TEST_ASSERT_EQUAL_INT(0, hidePartialSNP);
    TEST_ASSERT_EQUAL_INT(0, minimalTreeSeq);

    /* No demography events and no selection sweep, so no events created. */
    TEST_ASSERT_EQUAL_INT(0, eventNumber);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_smoke_minimum_viable_config);

    /* Simulation block */
    RUN_TEST(test_load_simulation_section);
    RUN_TEST(test_load_seed_array);
    RUN_TEST(test_simulation_rejects_zero_sample_size);
    RUN_TEST(test_simulation_rejects_negative_sample_size);
    RUN_TEST(test_simulation_rejects_nonpositive_num_replicates);
    RUN_TEST(test_simulation_rejects_nonpositive_num_sites);
    RUN_TEST(test_simulation_sets_single_pop_defaults);

    /* Genetics block */
    RUN_TEST(test_load_genetics_section);
    RUN_TEST(test_load_gene_conversion);
    RUN_TEST(test_load_gene_conversion_crossover_ratio);
    RUN_TEST(test_genetics_rejects_negative_mutation_rate);
    RUN_TEST(test_genetics_rejects_negative_recombination_rate);
    RUN_TEST(test_genetics_rejects_negative_gene_conversion_rate);
    RUN_TEST(test_genetics_rejects_negative_crossover_ratio);
    RUN_TEST(test_genetics_rejects_negative_tract_length);
    RUN_TEST(test_genetics_rejects_rate_and_crossover_ratio_both_set);

    /* Demography block */
    RUN_TEST(test_load_populations_section);
    RUN_TEST(test_load_population_size_change_event);
    RUN_TEST(test_load_population_split_event);
    RUN_TEST(test_demography_demes_filename_loads_events);
    RUN_TEST(test_demography_rejects_negative_deme_sample_size);
    RUN_TEST(test_demography_rejects_deme_sum_mismatch);
    RUN_TEST(test_demography_rejects_negative_effective_population_size);
    RUN_TEST(test_demography_rejects_nonintegral_effective_population_size);
    RUN_TEST(test_demography_rejects_effective_population_size_over_intmax);
    RUN_TEST(test_demography_population_size_change_scales_time_by_2);
    RUN_TEST(test_demography_population_size_change_rejects_bad_pop_index);
    RUN_TEST(test_demography_population_size_change_rejects_negative_time);
    RUN_TEST(test_demography_population_size_change_rejects_nonpositive_size);
    RUN_TEST(test_demography_population_split_rejects_same_derived_ancestral);
    RUN_TEST(test_demography_population_split_rejects_bad_pop_index);
    RUN_TEST(test_demography_population_split_rejects_nonpositive_time);
    RUN_TEST(test_demography_rejects_migration_rate_changes_not_implemented);
    RUN_TEST(test_demography_rejects_ancient_samples_not_implemented);
    RUN_TEST(test_demography_migration_matrix_requires_square);
    RUN_TEST(test_demography_migration_matrix_requires_square_rows);
    RUN_TEST(test_demography_migration_matrix_rejects_nonzero_diagonal);
    RUN_TEST(test_demography_migration_matrix_rejects_negative_rate);
    RUN_TEST(test_demography_migration_matrix_sets_mig_flag);
    RUN_TEST(test_demography_rejects_events_and_demes_filename_both_set);
    RUN_TEST(test_demography_rejects_migration_matrix_and_demes_filename_both_set);

    /* Selection block */
    RUN_TEST(test_load_selection_section);
    RUN_TEST(test_load_soft_sweep);
    RUN_TEST(test_selection_rejects_negative_selection_coefficient);
    RUN_TEST(test_selection_rejects_nonzero_coefficient_for_neutral_mode);
    RUN_TEST(test_selection_rejects_sweep_position_above_one);
    RUN_TEST(test_selection_rejects_sweep_position_below_zero);
    RUN_TEST(test_selection_rejects_initial_frequency_ge_1);
    RUN_TEST(test_selection_rejects_final_frequency_ge_1);
    RUN_TEST(test_selection_rejects_initial_ge_final_when_both_set);
    RUN_TEST(test_selection_rejects_negative_beneficial_mutation_rate);
    RUN_TEST(test_selection_rejects_negative_recurrent_sweep_rate);
    RUN_TEST(test_selection_rejects_fixation_and_recurrent_both_set);
    RUN_TEST(test_selection_requires_fixation_time_ago_when_not_recurrent);
    RUN_TEST(test_selection_emits_sweep_event_at_tau);
    RUN_TEST(test_selection_recurrent_sets_flags_without_event);

    /* Output block + load-level + end-to-end */
    RUN_TEST(test_load_output_section);
    RUN_TEST(test_load_tskit_output);
    RUN_TEST(test_output_tree_sequence_requires_filename);
    RUN_TEST(test_output_filename_without_tree_sequence_mode_rejected);
    RUN_TEST(test_output_unsimplified_tree_sequence_unsets_minimal);
    RUN_TEST(test_load_empty_yaml);
    RUN_TEST(test_load_nonexistent_file);
    RUN_TEST(test_apply_full_config);
    RUN_TEST(test_apply_zero_block_yaml);
    return UNITY_END();
}
#endif

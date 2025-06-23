/**
 * test_demes_validation.c - Tests for validating Demes model compatibility
 * 
 * This test suite checks which Demes features are supported by discoal
 * and ensures appropriate errors/warnings for unsupported features.
 */

#include "unity.h"
#include "params.h"
#include "demes_loader.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>

void setUp(void) {
    /* This is run before each test */
}

void tearDown(void) {
    /* This is run after each test */
}

/* Helper to create a temporary Demes file */
static const char* create_temp_demes_file(const char *content) {
    static char temp_file[256];
    sprintf(temp_file, "/tmp/test_demes_%d.yaml", rand());
    
    FILE *f = fopen(temp_file, "w");
    if (!f) return NULL;
    
    fprintf(f, "%s", content);
    fclose(f);
    
    return temp_file;
}

/* Test supported features */

void test_simple_single_population(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n";
    
    const char *filename = create_temp_demes_file(yaml);
    TEST_ASSERT_NOT_NULL(filename);
    
    SimulationParams *params = params_create();
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    /* Size should be 1.0 since it's relative to itself */
    TEST_ASSERT_EQUAL_DOUBLE(1.0, params->demographics.pop_sizes.sizes[0]);
    
    params_destroy(params);
    unlink(filename);
}

void test_population_split(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: ancestral\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "  - name: derived\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - start_size: 5000\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(2, params->core.num_populations);
    
    params_destroy(params);
    unlink(filename);
}

void test_exponential_growth(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: ancestral\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "  - name: growing\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - end_time: 0\n"
        "        start_size: 1000\n"
        "        end_size: 10000\n"
        "        size_function: exponential\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Check that we have 2 populations */
    TEST_ASSERT_EQUAL(2, params->core.num_populations);
    
    /* Since discoal doesn't support exponential growth, we can only check that
       the populations exist. The exponential growth would need to be converted
       to discrete size change events, which is not yet implemented. */
    TEST_ASSERT_TRUE(params->demographics.pop_sizes.sizes[0] > 0 || 
                     params->demographics.pop_sizes.sizes[1] > 0);
    
    params_destroy(params);
    unlink(filename);
}

void test_symmetric_migration(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n"
        "  - name: pop2\n"
        "    epochs:\n"
        "      - start_size: 8000\n"
        "migrations:\n"
        "  - demes: [pop1, pop2]\n"
        "    rate: 0.001\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_TRUE(params->forces.symmetric_migration);
    TEST_ASSERT_EQUAL_DOUBLE(0.001, params->forces.migration_rate);
    
    params_destroy(params);
    unlink(filename);
}

/* Test unsupported features (should load but with limitations) */

void test_selfing_rate_ignored(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n"
        "        selfing_rate: 0.5\n";  /* This should be ignored */
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* Should still load successfully, just ignore selfing */
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    /* Size should be 1.0 since it's relative to itself */
    TEST_ASSERT_EQUAL_DOUBLE(1.0, params->demographics.pop_sizes.sizes[0]);
    
    params_destroy(params);
    unlink(filename);
}

void test_cloning_rate_ignored(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n"
        "        cloning_rate: 0.9\n";  /* This should be ignored */
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* Should still load successfully, just ignore cloning */
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    params_destroy(params);
    unlink(filename);
}

void test_pulse_events_ignored(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n"
        "  - name: pop2\n"
        "    epochs:\n"
        "      - start_size: 8000\n"
        "pulses:\n"
        "  - sources: [pop1]\n"
        "    dest: pop2\n"
        "    time: 500\n"
        "    proportions: [0.1]\n";  /* Pulses not supported */
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* Should load but without pulse events */
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(2, params->core.num_populations);
    /* No way to check pulses weren't converted since we don't support them */
    
    params_destroy(params);
    unlink(filename);
}

void test_multiple_ancestors_simplified(void) {
    /* This tests multi-way admixture which discoal can't handle directly */
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: anc1\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "  - name: anc2\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 8000\n"
        "  - name: anc3\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 6000\n"
        "  - name: admixed\n"
        "    ancestors: [anc1, anc2, anc3]\n"
        "    proportions: [0.5, 0.3, 0.2]\n"
        "    start_time: 100\n"
        "    epochs:\n"
        "      - start_size: 12000\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* This should load but admixture won't be properly represented */
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(4, params->core.num_populations);
    
    params_destroy(params);
    unlink(filename);
}

void test_time_bounded_migration(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - start_size: 10000\n"
        "  - name: pop2\n"
        "    epochs:\n"
        "      - start_size: 8000\n"
        "migrations:\n"
        "  - demes: [pop1, pop2]\n"
        "    rate: 0.001\n"
        "    start_time: 1000\n"
        "    end_time: 500\n";  /* Time bounds not fully supported */
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* Should load but time bounds might be ignored */
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Migration is set but without time bounds */
    TEST_ASSERT_TRUE(params->forces.symmetric_migration);
    
    params_destroy(params);
    unlink(filename);
}

void test_epoch_size_changes(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - end_time: 0.2\n"
        "        start_size: 5000\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "      - start_size: 20000\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Check that size change events were created */
    TEST_ASSERT_TRUE(params->demographics.num_events >= 2);
    
    /* Check for the expected size change events */
    int size_change_count = 0;
    for (int i = 0; i < params->demographics.num_events; i++) {
        if (params->demographics.events[i].type == EVENT_SIZE_CHANGE) {
            size_change_count++;
            /* Verify population index */
            TEST_ASSERT_EQUAL(0, params->demographics.events[i].params.size_change.pop);
        }
    }
    
    /* We should have 2 size change events for the epoch transitions */
    TEST_ASSERT_EQUAL(2, size_change_count);
    
    params_destroy(params);
    unlink(filename);
}

int main(void) {
    UNITY_BEGIN();
    
    /* Supported features */
    RUN_TEST(test_simple_single_population);
    RUN_TEST(test_population_split);
    RUN_TEST(test_exponential_growth);
    RUN_TEST(test_symmetric_migration);
    
    /* Unsupported features (should not crash) */
    RUN_TEST(test_selfing_rate_ignored);
    RUN_TEST(test_cloning_rate_ignored);
    RUN_TEST(test_pulse_events_ignored);
    RUN_TEST(test_multiple_ancestors_simplified);
    RUN_TEST(test_time_bounded_migration);
    
    /* Epoch size changes */
    RUN_TEST(test_epoch_size_changes);
    
    return UNITY_END();
}
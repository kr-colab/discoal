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
#include <fcntl.h>

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

void test_population_split_events(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: ancestral\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "  - name: derived1\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - start_size: 5000\n"
        "  - name: derived2\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - start_size: 8000\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Should have 3 populations */
    TEST_ASSERT_EQUAL(3, params->core.num_populations);
    
    /* Check for join events (splits are joins backward in time) */
    int join_count = 0;
    for (int i = 0; i < params->demographics.num_events; i++) {
        if (params->demographics.events[i].type == EVENT_JOIN) {
            join_count++;
            /* Verify the join time is 0.1 */
            TEST_ASSERT_EQUAL_DOUBLE(0.1, params->demographics.events[i].time);
        }
    }
    
    /* We should have 2 join events for the two derived populations */
    TEST_ASSERT_EQUAL(2, join_count);
    
    params_destroy(params);
    unlink(filename);
}

void test_migration_events(void) {
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
        "  - source: pop1\n"
        "    dest: pop2\n"
        "    start_time: 0.3\n"
        "    end_time: 0.1\n"
        "    rate: 0.002\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Check for migration events */
    int migration_count = 0;
    double first_time = -1, second_time = -1;
    double first_rate = -1, second_rate = -1;
    
    for (int i = 0; i < params->demographics.num_events; i++) {
        if (params->demographics.events[i].type == EVENT_MIGRATION_CHANGE) {
            migration_count++;
            if (migration_count == 1) {
                first_time = params->demographics.events[i].time;
                first_rate = params->demographics.events[i].params.migration.rate;
            } else if (migration_count == 2) {
                second_time = params->demographics.events[i].time;
                second_rate = params->demographics.events[i].params.migration.rate;
            }
        }
    }
    
    /* Should have 2 migration events (start and stop) */
    TEST_ASSERT_EQUAL(2, migration_count);
    
    /* After sorting, the stop event (0.1) comes before start event (0.3) */
    TEST_ASSERT_EQUAL_DOUBLE(0.1, first_time);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, first_rate);    /* Migration stops */
    TEST_ASSERT_EQUAL_DOUBLE(0.3, second_time);
    TEST_ASSERT_EQUAL_DOUBLE(0.002, second_rate); /* Migration starts */
    
    params_destroy(params);
    unlink(filename);
}

void test_warnings_generated(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: pop1\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "        selfing_rate: 0.1\n"
        "        cloning_rate: 0.05\n"
        "      - start_size: 5000\n"
        "        end_size: 10000\n"
        "        size_function: exponential\n"
        "  - name: pop2\n"
        "    ancestors: [pop1]\n"
        "    start_time: 0.08\n"
        "    epochs:\n"
        "      - start_size: 8000\n"
        "pulses:\n"
        "  - sources: [pop1]\n"
        "    dest: pop2\n"
        "    time: 0.05\n"
        "    proportions: [0.1]\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    /* Capture stderr to verify warnings are printed */
    int stderr_backup = dup(STDERR_FILENO);
    int pipefd[2];
    pipe(pipefd);
    dup2(pipefd[1], STDERR_FILENO);
    close(pipefd[1]);
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Restore stderr */
    fflush(stderr);
    dup2(stderr_backup, STDERR_FILENO);
    close(stderr_backup);
    
    /* Read captured output */
    char buffer[4096] = {0};
    read(pipefd[0], buffer, sizeof(buffer) - 1);
    close(pipefd[0]);
    
    /* Verify warnings were generated */
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Selfing rate"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Cloning rate"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Exponential growth"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "pulse event"));
    
    params_destroy(params);
    unlink(filename);
}

void test_event_sorting(void) {
    const char *yaml = 
        "time_units: 4N\n"
        "generation_time: 1\n"
        "demes:\n"
        "  - name: ancestral\n"
        "    epochs:\n"
        "      - end_time: 0.2\n"
        "        start_size: 20000\n"
        "  - name: pop1\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - end_time: 0.1\n"
        "        start_size: 10000\n"
        "      - start_size: 15000\n"
        "  - name: pop2\n"
        "    ancestors: [ancestral]\n"
        "    epochs:\n"
        "      - start_size: 8000\n"
        "migrations:\n"
        "  - source: pop1\n"
        "    dest: pop2\n"
        "    start_time: 0.15\n"
        "    end_time: 0.05\n"
        "    rate: 0.002\n";
    
    const char *filename = create_temp_demes_file(yaml);
    SimulationParams *params = params_create();
    
    int ret = demes_load_demographics(params, filename);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Check that events are sorted by time */
    TEST_ASSERT_TRUE(params->demographics.num_events >= 4);
    
    /* Verify events are in ascending time order */
    for (int i = 1; i < params->demographics.num_events; i++) {
        TEST_ASSERT_TRUE(params->demographics.events[i].time >= 
                        params->demographics.events[i-1].time);
    }
    
    /* Check specific event order */
    /* Should have events at times: 0.05, 0.1, 0.15, 0.2, 0.2 */
    TEST_ASSERT_EQUAL_DOUBLE(0.05, params->demographics.events[0].time);
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->demographics.events[1].time);
    TEST_ASSERT_EQUAL_DOUBLE(0.15, params->demographics.events[2].time);
    TEST_ASSERT_EQUAL_DOUBLE(0.2, params->demographics.events[3].time);
    TEST_ASSERT_EQUAL_DOUBLE(0.2, params->demographics.events[4].time);
    
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
    
    /* Population splits */
    RUN_TEST(test_population_split_events);
    
    /* Migration events */
    RUN_TEST(test_migration_events);
    
    /* Event sorting */
    RUN_TEST(test_event_sorting);
    
    /* Warnings */
    RUN_TEST(test_warnings_generated);
    
    return UNITY_END();
}
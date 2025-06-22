/**
 * test_params.c - Unit tests for parameter parsing and validation
 */

#include "unity.h"
#include "params.h"
#include <string.h>
#include <stdlib.h>

void setUp(void) {
    /* This is run before each test */
}

void tearDown(void) {
    /* This is run after each test */
}

void test_params_create_destroy(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test default values */
    TEST_ASSERT_EQUAL(0, params->core.total_samples);
    TEST_ASSERT_EQUAL(1, params->core.num_replicates);
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, params->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, params->forces.rho);
    
    params_destroy(params);
}

void test_params_basic_parsing(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test basic command line */
    char *argv[] = {"discoal", "100", "1000", "10000"};
    int argc = 4;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL(100, params->core.sample_sizes[0]);
    
    params_destroy(params);
}

void test_params_mutation_recombination(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "50", "100", "5000", "-t", "100.5", "-r", "200.0"};
    int argc = 8;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.5, params->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(200.0, params->forces.rho);
    
    params_destroy(params);
}

void test_params_multiple_populations(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100", "10", "1000", "-p", "3", "30", "40", "30"};
    int argc = 9;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(3, params->core.num_populations);
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(30, params->core.sample_sizes[0]);
    TEST_ASSERT_EQUAL(40, params->core.sample_sizes[1]);
    TEST_ASSERT_EQUAL(30, params->core.sample_sizes[2]);
    
    params_destroy(params);
}

void test_params_selection(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100", "10", "1000", "-a", "100", "-x", "0.5", "-ws", "0.01"};
    int argc = 10;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->selection.alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, params->selection.sweep_position);
    TEST_ASSERT_EQUAL_DOUBLE(0.02, params->selection.tau);  /* 0.01 * 2 */
    TEST_ASSERT_EQUAL('s', params->selection.sweep_mode);
    
    params_destroy(params);
}

void test_params_partial_sweep(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100", "10", "1000", "-c", "0.8", "-a", "100"};
    int argc = 8;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_TRUE(params->selection.partial_sweep);
    TEST_ASSERT_EQUAL_DOUBLE(0.8, params->selection.partial_sweep_final_freq);
    TEST_ASSERT_EQUAL('s', params->selection.sweep_mode);  /* Should be stochastic */
    
    params_destroy(params);
}

void test_params_gene_conversion(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test -gr flag */
    char *argv[] = {"discoal", "100", "10", "1000", "-r", "100", "-gr", "0.1", "500"};
    int argc = 9;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.rho);
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->forces.gene_conversion_ratio);
    TEST_ASSERT_EQUAL(500, params->forces.gc_tract_mean);
    
    params_destroy(params);
}

void test_params_migration(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test symmetric migration */
    char *argv[] = {"discoal", "100", "10", "1000", "-p", "2", "50", "50", "-M", "1.5"};
    int argc = 10;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(1.5, params->forces.migration_matrix[0][1]);
    TEST_ASSERT_EQUAL_DOUBLE(1.5, params->forces.migration_matrix[1][0]);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, params->forces.migration_matrix[0][0]);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, params->forces.migration_matrix[1][1]);
    
    params_destroy(params);
}

void test_params_demographic_events(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100", "10", "1000", "-en", "0.05", "0", "2.0", "-en", "0.1", "0", "0.5"};
    int argc = 12;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(2, params->demographics.num_events);
    
    /* Events should be sorted by time */
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->demographics.events[0].time);  /* 0.05 * 2 */
    TEST_ASSERT_EQUAL(EVENT_SIZE_CHANGE, params->demographics.events[0].type);
    TEST_ASSERT_EQUAL(0, params->demographics.events[0].params.size_change.pop);
    TEST_ASSERT_EQUAL_DOUBLE(2.0, params->demographics.events[0].params.size_change.size);
    
    TEST_ASSERT_EQUAL_DOUBLE(0.2, params->demographics.events[1].time);  /* 0.1 * 2 */
    TEST_ASSERT_EQUAL(EVENT_SIZE_CHANGE, params->demographics.events[1].type);
    TEST_ASSERT_EQUAL(0, params->demographics.events[1].params.size_change.pop);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, params->demographics.events[1].params.size_change.size);
    
    params_destroy(params);
}

void test_params_validation_errors(void) {
    SimulationParams *params = params_create();
    if (!params) {
        TEST_FAIL_MESSAGE("Failed to create params");
        return;
    }
    
    char error_msg[256];
    
    /* Set up valid base parameters first */
    params->core.num_replicates = 10;
    params->core.num_sites = 1000;
    params->core.num_populations = 1;
    
    /* Test negative samples */
    params->core.total_samples = -10;
    error_msg[0] = '\0';  /* Clear error message */
    int ret = params_validate(params, error_msg, sizeof(error_msg));
    if (ret == 0) {
        TEST_FAIL_MESSAGE("Validation should fail for negative samples but returned success");
    }
    
    /* Test invalid gene conversion ratio */
    params->core.total_samples = 100;
    params->core.sample_sizes[0] = 100;  /* Make sure sample sizes match */
    params->forces.gene_conversion_ratio = 1.5;
    error_msg[0] = '\0';  /* Clear error message */
    ret = params_validate(params, error_msg, sizeof(error_msg));
    if (ret == 0) {
        TEST_FAIL_MESSAGE("Validation should fail for GC ratio > 1 but returned success");
    }
    
    params_destroy(params);
}

void test_params_tree_sequence_output(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100", "10", "1000", "-ts", "output.trees", "-F"};
    int argc = 7;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_TRUE(params->output.tree_sequence_output);
    TEST_ASSERT_EQUAL_STRING("output.trees", params->output.output_filename);
    TEST_ASSERT_FALSE(params->output.minimal_tree_seq);  /* -F flag disables minimal mode */
    
    params_destroy(params);
}

void test_params_copy(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Set up some values */
    char *argv[] = {"discoal", "100", "10", "1000", "-t", "50", "-en", "0.1", "0", "2.0"};
    int argc = 10;
    params_load_from_args(params, argc, argv);
    
    /* Copy the params */
    SimulationParams *copy = params_copy(params);
    TEST_ASSERT_NOT_NULL(copy);
    
    /* Verify copy */
    TEST_ASSERT_EQUAL(params->core.total_samples, copy->core.total_samples);
    TEST_ASSERT_EQUAL_DOUBLE(params->forces.theta, copy->forces.theta);
    TEST_ASSERT_EQUAL(params->demographics.num_events, copy->demographics.num_events);
    
    /* Verify deep copy of events */
    TEST_ASSERT_NOT_EQUAL(params->demographics.events, copy->demographics.events);
    TEST_ASSERT_EQUAL_DOUBLE(params->demographics.events[0].time, 
                            copy->demographics.events[0].time);
    
    params_destroy(params);
    params_destroy(copy);
}

int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_params_create_destroy);
    RUN_TEST(test_params_basic_parsing);
    RUN_TEST(test_params_mutation_recombination);
    RUN_TEST(test_params_multiple_populations);
    RUN_TEST(test_params_selection);
    RUN_TEST(test_params_partial_sweep);
    RUN_TEST(test_params_gene_conversion);
    RUN_TEST(test_params_migration);
    RUN_TEST(test_params_demographic_events);
    RUN_TEST(test_params_validation_errors);
    RUN_TEST(test_params_tree_sequence_output);
    RUN_TEST(test_params_copy);
    
    return UNITY_END();
}
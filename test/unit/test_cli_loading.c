/**
 * test_cli_loading.c - Unit tests for command-line parameter loading
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

void test_basic_positional_args(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "50", "1000", "10000"};
    int argc = 4;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(50, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL(50, params->core.sample_sizes[0]);
    
    params_destroy(params);
}

void test_theta_option(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "50", "1000", "10000", "-t", "100"};
    int argc = 6;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.theta);
    
    params_destroy(params);
}

void test_yaml_loading_flag(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "-yaml", "test/demes_examples/single_pop.yaml"};
    int argc = 3;
    
    /* This will fail because the YAML file doesn't have full discoal params */
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_NOT_EQUAL(0, ret);
    
    params_destroy(params);
}

void test_demes_loading_flag(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "-demes", "test/demes_examples/single_pop.yaml", 
                    "50", "1000", "10000"};
    int argc = 6;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(50, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL_DOUBLE(10000.0, params->demographics.pop_sizes.sizes[0]);
    
    params_destroy(params);
}

void test_demes_two_populations(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "-demes", "test/demes_examples/two_pop_migration.yaml", 
                    "50", "1000", "10000", "-t", "100"};
    int argc = 8;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL(2, params->core.num_populations);
    TEST_ASSERT_EQUAL(25, params->core.sample_sizes[0]);
    TEST_ASSERT_EQUAL(25, params->core.sample_sizes[1]);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.theta);
    
    /* Check that migration was set */
    TEST_ASSERT_TRUE(params->forces.symmetric_migration);
    TEST_ASSERT_EQUAL_DOUBLE(0.0001, params->forces.migration_rate);
    
    params_destroy(params);
}

void test_sample_size_limit(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    char *argv[] = {"discoal", "100000", "1000", "10000"};
    int argc = 4;
    
    int ret = params_load_from_args(params, argc, argv);
    TEST_ASSERT_NOT_EQUAL(0, ret);  /* Should fail due to sample size limit */
    
    params_destroy(params);
}

int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_basic_positional_args);
    RUN_TEST(test_theta_option);
    RUN_TEST(test_yaml_loading_flag);
    RUN_TEST(test_demes_loading_flag);
    RUN_TEST(test_demes_two_populations);
    RUN_TEST(test_sample_size_limit);
    
    return UNITY_END();
}
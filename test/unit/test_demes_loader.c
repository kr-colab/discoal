/**
 * test_demes_loader.c - Unit tests for Demes integration
 */

#include "unity.h"
#include "params.h"
#include "demes_loader.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

void setUp(void) {
    /* This is run before each test */
}

void tearDown(void) {
    /* This is run after each test */
}

void test_demes_single_population(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test loading single population Demes file */
    int ret = demes_load_demographics(params, "test/demes_examples/single_pop.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify population parameters */
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL_DOUBLE(10000.0, params->demographics.pop_sizes.sizes[0]);
    TEST_ASSERT_EQUAL_DOUBLE(10000.0, params->demographics.pop_sizes.initial_sizes[0]);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, params->demographics.pop_sizes.growth_rates[0]);
    
    params_destroy(params);
}

void test_demes_bottleneck(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test loading bottleneck Demes file */
    int ret = demes_load_demographics(params, "test/demes_examples/bottleneck.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify we have 2 populations (demes) */
    TEST_ASSERT_EQUAL(2, params->core.num_populations);
    
    /* Verify population sizes */
    /* Note: The exact mapping depends on how demes are ordered */
    TEST_ASSERT_TRUE(params->demographics.pop_sizes.sizes[0] > 0);
    TEST_ASSERT_TRUE(params->demographics.pop_sizes.sizes[1] > 0);
    
    params_destroy(params);
}

void test_demes_time_conversion(void) {
    struct demes_graph graph;
    
    /* Test generations */
    graph.time_units = (demes_char_t*)"generations";
    graph.generation_time = 1.0;
    TEST_ASSERT_EQUAL_DOUBLE(100.0, demes_time_to_generations(100.0, &graph));
    
    /* Test years */
    graph.time_units = (demes_char_t*)"years";
    graph.generation_time = 25.0;
    TEST_ASSERT_EQUAL_DOUBLE(4.0, demes_time_to_generations(100.0, &graph));
    TEST_ASSERT_EQUAL_DOUBLE(40.0, demes_time_to_generations(1000.0, &graph));
}

void test_yaml_with_demes_file(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Create a temporary YAML file with demes_file reference */
    const char *temp_file = "/tmp/test_demes_ref.yaml";
    FILE *f = fopen(temp_file, "w");
    TEST_ASSERT_NOT_NULL(f);
    
    fprintf(f, "simulation:\n");
    fprintf(f, "  samples: 100\n");
    fprintf(f, "  replicates: 1000\n");
    fprintf(f, "  sites: 10000\n");
    fprintf(f, "\n");
    fprintf(f, "evolution:\n");
    fprintf(f, "  mutation:\n");
    fprintf(f, "    theta: 100.0\n");
    fprintf(f, "\n");
    fprintf(f, "demographics:\n");
    fprintf(f, "  demes_file: test/demes_examples/single_pop.yaml\n");
    fclose(f);
    
    /* Load with Demes support */
    int ret = yaml_load_params_with_demes(params, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify standard parameters loaded */
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.theta);
    
    /* Verify demographics loaded from Demes file */
    TEST_ASSERT_EQUAL(1, params->core.num_populations);
    TEST_ASSERT_EQUAL_DOUBLE(10000.0, params->demographics.pop_sizes.sizes[0]);
    
    /* Clean up */
    unlink(temp_file);
    params_destroy(params);
}

int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_demes_single_population);
    RUN_TEST(test_demes_bottleneck);
    RUN_TEST(test_demes_time_conversion);
    RUN_TEST(test_yaml_with_demes_file);
    
    return UNITY_END();
}
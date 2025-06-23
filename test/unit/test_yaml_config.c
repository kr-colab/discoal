/**
 * test_yaml_config.c - Unit tests for YAML configuration loading
 */

#include "unity.h"
#include "params.h"
#include "version.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

void setUp(void) {
    /* This is run before each test */
}

void tearDown(void) {
    /* This is run after each test */
}

void test_yaml_basic_neutral(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Test loading basic neutral config */
    int ret = params_load_from_yaml(params, "test/configs/basic_neutral.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify loaded values */
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.rho);
    TEST_ASSERT_EQUAL(OUTPUT_MS, params->output.format);
    
    params_destroy(params);
}

void test_yaml_selection_sweep(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    int ret = params_load_from_yaml(params, "test/configs/selection_sweep.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify selection parameters */
    TEST_ASSERT_EQUAL(50, params->core.total_samples);
    TEST_ASSERT_EQUAL_DOUBLE(10000.0, params->core.Ne);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->selection.alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, params->selection.sweep_position);
    TEST_ASSERT_EQUAL(SWEEP_STOCHASTIC, params->selection.sweep_mode);
    TEST_ASSERT_EQUAL_DOUBLE(0.02, params->selection.tau);  /* 0.01 * 2 */
    
    /* Verify tree sequence output */
    TEST_ASSERT_TRUE(params->output.tree_sequence_output);
    TEST_ASSERT_EQUAL_STRING("sweep_simulation.trees", params->output.output_filename);
    TEST_ASSERT_TRUE(params->output.minimal_tree_seq);
    
    params_destroy(params);
}

void test_yaml_multi_population(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    int ret = params_load_from_yaml(params, "test/configs/multi_population.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify population structure */
    TEST_ASSERT_EQUAL(3, params->core.num_populations);
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(40, params->core.sample_sizes[0]);
    TEST_ASSERT_EQUAL(30, params->core.sample_sizes[1]);
    TEST_ASSERT_EQUAL(30, params->core.sample_sizes[2]);
    
    /* Verify migration */
    TEST_ASSERT_TRUE(params->forces.symmetric_migration);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, params->forces.migration_rate);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, params->forces.migration_matrix[0][1]);
    TEST_ASSERT_EQUAL_DOUBLE(1.0, params->forces.migration_matrix[1][0]);
    
    /* Verify demographic events */
    TEST_ASSERT_EQUAL_MESSAGE(2, params->demographics.num_events, "Should have 2 demographic events");
    
    /* First event - size change */
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->demographics.events[0].time);  /* 0.05 * 2 */
    TEST_ASSERT_EQUAL(EVENT_SIZE_CHANGE, params->demographics.events[0].type);
    TEST_ASSERT_EQUAL(0, params->demographics.events[0].params.size_change.pop);
    TEST_ASSERT_EQUAL_DOUBLE(2.0, params->demographics.events[0].params.size_change.size);
    
    /* Second event - join */
    TEST_ASSERT_EQUAL_DOUBLE(0.2, params->demographics.events[1].time);  /* 0.1 * 2 */
    TEST_ASSERT_EQUAL(EVENT_JOIN, params->demographics.events[1].type);
    
    params_destroy(params);
}

void test_yaml_gene_conversion(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    int ret = params_load_from_yaml(params, "test/configs/gene_conversion.yaml");
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify gene conversion parameters */
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.rho);
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->forces.gene_conversion_ratio);
    TEST_ASSERT_EQUAL(500, params->forces.gc_tract_mean);
    
    params_destroy(params);
}

void test_yaml_save_and_reload(void) {
    SimulationParams *params1 = params_create();
    TEST_ASSERT_NOT_NULL(params1);
    
    /* Set up some parameters */
    params1->core.total_samples = 75;
    params1->core.num_replicates = 250;
    params1->core.num_sites = 7500;
    params1->forces.theta = 55.5;
    params1->forces.rho = 77.7;
    params1->selection.alpha = 150.0;
    params1->selection.sweep_position = 0.3;
    params1->selection.sweep_mode = SWEEP_DETERMINISTIC;
    
    /* Add a demographic event */
    DemographicEvent event;
    memset(&event, 0, sizeof(event));
    event.time = 0.2;  /* 0.1 * 2 in coalescent units */
    event.type = EVENT_SIZE_CHANGE;
    event.params.size_change.pop = 0;
    event.params.size_change.size = 3.0;
    params_add_demographic_event(params1, &event);
    
    /* Save to YAML */
    const char *temp_file = "/tmp/test_params_save.yaml";
    int ret = params_save_to_yaml(params1, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Load into new params */
    SimulationParams *params2 = params_create();
    TEST_ASSERT_NOT_NULL(params2);
    
    ret = params_load_from_yaml(params2, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify values match */
    TEST_ASSERT_EQUAL(params1->core.total_samples, params2->core.total_samples);
    TEST_ASSERT_EQUAL(params1->core.num_replicates, params2->core.num_replicates);
    TEST_ASSERT_EQUAL(params1->core.num_sites, params2->core.num_sites);
    TEST_ASSERT_EQUAL_DOUBLE(params1->forces.theta, params2->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(params1->forces.rho, params2->forces.rho);
    TEST_ASSERT_EQUAL_DOUBLE(params1->selection.alpha, params2->selection.alpha);
    TEST_ASSERT_EQUAL_DOUBLE(params1->selection.sweep_position, params2->selection.sweep_position);
    TEST_ASSERT_EQUAL(params1->selection.sweep_mode, params2->selection.sweep_mode);
    
    /* Verify demographic event */
    TEST_ASSERT_EQUAL(1, params2->demographics.num_events);
    TEST_ASSERT_EQUAL_DOUBLE(0.2, params2->demographics.events[0].time);
    
    /* Clean up */
    unlink(temp_file);
    params_destroy(params1);
    params_destroy(params2);
}

void test_yaml_invalid_file(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Try to load non-existent file */
    int ret = params_load_from_yaml(params, "/tmp/this_file_does_not_exist.yaml");
    TEST_ASSERT_NOT_EQUAL(0, ret);
    
    params_destroy(params);
}

void test_yaml_validation_errors(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Create a YAML file with invalid parameters */
    const char *temp_file = "/tmp/test_invalid_params.yaml";
    FILE *f = fopen(temp_file, "w");
    TEST_ASSERT_NOT_NULL(f);
    
    fprintf(f, "simulation:\n");
    fprintf(f, "  samples: -50\n");  /* Invalid negative samples */
    fprintf(f, "  replicates: 100\n");
    fprintf(f, "  sites: 1000\n");
    fclose(f);
    
    /* Loading should fail due to validation */
    int ret = params_load_from_yaml(params, temp_file);
    TEST_ASSERT_NOT_EQUAL(0, ret);
    
    /* Clean up */
    unlink(temp_file);
    params_destroy(params);
}

int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_yaml_basic_neutral);
    RUN_TEST(test_yaml_selection_sweep);
    RUN_TEST(test_yaml_multi_population);
    RUN_TEST(test_yaml_gene_conversion);
    RUN_TEST(test_yaml_save_and_reload);
    RUN_TEST(test_yaml_invalid_file);
    RUN_TEST(test_yaml_validation_errors);
    
    return UNITY_END();
}
/**
 * test_yaml_loader.c - Unit tests for YAML configuration loading with libyaml
 */

#include "unity.h"
#include "params.h"
#include "yaml_loader.h"
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

void setUp(void) {
    /* This is run before each test */
}

void tearDown(void) {
    /* This is run after each test */
}

void test_yaml_basic_loading(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    /* Create a temporary YAML file */
    const char *temp_file = "/tmp/test_basic.yaml";
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
    fprintf(f, "  recombination:\n");
    fprintf(f, "    rho: 100.0\n");
    fprintf(f, "\n");
    fprintf(f, "output:\n");
    fprintf(f, "  format: ms\n");
    fclose(f);
    
    /* Load the YAML file */
    int ret = yaml_load_params(params, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify loaded values */
    TEST_ASSERT_EQUAL(100, params->core.total_samples);
    TEST_ASSERT_EQUAL(1000, params->core.num_replicates);
    TEST_ASSERT_EQUAL(10000, params->core.num_sites);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.rho);
    TEST_ASSERT_EQUAL(OUTPUT_MS, params->output.format);
    
    /* Clean up */
    unlink(temp_file);
    params_destroy(params);
}

void test_yaml_selection_params(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    const char *temp_file = "/tmp/test_selection.yaml";
    FILE *f = fopen(temp_file, "w");
    TEST_ASSERT_NOT_NULL(f);
    
    fprintf(f, "simulation:\n");
    fprintf(f, "  samples: 50\n");
    fprintf(f, "  replicates: 100\n");
    fprintf(f, "  sites: 5000\n");
    fprintf(f, "\n");
    fprintf(f, "selection:\n");
    fprintf(f, "  coefficient: 100.0\n");
    fprintf(f, "  position: 0.5\n");
    fprintf(f, "  mode: stochastic\n");
    fprintf(f, "  time: 0.01\n");
    fclose(f);
    
    int ret = yaml_load_params(params, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->selection.alpha);
    TEST_ASSERT_EQUAL_DOUBLE(0.5, params->selection.sweep_position);
    TEST_ASSERT_EQUAL(SWEEP_STOCHASTIC, params->selection.sweep_mode);
    TEST_ASSERT_EQUAL_DOUBLE(0.02, params->selection.tau);  /* 0.01 * 2 */
    
    unlink(temp_file);
    params_destroy(params);
}

void test_yaml_gene_conversion(void) {
    SimulationParams *params = params_create();
    TEST_ASSERT_NOT_NULL(params);
    
    const char *temp_file = "/tmp/test_gc.yaml";
    FILE *f = fopen(temp_file, "w");
    TEST_ASSERT_NOT_NULL(f);
    
    fprintf(f, "simulation:\n");
    fprintf(f, "  samples: 50\n");
    fprintf(f, "  replicates: 500\n");
    fprintf(f, "  sites: 5000\n");
    fprintf(f, "\n");
    fprintf(f, "evolution:\n");
    fprintf(f, "  recombination:\n");
    fprintf(f, "    rho: 100.0\n");
    fprintf(f, "    gene_conversion:\n");
    fprintf(f, "      ratio: 0.1\n");
    fprintf(f, "      tract_length: 500\n");
    fclose(f);
    
    int ret = yaml_load_params(params, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    TEST_ASSERT_EQUAL_DOUBLE(100.0, params->forces.rho);
    TEST_ASSERT_EQUAL_DOUBLE(0.1, params->forces.gene_conversion_ratio);
    TEST_ASSERT_EQUAL(500, params->forces.gc_tract_mean);
    
    unlink(temp_file);
    params_destroy(params);
}

void test_yaml_save_reload(void) {
    SimulationParams *params1 = params_create();
    TEST_ASSERT_NOT_NULL(params1);
    
    /* Set up some parameters */
    params1->core.total_samples = 75;
    params1->core.num_replicates = 250;
    params1->core.num_sites = 7500;
    params1->forces.theta = 55.5;
    params1->forces.rho = 77.7;
    
    /* Save to YAML */
    const char *temp_file = "/tmp/test_save_reload.yaml";
    int ret = yaml_save_params(params1, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Load into new params */
    SimulationParams *params2 = params_create();
    TEST_ASSERT_NOT_NULL(params2);
    
    ret = yaml_load_params(params2, temp_file);
    TEST_ASSERT_EQUAL(0, ret);
    
    /* Verify basic values match */
    TEST_ASSERT_EQUAL(params1->core.total_samples, params2->core.total_samples);
    TEST_ASSERT_EQUAL(params1->core.num_replicates, params2->core.num_replicates);
    TEST_ASSERT_EQUAL(params1->core.num_sites, params2->core.num_sites);
    TEST_ASSERT_EQUAL_DOUBLE(params1->forces.theta, params2->forces.theta);
    TEST_ASSERT_EQUAL_DOUBLE(params1->forces.rho, params2->forces.rho);
    
    /* Clean up */
    unlink(temp_file);
    params_destroy(params1);
    params_destroy(params2);
}

int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_yaml_basic_loading);
    RUN_TEST(test_yaml_selection_params);
    RUN_TEST(test_yaml_gene_conversion);
    RUN_TEST(test_yaml_save_reload);
    
    return UNITY_END();
}
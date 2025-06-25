#include "unity.h"
#include "configInterface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Test fixture
static SimulationConfig testConfig;
static char tempFilename[256];

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    // Initialize test config to default values
    initializeDefaultConfig(&testConfig);
    
    // Create temporary filename for test YAML files
    sprintf(tempFilename, "/tmp/test_config_%d.yaml", getpid());
}

void tearDown(void) {
    // Clean up temporary file if it exists
    unlink(tempFilename);
}
#endif

// Helper function to create a YAML file with specific content
static void createYAMLFile(const char *content) {
    FILE *fp = fopen(tempFilename, "w");
    TEST_ASSERT_NOT_NULL(fp);
    fprintf(fp, "%s", content);
    fclose(fp);
}

// Test initialization of default configuration
void test_initialize_default_config(void) {
    SimulationConfig config;
    initializeDefaultConfig(&config);
    
    // Test default values
    TEST_ASSERT_EQUAL(0, config.sample_size);
    TEST_ASSERT_EQUAL(0, config.num_replicates);
    TEST_ASSERT_EQUAL(0, config.num_sites);
    TEST_ASSERT_EQUAL(0, config.seed1);
    TEST_ASSERT_EQUAL(0, config.seed2);
    TEST_ASSERT_EQUAL(0, config.has_seed);
    
    TEST_ASSERT_EQUAL(1, config.npops);
    TEST_ASSERT_EQUAL(1000000, config.effective_popn_size);
    TEST_ASSERT_EQUAL(0, config.has_populations);
    
    TEST_ASSERT_EQUAL(0.0, config.theta);
    TEST_ASSERT_EQUAL(0.0, config.rho);
    TEST_ASSERT_EQUAL(0.0, config.gamma);
    TEST_ASSERT_EQUAL(0, config.gc_mean);
    TEST_ASSERT_EQUAL(0.0, config.gamma_co_ratio);
    TEST_ASSERT_EQUAL(0, config.gamma_co_ratio_mode);
    TEST_ASSERT_EQUAL(0, config.has_genetics);
    
    TEST_ASSERT_EQUAL(0.0, config.alpha);
    TEST_ASSERT_EQUAL(0.5, config.sweep_site);
    TEST_ASSERT_EQUAL('s', config.sweep_mode);
    TEST_ASSERT_EQUAL(0, config.has_selection);
    
    TEST_ASSERT_EQUAL('h', config.output_style);
    TEST_ASSERT_EQUAL(0, config.finite_output_flag);
    TEST_ASSERT_EQUAL(0, config.tskit_output_mode);
    TEST_ASSERT_EQUAL(1, config.minimal_tree_seq);
    TEST_ASSERT_EQUAL(0, config.has_output);
    
    TEST_ASSERT_EQUAL(0, config.use_demes);
    TEST_ASSERT_EQUAL(0, config.num_explicit_events);
}

// Test loading simulation section
void test_load_simulation_section(void) {
    const char *yaml_content = 
        "simulation:\n"
        "  sample_size: 20\n"
        "  num_replicates: 10\n"
        "  num_sites: 100000\n"
        "  effective_popn_size: 2000000\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(20, config.sample_size);
    TEST_ASSERT_EQUAL(10, config.num_replicates);
    TEST_ASSERT_EQUAL(100000, config.num_sites);
    TEST_ASSERT_EQUAL(2000000, config.effective_popn_size);
}

// Test loading seed as scalar
void test_load_seed_scalar(void) {
    const char *yaml_content = 
        "simulation:\n"
        "  seed: 12345\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(12345, config.seed1);
    TEST_ASSERT_EQUAL(1, config.has_seed);
}

// Test loading seed as array
void test_load_seed_array(void) {
    const char *yaml_content = 
        "simulation:\n"
        "  seed: [12345, 67890]\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(12345, config.seed1);
    TEST_ASSERT_EQUAL(67890, config.seed2);
    TEST_ASSERT_EQUAL(1, config.has_seed);
}

// Test loading genetics section
void test_load_genetics_section(void) {
    const char *yaml_content = 
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(10.0, config.theta);
    TEST_ASSERT_EQUAL(5.0, config.rho);
    TEST_ASSERT_EQUAL(1, config.has_genetics);
}

// Test loading gene conversion parameters
void test_load_gene_conversion(void) {
    const char *yaml_content = 
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  gene_conversion:\n"
        "    rate: 2.5\n"
        "    tract_length: 500\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(10.0, config.theta);
    TEST_ASSERT_EQUAL(5.0, config.rho);
    TEST_ASSERT_EQUAL(2.5, config.gamma);
    TEST_ASSERT_EQUAL(500, config.gc_mean);
    TEST_ASSERT_EQUAL(1, config.has_genetics);
}

// Test loading gene conversion with crossover ratio
void test_load_gene_conversion_crossover_ratio(void) {
    const char *yaml_content = 
        "genetics:\n"
        "  mutation_rate: 10.0\n"
        "  recombination_rate: 5.0\n"
        "  gene_conversion:\n"
        "    crossover_ratio: 1.5\n"
        "    tract_length: 500\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(10.0, config.theta);
    TEST_ASSERT_EQUAL(5.0, config.rho);
    TEST_ASSERT_EQUAL(1.5, config.gamma_co_ratio);
    TEST_ASSERT_EQUAL(1, config.gamma_co_ratio_mode);
    TEST_ASSERT_EQUAL(500, config.gc_mean);
    TEST_ASSERT_EQUAL(1, config.has_genetics);
}

// Test loading populations section
void test_load_populations_section(void) {
    const char *yaml_content = 
        "populations:\n"
        "  sample_sizes: [10, 5, 8]\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(3, config.npops);
    TEST_ASSERT_EQUAL(10, config.sample_sizes[0]);
    TEST_ASSERT_EQUAL(5, config.sample_sizes[1]);
    TEST_ASSERT_EQUAL(8, config.sample_sizes[2]);
    TEST_ASSERT_EQUAL(23, config.sample_size);  // Total sample size
    TEST_ASSERT_EQUAL(1, config.has_populations);
}

// Test loading selection section
void test_load_selection_section(void) {
    const char *yaml_content = 
        "selection:\n"
        "  alpha: 100.0\n"
        "  sweep_site: 0.6\n"
        "  sweep_mode: deterministic\n"
        "  tau: 0.05\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(100.0, config.alpha);
    TEST_ASSERT_EQUAL(0.6, config.sweep_site);
    TEST_ASSERT_EQUAL('d', config.sweep_mode);
    TEST_ASSERT_EQUAL(0.05, config.tau);
    TEST_ASSERT_EQUAL(1, config.has_selection);
}

// Test loading soft sweep parameters
void test_load_soft_sweep(void) {
    const char *yaml_content = 
        "selection:\n"
        "  alpha: 100.0\n"
        "  sweep_site: 0.5\n"
        "  soft_sweep:\n"
        "    initial_frequency: 0.01\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(100.0, config.alpha);
    TEST_ASSERT_EQUAL(0.5, config.sweep_site);
    TEST_ASSERT_EQUAL(0.01, config.f0);
    TEST_ASSERT_EQUAL(1, config.soft_sweep_mode);
    TEST_ASSERT_EQUAL(1, config.has_selection);
}

// Test loading output section
void test_load_output_section(void) {
    const char *yaml_content = 
        "output:\n"
        "  style: snp_matrix\n"
        "  finite_output: true\n"
        "  hide_partial_snp: true\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL('s', config.output_style);
    TEST_ASSERT_EQUAL(1, config.finite_output_flag);
    TEST_ASSERT_EQUAL(1, config.hide_partial_snp);
    TEST_ASSERT_EQUAL(1, config.has_output);
}

// Test loading tskit output parameters
void test_load_tskit_output(void) {
    const char *yaml_content = 
        "output:\n"
        "  tskit:\n"
        "    enabled: true\n"
        "    filename: test_output.trees\n"
        "    minimal: false\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(1, config.tskit_output_mode);
    TEST_ASSERT_EQUAL_STRING("test_output.trees", config.tskit_output_filename);
    TEST_ASSERT_EQUAL(0, config.minimal_tree_seq);
    TEST_ASSERT_EQUAL(1, config.has_output);
}

// Test loading events section with population size change
void test_load_population_size_change_event(void) {
    const char *yaml_content = 
        "events:\n"
        "  - type: population_size_change\n"
        "    time: 0.05\n"
        "    population: 0\n"
        "    new_size: 0.5\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(1, config.num_explicit_events);
    TEST_ASSERT_EQUAL('g', config.explicit_events[0].type);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.05, config.explicit_events[0].time);
    TEST_ASSERT_EQUAL(0, config.explicit_events[0].popID);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.5, config.explicit_events[0].popnSize);
}

// Test loading migration rate event
void test_load_migration_rate_event(void) {
    const char *yaml_content = 
        "events:\n"
        "  - type: migration_rate_change\n"
        "    time: 0.1\n"
        "    source: 0\n"
        "    destination: 1\n"
        "    rate: 2.5\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(1, config.num_explicit_events);
    TEST_ASSERT_EQUAL('m', config.explicit_events[0].type);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.1, config.explicit_events[0].time);
    TEST_ASSERT_EQUAL(0, config.explicit_events[0].popID);
    TEST_ASSERT_EQUAL(1, config.explicit_events[0].popID2);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 2.5, config.explicit_events[0].popnSize);
}

// Test loading population split event
void test_load_population_split_event(void) {
    const char *yaml_content = 
        "events:\n"
        "  - type: population_split\n"
        "    time: 0.2\n"
        "    derived: 1\n"
        "    ancestral: 0\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(1, config.num_explicit_events);
    TEST_ASSERT_EQUAL('p', config.explicit_events[0].type);
    TEST_ASSERT_FLOAT_WITHIN(0.0001, 0.2, config.explicit_events[0].time);
    TEST_ASSERT_EQUAL(1, config.explicit_events[0].popID);
    TEST_ASSERT_EQUAL(0, config.explicit_events[0].popID2);
}

// Test loading demes file specification
void test_load_demes_file(void) {
    const char *yaml_content = 
        "demes: test_demes.yaml\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_EQUAL(1, config.use_demes);
    TEST_ASSERT_EQUAL_STRING("test_demes.yaml", config.demes_file);
}

// Test loading empty YAML file
void test_load_empty_yaml(void) {
    const char *yaml_content = "";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    // All values should remain at defaults
    TEST_ASSERT_EQUAL(0, config.sample_size);
    TEST_ASSERT_EQUAL(0, config.has_seed);
    TEST_ASSERT_EQUAL(0, config.has_genetics);
}

// Test loading invalid YAML file
void test_load_invalid_yaml(void) {
    const char *yaml_content = 
        "simulation:\n"
        "  sample_size: [\n"  // Invalid - unclosed bracket
        "  num_replicates: 10\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_NOT_EQUAL(0, result);
}

// Test loading non-existent file
void test_load_nonexistent_file(void) {
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile("/tmp/this_file_does_not_exist.yaml", &config);
    
    TEST_ASSERT_NOT_EQUAL(0, result);
}

// Test loading complete configuration
void test_load_complete_configuration(void) {
    const char *yaml_content = 
        "simulation:\n"
        "  sample_size: 50\n"
        "  num_replicates: 100\n"
        "  num_sites: 200000\n"
        "  seed: [11111, 22222]\n"
        "  effective_popn_size: 5000000\n"
        "\n"
        "genetics:\n"
        "  mutation_rate: 20.0\n"
        "  recombination_rate: 10.0\n"
        "  gene_conversion:\n"
        "    rate: 5.0\n"
        "    tract_length: 300\n"
        "\n"
        "populations:\n"
        "  sample_sizes: [25, 25]\n"
        "\n"
        "selection:\n"
        "  alpha: 50.0\n"
        "  sweep_site: 0.4\n"
        "  sweep_mode: stochastic\n"
        "  tau: 0.02\n"
        "\n"
        "output:\n"
        "  style: haplotype\n"
        "  tskit:\n"
        "    enabled: true\n"
        "    filename: output.trees\n"
        "    minimal: true\n"
        "\n"
        "events:\n"
        "  - type: population_size_change\n"
        "    time: 0.01\n"
        "    population: 0\n"
        "    new_size: 0.8\n"
        "  - type: migration_rate_change\n"
        "    time: 0.02\n"
        "    source: 0\n"
        "    destination: 1\n"
        "    rate: 1.0\n";
    
    createYAMLFile(yaml_content);
    
    SimulationConfig config;
    initializeDefaultConfig(&config);
    int result = loadConfigFile(tempFilename, &config);
    
    TEST_ASSERT_EQUAL(0, result);
    
    // Verify all sections were loaded
    TEST_ASSERT_EQUAL(50, config.sample_size);
    TEST_ASSERT_EQUAL(100, config.num_replicates);
    TEST_ASSERT_EQUAL(200000, config.num_sites);
    TEST_ASSERT_EQUAL(11111, config.seed1);
    TEST_ASSERT_EQUAL(22222, config.seed2);
    TEST_ASSERT_EQUAL(1, config.has_seed);
    
    TEST_ASSERT_EQUAL(20.0, config.theta);
    TEST_ASSERT_EQUAL(10.0, config.rho);
    TEST_ASSERT_EQUAL(5.0, config.gamma);
    TEST_ASSERT_EQUAL(300, config.gc_mean);
    TEST_ASSERT_EQUAL(1, config.has_genetics);
    
    TEST_ASSERT_EQUAL(2, config.npops);
    TEST_ASSERT_EQUAL(25, config.sample_sizes[0]);
    TEST_ASSERT_EQUAL(25, config.sample_sizes[1]);
    TEST_ASSERT_EQUAL(1, config.has_populations);
    
    TEST_ASSERT_EQUAL(50.0, config.alpha);
    TEST_ASSERT_EQUAL(0.4, config.sweep_site);
    TEST_ASSERT_EQUAL('s', config.sweep_mode);
    TEST_ASSERT_EQUAL(0.02, config.tau);
    TEST_ASSERT_EQUAL(1, config.has_selection);
    
    TEST_ASSERT_EQUAL('h', config.output_style);
    TEST_ASSERT_EQUAL(1, config.tskit_output_mode);
    TEST_ASSERT_EQUAL_STRING("output.trees", config.tskit_output_filename);
    TEST_ASSERT_EQUAL(1, config.minimal_tree_seq);
    TEST_ASSERT_EQUAL(1, config.has_output);
    
    TEST_ASSERT_EQUAL(2, config.num_explicit_events);
}

// Main test runner - only if not using external runner
#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_initialize_default_config);
    RUN_TEST(test_load_simulation_section);
    RUN_TEST(test_load_seed_scalar);
    RUN_TEST(test_load_seed_array);
    RUN_TEST(test_load_genetics_section);
    RUN_TEST(test_load_gene_conversion);
    RUN_TEST(test_load_gene_conversion_crossover_ratio);
    RUN_TEST(test_load_populations_section);
    RUN_TEST(test_load_selection_section);
    RUN_TEST(test_load_soft_sweep);
    RUN_TEST(test_load_output_section);
    RUN_TEST(test_load_tskit_output);
    RUN_TEST(test_load_population_size_change_event);
    RUN_TEST(test_load_migration_rate_event);
    RUN_TEST(test_load_population_split_event);
    RUN_TEST(test_load_demes_file);
    RUN_TEST(test_load_empty_yaml);
    RUN_TEST(test_load_invalid_yaml);
    RUN_TEST(test_load_nonexistent_file);
    RUN_TEST(test_load_complete_configuration);
    
    return UNITY_END();
}
#endif
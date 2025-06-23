/**
 * test_phase1_yaml.c - Test minimal Phase 1 parameters via YAML loading
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../src/core/params.h"
#include "../src/core/yaml_loader.h"

void test_yaml_loading() {
    printf("Testing Phase 1 parameters via YAML loading...\n");
    
    SimulationParams *params = params_create();
    assert(params != NULL);
    
    // Load test YAML file
    int ret = yaml_load_params(params, "test/configs/test_phase1.yaml");
    assert(ret == 0);
    
    // Verify basic parameters
    assert(params->core.total_samples == 50);
    assert(params->core.num_replicates == 10);
    assert(params->core.num_sites == 10000);
    
    // Verify Ne parameter
    assert(params->core.Ne == 500000.0);
    printf("  ✓ Ne parameter loaded correctly: %.0f\n", params->core.Ne);
    
    // Verify recurrent sweep parameters
    assert(params->selection.recurrent_sweep == true);
    assert(params->selection.recurrent_sweep_rate == 0.001);
    printf("  ✓ Recurrent sweep parameters loaded correctly\n");
    
    // Verify mask parameter
    assert(params->output.mask == 1);
    printf("  ✓ Mask parameter loaded correctly: %d\n", params->output.mask);
    
    // Verify other selection parameters
    assert(params->selection.alpha == 100.0);
    assert(params->selection.sweep_position == 0.5);
    assert(params->selection.sweep_mode == SWEEP_STOCHASTIC);
    
    params_destroy(params);
}

int main() {
    printf("=== Testing Phase 1 YAML Loading ===\n\n");
    
    test_yaml_loading();
    
    printf("\n=== All YAML tests passed! ===\n");
    
    return 0;
}
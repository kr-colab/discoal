/**
 * test_minimal_phase1.c - Test minimal Phase 1 parameter additions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../src/core/params.h"

void test_recurrent_sweep_rate() {
    printf("Testing recurrent sweep rate parameter...\n");
    
    SimulationParams *params = params_create();
    assert(params != NULL);
    
    // Test default value
    assert(params->selection.recurrent_sweep_rate == 0.0);
    assert(params->selection.recurrent_sweep == false);
    
    // Simulate -R flag
    const char *argv[] = {"discoal", "10", "1", "1000", "-R", "0.001"};
    int argc = 6;
    
    int ret = params_load_from_args(params, argc, argv);
    assert(ret == 0);
    
    // Verify -R flag sets both fields
    assert(params->selection.recurrent_sweep == true);
    assert(params->selection.recurrent_sweep_rate == 0.001);
    assert(params->selection.sweep_mode == SWEEP_STOCHASTIC);
    
    printf("  ✓ Recurrent sweep rate correctly set\n");
    
    params_destroy(params);
}

void test_mask_parameter() {
    printf("Testing mask parameter...\n");
    
    SimulationParams *params = params_create();
    assert(params != NULL);
    
    // Test default value
    assert(params->output.mask == 0);
    
    printf("  ✓ Mask parameter initialized correctly\n");
    
    params_destroy(params);
}

void test_ne_parameter() {
    printf("Testing Ne parameter...\n");
    
    SimulationParams *params = params_create();
    assert(params != NULL);
    
    // Test default value
    assert(params->core.Ne == 1.0);
    
    // Simulate -N flag
    const char *argv[] = {"discoal", "10", "1", "1000", "-N", "500000"};
    int argc = 6;
    
    int ret = params_load_from_args(params, argc, argv);
    assert(ret == 0);
    
    assert(params->core.Ne == 500000.0);
    
    printf("  ✓ Ne parameter correctly set\n");
    
    params_destroy(params);
}

int main() {
    printf("=== Testing Minimal Phase 1 Parameter Additions ===\n\n");
    
    test_recurrent_sweep_rate();
    test_mask_parameter();
    test_ne_parameter();
    
    printf("\n=== All tests passed! ===\n");
    
    return 0;
}
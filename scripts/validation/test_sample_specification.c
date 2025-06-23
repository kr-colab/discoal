/**
 * test_sample_specification.c - Test sample specification with Demes files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../src/core/params.h"
#include "../../src/core/yaml_loader.h"
#include "../../src/core/demes_loader.h"

void print_sample_distribution(SimulationParams *params) {
    printf("Sample distribution across %d populations:\n", params->core.num_populations);
    for (int i = 0; i < params->core.num_populations; i++) {
        printf("  Population %d: %d samples\n", i, params->core.sample_sizes[i]);
    }
    printf("  Total samples: %d\n", params->core.total_samples);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <config.yaml>\n", argv[0]);
        return 1;
    }
    
    /* Create parameters */
    SimulationParams *params = params_create();
    if (!params) {
        fprintf(stderr, "Failed to create parameters\n");
        return 1;
    }
    
    /* Load configuration with Demes and sample specs */
    int ret = yaml_load_params_with_demes(params, argv[1]);
    if (ret != 0) {
        fprintf(stderr, "Failed to load configuration from '%s'\n", argv[1]);
        params_destroy(params);
        return 1;
    }
    
    /* Print results */
    printf("Configuration loaded successfully!\n\n");
    
    /* Print basic info */
    printf("Simulation parameters:\n");
    printf("  Replicates: %d\n", params->core.num_replicates);
    printf("  Sites: %d\n", params->core.num_sites);
    printf("\n");
    
    /* Print sample distribution */
    print_sample_distribution(params);
    printf("\n");
    
    /* Print sample specifications if present */
    if (params->demographics.sample_specs) {
        printf("Sample specifications:\n");
        for (int i = 0; i < params->demographics.num_sample_specs; i++) {
            SampleSpec *spec = &params->demographics.sample_specs[i];
            printf("  %s: %d samples at time %.4f\n", 
                   spec->population, spec->size, spec->time);
        }
    }
    
    /* Clean up */
    params_destroy(params);
    
    return 0;
}
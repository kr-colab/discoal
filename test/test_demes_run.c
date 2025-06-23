/**
 * test_demes_run.c - Test program to verify Demes loading functionality
 */

#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "demes_loader.h"

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <demes_file.yaml>\n", argv[0]);
        return 1;
    }
    
    SimulationParams *params = params_create();
    if (!params) {
        fprintf(stderr, "Failed to create parameters\n");
        return 1;
    }
    
    int ret = demes_load_demographics(params, argv[1]);
    if (ret != 0) {
        fprintf(stderr, "Failed to load Demes file: %s\n", argv[1]);
        params_destroy(params);
        return 1;
    }
    
    printf("Successfully loaded Demes file: %s\n", argv[1]);
    printf("Number of populations: %d\n", params->core.num_populations);
    printf("Population sizes:\n");
    for (int i = 0; i < params->core.num_populations; i++) {
        printf("  Population %d: size = %.6f (relative to reference)\n", 
               i, params->demographics.pop_sizes.sizes[i]);
    }
    
    printf("\nDemographic events: %d\n", params->demographics.num_events);
    for (int i = 0; i < params->demographics.num_events; i++) {
        DemographicEvent *event = &params->demographics.events[i];
        printf("  Event %d: time=%.6f, type=", i, event->time);
        switch (event->type) {
            case EVENT_SIZE_CHANGE:
                printf("SIZE_CHANGE, pop=%d, size=%.6f\n", 
                       event->params.size_change.pop,
                       event->params.size_change.size);
                break;
            case EVENT_SPLIT:
                printf("SPLIT, from=%d, to1=%d, to2=%d, prop=%.3f\n",
                       event->params.split.from,
                       event->params.split.to1,
                       event->params.split.to2,
                       event->params.split.prop);
                break;
            case EVENT_JOIN:
                printf("JOIN, pop1=%d, pop2=%d, dest=%d\n",
                       event->params.join.pop1,
                       event->params.join.pop2,
                       event->params.join.dest);
                break;
            case EVENT_ADMIX:
                printf("ADMIX, from=%d, to=%d, prop=%.3f\n",
                       event->params.admix.from,
                       event->params.admix.to,
                       event->params.admix.prop);
                break;
            case EVENT_MIGRATION_CHANGE:
                printf("MIGRATION_CHANGE, from=%d, to=%d, rate=%.6f\n",
                       event->params.migration.from,
                       event->params.migration.to,
                       event->params.migration.rate);
                break;
            default:
                printf("UNKNOWN\n");
        }
    }
    
    params_destroy(params);
    return 0;
}
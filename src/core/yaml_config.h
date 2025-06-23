/**
 * yaml_config.h - YAML configuration support for discoal parameters
 */

#ifndef YAML_CONFIG_H
#define YAML_CONFIG_H

#include "params.h"
#include <stdio.h>

/* Check if we have YAML support compiled in */
#ifdef WITH_YAML

/* Function to load parameters from a YAML file */
int params_load_from_yaml(SimulationParams *params, const char *filename);

/* Function to save parameters to a YAML file */
int params_save_to_yaml(const SimulationParams *params, const char *filename);

#else

/* Stub functions when YAML support is not compiled in */
static inline int params_load_from_yaml(SimulationParams *params, const char *filename) {
    fprintf(stderr, "Error: YAML support not compiled in. Rebuild with -DWITH_YAML and link with -lyaml\n");
    return -1;
}

static inline int params_save_to_yaml(const SimulationParams *params, const char *filename) {
    fprintf(stderr, "Error: YAML support not compiled in. Rebuild with -DWITH_YAML and link with -lyaml\n");
    return -1;
}

#endif /* WITH_YAML */

#endif /* YAML_CONFIG_H */
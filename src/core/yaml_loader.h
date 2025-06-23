/**
 * yaml_loader.h - YAML configuration support using libyaml
 * 
 * This module provides functions to load simulation parameters
 * from YAML files using the libyaml library.
 */

#ifndef YAML_LOADER_H
#define YAML_LOADER_H

#include "params.h"
#include <yaml.h>

/**
 * Load simulation parameters from a YAML file using libyaml
 * 
 * @param params The SimulationParams structure to populate
 * @param filename Path to the YAML configuration file
 * @return 0 on success, -1 on error
 */
int yaml_load_params(SimulationParams *params, const char *filename);

/**
 * Save simulation parameters to a YAML file
 * 
 * @param params The SimulationParams structure to save
 * @param filename Path to the output YAML file
 * @return 0 on success, -1 on error
 */
int yaml_save_params(const SimulationParams *params, const char *filename);

/**
 * Parse a YAML document into simulation parameters
 * 
 * @param parser The libyaml parser with loaded document
 * @param params The SimulationParams structure to populate
 * @return 0 on success, -1 on error
 */
int yaml_parse_document(yaml_parser_t *parser, SimulationParams *params);

/**
 * Extended parameter loading that supports Demes demographics
 * 
 * This function loads discoal-specific parameters from a YAML file
 * and can reference an external Demes file for demographics.
 * 
 * @param params The SimulationParams structure to populate
 * @param config_file Path to the main configuration file
 * @param demes_file Optional path to Demes file (can be NULL)
 * @return 0 on success, -1 on error
 */
int yaml_load_params_with_demes(SimulationParams *params, 
                               const char *config_file,
                               const char *demes_file);

#endif /* YAML_LOADER_H */
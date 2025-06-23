/**
 * demes_loader.h - Integration with demes-c for demographic model loading
 * 
 * This module provides functions to load demographic models from Demes
 * specification files and convert them to discoal's internal format.
 */

#ifndef DEMES_LOADER_H
#define DEMES_LOADER_H

#include "params.h"
#include "../../extern/demes-c/demes.h"

/**
 * Load demographics from a Demes file into SimulationParams
 * 
 * This function reads a Demes YAML file and populates the demographics
 * portion of the SimulationParams structure.
 * 
 * @param params The SimulationParams structure to populate
 * @param filename Path to the Demes YAML file
 * @return 0 on success, -1 on error
 */
int demes_load_demographics(SimulationParams *params, const char *filename);

/**
 * Convert a Demes graph to discoal demographics
 * 
 * This function takes a parsed Demes graph and converts it to discoal's
 * internal demographic representation.
 * 
 * @param params The SimulationParams structure to populate
 * @param graph The parsed Demes graph
 * @return 0 on success, -1 on error
 */
int demes_convert_graph(SimulationParams *params, struct demes_graph *graph);

/**
 * Extended YAML loading with Demes support
 * 
 * Load discoal parameters from YAML with optional Demes file reference
 * for demographics. The YAML can contain a "demographics" section with
 * a "demes_file" field pointing to an external Demes file.
 * 
 * @param params The SimulationParams structure to populate
 * @param config_file Path to the discoal YAML configuration
 * @return 0 on success, -1 on error
 */
int yaml_load_params_with_demes(SimulationParams *params, const char *config_file);

/**
 * Helper function to convert Demes time units to coalescent units
 * 
 * Note: discoal requires time_units="4N" in Demes files
 * 
 * @param time Time value in Demes units
 * @param graph The Demes graph
 * @return Time in coalescent units, or -1 if time_units is not "4N"
 */
double demes_time_to_coalescent(double time, struct demes_graph *graph);

/**
 * Helper function to map Demes deme names to population indices
 * 
 * @param params The SimulationParams with population info
 * @param deme_name The name of the deme
 * @return Population index, or -1 if not found
 */
int demes_get_population_index(SimulationParams *params, const char *deme_name);

#endif /* DEMES_LOADER_H */
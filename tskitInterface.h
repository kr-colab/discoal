// tskitInterface.h
// Interface for tskit tree sequence recording in discoal
//
// IMPORTANT: Time Scaling Convention
// ----------------------------------
// discoal internally uses time units where the expected coalescence time for two lineages
// in a population of size N is 2N generations. However, msprime and standard coalescent
// theory use time units where the expected coalescence time is 2N generations for a 
// population of diploid size N.
// 
// To maintain compatibility with msprime and other tools, all node times are automatically
// divided by 2.0 when writing to the tskit tree sequence. This scaling is applied
// internally in tskit_add_node(), so users don't need to manually adjust times.
//
// This means:
// - discoal internal time t -> tskit time t/2
// - Tree heights in tskit output will be ~0.5x the internal discoal values
// - This matches msprime's convention for ms-compatible parameters with Ne=0.5

#ifndef __TSKIT_INTERFACE_H__
#define __TSKIT_INTERFACE_H__

#include <tskit.h>
#include "discoal.h"

// Global tskit table collection
extern tsk_table_collection_t *tsk_tables;

// Mapping between discoal node IDs and tskit node IDs
extern tsk_id_t *node_id_map;
extern int node_id_map_capacity;

// Initialize tskit recording
int tskit_initialize(double sequence_length);

// Record a new node (called when creating nodes in discoal)
tsk_id_t tskit_add_node(double time, int population, int is_sample);

// Record edges during coalescence
int tskit_add_edges(tsk_id_t parent, tsk_id_t child, double left, double right);

// Record a site
tsk_id_t tskit_add_site(double position, const char *ancestral_state);

// Record a mutation
int tskit_add_mutation(tsk_id_t site, tsk_id_t node, const char *derived_state);

// Finalize and write tree sequence
int tskit_finalize(const char *filename);

// Clean up tskit resources
void tskit_cleanup(void);

// Map discoal node to tskit node ID
tsk_id_t get_tskit_node_id(rootedNode *node);

// Store mapping from discoal node to tskit node ID
void set_tskit_node_id(rootedNode *node, tsk_id_t tsk_id);
void set_tskit_node_id_at_index(int index, tsk_id_t tsk_id);

// Record all mutations after they've been placed on the tree
int tskit_record_mutations(void);

// Place mutations using node-based algorithm (RNG-compatible with original dropMutations)
int tskit_place_mutations_node_based(double theta);

// Place mutations using edge-based algorithm (similar to msprime's design)
int tskit_place_mutations_edge_based(double theta);

// Place mutations directly on tskit edges (defaults to node-based for compatibility)
int tskit_place_mutations_directly(double theta);

// Populate discoal mutation arrays from tskit data (for ms output compatibility)
int tskit_populate_discoal_mutations(void);

#endif
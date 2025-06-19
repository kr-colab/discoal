// tskitInterface.h
// Interface for tskit tree sequence recording in discoal
//
// IMPORTANT: Time Scaling Convention
// ----------------------------------
// discoal internally uses time units where the expected coalescence time for two lineages
// in a population of size N is 2N generations. This is the standard coalescent time scale.
// 
// Node times are recorded directly in these coalescent units without conversion.
// The time_units attribute is set to "coalescent units (2N generations)" to make this
// scaling explicit.
//
// This means:
// - discoal internal time t -> tskit time t (no conversion)
// - Tree heights represent coalescent time (scaled by 2N)
// - To convert to generations, multiply times by 2N
// - To match msprime conventions (which uses Ne=0.5 for haploid samples), 
//   times would need to be divided by 2

#ifndef __TSKIT_INTERFACE_H__
#define __TSKIT_INTERFACE_H__

#include <tskit.h>
#include "discoal.h"

// Global tskit table collection
extern tsk_table_collection_t *tsk_tables;

// Mapping between discoal node IDs and tskit node IDs
extern tsk_id_t *node_id_map;
extern int node_id_map_capacity;

// Edge buffering for squashing optimization
extern tsk_edge_t *buffered_edges;
extern int num_buffered_edges;
extern int max_buffered_edges;

// Initialize tskit recording
int tskit_initialize(double sequence_length);

// Record a new node (called when creating nodes in discoal)
tsk_id_t tskit_add_node(double time, int population, int is_sample);

// Record edges during coalescence (with buffering/squashing)
int tskit_add_edges(tsk_id_t parent, tsk_id_t child, double left, double right);

// Flush buffered edges (apply squashing and add to table)
int tskit_flush_edges(void);

// Special flush for periodic node sweeping (to track frequency)
int tskit_flush_edges_periodic(void);

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

// Place mutations using edge-based algorithm (similar to msprime's design)
int tskit_place_mutations_edge_based(double theta);

// Place mutations directly on tskit edges (uses edge-based algorithm)
int tskit_place_mutations_directly(double theta);

// Populate discoal mutation arrays from tskit data (for ms output compatibility)
int tskit_populate_discoal_mutations(void);

// Record sweep mutations for nodes that carry them
int tskit_record_sweep_mutations(double sweepSite);

// Store command line for provenance (called from main)
void tskit_store_command_line(int argc, const char **argv);

// Functions to enable/disable edge squashing (for performance comparison)
void tskit_enable_edge_squashing(void);
void tskit_disable_edge_squashing(void);

#endif
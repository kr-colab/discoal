// tskitInterface.h
// Interface for tskit tree sequence recording in discoal

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

#endif
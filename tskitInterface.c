// tskitInterface.c
// Implementation of tskit tree sequence recording for discoal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "tskitInterface.h"

// Global tskit table collection
tsk_table_collection_t *tsk_tables = NULL;

// Mapping between discoal node IDs and tskit node IDs
// We'll use the allNodes array index as the key
tsk_id_t *node_id_map = NULL;
int node_id_map_capacity = 0;

// Initialize tskit recording
int tskit_initialize(double sequence_length) {
    int ret;
    
    // Allocate table collection
    tsk_tables = malloc(sizeof(tsk_table_collection_t));
    if (tsk_tables == NULL) {
        return -1;
    }
    
    // Initialize table collection
    ret = tsk_table_collection_init(tsk_tables, 0);
    if (ret != 0) {
        free(tsk_tables);
        tsk_tables = NULL;
        return ret;
    }
    
    // Set sequence length
    tsk_tables->sequence_length = sequence_length;
    
    // Add population records for all populations in discoal
    // We need to add enough populations to cover the maximum population ID used
    extern int npops;  // From discoal globals
    for (int i = 0; i < npops; i++) {
        ret = tsk_population_table_add_row(&tsk_tables->populations, NULL, 0);
        if (ret < 0) {
            tsk_table_collection_free(tsk_tables);
            free(tsk_tables);
            tsk_tables = NULL;
            return ret;
        }
    }
    
    // Initialize node ID mapping
    node_id_map_capacity = 1000;  // Start with reasonable size
    node_id_map = malloc(sizeof(tsk_id_t) * node_id_map_capacity);
    if (node_id_map == NULL) {
        tsk_table_collection_free(tsk_tables);
        free(tsk_tables);
        tsk_tables = NULL;
        return -1;
    }
    
    // Initialize all mappings to TSK_NULL
    for (int i = 0; i < node_id_map_capacity; i++) {
        node_id_map[i] = TSK_NULL;
    }
    
    return 0;
}

// Ensure node_id_map has enough capacity
static void ensure_node_map_capacity(int required_size) {
    if (required_size >= node_id_map_capacity) {
        int new_capacity = node_id_map_capacity;
        while (new_capacity <= required_size) {
            new_capacity *= 2;
        }
        
        tsk_id_t *new_map = realloc(node_id_map, sizeof(tsk_id_t) * new_capacity);
        if (new_map == NULL) {
            fprintf(stderr, "Error: Failed to reallocate node_id_map\n");
            exit(1);
        }
        
        // Initialize new entries to TSK_NULL
        for (int i = node_id_map_capacity; i < new_capacity; i++) {
            new_map[i] = TSK_NULL;
        }
        
        node_id_map = new_map;
        node_id_map_capacity = new_capacity;
    }
}

// Record a new node
tsk_id_t tskit_add_node(double time, int population, int is_sample) {
    tsk_flags_t flags = is_sample ? TSK_NODE_IS_SAMPLE : 0;
    tsk_id_t node_id;
    
    if (tsk_tables == NULL) {
        return TSK_NULL;
    }
    
    // tskit doesn't allow negative population IDs, so map -1 to TSK_NULL
    tsk_id_t tsk_population = (population < 0) ? TSK_NULL : population;
    
    node_id = tsk_node_table_add_row(&tsk_tables->nodes, 
                                     flags, 
                                     time / 2.0,  // Scale time by 0.5 to match msprime convention
                                     tsk_population, 
                                     TSK_NULL,  // individual
                                     NULL, 0);  // metadata
    
    // Debug output
    if (is_sample) {
        fprintf(stderr, "Added sample node %lld: time=%f (scaled=%f), pop=%d\n", 
                (long long)node_id, time, time / 2.0, tsk_population);
    }
    
    return node_id;
}

// Record edges during coalescence
int tskit_add_edges(tsk_id_t parent, tsk_id_t child, double left, double right) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    ret = tsk_edge_table_add_row(&tsk_tables->edges,
                                 left, right,
                                 parent, child,
                                 NULL, 0);  // metadata
    
    // Debug output
    if (right - left < 0.0001) {
        fprintf(stderr, "WARNING: Very small edge interval: parent=%lld, child=%lld, left=%.10f, right=%.10f, length=%.10f\n", 
                (long long)parent, (long long)child, left, right, right - left);
    }
    
    return ret;
}

// Record a site
tsk_id_t tskit_add_site(double position, const char *ancestral_state) {
    tsk_id_t site_id;
    
    if (tsk_tables == NULL) {
        return TSK_NULL;
    }
    
    site_id = tsk_site_table_add_row(&tsk_tables->sites,
                                     position,
                                     ancestral_state,
                                     strlen(ancestral_state),
                                     NULL, 0);  // metadata
    
    return site_id;
}

// Record a mutation
int tskit_add_mutation(tsk_id_t site, tsk_id_t node, const char *derived_state) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    ret = tsk_mutation_table_add_row(&tsk_tables->mutations,
                                     site,
                                     node,
                                     TSK_NULL,  // parent mutation
                                     TSK_UNKNOWN_TIME,  // time
                                     derived_state,
                                     strlen(derived_state),
                                     NULL, 0);  // metadata
    
    return ret;
}

// Finalize and write tree sequence
int tskit_finalize(const char *filename) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    
    // Sort tables (required for valid tree sequences)
    ret = tsk_table_collection_sort(tsk_tables, NULL, 0);
    if (ret != 0) {
        fprintf(stderr, "Error sorting tables: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    // Note: We're not simplifying because with recombination, it's natural
    // to have multiple roots - different genomic regions may not fully coalesce
    
    // Build index after simplification
    ret = tsk_table_collection_build_index(tsk_tables, 0);
    if (ret != 0) {
        fprintf(stderr, "Error building index: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    // Write to file
    ret = tsk_table_collection_dump(tsk_tables, filename, 0);
    if (ret != 0) {
        fprintf(stderr, "Error writing tree sequence: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    return 0;
}

// Clean up tskit resources
void tskit_cleanup(void) {
    if (tsk_tables != NULL) {
        tsk_table_collection_free(tsk_tables);
        free(tsk_tables);
        tsk_tables = NULL;
    }
    
    if (node_id_map != NULL) {
        free(node_id_map);
        node_id_map = NULL;
        node_id_map_capacity = 0;
    }
}

// Map discoal node to tskit node ID
tsk_id_t get_tskit_node_id(rootedNode *node) {
    // Find the node's index in allNodes array
    extern int totNodeNumber;
    extern rootedNode **allNodes;
    for (int i = 0; i < totNodeNumber; i++) {
        if (allNodes[i] == node) {
            if (i < node_id_map_capacity) {
                // fprintf(stderr, "Found mapping: node %p at index %d -> tskit ID %lld\n", 
                //         (void*)node, i, (long long)node_id_map[i]);
                return node_id_map[i];
            }
            break;
        }
    }
    fprintf(stderr, "WARNING: Could not find mapping for node %p\n", (void*)node);
    return TSK_NULL;
}

// Store mapping from discoal node to tskit node ID using index
void set_tskit_node_id_at_index(int index, tsk_id_t tsk_id) {
    ensure_node_map_capacity(index);
    node_id_map[index] = tsk_id;
    // fprintf(stderr, "Stored mapping at index %d -> tskit ID %lld\n", 
    //         index, (long long)tsk_id);
}

// Store mapping from discoal node to tskit node ID
void set_tskit_node_id(rootedNode *node, tsk_id_t tsk_id) {
    // Find the node's index in allNodes array
    extern int totNodeNumber;
    extern rootedNode **allNodes;
    for (int i = 0; i < totNodeNumber; i++) {
        if (allNodes[i] == node) {
            ensure_node_map_capacity(i);
            node_id_map[i] = tsk_id;
            // fprintf(stderr, "Stored mapping: node %p at index %d -> tskit ID %lld\n", 
            //         (void*)node, i, (long long)tsk_id);
            return;
        }
    }
    fprintf(stderr, "WARNING: Could not find node %p in allNodes array\n", (void*)node);
}

// Record all mutations after they've been placed on the tree
int tskit_record_mutations(void) {
    if (tsk_tables == NULL) {
        return -1;
    }
    
    extern int totNodeNumber;
    extern rootedNode **allNodes;
    extern int nSites;
    
    // First, collect all unique mutation positions and create sites
    // We'll use a simple array to track sites we've already created
    double *site_positions = malloc(sizeof(double) * MAXMUTS);
    tsk_id_t *site_ids = malloc(sizeof(tsk_id_t) * MAXMUTS);
    int num_sites = 0;
    
    if (site_positions == NULL || site_ids == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for mutation recording\n");
        return -1;
    }
    
    // Iterate through all nodes to find unique mutation positions
    for (int i = 0; i < totNodeNumber; i++) {
        rootedNode *node = allNodes[i];
        if (node == NULL) continue;
        
        for (int j = 0; j < node->mutationNumber; j++) {
            double mut_pos = node->muts[j];
            
            // Convert from discoal's [0,1] to actual position [0, sequence_length)
            double actual_pos = mut_pos * tsk_tables->sequence_length;
            
            // Check if we've already created a site for this position
            int found = 0;
            tsk_id_t site_id = TSK_NULL;
            for (int k = 0; k < num_sites; k++) {
                if (site_positions[k] == actual_pos) {
                    found = 1;
                    site_id = site_ids[k];
                    break;
                }
            }
            
            // If not found, create a new site
            if (!found) {
                site_id = tskit_add_site(actual_pos, "0");  // ancestral state is "0"
                if (site_id < 0) {
                    fprintf(stderr, "Error: Failed to add site at position %f\n", actual_pos);
                    free(site_positions);
                    free(site_ids);
                    return -1;
                }
                site_positions[num_sites] = actual_pos;
                site_ids[num_sites] = site_id;
                num_sites++;
            }
            
            // Now add the mutation at this site on this node
            tsk_id_t tskit_node_id = get_tskit_node_id(node);
            if (tskit_node_id != TSK_NULL) {
                int ret = tskit_add_mutation(site_id, tskit_node_id, "1");  // derived state is "1"
                if (ret < 0) {
                    fprintf(stderr, "Error: Failed to add mutation at site %lld on node %lld\n", 
                            (long long)site_id, (long long)tskit_node_id);
                    free(site_positions);
                    free(site_ids);
                    return -1;
                }
            }
        }
    }
    
    free(site_positions);
    free(site_ids);
    
    // fprintf(stderr, "Successfully recorded %d sites with mutations in tskit\n", num_sites);
    return 0;
}